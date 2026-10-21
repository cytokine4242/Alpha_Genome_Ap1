"""Shared helpers for the AlphaGenome AP-1 perturbation jobs: sequence windows, edits, model wrapper, readouts.

Coordinates are 0-based half-open (BED-like), as in the TOBIAS table. Sequences are numpy uint8 codes
A=0 C=1 G=2 T=3 N=4 (N encodes to an all-zero one-hot row).
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LUT = np.full(256, 4, dtype=np.uint8)
for i, ch in enumerate("ACGT"):
    LUT[ord(ch)] = i
    LUT[ord(ch.lower())] = i
BASES = np.array(list("ACGTN"))
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
TRANSVERSION = np.array([1, 0, 3, 2, 4], dtype=np.uint8)  # A<->C, G<->T, used for the targeted core edit
CORE_RE = re.compile("TGA[CG]TCA")
ONEHOT = np.eye(5, 4, dtype=np.float32)  # row 4 (N) is all zeros


def load_config(path=None):
    path = Path(path) if path else ROOT / "config.json"
    with open(path) as f:
        cfg = json.load(f)
    if os.environ.get("AG_BATCH_SIZE"):  # set by submit_all.sh (BATCH=N)
        cfg["model"]["batch_size"] = int(os.environ["AG_BATCH_SIZE"])
    return cfg


def encode(seq: str) -> np.ndarray:
    return LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def decode(codes: np.ndarray) -> str:
    return "".join(BASES[codes])


def onehot(codes: np.ndarray) -> np.ndarray:
    return ONEHOT[codes]


# ---------------------------------------------------------------- genome / annotation
class Genome:
    def __init__(self, fasta):
        from pyfaidx import Fasta
        self.fa = Fasta(str(fasta), as_raw=True, sequence_always_upper=True)

    def length(self, chrom):
        return len(self.fa[chrom])

    def fetch(self, chrom, start, end):
        if chrom not in self.fa or start < 0 or end > len(self.fa[chrom]):
            return None
        return encode(self.fa[chrom][start:end])

    def window(self, chrom, mid, size):
        """Window of `size` bp centred on `mid`. Returns (codes, window_start) or (None, None) off the chromosome."""
        start = int(mid) - size // 2
        codes = self.fetch(chrom, start, start + size)
        return (codes, start) if codes is not None and len(codes) == size else (None, None)


class TssIndex:
    """Nearest annotated transcription start site (any gene, any strand), from a GTF with transcript or gene lines."""

    def __init__(self, gtf):
        import gzip
        import pandas as pd
        opener = gzip.open if str(gtf).endswith(".gz") else open
        rows = []
        with opener(gtf, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.split("\t")
                if len(p) < 9 or p[2] != "transcript":
                    continue
                tss = int(p[3]) - 1 if p[6] == "+" else int(p[4]) - 1
                rows.append((p[0], tss))
        df = pd.DataFrame(rows, columns=["chrom", "tss"]).drop_duplicates()
        self.by_chrom = {c: np.sort(g["tss"].to_numpy()) for c, g in df.groupby("chrom")}

    def nearest(self, chrom, pos):
        a = self.by_chrom.get(chrom)
        if a is None or len(a) == 0:
            return None
        i = np.searchsorted(a, pos)
        cand = a[max(i - 1, 0): i + 1]
        return int(cand[np.argmin(np.abs(cand - pos))])


# ---------------------------------------------------------------- sequence edits
def shuffle_region(codes, a, b, rng):
    out = codes.copy()
    seg = out[a:b].copy()
    rng.shuffle(seg)
    out[a:b] = seg
    return out


def random_mutate_region(codes, a, b, rng):
    """Replace every base in [a, b) by a different random base."""
    out = codes.copy()
    out[a:b] = (out[a:b] + rng.integers(1, 4, size=b - a)).astype(np.uint8) % 4
    return out


def dinuc_shuffle_codes(seq, rng):
    """Altschul-Erikson dinucleotide-preserving shuffle of a 0-3 coded sequence (same first/last base)."""
    s = [int(x) for x in seq]
    n = len(s)
    if n < 4:
        return np.array(s, dtype=np.uint8)
    edges = {v: [] for v in range(4)}
    for a, b in zip(s[:-1], s[1:]):
        edges[a].append(b)
    last = s[-1]
    while True:  # random last-exit edges must form a tree rooted at `last`
        last_edge = {v: edges[v][int(rng.integers(len(edges[v])))] for v in range(4) if v != last and edges[v]}
        ok = True
        for v in last_edge:
            cur, steps = v, 0
            while cur != last and steps <= 4:
                if cur not in last_edge:
                    ok = False
                    break
                cur = last_edge[cur]
                steps += 1
            if cur != last:
                ok = False
            if not ok:
                break
        if ok:
            break
    for v in range(4):
        lst = list(edges[v])
        if v in last_edge:
            lst.remove(last_edge[v])
            rng.shuffle(lst)
            lst.append(last_edge[v])
        else:
            rng.shuffle(lst)
        edges[v] = lst
    ptr = {v: 0 for v in range(4)}
    out = [s[0]]
    prev = s[0]
    for _ in range(n - 1):
        nxt = edges[prev][ptr[prev]]
        ptr[prev] += 1
        out.append(nxt)
        prev = nxt
    return np.array(out, dtype=np.uint8)


def dinuc_region(codes, a, b, rng):
    out = codes.copy()
    out[a:b] = dinuc_shuffle_codes(out[a:b], rng)
    return out


def core_interval(codes, motif, motif_len=11):
    """(core_start, core_end, found) of the TGA[CG]TCA core inside the motif interval (offsets in the window).
    Falls back to the central 7 bp if the motif holds no exact core."""
    a, b = motif
    m = CORE_RE.search(decode(codes[a:b]))
    if m:
        return a + m.start(), a + m.end(), True
    mid = (a + b) // 2
    return mid - 3, mid + 4, False


def core_edit(codes, motif):
    c0, c1, found = core_interval(codes, motif)
    out = codes.copy()
    out[c0:c1] = TRANSVERSION[out[c0:c1]]
    return out, found


def flank_edit(codes, regions, width, rng):
    """Mutate the flanks but keep the AP-1 core. motif flanks: motif minus core. 170 flanks: 170 bp window minus motif."""
    out = codes.copy()
    c0, c1, found = core_interval(codes, regions["motif"])
    if width == "motif":
        for a, b in ((regions["motif"][0], c0), (c1, regions["motif"][1])):
            if b > a:
                out = random_mutate_region(out, a, b, rng)
    else:
        m0, m1 = regions["motif"]
        a, b = regions["w170"]
        for lo, hi in ((a, m0), (m1, b)):
            if hi > lo:
                out = random_mutate_region(out, lo, hi, rng)
    return out, found


def has_n(codes, a, b):
    return bool((codes[a:b] == 4).any())


# conditions: (name, kind, width). width in motif / 170 / peak
def build_conditions(reps):
    conds = []
    for kind in ("shuffle", "random", "dinuc"):
        for width in ("motif", "170", "peak"):
            conds.append((f"{kind}@{width}", kind, width, reps[kind]))
    conds.append(("core@motif", "core", "motif", 1))
    conds.append(("flank@motif", "flank", "motif", reps["flank"]))
    conds.append(("flank@170", "flank", "170", reps["flank"]))
    return conds


def apply_condition(codes, regions, kind, width, rng):
    """Return edited window, or None if the edit region holds N. regions: dict of (a, b) offsets in the window."""
    if kind == "core":
        return core_edit(codes, regions["motif"])[0]
    if kind == "flank":
        return flank_edit(codes, regions, width, rng)[0]
    key = {"motif": "motif", "170": "w170", "peak": "peak"}[width]
    a, b = regions[key]
    if b - a < 2 or has_n(codes, a, b):
        return None
    if kind == "shuffle":
        return shuffle_region(codes, a, b, rng)
    if kind == "random":
        return random_mutate_region(codes, a, b, rng)
    if kind == "dinuc":
        return dinuc_region(codes, a, b, rng)
    raise ValueError(kind)


def site_regions(row, win_start, size, cfg):
    """Window-relative intervals for a site row (needs start, end, motif_mid, peak_start, peak_end)."""
    p = cfg["perturb"]
    clip = lambda a, b: (max(int(a) - win_start, 0), min(int(b) - win_start, size))
    mid = int(row["motif_mid"])
    half = p["w170"] // 2
    peak = clip(row["peak_start"], row["peak_end"])
    return {
        "motif": clip(row["start"], row["end"]),
        "w170": clip(mid - half, mid + half),
        "peak": peak,
        "mid": mid - win_start,
    }


# ---------------------------------------------------------------- model
def load_model(cfg):
    import torch
    from alphagenome_pytorch import AlphaGenome
    path = cfg["paths"]["weights"]
    device = cfg["model"]["device"] if torch.cuda.is_available() else "cpu"
    if str(path).endswith(".safetensors"):
        from safetensors.torch import load_file
        model = AlphaGenome()
        res = model.load_state_dict(load_file(path), strict=False)
        print(f"loaded {path}: {len(res.missing_keys)} missing, {len(res.unexpected_keys)} unexpected keys", flush=True)
        model = model.to(device).eval()
    else:
        model = AlphaGenome.from_pretrained(path, device=device)
    return model, device


class AG:
    """AlphaGenome mouse predictor returning averaged-track ATAC (1 bp) and RNA (128 bp) per trackset."""

    def __init__(self, cfg):
        import torch
        from alphagenome_pytorch.named_outputs import TrackMetadataCatalog
        self.torch = torch
        self.cfg = cfg
        self.size = cfg["model"]["window"]
        self.org = cfg["model"]["organism_index"]
        self.bs = cfg["model"]["batch_size"]
        self.model, self.device = load_model(cfg)
        self.model.set_track_metadata_catalog(TrackMetadataCatalog.load_builtin("mouse"))
        self.tracksets = self._resolve_tracks()

    def _named(self, x):
        with self.torch.no_grad():
            return self.model.predict(x, organism_index=self.org, named_outputs=True)

    def track_table(self):
        x = self.torch.zeros((1, self.size, 4), dtype=self.torch.float32, device=self.device)
        named = self._named(x)
        rows = []
        for head, res in (("atac", 1), ("rna_seq", 128), ("dnase", 128), ("cage", 128)):
            try:
                for i, t in enumerate(getattr(named, head)[res].tracks):
                    d = t.to_dict()
                    d.update(head=head, index=i)
                    rows.append(d)
            except Exception as e:  # head missing for mouse
                print(f"skip head {head}: {e}")
        return rows

    def _resolve_tracks(self):
        rows = self.track_table()
        out = {}
        for name, spec in self.cfg["tracksets"].items():
            out[name] = {}
            for head, pat in spec.items():
                rx = re.compile(pat, re.I)
                idx = [r["index"] for r in rows if r["head"] == head and rx.search(
                    " | ".join(str(r.get(k, "")) for k in ("biosample_name", "ontology_curie", "track_name", "assay_title")))]
                names = [r.get("track_name") for r in rows if r["head"] == head and r["index"] in idx]
                print(f"trackset {name} {head}: {len(idx)} tracks {names}", flush=True)
                if not idx:
                    raise SystemExit(f"trackset {name}/{head}: pattern '{pat}' matches no mouse track; edit config.json "
                                     "(see results/mouse_tracks.csv from 00_check_env)")
                out[name][head] = idx
        return out

    def predict(self, seqs):
        """seqs: list of window codes. Returns list of per-sequence dicts {trackset: {'atac': (L,), 'rna': (L/128,)}}."""
        torch = self.torch
        res = []
        for i in range(0, len(seqs), self.bs):
            x = torch.from_numpy(np.stack([onehot(s) for s in seqs[i:i + self.bs]])).to(self.device)
            named = self._named(x)
            atac = named.atac[1].tensor
            rna = named.rna_seq[128].tensor
            per = {}
            for ts, spec in self.tracksets.items():
                per[ts] = (atac[..., spec["atac"]].mean(-1).float().cpu().numpy(),
                           rna[..., spec["rna_seq"]].mean(-1).float().cpu().numpy())
            for b in range(x.shape[0]):
                res.append({ts: {"atac": per[ts][0][b], "rna": per[ts][1][b]} for ts in per})
        return res


def readout(pred, regions, tss_rel, cfg):
    """Readouts for one prediction (one trackset). regions: motif / w170 / peak (a, b) and mid, window offsets."""
    p = cfg["perturb"]
    atac, rna = pred["atac"], pred["rna"]
    out = {}
    for k, key in (("atac_motif", "motif"), ("atac_170", "w170"), ("atac_peak", "peak")):
        a, b = regions[key]
        out[k] = float(atac[a:b].mean()) if b > a else np.nan
    mid, loc = regions["mid"], p["local_bp"]
    lo, hi = max(mid - loc, 0) // 128, min(mid + loc, len(rna) * 128) // 128 + 1
    out["rna_mean50kb"] = float(rna[lo:hi].mean())
    out["rna_max50kb"] = float(rna[lo:hi].max())
    if tss_rel is not None and 0 <= tss_rel < len(rna) * 128:
        t, pr = tss_rel, p["tss_prox_bp"]
        out["rna_prox2kb"] = float(rna[max(t - pr, 0) // 128: min(t + pr, len(rna) * 128) // 128 + 1].mean())
    else:
        out["rna_prox2kb"] = np.nan
    return out


def append_rows(path, rows):
    import pandas as pd
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, mode="a", header=not path.exists(), index=False)


def done_ids(path, col="site_id"):
    import pandas as pd
    path = Path(path)
    return set(pd.read_csv(path, usecols=[col])[col]) if path.exists() else set()
