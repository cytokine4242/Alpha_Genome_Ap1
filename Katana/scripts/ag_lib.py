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
CORE_RE = re.compile("TGA[CG]TCA")
ONEHOT = np.eye(5, 4, dtype=np.float32)  # row 4 (N) is all zeros


def load_config(path=None):
    path = Path(path or os.environ.get("AG_CONFIG") or ROOT / "config.json")
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


HALF_SITES = ("TGAC", "GTCA", "TGAG", "CTCA")  # AP-1 half-sites TGA[CG] and their reverse complements


def core_edit(codes, motif, rng):
    """Disrupt the AP-1 core: every core base is replaced by a different random base, redrawn until no AP-1 half-site
    (TGAC/TGAG, or their reverse complements GTCA/CTCA) is left in the core or at its junction with the flanks.
    (A fixed transversion A<->C, G<->T would turn TGAGTCA into GTCTGAC, which still holds the half-site TGAC.)"""
    c0, c1, found = core_interval(codes, motif)
    for _ in range(200):
        out = codes.copy()
        out[c0:c1] = (out[c0:c1] + rng.integers(1, 4, size=c1 - c0)).astype(np.uint8) % 4
        ctx = decode(out[max(c0 - 3, 0):c1 + 3])
        if not any(h in ctx for h in HALF_SITES):
            return out, found
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


# conditions: (name, kind, width, reps). width in motif / 170 / peak. reps = min(reps[kind], caps[width]): large edits have
# large, stable effects and need fewer replicates than motif-width edits.
def build_conditions(reps, caps=None, kinds=None):
    caps = caps or {}
    conds = []
    for kind in ("shuffle", "random", "dinuc"):
        for width in ("motif", "170", "peak"):
            conds.append((f"{kind}@{width}", kind, width, min(reps[kind], caps.get(width, 10**6))))
    conds.append(("core@motif", "core", "motif", reps.get("core", 3)))
    conds.append(("flank@motif", "flank", "motif", min(reps["flank"], caps.get("motif", 10**6))))
    conds.append(("flank@170", "flank", "170", min(reps["flank"], caps.get("170", 10**6))))
    return [c for c in conds if not kinds or c[1] in kinds]


def apply_condition(codes, regions, kind, width, rng):
    """Return edited window, or None if the edit region holds N. regions: dict of (a, b) offsets in the window."""
    if kind == "core":
        return core_edit(codes, regions["motif"], rng)[0]
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
    if hasattr(AlphaGenome, "from_pretrained"):
        try:
            model = AlphaGenome.from_pretrained(path, device=device)
            print(f"loaded {path} with AlphaGenome.from_pretrained", flush=True)
            return model, device
        except Exception as e:  # e.g. a .safetensors file the package cannot read directly
            print(f"from_pretrained failed ({type(e).__name__}: {str(e)[:300]}); trying a raw state dict", flush=True)
    if str(path).endswith(".safetensors"):
        from safetensors.torch import load_file
        state = load_file(path)
    else:
        state = torch.load(path, map_location="cpu", weights_only=True)
    model = AlphaGenome()
    res = model.load_state_dict(state, strict=False)
    print(f"loaded {path}: {len(res.missing_keys)} missing, {len(res.unexpected_keys)} unexpected keys", flush=True)
    if res.missing_keys or res.unexpected_keys:
        print("missing (first 8):", res.missing_keys[:8], "\nunexpected (first 8):", res.unexpected_keys[:8], flush=True)
        raise SystemExit("weights do not match the installed alphagenome_pytorch model (key names differ). "
                         "Install the version the weights were made for; run scripts/00_inspect_weights.py for details.")
    return model.to(device).eval(), device


ACC_HEADS = {"atac": "atac", "dnase": "dnase"}     # accessibility, read at 1 bp; columns <prefix>_motif/_170/_peak
EXPR_HEADS = {"rna_seq": "rna", "cage": "cage"}      # expression, read at 128 bp; columns <prefix>_prox2kb/_mean50kb/_max50kb
TRACK_FIELDS = ("biosample_name", "ontology_curie", "track_name", "assay_title", "biosample_life_stage", "strand", "biosample_type")


def track_text(row):
    return " | ".join(str(row.get(k, "")) for k in TRACK_FIELDS)


class AG:
    """AlphaGenome mouse predictor: per-trackset averaged tracks. ATAC/DNase at 1 bp, RNA-seq/CAGE at 128 bp.

    config tracksets: {name: {head: regex}}, head in atac / dnase / rna_seq / cage. The regex (case-insensitive) is searched in
    "biosample_name | ontology_curie | track_name | assay_title | life_stage | strand | biosample_type".
    Only the heads used by the tracksets are computed. strict=False reports unmatched tracksets instead of exiting."""

    def __init__(self, cfg, strict=True):
        import torch
        from alphagenome_pytorch.named_outputs import TrackMetadataCatalog
        self.torch = torch
        self.cfg = cfg
        self.size = cfg["model"]["window"]
        self.org = cfg["model"]["organism_index"]
        self.bs = cfg["model"]["batch_size"]
        self.model, self.device = load_model(cfg)
        self.model.set_track_metadata_catalog(TrackMetadataCatalog.load_builtin("mouse"))
        self.heads = None
        self.tracksets = self._resolve_tracks(strict)
        self.heads = sorted({h for spec in self.tracksets.values() for h in spec}) or None

    def _named(self, x, all_heads=False):
        kw = {}
        if not all_heads and self.heads:
            kw["heads"] = tuple(self.heads)
        with self.torch.no_grad():
            return self.model.predict(x, organism_index=self.org, named_outputs=True, **kw)

    def track_table(self):
        x = self.torch.zeros((1, self.size, 4), dtype=self.torch.float32, device=self.device)
        named = self._named(x, all_heads=True)
        rows = []
        for head, res in (("atac", 1), ("rna_seq", 128), ("dnase", 128), ("cage", 128), ("chip_tf", 128), ("chip_histone", 128)):
            try:
                for i, t in enumerate(getattr(named, head)[res].tracks):
                    d = t.to_dict()
                    d.update(head=head, index=i)
                    rows.append(d)
            except Exception as e:  # head missing for mouse
                print(f"skip head {head}: {e}")
        return rows

    def _resolve_tracks(self, strict=True):
        rows = self.track_table()
        out = {}
        for name, spec in self.cfg["tracksets"].items():
            out[name] = {}
            for head, pat in spec.items():
                if head not in ACC_HEADS and head not in EXPR_HEADS:
                    raise SystemExit(f"trackset {name}: unsupported head '{head}' (use atac, dnase, rna_seq, cage)")
                rx = re.compile(pat, re.I)
                hit = [r for r in rows if r["head"] == head and rx.search(track_text(r))]
                print(f"trackset {name} {head} /{pat}/: {len(hit)} tracks", flush=True)
                for r in hit:
                    print(f"    [{r['index']}] {r.get('track_name')} | {r.get('biosample_name')} | {r.get('biosample_life_stage')} | strand {r.get('strand')}", flush=True)
                if not hit:
                    msg = (f"trackset {name}/{head}: pattern '{pat}' matches no mouse track; edit config.json "
                           "(see results/mouse_tracks.csv from 00_check_env)")
                    if strict:
                        raise SystemExit(msg)
                    print("TRACKSET PROBLEM:", msg, flush=True)
                    continue
                out[name][head] = [r["index"] for r in hit]
        return out

    def predict(self, seqs):
        """seqs: list of window codes -> list of {trackset: {head: 1-D array (L at 1 bp, L/128 at 128 bp)}}."""
        torch = self.torch
        res = []
        for i in range(0, len(seqs), self.bs):
            x = torch.from_numpy(np.stack([onehot(s) for s in seqs[i:i + self.bs]])).to(self.device)
            named = self._named(x)
            per = {}
            for ts, spec in self.tracksets.items():
                per[ts] = {}
                for head, idx in spec.items():
                    tens = getattr(named, head)[1 if head in ACC_HEADS else 128].tensor
                    per[ts][head] = tens[..., idx].mean(-1).float().cpu().numpy()
            for b in range(x.shape[0]):
                res.append({ts: {h: arr[b] for h, arr in per[ts].items()} for ts in per})
        return res


def readout(pred, regions, tss_rel, cfg):
    """Readouts for one prediction of one trackset. regions: motif / w170 / peak (a, b) and mid, window offsets."""
    p = cfg["perturb"]
    out = {}
    for head, arr in pred.items():
        if head in ACC_HEADS:
            for suffix, key in (("motif", "motif"), ("170", "w170"), ("peak", "peak")):
                a, b = regions[key]
                out[f"{ACC_HEADS[head]}_{suffix}"] = float(arr[a:b].mean()) if b > a else np.nan
        else:
            pre, n = EXPR_HEADS[head], len(arr)
            mid, loc = regions["mid"], p["local_bp"]
            lo, hi = max(mid - loc, 0) // 128, min(mid + loc, n * 128) // 128 + 1
            out[f"{pre}_mean50kb"] = float(arr[lo:hi].mean())
            out[f"{pre}_max50kb"] = float(arr[lo:hi].max())
            if tss_rel is not None and 0 <= tss_rel < n * 128:
                t, pr = tss_rel, p["tss_prox_bp"]
                out[f"{pre}_prox2kb"] = float(arr[max(t - pr, 0) // 128: min(t + pr, n * 128) // 128 + 1].mean())
            else:
                out[f"{pre}_prox2kb"] = np.nan
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
