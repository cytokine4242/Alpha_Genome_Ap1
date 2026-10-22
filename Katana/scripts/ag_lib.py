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
    if os.environ.get("AG_WINDOW"):  # context-length sensitivity runs
        cfg["model"]["window"] = int(os.environ["AG_WINDOW"])
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


# conditions: (name, kind, width, reps). width in motif / 170 / peak for the standard kinds. reps = min(reps[kind], caps[width]): large edits have
# large, stable effects and need fewer replicates than motif-width edits.
# Extra kinds, only built when named in `kinds`:
#   offset   shuffle of a 170 bp region centred OFFSETS bp from the motif midpoint (matched control for the 170 bp shuffle)
#   titrate  shuffle of a region centred on the motif of width TITRATE_WIDTHS (dose-response; also tests the motif-interval definition)
#   sat      single-base substitution of every base of the 11 bp motif to each alternative base (saturation mutagenesis)
OFFSETS = (-600, -300, 300, 600)
TITRATE_WIDTHS = (7, 11, 21, 41, 85, 170, 340, 680)
EXTRA_KINDS = {"offset", "titrate", "sat"}


def build_conditions(reps, caps=None, kinds=None):
    caps = caps or {}
    conds = []
    for kind in ("shuffle", "random", "dinuc"):
        for width in ("motif", "170", "peak"):
            conds.append((f"{kind}@{width}", kind, width, min(reps[kind], caps.get(width, 10**6))))
    conds.append(("core@motif", "core", "motif", reps.get("core", 3)))
    conds.append(("flank@motif", "flank", "motif", min(reps["flank"], caps.get("motif", 10**6))))
    conds.append(("flank@170", "flank", "170", min(reps["flank"], caps.get("170", 10**6))))
    if kinds and (set(kinds) & EXTRA_KINDS):
        for off in OFFSETS:
            conds.append((f"offset@{off:+d}", "offset", str(off), reps.get("offset", 5)))
        for w in TITRATE_WIDTHS:
            conds.append((f"titrate@{w}", "titrate", str(w), reps.get("titrate", 5)))
        for pos in range(11):
            for alt in "ACGT":
                conds.append((f"sat@{pos}{alt}", "sat", f"{pos}:{alt}", 1))
    return [c for c in conds if not kinds or c[1] in kinds]


def apply_condition(codes, regions, kind, width, rng):
    """Return edited window, or None if the edit region holds N. regions: dict of (a, b) offsets in the window."""
    if kind == "core":
        return core_edit(codes, regions["motif"], rng)[0]
    if kind == "flank":
        return flank_edit(codes, regions, width, rng)[0]
    if kind in ("offset", "titrate"):
        mid = regions["mid"]
        if kind == "offset":
            a = mid + int(width) - 85
            b = a + 170
        else:
            w = int(width)
            a = mid - w // 2
            b = a + w
        if a < 0 or b > len(codes) or has_n(codes, a, b):
            return None
        return shuffle_region(codes, a, b, rng)
    if kind == "sat":
        pos, alt = width.split(":")
        m0, m1 = regions["motif"]
        i = m0 + int(pos)
        alt_code = "ACGT".index(alt)
        if i >= m1 or codes[i] == alt_code or codes[i] == 4:
            return None
        out = codes.copy()
        out[i] = alt_code
        return out
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
FT_HEAD = "ft_atac"   # fine-tuned ATAC head (config "finetuned"), returned under "atac" so readouts and columns stay the same; its 1 bp output if it has one
EXPR_HEADS = {"rna_seq": "rna", "cage": "cage"}      # expression, read at 128 bp; columns <prefix>_prox2kb/_mean50kb/_max50kb
TRACK_FIELDS = ("biosample_name", "ontology_curie", "track_name", "assay_title", "biosample_life_stage", "strand", "biosample_type")


def catalog_tracks(organism, heads=("atac", "dnase", "rna_seq", "cage", "chip_tf", "chip_histone")):
    """Track metadata straight from the package's built-in catalog (CPU only, no model). organism: 0 human, 1 mouse."""
    from alphagenome_pytorch.named_outputs import TrackMetadataCatalog
    cat = TrackMetadataCatalog.load_builtin(organism)
    rows = []
    for head in heads:
        try:
            tracks = cat.get_tracks(head, organism=organism)
        except Exception as e:
            print(f"skip head {head} (organism {organism}): {type(e).__name__}: {str(e)[:120]}")
            continue
        for tr in tracks:
            if getattr(tr, "is_padding", False):
                continue
            d = {"head": head, "organism": organism, "index": tr.track_index, "track_name": tr.track_name}
            d.update(getattr(tr, "extras", {}) or {})
            rows.append(d)
    return rows


def track_text(row):
    return " | ".join(str(row.get(k, "")) for k in TRACK_FIELDS)


ORG_NAME = {0: "human", 1: "mouse"}


class AG:
    """AlphaGenome predictor with per-trackset averaged tracks. ATAC/DNase at 1 bp, RNA-seq/CAGE at 128 bp.

    config tracksets: {name: {"organism": 0|1 (optional, default model.organism_index), head: regex, ...}}, head in
    atac / dnase / rna_seq / cage. The regex (case-insensitive) is searched in
    "biosample_name | ontology_curie | track_name | assay_title | life_stage | strand | biosample_type" of that organism's catalog.
    A trackset with organism 0 runs the human output heads on the mouse sequence (an extra forward pass per sequence).
    Only the heads used by the tracksets are computed. strict=False reports unmatched tracksets instead of exiting.

    Fine-tuned head: with config "finetuned": {"checkpoint": <delta .pth>, "label_regex": "_(YS|MAS|MAE|YE)_"} and paths.weights = the base
    weights it was trained on, a trackset may use "ft_atac": regex over the fine-tuned track names (averaged, returned as "atac") and
    "ft_split": true to also return every matched track as "atac@<label>" (columns atac_motif_<label>, ...). The new head runs on the
    embeddings of the same forward pass as the built-in heads, so it costs no extra trunk pass. A head with a 1 bp output (trained on 1 bp
    insertion tracks) is read at 1 bp like the built-in ATAC. A 128 bp-only head (e.g. trained on binned coverage) is refused unless
    finetuned.allow_128bp is true, because its motif and 170 bp readouts would just be the one or two 128 bp bins they overlap."""

    def __init__(self, cfg, strict=True):
        import torch
        self.torch = torch
        self.cfg = cfg
        self.size = cfg["model"]["window"]
        self.default_org = cfg["model"]["organism_index"]
        self.bs = cfg["model"]["batch_size"]
        self.model, self.device = load_model(cfg)
        self._cats, self._cat_org = {}, None
        self.ft_split = {}                                # trackset -> [(fine-tuned track index, label)] for ft_split
        self.ft_head = self._load_finetuned(cfg.get("finetuned"))
        self.tracksets, self.org_of = self._resolve_tracks(strict)
        self.heads = {}  # organism -> sorted built-in heads to compute
        for ts, spec in self.tracksets.items():
            self.heads.setdefault(self.org_of[ts], set()).update(h for h in spec if h != FT_HEAD)
        self.heads = {o: sorted(h) for o, h in self.heads.items() if h}
        self.need_ft = any(FT_HEAD in spec for spec in self.tracksets.values())

    def _load_finetuned(self, spec):
        """The fine-tuned head module (on the GPU) from a delta checkpoint, after checking it was trained on these base weights."""
        if not spec:
            return None
        from alphagenome_pytorch.extensions.finetuning.checkpointing import compute_base_model_weights_hash_from_file, load_finetuned_model
        base = self.cfg["paths"]["weights"]
        spec["checkpoint"] = os.environ.get("AG_FT_CHECKPOINT") or spec["checkpoint"]   # e.g. a pilot checkpoint for a smoke test
        ft, meta = load_finetuned_model(spec["checkpoint"], base, device="cpu")
        want = meta.get("base_model_weights_hash")
        have = compute_base_model_weights_hash_from_file(base)
        if want and want != have:
            raise SystemExit(f"{spec['checkpoint']} was trained on other base weights than {base} ({want} vs {have}); set paths.weights to its base model")
        modality = spec.get("modality", "atac")
        self._ft_names = list(meta["track_names"][modality])
        self._ft_res = tuple(meta["resolutions"][modality])
        if 1 not in self._ft_res and not spec.get("allow_128bp"):
            raise SystemExit(f"fine-tuned head resolutions {self._ft_res}: no 1 bp output, so motif-level readouts would be 128 bp bins; "
                             "use a head trained on 1 bp targets (or set finetuned.allow_128bp to accept that)")
        self._ft_label = re.compile(spec.get("label_regex", "_([^_]+)_"))
        head = ft.heads[modality].to(self.device).eval()
        del ft
        print(f"fine-tuned head {spec['checkpoint']}: epoch {meta.get('epoch')}, val_loss {meta.get('val_loss')}, base weights hash matches, "
              f"{len(self._ft_names)} tracks at {self._ft_res} bp: {self._ft_names}", flush=True)
        return head

    def _set_catalog(self, org):
        from alphagenome_pytorch.named_outputs import TrackMetadataCatalog
        if org not in self._cats:
            self._cats[org] = TrackMetadataCatalog.load_builtin(ORG_NAME[org])
        if self._cat_org != org:
            self.model.set_track_metadata_catalog(self._cats[org])
            self._cat_org = org

    def _named(self, x, org=None, heads=None):
        org = self.default_org if org is None else org
        self._set_catalog(org)
        kw = {"heads": tuple(heads)} if heads else {}
        with self.torch.no_grad():
            return self.model.predict(x, organism_index=org, named_outputs=True, **kw)

    def forward_all(self, x):
        """One forward pass per organism used by the tracksets (only the heads needed). Returns ({organism: named outputs}, fine-tuned
        head output at 1 bp (B, L, tracks) or None). The fine-tuned head reads the embeddings of the default organism's pass."""
        if not self.need_ft:
            return {org: self._named(x, org, heads) for org, heads in self.heads.items()}, None
        torch = self.torch
        named = {org: self._named(x, org, heads) for org, heads in self.heads.items() if org != self.default_org}
        org = self.default_org
        heads = self.heads.get(org)
        self._set_catalog(org)
        oi = self.model._normalize_organism_index(org, x)
        with torch.no_grad():
            kw = {"heads": tuple(heads)} if heads else {"embeddings_only": True}
            raw = self.model.predict(x, organism_index=org, named_outputs=False, return_embeddings=True, resolutions=(1, 128), **kw)
            emb = {}
            for res in self._ft_res:
                e, n = raw[f"embeddings_{res}bp"], x.shape[1] // res
                emb[res] = (e.transpose(1, 2) if e.shape[1] == n else e).float().contiguous()   # the head wants (B, C, L/res)
            out = self.ft_head(emb, oi, return_scaled=False, channels_last=True)
            ft = out[1].float() if 1 in out else torch.repeat_interleave(out[128].float(), 128, dim=1)   # (B, L, tracks) at 1 bp
            if heads:
                raw = {k: v for k, v in raw.items() if not k.startswith("embeddings")}
                named[org] = self.model.named_outputs(raw, organism=oi, strict_metadata=False, channels_last=True, include_padding=False)
        return named, ft

    def track_table(self, organism=None):
        """All tracks of an organism's catalog (CPU only)."""
        return catalog_tracks(self.default_org if organism is None else organism)

    def _resolve_tracks(self, strict=True):
        out, org_of, cache = {}, {}, {}
        for name, spec in self.cfg["tracksets"].items():
            org = spec.get("organism", self.default_org)
            org_of[name] = org
            rows = cache.setdefault(org, catalog_tracks(org))
            if FT_HEAD in spec:
                if "atac" in spec:
                    raise SystemExit(f"trackset {name}: use atac or {FT_HEAD}, not both (both write the atac columns)")
                if self.ft_head is None:
                    raise SystemExit(f"trackset {name} uses {FT_HEAD} but the config has no 'finetuned' checkpoint")
                if org != self.default_org:
                    raise SystemExit(f"trackset {name}: {FT_HEAD} runs on the model's default organism only")
                rx = re.compile(spec[FT_HEAD], re.I)
                idx = [i for i, t in enumerate(self._ft_names) if rx.search(t)]
                print(f"trackset {name} (fine-tuned) /{spec[FT_HEAD]}/: {len(idx)} tracks: {[self._ft_names[i] for i in idx]}", flush=True)
                if not idx:
                    raise SystemExit(f"trackset {name}: {FT_HEAD} pattern matches none of {self._ft_names}")
                out[name] = {FT_HEAD: idx}
                if spec.get("ft_split"):
                    labels = [(m.group(1) if (m := self._ft_label.search(self._ft_names[i])) else str(i)) for i in idx]
                    if len(set(labels)) != len(labels):
                        raise SystemExit(f"trackset {name}: label_regex gives duplicate labels {labels}")
                    self.ft_split[name] = list(zip(idx, labels))
            else:
                out[name] = {}
            for head, pat in spec.items():
                if head in ("organism", FT_HEAD, "ft_split"):
                    continue
                if head not in ACC_HEADS and head not in EXPR_HEADS:
                    raise SystemExit(f"trackset {name}: unsupported head '{head}' (use atac, dnase, rna_seq, cage, {FT_HEAD})")
                rx = re.compile(pat, re.I)
                hit = [r for r in rows if r["head"] == head and rx.search(track_text(r))]
                print(f"trackset {name} ({ORG_NAME[org]}) {head} /{pat}/: {len(hit)} tracks", flush=True)
                for r in hit:
                    print(f"    [{r['index']}] {r.get('track_name')} | {r.get('biosample_name')} | {r.get('biosample_life_stage')} | strand {r.get('strand')}", flush=True)
                if not hit:
                    msg = (f"trackset {name}/{head}: pattern '{pat}' matches no {ORG_NAME[org]} track; edit config.json "
                           "(see results/*tracks.csv, 00_list_tracks.py)")
                    if strict:
                        raise SystemExit(msg)
                    print("TRACKSET PROBLEM:", msg, flush=True)
                    continue
                out[name][head] = [r["index"] for r in hit]
        return out, org_of

    def predict(self, seqs):
        """seqs: list of window codes -> list of {trackset: {head: 1-D array (L at 1 bp, L/128 at 128 bp)}}."""
        torch = self.torch
        res = []
        for i in range(0, len(seqs), self.bs):
            x = torch.from_numpy(np.stack([onehot(s) for s in seqs[i:i + self.bs]])).to(self.device)
            named_by_org, ft = self.forward_all(x)
            per = {}
            for ts, spec in self.tracksets.items():
                per[ts] = {}
                for head, idx in spec.items():
                    if head == FT_HEAD:
                        per[ts]["atac"] = ft[..., idx].mean(-1).cpu().numpy()
                        for j, lab in self.ft_split.get(ts, []):
                            per[ts][f"atac@{lab}"] = ft[..., j].cpu().numpy()
                        continue
                    tens = getattr(named_by_org[self.org_of[ts]], head)[1 if head in ACC_HEADS else 128].tensor
                    per[ts][head] = tens[..., idx].mean(-1).float().cpu().numpy()
            for b in range(x.shape[0]):
                res.append({ts: {h: arr[b] for h, arr in per[ts].items()} for ts in per})
        return res


def readout(pred, regions, tss_rel, cfg):
    """Readouts for one prediction of one trackset. regions: motif / w170 / peak (a, b) and mid, window offsets."""
    p = cfg["perturb"]
    out = {}
    for head, arr in pred.items():
        base, _, label = head.partition("@")            # "atac@YS": one fine-tuned track -> atac_motif_YS, ...
        if base in ACC_HEADS:
            for suffix, key in (("motif", "motif"), ("170", "w170"), ("peak", "peak")):
                a, b = regions[key]
                out[f"{ACC_HEADS[base]}_{suffix}{'_' + label if label else ''}"] = float(arr[a:b].mean()) if b > a else np.nan
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
