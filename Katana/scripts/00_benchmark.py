#!/usr/bin/env python3
"""Where does the forward-pass time go? Times variants of the real prediction call on the current GPU and prints seconds per
sequence and peak GPU memory: batch size, with/without the 1 bp decoder (resolutions), heads (atac / atac+rna_seq). Also prints
the parameter dtype and whether TF32 / autocast are on. Run on the GPU type you want to know about (GPU=L40S bash submit_all.sh benchmark).

  python 00_benchmark.py
"""
from __future__ import annotations

import time

import torch

import ag_lib as L


def timed(fn, n=3):
    fn()  # warm-up
    torch.cuda.synchronize()
    ts = []
    for _ in range(n):
        t0 = time.time()
        fn()
        torch.cuda.synchronize()
        ts.append(time.time() - t0)
    return min(ts)


def main():
    cfg = L.load_config()
    ag = L.AG(cfg, strict=False)
    m = ag.model
    print("GPU:", torch.cuda.get_device_name(0), "| params dtype:", next(m.parameters()).dtype,
          "| tf32 matmul:", torch.backends.cuda.matmul.allow_tf32, "| tf32 cudnn:", torch.backends.cudnn.allow_tf32)
    org = ag.default_org

    def run(batch, heads, res):
        x = torch.rand(batch, ag.size, 4, device=ag.device)
        x = (x == x.max(-1, keepdim=True).values).float()  # one-hot-like
        kw = {"heads": tuple(heads)} if heads else {}
        if res:
            kw["resolutions"] = tuple(res)
        with torch.no_grad():
            m.predict(x, organism_index=org, named_outputs=False, **kw)

    rows = []
    for batch, heads, res, label in (
        (4, ["atac", "rna_seq"], None, "baseline: atac+rna_seq, 1bp+128bp"),
        (1, ["atac", "rna_seq"], None, "batch 1"), (2, ["atac", "rna_seq"], None, "batch 2"),
        (8, ["atac", "rna_seq"], None, "batch 8"), (16, ["atac", "rna_seq"], None, "batch 16"),
        (4, ["atac"], None, "atac only (1bp+128bp)"),
        (4, ["atac", "rna_seq"], [128], "no 1bp decoder (resolutions=128)"),
        (8, ["atac", "rna_seq"], [128], "no 1bp decoder, batch 8"),
        (4, ["rna_seq"], [128], "rna_seq only, 128bp"),
    ):
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        try:
            dt = timed(lambda: run(batch, heads, res))
            rows.append((label, batch, dt / batch, torch.cuda.max_memory_allocated() / 2**30))
        except Exception as e:
            rows.append((label, batch, float("nan"), float("nan")))
            print(f"  {label}: failed {type(e).__name__}: {str(e)[:80]}")
            torch.cuda.empty_cache()
    base = rows[0][2]
    print(f"\n{'variant':42s}{'batch':>6s}{'s/seq':>9s}{'vs base':>9s}{'peak GB':>9s}")
    for label, b, s, g in rows:
        print(f"{label:42s}{b:6d}{s:9.3f}{s / base:9.2f}{g:9.1f}")


if __name__ == "__main__":
    main()
