#!/usr/bin/env python3
"""CPU-only diagnostics (fine on a login node or interactive job): which alphagenome_pytorch is installed, whether it has
named_outputs, how the weights file's key names compare with the model's, and how the package itself loads weights.

  python 00_inspect_weights.py [weights.safetensors]
"""
from __future__ import annotations

import importlib
import inspect
import sys
from pathlib import Path

import ag_lib as L


def main():
    cfg = L.load_config()
    wpath = sys.argv[1] if len(sys.argv) > 1 else cfg["paths"]["weights"]
    import alphagenome_pytorch as ap
    print("package file:", ap.__file__, "| version:", getattr(ap, "__version__", "n/a"))
    try:
        from importlib.metadata import version
        print("pip version:", version("alphagenome-pytorch"))
    except Exception as e:
        print("pip version unknown:", e)
    print("submodules:", sorted(p.stem for p in Path(ap.__file__).parent.glob("*.py")),
          sorted(p.name for p in Path(ap.__file__).parent.iterdir() if p.is_dir()))
    try:
        importlib.import_module("alphagenome_pytorch.named_outputs")
        print("named_outputs: present")
    except Exception as e:
        print("named_outputs: MISSING:", e)
    print("exported names:", [n for n in dir(ap) if not n.startswith("_")][:60])
    from alphagenome_pytorch import AlphaGenome
    print("AlphaGenome.from_pretrained signature:", inspect.signature(AlphaGenome.from_pretrained))

    from safetensors import safe_open
    with safe_open(wpath, "pt") as f:
        wk = list(f.keys())
        md = f.metadata()
    model = AlphaGenome()
    mk = list(model.state_dict().keys())
    print(f"\nweights: {len(wk)} tensors | metadata: {md}\nmodel:   {len(mk)} tensors")
    print("weight keys (first 12):", *wk[:12], sep="\n  ")
    print("model keys (first 12):", *mk[:12], sep="\n  ")
    sw, sm = set(wk), set(mk)
    print(f"exact name overlap: {len(sw & sm)}; only in weights: {len(sw - sm)}; only in model: {len(sm - sw)}")
    for pre in ("model.", "module.", "alphagenome."):
        stripped = {k[len(pre):] for k in wk if k.startswith(pre)}
        if stripped:
            print(f"after stripping '{pre}': overlap {len(stripped & sm)}")


if __name__ == "__main__":
    main()
