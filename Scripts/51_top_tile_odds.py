#!/usr/bin/env python3
"""Does an element's MOST important window carry a given motif family?

Comparing median tile log2FC between classes is confounded: responsive sites have less total dependence to begin
with, so every family looks weaker in them. Ranking within each element removes that. For each site we take the
k most damaging 10 bp windows and ask whether any of them overlaps a footprint of each family; the class comparison
is then an odds ratio, computed twice:
  unconditional  -- P(family is in the top k)                      : mixes 'is it there' with 'does it matter'
  conditional    -- P(family is in the top k | present at all)     : pure 'when it is there, is it what matters most'
Outputs Analysis/top_tile_odds/.   python3 51_top_tile_odds.py [--topk 3]
"""
from __future__ import annotations

import argparse
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "top_tile_odds"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument("--topk", type=int, default=3)
a = ap.parse_args()
rng = np.random.default_rng(51)
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]

sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
S = sl.join(meas.rename("measured"))
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))


def pairs(frame, col, case, ctrl):
    m_ = frame[frame[col].isin([case, ctrl])].dropna(subset=MCOLS).copy()
    Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
    avail, out = set(m_.index[m_[col] == ctrl]), []
    for cid in rng.permutation(m_.index[m_[col] == case].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d))
        if d[j] < 0.5:
            out.append((cid, cand[j])); avail.discard(cand[j])
    return pd.DataFrame(out, columns=["case", "control"])


# ---------------------------------------------------------------- tiles, ranked within each element
sc = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "pos", "atac_170", "trackset"])
                for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[sc.trackset == "CM"]
nat = sc[sc.kind == "native"].groupby("site_id").atac_170.median()
pc = 0.05 * nat.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
T = t.groupby(["site_id", "pos"], as_index=False).eff.mean()
T["rank"] = T.groupby("site_id").eff.rank(method="first")          # 1 = most damaging
TOP = T[T["rank"] <= a.topk][["site_id", "pos"]]

FP = pd.read_csv(ROOT / "Analysis/partners/tables/footprints_near_sites.csv")
FP = FP[~FP.own_ap1].copy()
fams = [f for f, n in FP.family.value_counts().items() if n >= 300 and f != "Other"]

present, intop = {}, {}
for fam in fams:
    f = FP[FP.family == fam][["site_id", "rel_start", "rel_end"]]
    present[fam] = set(f.site_id)
    j = TOP.merge(f, on="site_id", how="inner")
    intop[fam] = set(j.site_id[(j.pos < j.rel_end) & (j.pos + 10 > j.rel_start)])

GROUPS = [("Dynamic (all)", "stable_dynamic", "Dynamic"), ("Age Up", "age_class", "Age Up"),
          ("Age Down", "age_class", "Age Down"), ("Exercise Up", "ex_class", "Exercise Up"),
          ("Exercise Down", "ex_class", "Exercise Down")]
stable = set(S.index[S.stable_dynamic == "Stable"])
rows = []
for nm, col, val in GROUPS:
    f = S.copy()
    f["grp"] = np.where(f.index.isin(stable), "Stable", np.where(f[col] == val, nm, ""))
    P = pairs(f[f.grp != ""], "grp", nm, "Stable")
    for fam in fams:
        for mode, A, B in (("all sites", list(S.index[S[col] == val]), list(stable)),
                           ("matched", list(P.case), list(P.control))):
            for cond, Ax, Bx in (("unconditional", A, B),
                                 ("given present", [i for i in A if i in present[fam]], [i for i in B if i in present[fam]])):
                if min(len(Ax), len(Bx)) < 20:
                    continue
                a_ = sum(i in intop[fam] for i in Ax); b_ = sum(i in intop[fam] for i in Bx)
                tab = [[a_, len(Ax) - a_], [b_, len(Bx) - b_]]
                orr, p = fisher_exact(tab)
                rows.append({"class": nm, "family": fam, "mode": mode, "conditioning": cond,
                             "n_class": len(Ax), "n_ctrl": len(Bx), "pct_class": 100 * a_ / len(Ax),
                             "pct_ctrl": 100 * b_ / len(Bx), "odds_ratio": orr, "p": p})
D = pd.DataFrame(rows)
for (m_, c_), g in D.groupby(["mode", "conditioning"]):
    o = np.argsort(g.p.values); r = np.empty(len(g)); r[o] = np.arange(1, len(g) + 1)
    D.loc[g.index, "fdr"] = np.minimum(g.p.values * len(g) / r, 1)
D.to_csv(OUT / "tables" / f"top{a.topk}_tile_odds_by_class.csv", index=False)
pd.set_option("display.width", 220)
print(f"== top {a.topk} windows per element; {len(fams)} families, {len(D)} tests\n")
for key in (("matched", "given present"), ("matched", "unconditional"), ("all sites", "given present")):
    g = D[(D["mode"] == key[0]) & (D.conditioning == key[1])].sort_values("p").head(6)
    print(f"-- {key[0]}, {key[1]}")
    print(g[["class", "family", "n_class", "n_ctrl", "pct_class", "pct_ctrl", "odds_ratio", "p", "fdr"]].round(3).to_string(index=False), "\n")
