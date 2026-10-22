#!/usr/bin/env python3
"""The AP-1 affinity titration, analysed per site and compared across every class definition (all 3,678 cardiomyocyte elements).
Conditions (15_scan_titrate.py): core_mm1..7 weaken the 7 bp core by 1-7 bases; core_upgrade makes an imperfect core perfect (only where there is none);
motif_rc reverses the 11 bp motif; core_halfsite keeps only the first half-site; dup/move place a second or relocated copy at +-20/40/80/150 bp,
each against its own shuffled control (ctrl, ctrlmove).
Per site this gives a dose-response slope, a saturation level, an orientation effect, and the net value of adding or moving a motif.
Classes are compared with 1:1 matching on predicted accessibility (motif and peak), peak width, TSS distance and MEASURED accessibility; BH within contrast.
Outputs Analysis/titration/{tables,figures}.   python3 47_titration_by_class.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "titration"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(47)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY, RED = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993", "#b3261e"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
MCOLS = ["l_170", "l_peak", "l_pw", "l_tss", "l_meas"]
OFFSETS = [-150, -80, -40, -20, 20, 40, 80, 150]


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_med(v, n=400):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5]) if len(v) > 5 else (np.nan, np.nan)


# ---------------------------------------------------------------- effects per site and condition
tt = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/titrate*_cm_chunk*.csv")))], ignore_index=True)
tt = tt[tt.trackset == "CM"]
E = {}
for rd, col in (("170", "atac_170"), ("motif", "atac_motif")):
    nat = tt[tt.kind == "native"].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = tt[tt.kind != "native"].copy()
    e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    E[rd] = e.groupby(["site_id", "kind", "pos"], as_index=False).eff.mean()
core_found = tt.drop_duplicates("site_id").set_index("site_id").core_found
print(f"titration: {tt.site_id.nunique():,} sites; {int(core_found.sum()):,} with an exact core, {int((~core_found).sum()):,} without")

feat = {}
for rd in ("170", "motif"):
    e = E[rd].set_index(["site_id", "kind", "pos"]).eff
    mm = E[rd][E[rd].kind == "core_mm"]
    rows = []
    for sid, w in mm.groupby("site_id"):
        w = w.sort_values("pos")
        k, y = w.pos.values.astype(float), w.eff.values
        slope = float(np.polyfit(k, y, 1)[0]) if len(k) > 2 else np.nan
        rows.append({"site_id": sid, f"dose_slope_{rd}": slope, f"mm1_{rd}": y[0] if len(y) else np.nan,
                     f"mm7_{rd}": y[-1] if len(y) else np.nan,
                     f"dose_curvature_{rd}": float(y[-1] - 7 * slope) if len(y) else np.nan})
    f = pd.DataFrame(rows).set_index("site_id")
    for cond, name in (("motif_rc", "reverse_complement"), ("core_halfsite", "halfsite"), ("core_upgrade", "upgrade")):
        s = E[rd][E[rd].kind == cond].set_index("site_id").eff
        f[f"{name}_{rd}"] = s.reindex(f.index)
    # net value of a second copy, and of moving the motif, each against its own shuffled control
    for a_, b_, name in (("dup", "ctrl", "add_copy"), ("move", "ctrlmove", "relocate")):
        A = E[rd][E[rd].kind == a_].pivot_table(index="site_id", columns="pos", values="eff")
        B = E[rd][E[rd].kind == b_].pivot_table(index="site_id", columns="pos", values="eff")
        net = (A - B).reindex(f.index)
        f[f"{name}_net_{rd}"] = net.mean(1)
        near = [o for o in OFFSETS if abs(o) <= 40 and o in net.columns]
        far = [o for o in OFFSETS if abs(o) >= 80 and o in net.columns]
        f[f"{name}_near_{rd}"] = net[near].mean(1) if near else np.nan
        f[f"{name}_far_{rd}"] = net[far].mean(1) if far else np.nan
        f[f"{name}_distance_decay_{rd}"] = f[f"{name}_near_{rd}"] - f[f"{name}_far_{rd}"]
    feat[rd] = f
TF = feat["170"].join(feat["motif"]).join(core_found.rename("exact_core_scan"))
TF.to_csv(OUT / "tables" / "titration_features_per_site.csv")

# ---------------------------------------------------------------- site table
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
core = pd.read_csv(ROOT / "Data/derived/core_quality.csv").set_index("site_id")
S = sl.join(scl[["age_score_class", "score_dynamic"]]).join(meas.rename("measured")).join(core[["exact_core", "core_mismatches"]]).join(TF, how="inner")
S["l_170"], S["l_peak"] = np.log10(S.native_170.clip(lower=1e-3)), np.log10(S.native_peak.clip(lower=1e-3))
S["l_pw"], S["l_tss"], S["l_meas"] = np.log10(S.peak_width), np.log10(S.abs_tss + 1), np.log10(S.measured.clip(lower=1e-3))
S["age_class"] = S.age_class.replace({"Neither YS nor MAS (= MAE only)": np.nan})
S["ex_class"] = S.ex_class.replace({"Neither MAS nor MAE (= YS only)": np.nan})

FEATURES = [("dose_slope_170", "dose-response slope (log2FC per core base changed)", "core"),
            ("mm1_170", "effect of changing 1 core base", "core"),
            ("mm7_170", "effect of changing all 7 core bases", "core"),
            ("halfsite_170", "effect of keeping only the first half-site", "core"),
            ("reverse_complement_170", "effect of reversing the 11 bp motif", "orientation"),
            ("upgrade_170", "gain from upgrading to a perfect core", "upgrade"),
            ("add_copy_net_170", "net value of a second motif copy (vs shuffled control)", "copy"),
            ("add_copy_distance_decay_170", "does a second copy help more when near?", "copy"),
            ("relocate_net_170", "net effect of relocating the motif (vs control)", "copy"),
            ("relocate_distance_decay_170", "does relocation cost less when near?", "copy")]
CONTRASTS = [("Dynamic vs Stable", "stable_dynamic", "Dynamic", "Stable"),
             ("Age Up vs Stable (threshold)", "age_class", "Age Up", "Stable"),
             ("Age Down vs Stable (threshold)", "age_class", "Age Down", "Stable"),
             ("Ex Up vs Stable (threshold)", "ex_class", "Exercise Up", "Stable"),
             ("Ex Down vs Stable (threshold)", "ex_class", "Exercise Down", "Stable"),
             ("Age Up vs Unchanged (Δscore)", "age_score_class", "Age Up (score)", "Unchanged"),
             ("Age Down vs Unchanged (Δscore)", "age_score_class", "Age Down (score)", "Unchanged"),
             ("Changing vs Steady (Δscore)", "score_dynamic", "Changing (score)", "Steady (score)")]


def make_pairs(frame, col, a_, b_):
    m_ = frame[frame[col].isin([a_, b_])].dropna(subset=MCOLS).copy()
    Z = (m_[MCOLS] - m_[MCOLS].mean()) / m_[MCOLS].std()
    avail, pairs = set(m_.index[m_[col] == b_]), []
    for cid in rng.permutation(m_.index[m_[col] == a_].tolist()):
        cand = np.array(sorted(avail))
        if not len(cand):
            break
        d_ = np.sqrt(((Z.loc[cand].values - Z.loc[cid].values) ** 2).sum(1))
        j = int(np.argmin(d_))
        if d_[j] < 0.5:
            pairs.append((cid, cand[j]))
            avail.discard(cand[j])
    return pd.DataFrame(pairs, columns=["case", "control"])


rows = []
for lab, col, a_, b_ in CONTRASTS:
    for v, vlab, fam in FEATURES:
        frame = S[S.exact_core_scan] if fam in ("core", "orientation") else (S[~S.exact_core_scan] if fam == "upgrade" else S)
        P = make_pairs(frame, col, a_, b_)
        if len(P) < 25:
            continue
        x, y = frame[v].reindex(P.case).values.astype(float), frame[v].reindex(P.control).values.astype(float)
        ok = ~(np.isnan(x) | np.isnan(y))
        if ok.sum() < 25 or not np.any(x[ok] != y[ok]):
            continue
        d = x[ok] - y[ok]
        sd = np.nanstd(frame[v].astype(float))
        rows.append({"contrast": lab, "family": fam, "feature": vlab, "n_pairs": int(ok.sum()),
                     "median_class": float(np.median(x[ok])), "median_reference": float(np.median(y[ok])),
                     "mean_difference": float(np.mean(d)), "standardised": float(np.mean(d) / sd) if sd else np.nan,
                     "p": float(wilcoxon(d).pvalue)})
R = pd.DataFrame(rows)
R["fdr"] = R.groupby("contrast").p.transform(lambda p: bh(p))
R["detects"] = np.where(R.fdr < 0.05, "YES", "no")
R.to_csv(OUT / "tables" / "titration_by_class.csv", index=False)

# does the dose-response track the core quality we measured in the sequence? (a validity check)
chk = S[S.exact_core_scan]
rho = spearmanr(S.core_mismatches, S.dose_slope_170, nan_policy="omit")
rho2 = spearmanr(S[~S.exact_core_scan].upgrade_170, S[~S.exact_core_scan].core_mismatches, nan_policy="omit")
pd.DataFrame([{"check": "dose slope vs core mismatches (all sites)", "rho": rho.statistic, "p": rho.pvalue, "n": int(S[["core_mismatches", "dose_slope_170"]].notna().all(1).sum())},
              {"check": "upgrade gain vs core mismatches (imperfect-core sites)", "rho": rho2.statistic, "p": rho2.pvalue, "n": int((~S.exact_core_scan).sum())}]).to_csv(OUT / "tables" / "validity_checks.csv", index=False)

# ---------------------------------------------------------------- figures
GRPS = [("Stable", "stable_dynamic", "Stable", BLUE), ("Dynamic", "stable_dynamic", "Dynamic", ORANGE)]
fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2))
ax = axes[0]
mm = E["170"][E["170"].kind == "core_mm"].join(S[["stable_dynamic", "age_score_class"]], on="site_id")
mm = mm[mm.site_id.isin(S.index[S.exact_core_scan])]
dose_rows = []
for lab, col, val, c in GRPS:
    w = mm[mm[col] == val]
    med = w.groupby("pos").eff.median()
    lo = w.groupby("pos").eff.apply(lambda v: boot_med(v.values)[0])
    hi = w.groupby("pos").eff.apply(lambda v: boot_med(v.values)[1])
    ax.plot(med.index, med.values, color=c, marker="o", ms=4.5, lw=1.8, label=f"{lab} (n = {w.site_id.nunique():,})")
    ax.fill_between(med.index, lo.values, hi.values, color=c, alpha=0.15, lw=0)
    dose_rows += [{"class": lab, "core_bases_changed": int(p_), "median_effect": v_} for p_, v_ in med.items()]
pd.DataFrame(dose_rows).to_csv(OUT / "tables" / "dose_response_curve.csv", index=False)
ax.set_xlabel("Core bases changed (of 7)")
ax.set_ylabel("Median effect on the 170 bp readout (log2FC)")
ax.set_title("Weakening the AP-1 core: dose-response", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)

ax = axes[1]
conds = [("upgrade_170", "upgrade an\nimperfect core", "upgrade"), ("reverse_complement_170", "reverse the\n11 bp motif", "orientation"),
         ("halfsite_170", "keep one\nhalf-site", "core"), ("add_copy_net_170", "add a second\ncopy (net)", "copy"), ("relocate_net_170", "relocate the\nmotif (net)", "copy")]
for j, (v, lab, fam) in enumerate(conds):
    frame = S[S.exact_core_scan] if fam in ("core", "orientation") else (S[~S.exact_core_scan] if fam == "upgrade" else S)
    for i, (glab, col, val, c) in enumerate(GRPS):
        w = frame[frame[col] == val][v].dropna()
        if len(w) < 10:
            continue
        lo, hi = boot_med(w.values)
        ax.errorbar(j + (i - 0.5) * 0.3, w.median(), yerr=[[w.median() - lo], [hi - w.median()]], fmt="o", color=c, ms=6, lw=1.5, capsize=0,
                    label=glab if j == 0 else None)
ax.axhline(0, color=MUTED, lw=0.9)
ax.set_xticks(range(len(conds)), [c[1] for c in conds], fontsize=8)
ax.set_ylabel("Median effect, 170 bp readout (log2FC)")
ax.set_title("Single edits, with 95% bootstrap CI", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)

ax = axes[2]
for name, lab, c in (("add_copy", "second copy (dup − ctrl)", AQUA), ("relocate", "relocated motif (move − ctrlmove)", ORANGE)):
    A = E["170"][E["170"].kind == ("dup" if name == "add_copy" else "move")].pivot_table(index="site_id", columns="pos", values="eff")
    B = E["170"][E["170"].kind == ("ctrl" if name == "add_copy" else "ctrlmove")].pivot_table(index="site_id", columns="pos", values="eff")
    net = (A - B)
    med = net.median()
    ax.plot(med.index, med.values, color=c, marker="o", ms=4.5, lw=1.8, label=lab)
ax.axhline(0, color=MUTED, lw=0.9)
ax.axvline(0, color=GRID, lw=6)
ax.set_xlabel("Offset of the second or relocated copy (bp)")
ax.set_ylabel("Net effect vs its own shuffled control (log2FC)")
ax.set_title("Does a motif work anywhere in the peak?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "figures" / "z1_titration_overview.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(11.0, 5.2))
labels = [f[1] for f in FEATURES]
cs = [c[0] for c in CONTRASTS]
M = R.pivot_table(index="feature", columns="contrast", values="standardised").reindex(index=labels, columns=cs)
Q = R.pivot_table(index="feature", columns="contrast", values="fdr").reindex(index=labels, columns=cs)
m = M.values
lim = np.nanpercentile(np.abs(m), 97) if np.isfinite(m).any() else 1
im = ax.imshow(m, cmap="PuOr_r", vmin=-lim, vmax=lim, aspect="auto")
ax.grid(False)
ax.set_xticks(range(len(cs)), [c.replace(" vs ", "\nvs ").replace(" (threshold)", "").replace(" (Δscore)", "\n(Δscore)") for c in cs], fontsize=7.5, rotation=30, ha="right")
ax.set_yticks(range(len(labels)), labels, fontsize=8)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        if not np.isfinite(m[i, j]):
            continue
        sig = Q.values[i, j] < 0.05
        ax.text(j, i, f"{m[i, j]:+.2f}" + ("*" if sig else ""), ha="center", va="center", fontsize=6.8,
                fontweight="bold" if sig else "normal", color="white" if abs(m[i, j]) > 0.6 * lim else INK)
        if sig:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.3))
ax.axvline(4.5, color=INK, lw=1.6)
ax.set_title("AP-1 affinity titration, by class\nstandardised difference from the reference; * and box = FDR < 0.05", fontsize=9.5, loc="left")
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("class minus reference (SD units)", fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "z2_titration_by_class.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 54)
print("\n== overall titration behaviour (all sites with an exact core)")
ov = []
for v, vlab, fam in FEATURES:
    frame = S[S.exact_core_scan] if fam in ("core", "orientation") else (S[~S.exact_core_scan] if fam == "upgrade" else S)
    w = frame[v].dropna()
    lo, hi = boot_med(w.values)
    ov.append({"feature": vlab, "n": len(w), "median": w.median(), "ci_lo": lo, "ci_hi": hi})
OV = pd.DataFrame(ov)
OV.to_csv(OUT / "tables" / "titration_overall.csv", index=False)
print(OV.round(3).to_string(index=False))
print("\n== validity checks\n", pd.read_csv(OUT / "tables" / "validity_checks.csv").round(3).to_string(index=False))
print("\n== significant class differences (FDR < 0.05)")
sig = R[R.fdr < 0.05].sort_values(["contrast", "fdr"])
print(sig[["contrast", "feature", "n_pairs", "median_class", "median_reference", "mean_difference", "standardised", "fdr"]].round(3).to_string(index=False) if len(sig) else "  none")
print("\n== detected per contrast\n", R.groupby("contrast").detects.apply(lambda v: f"{(v == 'YES').sum()} of {len(v)}").to_string())
