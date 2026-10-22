#!/usr/bin/env python3
"""Is there a spacing rule between an AP-1 site and its partner motifs, or helical (about 10.5 bp) phasing that would point to direct cooperation?
Uses every footprinted cardiomyocyte TFBS (Data/derived/tobias_cm_all_motifs.csv.gz) around all 3,679 cardiomyocyte AP-1 sites, not only the scanned subset.
Overlapping instances of the same family at one locus are merged first, so redundant PWMs of one family are not counted several times.
 A  distance distribution per family, against the pooled distribution of all footprints (which carries the peak shape)
 B  fine structure within 150 bp: is there a preferred spacing, or a steric hole close in?
 C  periodicity: power spectrum of the 1 bp distance histogram, tested against jittered nulls. NOTE Tn5 itself cuts with about 10.4 bp
    periodicity on nucleosomal DNA, so a peak there is not by itself evidence of protein-protein contact.
 D  does importance depend on spacing? (scanned sites only)   E  do the classes differ in partner spacing?
Outputs Analysis/distance_rules/{tables,figures}.   python3 40_partner_distance_rules.py"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "distance_rules"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(40)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
HALF = 500
KEY = ["MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD", "CTCF / cohesin", "AP-1 (Fos/Jun/Atf/Batf)", "NKX / homeobox cardiac", "TBX", "KLF / SP"]

def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


# family map: copied from 38_tfbs_partner_mapping.py so both scripts label families identically
import re

FAM = [("AP-1 (Fos/Jun/Atf/Batf)", r"^(Fos|Jun|Batf|Atf3|Bach|Nfe2|Maf)"), ("MEF2", r"^Mef2|^Borcs8"), ("GATA", r"^Gata|^Trps1"), ("TBX", r"^Tbx|^Tbr|^Eomes"),
       ("NKX / homeobox cardiac", r"^Nkx|^Hand|^Msx|^Pitx|^Isl"), ("NFAT", r"^Nfat"), ("NFI (Nfia/b/c/x)", r"^Nfi[abcx]"), ("TALE homeodomain (Meis/Pknox/Tgif/Pbx)", r"^Meis|^Pknox|^Tgif|^Pbx"),
       ("TEAD", r"^Tead"), ("SRF / MADS", r"^Srf"), ("ETS", r"^Ets|^Elk|^Erg|^Fli1|^Gabp|^Elf|^Ehf|^Spdef|^Etv"), ("RUNX", r"^Runx"),
       ("KLF / SP", r"^Klf|^Sp[0-9]|^Patz|^Zfp(281|341|383|467|740)|^Wt1|^Egr"), ("CTCF / cohesin", r"^Ctcf|^Ctcfl"), ("E-box bHLH (Myod/Tcf/Neuro)", r"^Myod|^Myog|^Myf|^Tcf(3|4|12)|^Neuro|^Ascl|^Ptf1|^Bhlh"),
       ("Nuclear receptor (Esrr/Nr/Ppar/Thr)", r"^Esrr|^Nr[0-9]|^Ppar|^Rxr|^Rar|^Thr|^Vdr|^Ar$|^Nr3|^Esr|^Hnf4"), ("MYC / MAX", r"^Myc|^Max|^Mxd|^Mnt|^Usf|^Srebf|^Tfe|^Mitf"),
       ("FOX", r"^Fox"), ("SOX", r"^Sox"), ("STAT / IRF", r"^Stat|^Irf"), ("NF-kB / Rel", r"^Rel|^Nfkb"), ("CEBP", r"^Cebp|^Ddit3"), ("Zinc finger, other", r"^Zfp|^Znf|^Zbtb|^Zic|^Zeb|^Snai"),
       ("Smad / Tgfb", r"^Smad"), ("RFX", r"^Rfx"), ("YY / Gli / Pax", r"^Yy[0-9]|^Gli|^Pax"), ("Heat shock / stress", r"^Hsf|^Xbp|^Creb|^Atf[0-9]")]


def family(name):
    for lab, rx in FAM:
        if re.match(rx, str(name), re.I):
            return lab
    return "Other"


# ---------------------------------------------------------------- footprints around every AP-1 site
tob = pd.read_csv(ROOT / "Data/derived/tobias_cm_all_motifs.csv.gz", usecols=["seqnames", "start", "end", "name"], dtype={"name": "category"})
tob["family"] = tob.name.map(family).astype("category")
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv")
sl = pd.read_csv(ROOT / "Analysis/cm/qc/tables/site_level_table_cm.csv").set_index("site_id")
scl = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
sites = sites.drop(columns=[c for c in ("stable_dynamic", "age_score_class") if c in sites]).join(scl[["age_score_class"]], on="site_id").join(sl[["stable_dynamic"]], on="site_id")

idx = {c: g.sort_values("start").reset_index(drop=True) for c, g in tob.groupby("seqnames", observed=True)}
rows = []
for r in sites.itertuples():
    g = idx.get(r.seqnames)
    if g is None:
        continue
    lo, hi = r.motif_mid - HALF, r.motif_mid + HALF
    s = g.start.values
    w = g.iloc[max(np.searchsorted(s, lo - 60), 0):np.searchsorted(s, hi)]
    w = w[w.end.values > lo]
    if len(w):
        rows.append(pd.DataFrame({"site_id": r.site_id, "family": w.family.values, "rs": w.start.values - r.motif_mid, "re": w.end.values - r.motif_mid}))
H = pd.concat(rows, ignore_index=True)
H = H[~((H.family == "AP-1 (Fos/Jun/Atf/Batf)") & (H.rs < 6) & (H.re > -5))]          # drop the site's own motif

# merge overlapping instances of one family at one locus, so redundant PWMs count once
merged = []
for (sid, fam_), w in H.sort_values(["site_id", "family", "rs"]).groupby(["site_id", "family"], observed=True):
    cs, ce = None, None
    for a, b in zip(w.rs.values, w.re.values):
        if cs is None:
            cs, ce = a, b
        elif a <= ce:
            ce = max(ce, b)
        else:
            merged.append((sid, fam_, cs, ce))
            cs, ce = a, b
    if cs is not None:
        merged.append((sid, fam_, cs, ce))
M = pd.DataFrame(merged, columns=["site_id", "family", "rs", "re"])
M["centre"] = (M.rs + M.re) / 2
M["dist"] = M.centre.abs()
M = M.join(sites.set_index("site_id")[["stable_dynamic", "age_score_class"]], on="site_id")
M.to_csv(OUT / "tables" / "merged_partner_instances.csv.gz", index=False, compression="gzip")
print(f"{len(sites):,} AP-1 sites; {len(H):,} footprint calls within +/-{HALF} bp -> {len(M):,} merged instances; "
      f"MEF2 {int((M.family == 'MEF2').sum()):,}, GATA {int((M.family == 'GATA').sum()):,}")

# ---------------------------------------------------------------- A/B. distance profile against the pooled background
BINS = np.arange(0, HALF + 1, 5)
ctr = (BINS[:-1] + BINS[1:]) / 2
all_h = np.histogram(M.dist, bins=BINS)[0].astype(float)
all_h /= all_h.sum()
prof_rows = []
for fam_ in KEY:
    w = M[M.family == fam_]
    if len(w) < 200:
        continue
    h = np.histogram(w.dist, bins=BINS)[0].astype(float)
    n = h.sum()
    hp = h / n
    enr = np.where(all_h > 0, hp / all_h, np.nan)
    for c_, e_, k_ in zip(ctr, enr, h):
        prof_rows.append({"family": fam_, "distance": c_, "n": int(k_), "share": k_ / n, "enrichment_vs_all_footprints": e_})
PROF = pd.DataFrame(prof_rows)
PROF.to_csv(OUT / "tables" / "a_distance_profile.csv", index=False)

near = []
for fam_ in KEY:
    w = M[M.family == fam_]
    if len(w) < 200:
        continue
    d = w.dist.values
    ref = M[M.family != fam_].dist.values
    near.append({"family": fam_, "n_instances": len(d), "median_distance": np.median(d), "share_within_25bp": (d <= 25).mean(), "share_within_50bp": (d <= 50).mean(),
                 "share_within_100bp": (d <= 100).mean(), "share_beyond_200bp": (d > 200).mean(),
                 "p_vs_other_families": mannwhitneyu(d, ref).pvalue})
NEAR = pd.DataFrame(near).sort_values("median_distance")
NEAR["fdr"] = bh(NEAR.p_vs_other_families)
NEAR.to_csv(OUT / "tables" / "b_distance_summary.csv", index=False)

# ---------------------------------------------------------------- C. periodicity
def spectrum(d, lo=15, hi=200, smooth=21):
    """Detrended 1 bp histogram over [lo, hi) -> (periods, power). The smooth window removes the broad decay so only fine structure is left."""
    h = np.histogram(d[(d >= lo) & (d < hi)], bins=np.arange(lo, hi + 1))[0].astype(float)
    k = np.ones(smooth) / smooth
    base = np.convolve(np.pad(h, smooth // 2, mode="edge"), k, mode="valid")[:len(h)]
    r = h - base
    r -= r.mean()
    f = np.fft.rfft(r * np.hanning(len(r)))
    power = np.abs(f) ** 2
    freq = np.fft.rfftfreq(len(r))
    with np.errstate(divide="ignore"):
        per = np.where(freq > 0, 1 / freq, np.inf)
    return per, power, h


per_rows, SPEC = [], {}
for fam_, label in [(f, f) for f in KEY] + [("important four", "MEF2 + GATA + NFI + TEAD")]:
    d = M[M.family.isin(["MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD"])].dist.values if fam_ == "important four" else M[M.family == fam_].dist.values
    if len(d) < 400:
        continue
    per, pw, h = spectrum(d)
    band = (per >= 9.5) & (per <= 11.5)                       # helical pitch of B-DNA
    other = (per >= 5) & (per <= 30) & ~band
    obs = pw[band].max() if band.any() else np.nan
    null = []
    for k in range(300):                                      # null: jitter each distance by up to +/-7 bp, which destroys fine phasing but keeps the broad shape
        dj = d + rng.integers(-7, 8, len(d))
        null.append(spectrum(dj[dj >= 0])[1][band].max())
    SPEC[label] = (per, pw, band, np.array(null))
    per_rows.append({"family": label, "n_instances": len(d), "peak_power_10_11bp": obs, "median_null_power": float(np.median(null)),
                     "ratio_to_null": obs / np.median(null), "p_vs_jitter_null": (np.sum(np.array(null) >= obs) + 1) / 301,
                     "peak_power_elsewhere_5_30bp": pw[other].max(), "obs_over_background_band": obs / pw[other].max()})
PER = pd.DataFrame(per_rows)
PER["fdr"] = bh(PER.p_vs_jitter_null)
PER.to_csv(OUT / "tables" / "c_periodicity.csv", index=False)

# ---------------------------------------------------------------- D. does importance depend on spacing? (scanned sites)
IMP = None
fs = sorted(glob.glob(str(ROOT / "Katana/results/scan_cm_chunk*.csv")))
if fs:
    site_fam = pd.read_csv(ROOT / "Analysis/partners/tables/i_site_level_family_importance.csv")
    nearest = M[M.family.isin(KEY)].groupby(["site_id", "family"], observed=True, as_index=False).dist.min()
    j = site_fam.merge(nearest, on=["site_id", "family"], how="inner")
    j["dbin"] = pd.cut(j.dist, [0, 25, 50, 85, 130, 200, 500], labels=["0-25", "26-50", "51-85", "86-130", "131-200", "201-500"])
    rows = []
    for (fam_, b_), w in j.groupby(["family", "dbin"], observed=True):
        if len(w) < 10:
            continue
        rows.append({"family": fam_, "distance_bin": b_, "n_sites": len(w), "median_importance_vs_free": w.vs_free.median()})
    IMP = pd.DataFrame(rows)
    IMP.to_csv(OUT / "tables" / "d_importance_by_distance.csv", index=False)
    rho_rows = []
    from scipy.stats import spearmanr
    for fam_, w in j.groupby("family", observed=True):
        if len(w) < 40:
            continue
        r = spearmanr(w.dist, w.vs_free)
        rho_rows.append({"family": fam_, "n_sites": len(w), "rho_distance_vs_importance": r.statistic, "p": r.pvalue})
    RHO = pd.DataFrame(rho_rows)
    RHO["fdr"] = bh(RHO.p)
    RHO.to_csv(OUT / "tables" / "d_distance_importance_correlation.csv", index=False)

# ---------------------------------------------------------------- E. class differences in spacing
rows = []
for fam_ in KEY:
    w = M[M.family == fam_]
    for col, a_, b_ in (("stable_dynamic", "Dynamic", "Stable"), ("age_score_class", "Age Up (score)", "Unchanged")):
        x, y = w[w[col] == a_].dist, w[w[col] == b_].dist
        if min(len(x), len(y)) < 60:
            continue
        rows.append({"family": fam_, "comparison": f"{a_} vs {b_}", "n_a": len(x), "n_b": len(y), "median_a": x.median(), "median_b": y.median(),
                     "share_within_50bp_a": (x <= 50).mean(), "share_within_50bp_b": (y <= 50).mean(), "p": mannwhitneyu(x, y).pvalue})
E = pd.DataFrame(rows)
E["fdr"] = E.groupby("comparison").p.transform(lambda v: bh(v))
E.to_csv(OUT / "tables" / "e_class_differences_in_spacing.csv", index=False)

# ---------------------------------------------------------------- F. the close-range peak: is there a preferred abutting spacing?
# The AP-1 motif is 11 bp (-5..+5), so a partner centred 11-25 bp away sits just outside it, edge to edge.
NEARW, FARW = (11, 26), (50, 151)
rows = []
ref = M[M.family == "KLF / SP"].dist.values                      # the abundant, unimportant family: the shape expected with no preference
ref_ratio = ((ref >= NEARW[0]) & (ref < NEARW[1])).sum() / max(((ref >= FARW[0]) & (ref < FARW[1])).sum(), 1) * (FARW[1] - FARW[0]) / (NEARW[1] - NEARW[0])
for fam_ in KEY:
    d = M[M.family == fam_].dist.values
    if len(d) < 200:
        continue
    n_near = int(((d >= NEARW[0]) & (d < NEARW[1])).sum())
    n_far = int(((d >= FARW[0]) & (d < FARW[1])).sum())
    dens_near = n_near / (NEARW[1] - NEARW[0])
    dens_far = n_far / (FARW[1] - FARW[0])
    # permutation: shuffle each instance's distance within [11,151) keeping the family's own far-range shape
    null = []
    pool = d[(d >= NEARW[0]) & (d < FARW[1])]
    for k in range(2000):
        s_ = rng.choice(pool, len(pool), replace=True)
        null.append(((s_ >= NEARW[0]) & (s_ < NEARW[1])).sum() / (NEARW[1] - NEARW[0]) / max(((s_ >= FARW[0]) & (s_ < FARW[1])).sum() / (FARW[1] - FARW[0]), 1e-9))
    rows.append({"family": fam_, "n_instances": len(d), "n_11_25bp": n_near, "n_50_150bp": n_far, "density_ratio_near_over_far": dens_near / dens_far if dens_far else np.nan,
                 "ratio_for_KLF_SP_reference": ref_ratio, "ci_lo": np.percentile(null, 2.5), "ci_hi": np.percentile(null, 97.5),
                 "p_vs_KLF_SP": (np.sum(np.array(null) <= ref_ratio) + 1) / 2001})
F = pd.DataFrame(rows).sort_values("density_ratio_near_over_far", ascending=False)
F["fdr"] = bh(F.p_vs_KLF_SP)
F.to_csv(OUT / "tables" / "f_close_range_peak.csv", index=False)

# 1 bp profile of the close range, for the figure and to locate the mode
fine = []
for fam_ in ("MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD", "KLF / SP"):
    d = M[M.family == fam_].dist.values
    if len(d) < 200:
        continue
    h = np.histogram(d, bins=np.arange(0, 61))[0].astype(float)
    for b_, k_ in zip(np.arange(0.5, 60), h):
        fine.append({"family": fam_, "distance": b_, "n": int(k_), "share": k_ / len(d)})
FINE = pd.DataFrame(fine)
FINE.to_csv(OUT / "tables" / "f_close_range_profile.csv", index=False)
mode_rows = []
for fam_, w in FINE.groupby("family"):
    w2 = w[(w.distance > 8) & (w.distance < 40)]
    mode_rows.append({"family": fam_, "modal_distance_9_40bp": float(w2.loc[w2.n.idxmax(), "distance"]), "n_at_mode": int(w2.n.max())})
pd.DataFrame(mode_rows).to_csv(OUT / "tables" / "f_modal_distance.csv", index=False)
print("== F close-range peak (density 11-25 bp relative to 50-150 bp)\n", F.round(3).to_string())
print("== F modal distance\n", pd.DataFrame(mode_rows).to_string())

# ---------------------------------------------------------------- figures
fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.2))
ax = axes[0]
for fam_, c in zip(["MEF2", "GATA", "NFI (Nfia/b/c/x)", "TEAD", "KLF / SP"], [ORANGE, BLUE, AQUA, "#8e44ad", GREY]):
    w = PROF[PROF.family == fam_]
    if len(w):
        ax.plot(w.distance, w.enrichment_vs_all_footprints, color=c, lw=1.7, label=fam_.split(" (")[0])
ax.axhline(1, color=MUTED, lw=0.9, ls="--")
ax.set_xlabel("Distance from the AP-1 motif centre (bp)")
ax.set_ylabel("Enrichment against all footprints\n(1 = same distance profile)")
ax.set_title("Where do the partners sit?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
ax = axes[1]
for fam_, c in zip(["MEF2", "GATA", "NFI (Nfia/b/c/x)", "KLF / SP"], [ORANGE, BLUE, AQUA, GREY]):
    d = M[M.family == fam_].dist.values
    if len(d) < 200:
        continue
    h = np.histogram(d[d < 150], bins=np.arange(0, 151, 3))[0].astype(float)
    ax.plot(np.arange(1.5, 150, 3), h / h.sum(), color=c, lw=1.6, label=fam_.split(" (")[0])
ax.set_xlabel("Distance from the AP-1 motif centre (bp)")
ax.set_ylabel("Share of that family's instances")
ax.set_title("Fine structure within 150 bp (3 bp bins)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "figures" / "r1_distance_profile.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.6, 4.0))
for fam_, c in (("MEF2", ORANGE), ("GATA", BLUE), ("NFI (Nfia/b/c/x)", AQUA), ("KLF / SP", GREY)):
    w = FINE[FINE.family == fam_]
    if len(w):
        ax.plot(w.distance, w.share.rolling(3, center=True, min_periods=1).mean(), color=c, lw=1.8, label=f"{fam_.split(' (')[0]} (n = {int(w.n.sum()):,})")
ax.axvspan(11, 25, color=ORANGE, alpha=0.1, lw=0)
ax.text(18, ax.get_ylim()[1] * 0.95, "abutting", fontsize=8, color=MUTED, ha="center")
ax.axvspan(0, 5.5, color=GRID, alpha=0.8, lw=0)
ax.text(2.7, ax.get_ylim()[1] * 0.95, "overlap", fontsize=8, color=MUTED, ha="center")
ax.set_xlabel("Distance between motif centres (bp), 1 bp bins, 3 bp rolling mean")
ax.set_ylabel("Share of that family's instances")
ax.set_title("Close range: a preferred abutting spacing for the important partners", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "figures" / "r3_close_range.png", bbox_inches="tight")
plt.close(fig)

if SPEC:
    labs = [l for l in ("MEF2 + GATA + NFI + TEAD", "MEF2", "KLF / SP") if l in SPEC]
    fig, axes = plt.subplots(1, len(labs), figsize=(4.4 * len(labs), 3.8), squeeze=False)
    for ax, lab in zip(axes[0], labs):
        per, pw, band, null = SPEC[lab]
        sel = (per >= 5) & (per <= 30)
        ax.plot(per[sel], pw[sel], color=INK, lw=1.5)
        ax.axvspan(9.5, 11.5, color=ORANGE, alpha=0.15, lw=0)
        ax.axhline(np.percentile(null, 95), color=ORANGE, lw=1.2, ls=":", label="95th percentile of jittered nulls")
        ax.set_xlabel("Period (bp)")
        ax.set_ylabel("Spectral power")
        r = PER[PER.family == lab].iloc[0]
        ax.set_title(f"{lab}\nn = {int(r.n_instances):,}, p = {r.p_vs_jitter_null:.2f}", fontsize=9.5, loc="left")
        ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "r2_periodicity.png", bbox_inches="tight")
    plt.close(fig)

pd.set_option("display.width", 240)
print("== B distance summary\n", NEAR.round(3).to_string())
print("== C periodicity\n", PER.round(3).to_string())
if IMP is not None:
    print("== D importance by distance\n", IMP.round(3).to_string())
    print("== D distance-importance correlation\n", RHO.round(3).to_string())
print("== E class differences in spacing\n", E.round(3).to_string())
