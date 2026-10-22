#!/usr/bin/env python3
"""Which transcription factor binding sites does AlphaGenome treat as important around an AP-1 site, and do different factors matter to different classes?
Maps every footprinted cardiomyocyte TFBS of the FULL TOBIAS table onto the tile scan and saturation mutagenesis of 15_scan_titrate.py.
 A  tile importance per TF / family, adjusted for the site's overall sensitivity and for distance from the AP-1 motif (two-way median adjustment)
 B  base-resolution check: |effect| of single-base changes inside vs outside each factor's footprint (+/-30 bp window)
 C  do different factors matter to Dynamic and to context-matched Stable sites? (importance and prevalence)
TOBIAS is used here for the LOCATION of bound motifs only; footprint scores are not used.   Inputs: Data/derived/tobias_cm_all_motifs.csv.gz (made by this script if missing), Katana/results/scan*_cm_chunk*.csv.
Outputs Analysis/partners/{tables,figures}.   python3 38_tfbs_partner_mapping.py"""
from __future__ import annotations

import glob
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Analysis" / "partners"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(38)
INK, MUTED, SURF, GRID, BLUE, ORANGE, AQUA, GREY = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df", "#2a78d6", "#eb6834", "#1baf7a", "#9a9993"
plt.rcParams.update({"font.size": 9.5, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "savefig.dpi": 200})
HALF, TILE = 300, 10


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def boot_med(v, n=400):
    v = np.asarray(v, float)
    return np.percentile([np.median(rng.choice(v, len(v))) for _ in range(n)], [2.5, 97.5])


# TF families: the name is the TOBIAS motif's gene symbol. Patterns are ordered, first match wins.
# Borcs8 is grouped with MEF2: it is the BORCS8-MEF2B readthrough, so the motif is MEF2B's under a non-TF gene symbol.
FAM = [("AP-1 (Fos/Jun/Atf/Batf)", r"^(Fos|Jun|Batf|Atf3|Bach|Nfe2|Maf)"), ("MEF2", r"^Mef2|^Borcs8"), ("GATA", r"^Gata|^Trps1"), ("TBX", r"^Tbx|^Tbr|^Eomes"),
       ("NKX / homeobox cardiac", r"^Nkx|^Hand|^Msx|^Pitx|^Isl"), ("NFAT", r"^Nfat"), ("NFI (Nfia/b/c/x)", r"^Nfi[abcx]"), ("TALE homeodomain (Meis/Pknox/Tgif/Pbx)", r"^Meis|^Pknox|^Tgif|^Pbx"), ("TEAD", r"^Tead"), ("SRF / MADS", r"^Srf"), ("ETS", r"^Ets|^Elk|^Erg|^Fli1|^Gabp|^Elf|^Ehf|^Spdef|^Etv"),
       ("RUNX", r"^Runx"), ("KLF / SP", r"^Klf|^Sp[0-9]|^Patz|^Zfp(281|341|383|467|740)|^Wt1|^Egr"), ("CTCF / cohesin", r"^Ctcf|^Ctcfl"), ("E-box bHLH (Myod/Tcf/Neuro)", r"^Myod|^Myog|^Myf|^Tcf(3|4|12)|^Neuro|^Ascl|^Ptf1|^Bhlh"),
       ("Nuclear receptor (Esrr/Nr/Ppar/Thr)", r"^Esrr|^Nr[0-9]|^Ppar|^Rxr|^Rar|^Thr|^Vdr|^Ar$|^Nr3|^Esr|^Hnf4"), ("MYC / MAX", r"^Myc|^Max|^Mxd|^Mnt|^Usf|^Srebf|^Tfe|^Mitf"),
       ("FOX", r"^Fox"), ("SOX", r"^Sox"), ("STAT / IRF", r"^Stat|^Irf"), ("NF-kB / Rel", r"^Rel|^Nfkb"), ("CEBP", r"^Cebp|^Ddit3"), ("Zinc finger, other", r"^Zfp|^Znf|^Zbtb|^Zic|^Zeb|^Snai"),
       ("Smad / Tgfb", r"^Smad"), ("RFX", r"^Rfx"), ("YY / Gli / Pax", r"^Yy[0-9]|^Gli|^Pax"), ("Heat shock / stress", r"^Hsf|^Xbp|^Creb|^Atf[0-9]")]


def family(name):
    for lab, rx in FAM:
        if re.match(rx, str(name), re.I):
            return lab
    return "Other"


# ---------------------------------------------------------------- inputs
TOB = ROOT / "Data/derived/tobias_cm_all_motifs.csv.gz"
if not TOB.exists():
    cols = ["seqnames", "start", "end", "in_MAE", "in_MAS", "in_YS", "MAS_score", "MAE_score", "YS_score", "category", "CellType", "name", "peak_coord"]
    it = pd.read_csv(ROOT / "Data/Tobias_df_Classified_Full.csv.gz", usecols=cols, chunksize=1_000_000, dtype={"name": "category", "CellType": "category", "category": "category"})
    pd.concat([c[c.CellType == "Cardiomyocyte"] for c in it], ignore_index=True).drop(columns=["CellType"]).to_csv(TOB, index=False, compression="gzip")
tob = pd.read_csv(TOB, usecols=["seqnames", "start", "end", "in_MAE", "in_MAS", "in_YS", "name"], dtype={"name": "category"})
tob["family"] = tob.name.map(family).astype("category")

sc = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(ROOT / "Katana/results/scan*_cm_chunk*.csv")))], ignore_index=True)
sc = sc[sc.trackset == "CM"]
n_pass = sc.groupby("site_id").size()
sc = sc[sc.site_id.isin(n_pass.index[n_pass == n_pass.max()])]                    # complete sites only
sites = pd.read_csv(ROOT / "Data/derived/cm_fos_sites.csv").set_index("site_id")
sc = sc.join(sites[["seqnames", "motif_mid"]], on="site_id")

EFF = {}
for rd, col in (("170", "atac_170"), ("motif", "atac_motif")):
    nat = sc[sc.kind == "native"].set_index("site_id")[col]
    pc = 0.05 * nat.median()
    e = sc[sc.kind != "native"].copy()
    e["eff"] = np.log2((e[col] + pc) / (e.site_id.map(nat) + pc))
    EFF[rd] = e
print(f"scan: {sc.site_id.nunique()} complete sites; TOBIAS cardiomyocyte footprints: {len(tob):,} in {tob.name.nunique()} motifs, {tob.family.nunique()} families")

# ---------------------------------------------------------------- overlap: every footprint within +/-HALF bp of each site
idx = {c: g.sort_values("start").reset_index(drop=True) for c, g in tob.groupby("seqnames", observed=True)}
meta = sc.drop_duplicates("site_id").set_index("site_id")[["group", "stable_dynamic", "seqnames", "motif_mid"]]
hits = []
for sid, r in meta.iterrows():
    g = idx.get(r.seqnames)
    if g is None:
        continue
    lo, hi = r.motif_mid - HALF - 30, r.motif_mid + HALF + 30
    a = np.searchsorted(g.start.values, hi)
    w = g.iloc[max(np.searchsorted(g.start.values, lo - 60), 0):a]
    w = w[w.end.values > lo]
    if len(w):
        hits.append(pd.DataFrame({"site_id": sid, "name": w.name.values, "family": w.family.values,
                                  "rel_start": w.start.values - r.motif_mid, "rel_end": w.end.values - r.motif_mid}))
H = pd.concat(hits, ignore_index=True)
H["own_ap1"] = (H.rel_start < 6) & (H.rel_end > -5) & (H.family == "AP-1 (Fos/Jun/Atf/Batf)")      # the site's own motif
H.to_csv(OUT / "tables" / "footprints_near_sites.csv", index=False)
per_site = H.groupby("site_id").size()
print(f"footprints within +/-{HALF} bp: {len(H):,}; median {per_site.median():.0f} per site")

# ---------------------------------------------------------------- A. tile importance, adjusted for site and distance
t = EFF["170"][EFF["170"].kind == "tile"].groupby(["site_id", "group", "pos"], as_index=False).eff.mean()
t["dist"] = (t.pos + 5).abs()
t["dbin"] = pd.cut(t.dist, [-1, 10, 25, 50, 85, 130, 200, 305], labels=["0-10", "11-25", "26-50", "51-85", "86-130", "131-200", "201-300"])
t["adj"] = t.eff - t.groupby("site_id").eff.transform("median") - t.groupby("dbin", observed=True).eff.transform("median") + t.eff.median()

tiles = t.set_index(["site_id", "pos"])
pairs = []                                                      # tile x family membership (a tile can carry several families)
for (sid, fam_), w in H[~H.own_ap1].groupby(["site_id", "family"], observed=True):
    if sid not in t.site_id.values:
        continue
    cov = np.zeros(2 * HALF // TILE, bool)
    for s_, e_ in zip(w.rel_start, w.rel_end):
        a, b = max((s_ + HALF) // TILE, 0), min((e_ - 1 + HALF) // TILE, len(cov) - 1)
        if b >= 0 and a < len(cov):
            cov[max(a, 0):b + 1] = True
    for k in np.flatnonzero(cov):
        pairs.append((sid, k * TILE - HALF, fam_))
P = pd.DataFrame(pairs, columns=["site_id", "pos", "family"])
TP = P.merge(t, on=["site_id", "pos"], how="inner")
own = P.site_id.nunique()
covered = TP.set_index(["site_id", "pos"]).index.unique()
t_free = t.set_index(["site_id", "pos"]).drop(index=covered, errors="ignore").reset_index()     # tiles with no footprint of any factor
print(f"tiles: {len(t):,}; with >=1 footprint {len(covered):,}; footprint-free {len(t_free):,}")

rows = []
for fam_, w in TP.groupby("family", observed=True):
    if len(w) < 40:
        continue
    lo, hi = boot_med(w.adj.values)
    rows.append({"family": fam_, "n_tiles": len(w), "n_sites": w.site_id.nunique(), "median_adj_effect": w.adj.median(), "ci_lo": lo, "ci_hi": hi,
                 "median_raw_effect": w.eff.median(), "diff_vs_footprint_free": w.adj.median() - t_free.adj.median(), "p": mannwhitneyu(w.adj, t_free.adj).pvalue})
A = pd.DataFrame(rows).sort_values("median_adj_effect")
A["fdr"] = bh(A.p)
A.to_csv(OUT / "tables" / "a_tile_importance_by_family.csv", index=False)

rows = []                                                        # same, per individual TF
TFP = P.copy()
for nm, w in H[~H.own_ap1].groupby("name", observed=True):
    pass
tf_pairs = []
for (sid, nm), w in H[~H.own_ap1].groupby(["site_id", "name"], observed=True):
    cov = np.zeros(2 * HALF // TILE, bool)
    for s_, e_ in zip(w.rel_start, w.rel_end):
        a, b = max((s_ + HALF) // TILE, 0), min((e_ - 1 + HALF) // TILE, len(cov) - 1)
        if b >= 0 and a < len(cov):
            cov[max(a, 0):b + 1] = True
    for k in np.flatnonzero(cov):
        tf_pairs.append((sid, k * TILE - HALF, nm))
TF = pd.DataFrame(tf_pairs, columns=["site_id", "pos", "name"]).merge(t, on=["site_id", "pos"])
for nm, w in TF.groupby("name", observed=True):
    if len(w) < 60:
        continue
    rows.append({"TF": nm, "family": family(nm), "n_tiles": len(w), "n_sites": w.site_id.nunique(), "median_adj_effect": w.adj.median(),
                 "diff_vs_footprint_free": w.adj.median() - t_free.adj.median(), "p": mannwhitneyu(w.adj, t_free.adj).pvalue})
ATF = pd.DataFrame(rows).sort_values("median_adj_effect")
ATF["fdr"] = bh(ATF.p)
ATF.to_csv(OUT / "tables" / "a_tile_importance_by_tf.csv", index=False)

# ---------------------------------------------------------------- B. base resolution: inside vs outside a footprint
im = EFF["170"][EFF["170"].kind == "ism"].groupby(["site_id", "group", "pos"], as_index=False).eff.agg(imp=lambda v: np.abs(v).mean())
im["dist"] = im.pos.abs()
im["adj"] = im.imp - im.groupby("site_id").imp.transform("median") - im.groupby("dist").imp.transform("median") + im.imp.median()
cover = {}
for sid, w in H[~H.own_ap1].groupby("site_id", observed=True):
    for s_, e_, fam_ in zip(w.rel_start, w.rel_end, w.family):
        for p_ in range(max(s_, -30), min(e_, 31)):
            cover.setdefault((sid, p_), set()).add(fam_)
im["fams"] = [cover.get((s_, p_), set()) for s_, p_ in zip(im.site_id, im.pos)]
free = im[im.fams.map(len) == 0]
rows = []
for fam_ in A.family:
    w = im[im.fams.map(lambda s, f=fam_: f in s)]
    if len(w) < 40:
        continue
    rows.append({"family": fam_, "n_positions": len(w), "n_sites": w.site_id.nunique(), "median_importance": w.imp.median(), "median_adj": w.adj.median(),
                 "diff_vs_footprint_free": w.adj.median() - free.adj.median(), "p": mannwhitneyu(w.adj, free.adj).pvalue})
B = pd.DataFrame(rows).sort_values("median_adj", ascending=False)
B["fdr"] = bh(B.p)
B.to_csv(OUT / "tables" / "b_base_importance_by_family.csv", index=False)

# ---------------------------------------------------------------- C. Dynamic vs context-matched Stable
G = {"dynamic": "Dynamic", "matched_stable": "Stable, context-matched"}
rows = []
for fam_, w in TP[TP.group.isin(G)].groupby("family", observed=True):
    a, b_ = w[w.group == "dynamic"], w[w.group == "matched_stable"]
    if min(len(a), len(b_)) < 25:
        continue
    sa = set(H[(H.family == fam_) & ~H.own_ap1].site_id)
    pa = meta[meta.group == "dynamic"].index.isin(list(sa)), meta[meta.group == "matched_stable"].index.isin(list(sa))
    tab = [[pa[0].sum(), (~pa[0]).sum()], [pa[1].sum(), (~pa[1]).sum()]]
    rows.append({"family": fam_, "n_tiles_dynamic": len(a), "n_tiles_stable": len(b_), "median_adj_dynamic": a.adj.median(), "median_adj_stable": b_.adj.median(),
                 "importance_diff": a.adj.median() - b_.adj.median(), "p_importance": mannwhitneyu(a.adj, b_.adj).pvalue,
                 "share_sites_dynamic": pa[0].mean(), "share_sites_stable": pa[1].mean(), "p_prevalence": fisher_exact(tab)[1]})
C = pd.DataFrame(rows)
C["fdr_importance"], C["fdr_prevalence"] = bh(C.p_importance), bh(C.p_prevalence)
C = C.sort_values("importance_diff")
C.to_csv(OUT / "tables" / "c_dynamic_vs_stable_by_family.csv", index=False)

# how many distinct factors / families are footprinted near a site, by group
rows = []
for g, lab in G.items():
    s_ = meta[meta.group == g].index
    h = H[H.site_id.isin(s_) & ~H.own_ap1]
    rows.append({"group": lab, "n_sites": len(s_), "median_footprints_per_site": h.groupby("site_id").size().reindex(s_).fillna(0).median(),
                 "median_families_per_site": h.groupby("site_id").family.nunique().reindex(s_).fillna(0).median(),
                 "median_tfs_per_site": h.groupby("site_id").name.nunique().reindex(s_).fillna(0).median()})
D = pd.DataFrame(rows)
_d = H[H.site_id.isin(meta[meta.group == "dynamic"].index) & ~H.own_ap1].groupby("site_id").family.nunique().reindex(meta[meta.group == "dynamic"].index).fillna(0)
_s = H[H.site_id.isin(meta[meta.group == "matched_stable"].index) & ~H.own_ap1].groupby("site_id").family.nunique().reindex(meta[meta.group == "matched_stable"].index).fillna(0)
D["p_families_vs_dynamic"] = [np.nan, mannwhitneyu(_d, _s).pvalue]
D.to_csv(OUT / "tables" / "c_partner_counts_by_group.csv", index=False)

# ---------------------------------------------------------------- E. is the partner-count difference just detection? (footprints need signal)
meas = pd.read_csv(ROOT / "Katana/data/measured_atac_sites.csv")
meas = meas[meas.site_table == "cm"].set_index("site_id").cm_cpm_per_kb
fams_per = H[~H.own_ap1].groupby("site_id").family.nunique()
S = meta[meta.group.isin(G)].copy()
S["n_families"] = fams_per.reindex(S.index).fillna(0)
S["measured"] = meas.reindex(S.index)
S["q"] = pd.qcut(S.measured, 4, labels=["Q1 lowest", "Q2", "Q3", "Q4 highest"])
rows = [{"stratum": "all sites", "n_dynamic": int((S.group == "dynamic").sum()), "n_stable": int((S.group == "matched_stable").sum()),
         "median_families_dynamic": S[S.group == "dynamic"].n_families.median(), "median_families_stable": S[S.group == "matched_stable"].n_families.median(),
         "p": mannwhitneyu(S[S.group == "dynamic"].n_families, S[S.group == "matched_stable"].n_families).pvalue,
         "median_measured_dynamic": S[S.group == "dynamic"].measured.median(), "median_measured_stable": S[S.group == "matched_stable"].measured.median(),
         "p_measured": mannwhitneyu(S[S.group == "dynamic"].measured.dropna(), S[S.group == "matched_stable"].measured.dropna()).pvalue}]
for q, w in S.groupby("q", observed=True):
    a, b_ = w[w.group == "dynamic"], w[w.group == "matched_stable"]
    rows.append({"stratum": f"measured accessibility {q}", "n_dynamic": len(a), "n_stable": len(b_), "median_families_dynamic": a.n_families.median(), "median_families_stable": b_.n_families.median(),
                 "p": mannwhitneyu(a.n_families, b_.n_families).pvalue if min(len(a), len(b_)) > 5 else np.nan,
                 "median_measured_dynamic": a.measured.median(), "median_measured_stable": b_.measured.median(), "p_measured": np.nan})
E = pd.DataFrame(rows)
E["rho_families_vs_measured"] = [float(pd.Series(S.n_families).corr(S.measured, method="spearman"))] + [np.nan] * (len(E) - 1)
E.to_csv(OUT / "tables" / "e_partner_count_vs_detection.csv", index=False)

# ---------------------------------------------------------------- F. does a MEF2 partner change how much the site needs its own AP-1 motif?
MOTIF_TILES = {-10, 0}
own_dep = t[t.pos.isin(MOTIF_TILES)].groupby(["site_id", "group"]).eff.min().rename("own_motif_tile").reset_index()
tm = EFF["motif"][EFF["motif"].kind == "tile"].groupby(["site_id", "pos"], as_index=False).eff.mean()
own_dep["own_motif_tile_read_motif"] = own_dep.site_id.map(tm[tm.pos.isin(MOTIF_TILES)].groupby("site_id").eff.min())
rows = []
for fam_ in ("MEF2", "GATA", "TEAD"):
    has = set(H[(H.family == fam_) & ~H.own_ap1].site_id)
    for grp, lab in [("all", "all sites")] + [(g, l) for g, l in G.items()]:
        w = own_dep if grp == "all" else own_dep[own_dep.group == grp]
        a, b_ = w[w.site_id.isin(has)], w[~w.site_id.isin(has)]
        if min(len(a), len(b_)) < 15:
            continue
        rows.append({"partner": fam_, "sites": lab, "n_with": len(a), "n_without": len(b_), "own_dep_with_partner": a.own_motif_tile.median(), "own_dep_without": b_.own_motif_tile.median(),
                     "difference": a.own_motif_tile.median() - b_.own_motif_tile.median(), "p": mannwhitneyu(a.own_motif_tile, b_.own_motif_tile).pvalue,
                     "own_dep_with_read_motif": a.own_motif_tile_read_motif.median(), "own_dep_without_read_motif": b_.own_motif_tile_read_motif.median(),
                     "p_read_motif": mannwhitneyu(a.own_motif_tile_read_motif.dropna(), b_.own_motif_tile_read_motif.dropna()).pvalue})
F = pd.DataFrame(rows)
F["fdr"] = bh(F.p)
F.to_csv(OUT / "tables" / "f_own_dependence_with_partner.csv", index=False)

# where do the important partners sit?
pos_rows = []
for fam_ in ("MEF2", "GATA", "TEAD", "CTCF / cohesin", "AP-1 (Fos/Jun/Atf/Batf)"):
    w = H[(H.family == fam_) & ~H.own_ap1]
    c = (w.rel_start + w.rel_end) / 2
    pos_rows.append({"family": fam_, "n_instances": len(w), "median_abs_distance": c.abs().median(), "share_within_50bp": (c.abs() <= 50).mean(), "share_within_100bp": (c.abs() <= 100).mean()})
pd.DataFrame(pos_rows).to_csv(OUT / "tables" / "f_partner_positions.csv", index=False)

# ---------------------------------------------------------------- H. how much does the adjustment matter? four ways of asking the same question
#   raw: no adjustment.  site: the site's own median removed.  site+distance (used above).
#   distance-matched: each footprint tile against footprint-free tiles of the SAME site at a similar distance (no modelling, no pooling across sites).
#   far only: tiles beyond +/-85 bp, where the AP-1 motif's own dependence has died away.
free_by_site = {}
for (sid), w in t_free.groupby("site_id"):
    free_by_site[sid] = (w.pos.values + 5, w.eff.values)
dm = []
for r in TP.itertuples():
    f = free_by_site.get(r.site_id)
    if f is None:
        continue
    d_ = abs(r.pos + 5)
    sel = np.abs(np.abs(f[0]) - d_) <= 25
    if sel.sum() >= 3:
        dm.append((r.family, r.eff - f[1][sel].mean()))
DM = pd.DataFrame(dm, columns=["family", "d"])
t["site_adj"] = t.eff - t.groupby("site_id").eff.transform("median")
TPs = P.merge(t[["site_id", "pos", "eff", "site_adj", "adj", "dist"]], on=["site_id", "pos"])
rows = []
for fam_, w in TPs.groupby("family", observed=True):
    if len(w) < 40:
        continue
    far = w[w.dist > 85]
    d_ = DM[DM.family == fam_].d
    rows.append({"family": fam_, "n_tiles": len(w), "raw": w.eff.median(), "site_adjusted": w.site_adj.median(), "site_and_distance_adjusted": w.adj.median(),
                 "distance_matched_within_site": d_.median() if len(d_) >= 30 else np.nan, "n_matched": len(d_),
                 "far_tiles_only_raw": far.eff.median() if len(far) >= 25 else np.nan, "n_far": len(far)})
Hs = pd.DataFrame(rows).sort_values("site_and_distance_adjusted")
base = {"raw": t.eff.median(), "site_adjusted": t_free.eff.median() - t.groupby("site_id").eff.transform("median").median(),
        "site_and_distance_adjusted": t_free.adj.median(), "distance_matched_within_site": 0.0, "far_tiles_only_raw": t_free[t_free.dist > 85].eff.median()}
for c in ("raw", "site_adjusted", "site_and_distance_adjusted", "far_tiles_only_raw"):
    Hs[f"{c}_vs_free"] = Hs[c] - (t_free[c if c != "far_tiles_only_raw" else "eff"][t_free.dist > 85].median() if c == "far_tiles_only_raw" else t_free[c if c in t_free else "eff"].median())
Hs.to_csv(OUT / "tables" / "h_adjustment_sensitivity.csv", index=False)
METH = ["raw", "site_adjusted", "site_and_distance_adjusted", "distance_matched_within_site", "far_tiles_only_raw"]
rk = Hs[METH].rank()
corr = rk.corr(method="spearman")
corr.to_csv(OUT / "tables" / "h_rank_agreement.csv")
print("== H sensitivity of the family ranking to the adjustment\n", Hs[["family", "n_tiles"] + METH + ["n_matched", "n_far"]].round(3).to_string())
print("== H rank agreement between methods (Spearman)\n", corr.round(3).to_string())
print("== distance profile: share of each family's tiles within 85 bp of the motif\n",
      TPs.groupby("family", observed=True).apply(lambda w: pd.Series({"n": len(w), "share_within_85bp": (w.dist <= 85).mean()})).round(3).sort_values("share_within_85bp", ascending=False).head(12).to_string(),
      "\nall tiles:", round((t.dist <= 85).mean(), 3), " footprint-free tiles:", round((t_free.dist <= 85).mean(), 3))

# ---------------------------------------------------------------- I. class differences with the SITE as the unit (tiles within a site are not independent)
sc_cls = pd.read_csv(ROOT / "Analysis/score_classes/tables/score_based_classes.csv").set_index("site_id")
meta2 = meta.join(sc_cls[["age_score_class"]])
site_fam = TP.groupby(["site_id", "family"], observed=True).adj.median().rename("site_med").reset_index()     # one value per site per family
site_free = t_free.groupby("site_id").adj.median().rename("free_med")
site_fam = site_fam.join(site_free, on="site_id").join(meta2[["group", "age_score_class"]], on="site_id")
site_fam["vs_free"] = site_fam.site_med - site_fam.free_med                                                    # within-site contrast: this factor against that site's footprint-free tiles
site_fam.to_csv(OUT / "tables" / "i_site_level_family_importance.csv", index=False)
COMPARISONS = [("Dynamic vs context-matched Stable", "group", "dynamic", "matched_stable"),
               ("Dynamic vs high-accessibility Stable", "group", "dynamic", "high_stable"),
               ("Age Up (score) vs Unchanged", "age_score_class", "Age Up (score)", "Unchanged")]
rows = []
for lab, col, a_, b_ in COMPARISONS:
    for fam_, w in site_fam.groupby("family", observed=True):
        x, y = w[w[col] == a_].vs_free.dropna(), w[w[col] == b_].vs_free.dropna()
        if min(len(x), len(y)) < 8:
            continue
        pooled = np.std(pd.concat([x, y]))
        rows.append({"comparison": lab, "family": fam_, "n_sites_a": len(x), "n_sites_b": len(y), "median_a": x.median(), "median_b": y.median(),
                     "difference": x.median() - y.median(), "p": mannwhitneyu(x, y).pvalue,
                     "min_detectable_diff_80pct": 2.8 * pooled * np.sqrt(1 / len(x) + 1 / len(y))})
I = pd.DataFrame(rows)
I["fdr"] = I.groupby("comparison").p.transform(lambda v: bh(v))
I = I.sort_values(["comparison", "p"])
I.to_csv(OUT / "tables" / "i_class_differences_site_level.csv", index=False)

# CTCF: what drives the tile-level result?
ct = site_fam[site_fam.family == "CTCF / cohesin"]
ctr = []
for g, lab in (("dynamic", "Dynamic"), ("matched_stable", "Stable, context-matched"), ("high_stable", "Stable, high accessibility")):
    w = ct[ct.group == g].vs_free.dropna()
    if len(w):
        ctr.append({"group": lab, "n_sites": len(w), "median_vs_free": w.median(), "q25": w.quantile(0.25), "q75": w.quantile(0.75),
                    "share_sites_below_-0.3": (w < -0.3).mean()})
CT = pd.DataFrame(ctr)
CT.to_csv(OUT / "tables" / "i_ctcf_site_level.csv", index=False)

# prevalence, stratified by measured accessibility (the same detection confound as the partner count)
S2 = meta[meta.group.isin(G)].copy()
S2["measured"] = meas.reindex(S2.index)
S2["q"] = pd.qcut(S2.measured, 4, labels=["Q1", "Q2", "Q3", "Q4"])
prv = []
for fam_ in ("MEF2", "GATA", "TEAD", "NFI (Nfia/b/c/x)", "CTCF / cohesin", "AP-1 (Fos/Jun/Atf/Batf)", "MYC / MAX", "Nuclear receptor (Esrr/Nr/Ppar/Thr)"):
    has = set(H[(H.family == fam_) & ~H.own_ap1].site_id)
    S2["has"] = S2.index.isin(has)
    a_, b_ = S2[S2.group == "dynamic"], S2[S2.group == "matched_stable"]
    row = {"family": fam_, "share_dynamic": a_.has.mean(), "share_stable": b_.has.mean(),
           "p_unstratified": fisher_exact([[a_.has.sum(), (~a_.has).sum()], [b_.has.sum(), (~b_.has).sum()]])[1]}
    tabs = [np.array([[w[w.group == "dynamic"].has.sum(), (~w[w.group == "dynamic"].has).sum()],
                      [w[w.group == "matched_stable"].has.sum(), (~w[w.group == "matched_stable"].has).sum()]]) for _, w in S2.groupby("q", observed=True)]
    num = sum(tb[0, 0] * tb[1, 1] / tb.sum() for tb in tabs if tb.sum())          # Mantel-Haenszel, stratified by measured accessibility
    den = sum(tb[0, 1] * tb[1, 0] / tb.sum() for tb in tabs if tb.sum())
    row["odds_ratio_stratified"] = num / den if den else np.nan
    # Mantel-Haenszel chi-square (continuity-corrected) and a site-bootstrap CI for the stratified odds ratio
    obs = sum(tb[0, 0] for tb in tabs)
    exp = sum(tb[0].sum() * tb[:, 0].sum() / tb.sum() for tb in tabs if tb.sum())
    var = sum(tb[0].sum() * tb[1].sum() * tb[:, 0].sum() * tb[:, 1].sum() / (tb.sum() ** 2 * (tb.sum() - 1)) for tb in tabs if tb.sum() > 1)
    row["p_mantel_haenszel"] = float(chi2.sf(max(abs(obs - exp) - 0.5, 0) ** 2 / var, 1)) if var > 0 else np.nan
    bs = []
    for k in range(400):
        s_ = S2.sample(len(S2), replace=True, random_state=k)
        tb_ = [np.array([[w[w.group == "dynamic"].has.sum(), (~w[w.group == "dynamic"].has).sum()],
                         [w[w.group == "matched_stable"].has.sum(), (~w[w.group == "matched_stable"].has).sum()]]) for _, w in s_.groupby("q", observed=True)]
        n_ = sum(tb[0, 0] * tb[1, 1] / tb.sum() for tb in tb_ if tb.sum())
        d_ = sum(tb[0, 1] * tb[1, 0] / tb.sum() for tb in tb_ if tb.sum())
        if d_:
            bs.append(n_ / d_)
    row["or_ci_lo"], row["or_ci_hi"] = np.percentile(bs, [2.5, 97.5]) if len(bs) > 50 else (np.nan, np.nan)
    prv.append(row)
PRV = pd.DataFrame(prv)
PRV["fdr_unstratified"], PRV["fdr_stratified"] = bh(PRV.p_unstratified), bh(PRV.p_mantel_haenszel)
PRV.to_csv(OUT / "tables" / "i_prevalence_stratified.csv", index=False)

print("== I class differences, site as the unit\n", I.round(3).to_string())
print("== I CTCF per site\n", CT.round(3).to_string())
print("== I prevalence, stratified by measured accessibility (OR 1 = no difference)\n", PRV.round(3).to_string())

# ---------------------------------------------------------------- figures
top = A.head(14)
fig, ax = plt.subplots(figsize=(8.6, 5.4))
ys = np.arange(len(top))[::-1]
sig = (top.fdr < 0.05).values
ax.errorbar(top.median_adj_effect, ys, xerr=[top.median_adj_effect - top.ci_lo, top.ci_hi - top.median_adj_effect], fmt="none", ecolor=GREY, elinewidth=1.3)
ax.scatter(top.median_adj_effect[~sig], ys[~sig], s=38, color=GREY, zorder=3, label="not significant (FDR &ge; 0.05)".replace("&ge;", "≥"))
ax.scatter(top.median_adj_effect[sig], ys[sig], s=48, color=ORANGE, zorder=4, label="FDR < 0.05 vs footprint-free tiles")
ax.axvline(t_free.adj.median(), color=BLUE, lw=1.2, ls="--", label="tiles with no footprint")
ax.set_yticks(ys, [f"{f}  (n = {n:,})" for f, n in zip(top.family, top.n_tiles)], fontsize=8.5)
ax.set_xlabel("Median effect of replacing a 10 bp tile that overlaps this factor's footprint (log2FC)\nadjusted for the site's overall sensitivity and for distance from the AP-1 motif")
ax.set_title("Which bound factors does the model treat as important near an AP-1 site?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1.0, -0.02))
ax.set_ylim(-1.6, len(top) - 0.4)
fig.tight_layout()
fig.savefig(OUT / "figures" / "p1_tile_importance_by_family.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.6, 5.0))
cc = C.head(16)
ys = np.arange(len(cc))[::-1]
ax.scatter(cc.median_adj_stable, ys, s=44, color=BLUE, label="Stable, context-matched", zorder=3)
ax.scatter(cc.median_adj_dynamic, ys, s=44, color=ORANGE, label="Dynamic", zorder=3)
for y_, a_, b_ in zip(ys, cc.median_adj_dynamic, cc.median_adj_stable):
    ax.plot([a_, b_], [y_, y_], color=GREY, lw=1.1, zorder=2)
for y_, q in zip(ys, cc.fdr_importance):
    if q < 0.05:
        ax.text(max(cc.median_adj_dynamic.max(), cc.median_adj_stable.max()) + 0.03, y_, "FDR < 0.05", fontsize=7.5, va="center", color=INK)
ax.set_yticks(ys, cc.family, fontsize=8.5)
ax.axvline(t_free.adj.median(), color=MUTED, lw=1, ls=":")
ax.set_xlabel("Median adjusted tile effect (log2FC); dotted line = tiles with no footprint")
ax.set_title("Do different factors matter to Dynamic and to Stable sites?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "p2_family_by_class.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axes[0]
bb = B.head(12)[::-1]
ax.barh(np.arange(len(bb)), bb.median_adj - free.adj.median(), color=AQUA, alpha=0.85, height=0.62)
ax.set_yticks(np.arange(len(bb)), bb.family, fontsize=8)
ax.axvline(0, color=MUTED, lw=0.9)
ax.set_xlabel("Extra |effect| of a single-base change inside this factor's footprint\n(log2FC, vs bases with no footprint, adjusted)")
ax.set_title("Base resolution: saturation mutagenesis, ±30 bp", fontsize=10, loc="left")
ax = axes[1]
cp = C.sort_values("share_sites_stable", ascending=False).head(12)[::-1]
ys = np.arange(len(cp))
ax.barh(ys + 0.19, cp.share_sites_dynamic, height=0.36, color=ORANGE, alpha=0.85, label="Dynamic")
ax.barh(ys - 0.19, cp.share_sites_stable, height=0.36, color=BLUE, alpha=0.85, label="Stable, context-matched")
ax.set_yticks(ys, cp.family, fontsize=8)
ax.set_xlabel("Share of sites with this factor footprinted within ±300 bp")
ax.set_title("Which partners are present?", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.tight_layout()
fig.savefig(OUT / "figures" / "p3_base_and_prevalence.png", bbox_inches="tight")
plt.close(fig)

pd.set_option("display.width", 250)
print("== A tile importance by family\n", A.round(3).to_string())
print("== A top / bottom individual TFs\n", pd.concat([ATF.head(12), ATF.tail(5)]).round(3).to_string())
print("== B base-resolution\n", B.round(3).to_string())
print("== C by class\n", C.round(3).to_string())
print("== partner counts\n", D.round(3).to_string())
print("== E detection control\n", E.round(3).to_string())
print("== F own dependence with a partner\n", F.round(3).to_string())
print("== partner positions\n", pd.read_csv(OUT / "tables" / "f_partner_positions.csv").round(3).to_string())
