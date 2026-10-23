"""Importance of a partner footprint as a function of its distance from the AP-1 motif.
For each 10 bp tile we know its position; a tile is 'on' a family if it overlaps that family's footprint.
Within each distance bin we compare tiles on the family against footprint-free tiles at the same distance,
so the overall fall-off of importance with distance is removed."""
import glob, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
R = "/Users/z5162865/WongLab/AlphaGenome_Ap1"
sc = pd.concat([pd.read_csv(f, usecols=["site_id", "kind", "pos", "atac_170", "trackset"])
                for f in sorted(glob.glob(f"{R}/Katana/results/scan*_cm_chunk*.csv"))], ignore_index=True)
sc = sc[sc.trackset == "CM"]
nat = sc[sc.kind == "native"].groupby("site_id").atac_170.median()
pc = 0.05 * nat.median()
t = sc[sc.kind == "tile"].copy()
t["eff"] = np.log2((t.atac_170 + pc) / (t.site_id.map(nat) + pc))
T = t.groupby(["site_id", "pos"], as_index=False).eff.mean()
T["centre"] = T.pos + 5
FP = pd.read_csv(f"{R}/Analysis/partners/tables/footprints_near_sites.csv")
FP = FP[~FP.own_ap1]
rows = []
BINS = [(0, 15), (15, 30), (30, 50), (50, 85), (85, 130), (130, 200), (200, 305)]
# tile -> set of families it overlaps
for fam in ("MEF2", "KLF / SP", "GATA", "NFI (Nfia/b/c/x)", "TEAD"):
    f = FP[FP.family == fam]
    key = T.merge(f[["site_id", "rel_start", "rel_end"]], on="site_id", how="left")
    key["hit"] = (key.pos < key.rel_end) & (key.pos + 10 > key.rel_start)
    on = key[key.hit.fillna(False)][["site_id", "pos"]].drop_duplicates()
    T2 = T.merge(on.assign(on=True), on=["site_id", "pos"], how="left")
    T2["on"] = T2.on.fillna(False)
    # footprint-free tiles: no footprint of ANY family
    anyfp = T.merge(FP[["site_id", "rel_start", "rel_end"]], on="site_id", how="left")
    anyfp["hit"] = (anyfp.pos < anyfp.rel_end) & (anyfp.pos + 10 > anyfp.rel_start)
    free = anyfp.groupby(["site_id", "pos"]).hit.max().rename("anyhit").reset_index()
    T2 = T2.merge(free, on=["site_id", "pos"], how="left")
    T2["d"] = T2.centre.abs()
    for lo, hi in BINS:
        w = T2[(T2.d >= lo) & (T2.d < hi)]
        a = w[w.on].eff
        b = w[(~w.on) & (~w.anyhit.fillna(False))].eff
        if len(a) < 40 or len(b) < 40:
            continue
        rows.append({"family": fam, "bin": f"{lo}–{hi}", "lo": lo, "hi": hi, "mid": (lo + hi) / 2,
                     "n_on": len(a), "n_free": len(b), "med_on": a.median(), "med_free": b.median(),
                     "importance": a.median() - b.median()})
D = pd.DataFrame(rows)
D.to_csv(f"{R}/Analysis/distance_rules/tables/g_importance_by_distance.csv", index=False)
pd.set_option("display.width", 200)
print(D.round(3).to_string(index=False))
