#!/usr/bin/env python3
"""Does the binding-pattern result of 33_dynamic_vs_stable_deep_dive.py replicate in the fibroblast sites (heart and Patski tracks)? Same full context adjustment (splines on native accessibility at the readout, peak accessibility, width, TSS distance). Writes Analysis/dynamic/tables/fibroblast_replication_by_pattern.csv"""
import glob, numpy as np, pandas as pd, patsy, statsmodels.formula.api as smf
R="Katana/results/"
def load(prefix, ts):
    d=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(R+prefix+"_chunk*.csv"))]); return d[d.trackset==ts]
def eff(d,col):
    nat=d[d.condition=="native"].set_index("site_id")[col]; pc=0.05*nat.median()
    e=d[d.condition!="native"].copy(); e["l"]=np.log2((e[col]+pc)/(e.site_id.map(nat)+pc))
    return e.groupby(["site_id","condition"]).l.mean().unstack(), nat
bs=patsy.bs
def run(site_csv, shuf, mut, ts, tag):
    s=pd.read_csv(site_csv).set_index("site_id")
    out={}
    for rd,col in (("motif","atac_motif"),("170","atac_170"),("peak","atac_peak")):
        t,_=eff(shuf,col); m,_=eff(mut,col)
        out[f"shuffle {rd}→{rd}"]=t[f"shuffle@{rd}"]; out[f"replace {rd}→{rd}"]=m[f"random@{rd}"]; out[f"dinuc {rd}→{rd}"]=m[f"dinuc@{rd}"]
    out["replace core→motif"]=eff(mut,"atac_motif")[0]["core@motif"]; out["replace surroundings→170"]=eff(mut,"atac_170")[0]["flank@170"]
    F=pd.DataFrame(out)
    nat=shuf[shuf.condition=="native"].set_index("site_id")[["atac_motif","atac_170","atac_peak"]]
    X=s[["binding_category","stable_dynamic","peak_width","motif_mid","seqnames"]].join(nat).join(F,how="inner")
    X["l_motif"],X["l_170"],X["l_peak"]=[np.log10(X[c].clip(lower=1e-3)) for c in ("atac_motif","atac_170","atac_peak")]
    X["l_pw"]=np.log10(X.peak_width)
    tss=shuf.drop_duplicates("site_id").set_index("site_id").tss_dist.abs(); X["l_tss"]=np.log10(tss.reindex(X.index)+1)
    rows=[]
    for f in F.columns:
        cov="l_"+f.split("→")[-1] if f.split("→")[-1] in ("motif","170","peak") else "l_170"
        cov={"motif":"l_motif","170":"l_170","peak":"l_peak"}[f.split("→")[-1]]
        base=X[X.stable_dynamic=="Stable"]
        for pat in sorted(X[X.stable_dynamic=="Dynamic"].binding_category.unique())+["All Dynamic"]:
            sub=X[X.stable_dynamic=="Dynamic"] if pat=="All Dynamic" else X[X.binding_category==pat]
            x=pd.concat([base.assign(g=0),sub.assign(g=1)]).dropna(subset=[f,cov,"l_peak","l_pw","l_tss"])
            terms=[(cov,5)]+[(c,4) for c in ("l_peak","l_pw","l_tss") if c!=cov]
            r=smf.ols(f"Q('{f}') ~ g + "+" + ".join(f"bs({c},df={d})" for c,d in terms),x).fit(cov_type="HC3")
            rows.append({"track":tag,"mutation":f,"pattern":pat,"n":int(x.g.sum()),"adj_diff":r.params["g"],"lo":r.conf_int().loc["g",0],"hi":r.conf_int().loc["g",1],"p":r.pvalues["g"]})
    return pd.DataFrame(rows)
fb=run("Data/derived/fb_fos_sites.csv",load("shuffle_fb","FB"),load("mut_fb","FB"),"FB","FB sites, Patski")
fbh=run("Data/derived/fb_fos_sites.csv",load("shuffle_fb","CM"),load("mut_fb","CM"),"CM","FB sites, heart")
R_=pd.concat([fb,fbh]); R_.to_csv("Analysis/dynamic/tables/fibroblast_replication_by_pattern.csv",index=False)
pd.set_option("display.width",200)
for tag,g in R_.groupby("track"):
    print("\n",tag); print(g[g.pattern.isin(["MAS+MAE","All Dynamic"])&g.mutation.isin(["shuffle 170→170","replace 170→170","dinuc 170→170","replace surroundings→170","shuffle peak→peak","dinuc peak→peak","shuffle motif→motif","replace core→motif"])][["mutation","pattern","n","adj_diff","lo","hi","p"]].round(3).to_string(index=False))
    print("n significant (p<0.01):", (g.p<0.01).sum(), "of", len(g), "| MAS+MAE:", (g[g.pattern=="MAS+MAE"].p<0.01).sum(), "of", (g.pattern=="MAS+MAE").sum())
