#!/usr/bin/env python3
"""Fibroblast site sets for the fine-tuned-head experiments (labels used exactly as given in the site table; seed 46):
  Katana/data/fb_subset_scan.csv   tile scan: every Dynamic site (group 'dynamic') + the same number of Stable sites drawn at random (group 'random_stable')
  Katana/data/fb_swap_sites.csv    peak swap: 50 sites per binding class (column swap_class): Stable (bound in all three conditions, stable_dynamic),
                                   Age Up / Age Down (age_class), Exercise Up / Exercise Down (ex_class). A site can carry both an age and an exercise
                                   class (e.g. MAS only = Age Up and Exercise Down), so each site is drawn once, rarest class first.
Rows are shuffled so every PBS chunk gets every group.   python3 46_make_fb_subsets.py"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED, N_SWAP = 46, 50
s = pd.read_csv(ROOT / "Katana" / "data" / "fb_fos_sites.csv")
s = s.drop(columns=[c for c in ("YS_score", "MAS_score", "MAE_score") if c in s])   # no TOBIAS footprint scores in these experiments

dyn = s[s.stable_dynamic == "Dynamic"]
stab = s[s.stable_dynamic == "Stable"].sample(n=len(dyn), random_state=SEED)
scan = pd.concat([dyn.assign(group="dynamic"), stab.assign(group="random_stable")]).sample(frac=1, random_state=SEED).reset_index(drop=True)
scan.to_csv(ROOT / "Katana" / "data" / "fb_subset_scan.csv", index=False)
print("scan subset:", scan.group.value_counts().to_dict())

classes = {"Age Down": s.age_class == "Age Down", "Exercise Up": s.ex_class == "Exercise Up", "Exercise Down": s.ex_class == "Exercise Down",
           "Age Up": s.age_class == "Age Up", "Stable": s.stable_dynamic == "Stable"}
order = sorted(classes, key=lambda c: int(classes[c].sum()))                    # rarest first, so overlaps do not starve a small class
taken, parts = set(), []
for c in order:
    pool = s[classes[c] & ~s.site_id.isin(taken)]
    pick = pool.sample(n=N_SWAP, random_state=SEED)
    taken |= set(pick.site_id)
    parts.append(pick.assign(swap_class=c))
    print(f"{c}: {int(classes[c].sum())} sites, {len(pool)} not yet taken, drew {len(pick)}")
swap = pd.concat(parts).sample(frac=1, random_state=SEED).reset_index(drop=True)
swap.to_csv(ROOT / "Katana" / "data" / "fb_swap_sites.csv", index=False)
print("swap sites:", swap.swap_class.value_counts().to_dict())
print(pd.crosstab(swap.swap_class, swap.binding_category))
