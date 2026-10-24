#!/usr/bin/env python3
"""One message, one panel: the projector versions of the fibroblast results for the lab meeting.

The analysis figures keep the caveats; these keep only the shape of the result, so it reads from the back of the room.
Every slide's wording, n and p live in the deck (build_fb_deck.js), not in the image.

  s_fb_transplant   elements moved into closed chromatin, Stable vs Dynamic        Analysis/fb_swap
  s_fb_scale        what sets a swapped site's level: the element or the place     Analysis/fb_swap
  s_fb_specificity  own AP-1 vs partner motifs, by age class                       Analysis/fb_partner_shuffle
  s_fb_share        how much of a site's motif dependence sits on AP-1             Analysis/fb_partner_shuffle

Writes into Analysis/slide_figures/.   python3 Scripts/slides/make_fb_simple_figures.py
"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Analysis" / "slide_figures"
OUT.mkdir(parents=True, exist_ok=True)
INK, MUTED, SURF, GRID = "#111111", "#595959", "#ffffff", "#e4e3df"
STABLE, ORANGE, BLUE, CTRL = "#4F4F4F", "#C2551C", "#1F6FB2", "#B3B2AE"
plt.rcParams.update({"font.size": 16, "axes.facecolor": SURF, "figure.facecolor": SURF, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
                     "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor=SURF)
    plt.close(fig)
    print("  wrote", OUT / f"{name}.png")


def lg(x):
    return np.log2(np.asarray(x, float).clip(min=1e-3))


# ================================================================ 1. transplanted into closed chromatin
files = sorted(glob.glob(str(ROOT / "Katana" / "results_ft" / "bg_fb_chunk*.csv")))
if files:
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d = d[d.trackset == "FB"]
    piv = d.pivot_table(index=["donor_id", "donor_class", "host_id", "width"], columns="kind",
                        values="atac_insert").reset_index()
    piv["lift"] = lg(piv.element) - lg(piv.shuffled)
    # the donor is the replicate: one median per donor over its 20 hosts, then one median per class
    per = piv[piv.width == "peak"].groupby(["donor_id", "donor_class"], as_index=False).lift.median()
    per["grp"] = np.where(per.donor_class == "Stable", "Stable", "Dynamic")
    per.to_csv(OUT / "s_fb_transplant_per_donor.csv", index=False)
else:
    print("  no bg_fb_chunk*.csv locally - transplant figure skipped")
    per = None

if per is not None:
    fig, ax = plt.subplots(figsize=(7.6, 3.5))
    for yi, g, col in ((1, "Stable", STABLE), (0, "Dynamic", ORANGE)):
        v = per[per.grp == g].lift
        ax.barh(yi, v.median(), height=0.52, color=col, zorder=3)
        ax.text(v.median() - 0.18, yi, f"{2 ** v.median():.0f}$\\times$", va="center", ha="right",
                fontsize=27, fontweight="bold", color="#ffffff", zorder=5)
    ax.set_yticks([1, 0], ["Stable", "Dynamic"], fontsize=19)
    ax.set_xticks(np.arange(0, 5.1, 1), ["1", "2", "4", "8", "16", "32"])
    ax.set_xlim(0, 5.1)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xlabel("fold opening of a closed site, after the element is moved in", fontsize=15)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    save(fig, "s_fb_transplant")

# ================================================================ 2. the element or the place?
sl = pd.read_csv(ROOT / "Analysis" / "fb_swap" / "tables" / "sequence_vs_location.csv")
fig, ax = plt.subplots(figsize=(7.6, 4.0))
w, xs = 0.34, np.arange(2)
seq = [sl[sl.width == "170"].beta_donor_sequence.iloc[0], sl[sl.width == "peak"].beta_donor_sequence.iloc[0]]
loc = [sl[sl.width == "170"].beta_recipient_location.iloc[0], sl[sl.width == "peak"].beta_recipient_location.iloc[0]]
ax.bar(xs - w / 2, seq, w, color=BLUE, zorder=3, label="what the element IS")
ax.bar(xs + w / 2, loc, w, color=CTRL, zorder=3, label="WHERE it sits")
for x, v, c in list(zip(xs - w / 2, seq, [INK] * 2)) + list(zip(xs + w / 2, loc, [MUTED] * 2)):
    ax.text(x, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=17, color=c)
ax.set_xticks(xs, ["a 170 bp element\n(the MPRA window)", "the whole peak\n(median 0.8 kb)"], fontsize=16)
ax.set_ylabel("sets the swapped site's\naccessibility", fontsize=14.5)
ax.set_ylim(0, 0.92)
ax.set_yticks([0, 0.25, 0.5, 0.75])
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.14), ncol=2, fontsize=14.5, handlelength=1.1)
ax.grid(axis="y", color=GRID, lw=0.8)
ax.set_axisbelow(True)
save(fig, "s_fb_scale")

# ================================================================ 3. own AP-1 vs partners, by class
PS = ROOT / "Analysis" / "fb_partner_shuffle" / "tables"
J = pd.read_csv(PS / "own_vs_partners_per_site.csv")
GR = [("Stable", "Stable"), ("AP-1 gain", "AP-1 gained\nwith age"), ("AP-1 loss", "AP-1 lost\nwith age")]
fig, ax = plt.subplots(figsize=(7.6, 4.3))
x, w = np.arange(len(GR)), 0.34
own = [J[J.grp == g].own_ap1.median() for g, _ in GR]
par = [J[J.grp == g].partners_mean_net.median() for g, _ in GR]
ax.bar(x - w / 2, own, w, color=ORANGE, zorder=3, label="its own AP-1 motif")
ax.bar(x + w / 2, par, w, color=CTRL, zorder=3, label="its partner motifs")
for xi, v in zip(x - w / 2, own):
    ax.text(xi, v + 0.07, f"{v:.2f}", ha="center", va="bottom", fontsize=16, color="#ffffff", zorder=5)
for xi, v in zip(x + w / 2, par):
    ax.text(xi, v - 0.06, f"{v:.2f}", ha="center", va="top", fontsize=15, color=MUTED)
ax.axhline(0, color=MUTED, lw=1.2)
ax.set_xticks(x, [lab for _, lab in GR], fontsize=16)
ax.set_ylabel("effect of shuffling it  (log$_2$)", fontsize=14.5)
ax.set_ylim(-1.42, 0.42)
ax.set_yticks([-1.2, -0.8, -0.4, 0])
ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.17), ncol=2, fontsize=14.5, handlelength=1.1)
ax.grid(axis="y", color=GRID, lw=0.8)
ax.set_axisbelow(True)
save(fig, "s_fb_specificity")

# ================================================================ 4. how AP-1-centric is a site?
J["ap1_share"] = J.own_ap1.abs() / (J.own_ap1.abs() + J.partners_mean_net.abs())
J["grp2"] = np.where(J.grp == "Stable", "Stable", "Dynamic")
fig, ax = plt.subplots(figsize=(6.4, 4.3))
data = [J[J.grp2 == g].ap1_share.dropna().values for g in ("Stable", "Dynamic")]
bp = ax.boxplot(data, positions=[0, 1], widths=0.5, patch_artist=True, showfliers=False,
                medianprops=dict(color="#ffffff", lw=2.6), whiskerprops=dict(color=MUTED, lw=1.2),
                capprops=dict(color=MUTED, lw=1.2), boxprops=dict(lw=0))
for patch, col in zip(bp["boxes"], (STABLE, ORANGE)):
    patch.set_facecolor(col)
for xi, v in enumerate(data):
    ax.text(xi, np.median(v) + 0.035, f"{100 * np.median(v):.0f}%", ha="center", va="bottom",
            fontsize=21, fontweight="bold", color=INK)
ax.set_xticks([0, 1], ["Stable", "Dynamic"], fontsize=19)
ax.set_ylabel("share of the site's motif\ndependence that is AP-1", fontsize=14.5)
ax.set_ylim(0.25, 1.06)
ax.set_yticks([0.4, 0.6, 0.8, 1.0], ["40%", "60%", "80%", "100%"])
ax.grid(axis="y", color=GRID, lw=0.8)
ax.set_axisbelow(True)
save(fig, "s_fb_share")
print("done")
