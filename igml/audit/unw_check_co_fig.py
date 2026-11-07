"""The rebuild's unwrapped pairs across Ridgecrest and of one batch against the lab's (2026-11-02; unw_check_co.py's npz):
top, the shortest coseismic pair (ours, the lab's, ours - lab in cycles; the north-west quarter boxed); bottom, the
coseismic triplet where the lab's closure is worst (ours, the lab's) and every sampled triplet's share of pixels a cycle
off, ours against the lab's, coseismic and batch triplets apart. Radar coordinates of the merged grid (every 4th cell;
descending: north up, west right).
usage (Mac, insar_sim env): python unw_check_co_fig.py CHECK.npz OUT.png [BATCH_LABEL]"""
import sys
from datetime import date, timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

INK, GRID = "#1f1f1f", "#e6e6e6"
CO, BAT = "#2a78d6", "#eb6834"                  # categorical slots 1 / 2 of the dataviz reference palette
MM = 55.465763 / (4 * np.pi)                    # mm of LOS change per radian (C-band)
HB = dict(orientation="horizontal", fraction=0.045, pad=0.02, aspect=40)

z = np.load(sys.argv[1])
out = sys.argv[2]
blabel = sys.argv[3] if len(sys.argv) > 3 else "batch"


def ymd(d):
    return (date(int(d[:4]), 1, 1) + timedelta(days=int(d[4:]))).isoformat()


def los(u):
    v = -u * MM
    return v - np.nanmedian(v)


fig, ax = plt.subplots(2, 3, figsize=(15, 12.5), layout="constrained")
pa, pb = (ymd(x) for x in str(z["fig_pair"]).split("_"))
o, l = los(z["fig_ours"]), los(z["fig_lab"])
lim = np.nanpercentile(np.abs(np.r_[o[np.isfinite(o)], l[np.isfinite(l)]]), 99.5)
for a, v, name in ((ax[0, 0], o, "the rebuild"), (ax[0, 1], l, "the lab's")):
    im = a.imshow(v, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
    a.set_title(f"{name}: unwrapped {pa} x {pb[5:]}", fontsize=10, color=INK)
    a.add_patch(Rectangle((v.shape[1] / 2, 0), v.shape[1] / 2, v.shape[0] / 4, fill=False, ec=INK, lw=1, ls="--"))
fig.colorbar(im, ax=list(ax[0, :2]), **HB).set_label("LOS change [mm] (median removed)", fontsize=8.5)
d = z["fig_ours"] - z["fig_lab"]
d = (d - np.nanmedian(d)) / (2 * np.pi)
im = ax[0, 2].imshow(d, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto", interpolation="nearest")
ax[0, 2].add_patch(Rectangle((d.shape[1] / 2, 0), d.shape[1] / 2, d.shape[0] / 4, fill=False, ec=INK, lw=1, ls="--"))
cyc = np.mean(np.abs(np.rint(d[np.isfinite(d)])) >= 1) * 100
ax[0, 2].set_title(f"rebuild - lab: {cyc:.2f} % of pixels a cycle apart", fontsize=10, color=INK)
fig.colorbar(im, ax=ax[0, 2], **HB).set_label("rebuild - lab [cycles]", fontsize=8.5)
a_, b_, c_ = str(z["fig_trip"]).split("_")
for a, key, name in ((ax[1, 0], "fig_clo_ours", "the rebuild"), (ax[1, 1], "fig_clo_lab", "the lab's")):
    c = z[key] / (2 * np.pi)
    im = a.imshow(c, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto", interpolation="nearest")
    frac = np.mean(np.abs(c[np.isfinite(c)]) > 0.5) * 100
    a.set_title(f"{name}: closure of {ymd(a_)} / {ymd(b_)[5:]} / {ymd(c_)[5:]}\n{frac:.2f} % of pixels a cycle off",
                fontsize=9.5, color=INK)
fig.colorbar(im, ax=list(ax[1, :2]), **HB).set_label("closure [cycles]", fontsize=8.5)
for a in ax.flat[:5]:
    a.set_xlabel("range cell (west ->)", fontsize=8); a.set_ylabel("azimuth cell (south ->)", fontsize=8)
    a.tick_params(labelsize=7)
s = ax[1, 2]
top = 0.0
for k_o, k_l, color, label in (("cot_ours", "cot_lab", CO, "triplets across Ridgecrest"),
                               ("bt_ours", "bt_lab", BAT, f"{blabel} triplets")):
    to, tl = z[k_o], z[k_l]
    m = np.isfinite(to) & np.isfinite(tl)
    if not m.any():
        continue
    s.scatter(tl[m], to[m], s=30, color=color, edgecolor="white", linewidth=0.6, zorder=3,
              label=f"{label}: rebuild lower on {int(np.sum(to[m] < tl[m]))} / {int(m.sum())}")
    top = max(top, np.nanpercentile(np.r_[tl[m], to[m]], 97) * 1.25)
top = max(top, 0.5)
s.plot([0, top], [0, top], color="#9a9a9a", lw=1, ls="--", zorder=2)
for k_o, k_l, color in (("cot_ours", "cot_lab", CO), ("bt_ours", "bt_lab", BAT)):
    to, tl = z[k_o], z[k_l]
    for x, y in zip(tl[np.isfinite(to) & np.isfinite(tl) & (tl > top)], to[np.isfinite(to) & np.isfinite(tl) & (tl > top)]):
        s.annotate(f"lab {x:.1f} %", xy=(top * 0.98, min(y, top * 0.97)), fontsize=7.5, color=INK, ha="right", va="bottom")
        s.scatter([top * 0.99], [min(y, top * 0.97)], marker=">", s=30, color=color, zorder=3)
s.set_xlim(0, top); s.set_ylim(0, top)
s.set_xlabel("the lab's closure: pixels a cycle off [%]", fontsize=8.5, color=INK)
s.set_ylabel("the rebuild's [%]", fontsize=8.5, color=INK)
s.set_title("triplet closure, rebuild against lab", fontsize=10, color=INK)
s.legend(fontsize=7.5, loc="upper left", framealpha=0.95)
s.grid(color=GRID, lw=0.6); s.set_axisbelow(True); s.tick_params(labelsize=7.5)
fig.suptitle(f"T173 unwrapped pairs across the Ridgecrest earthquake (SNAPHU DEFO) and {blabel}, rebuild against the lab's\n"
             f"{len(z['co_cyc'])} coseismic pairs: a whole cycle apart at {np.median(z['co_cyc']):.2f} % of pixels "
             f"(median; NW quarter {np.nanmedian(z['co_cyc_nw']):.2f} %, max {np.nanmax(z['co_cyc_nw']):.2f} %); "
             f"{len(z['b_cyc'])} {blabel} pairs: {np.median(z['b_cyc']):.2f} % (max {np.max(z['b_cyc']):.2f} %)",
             fontsize=11.5, color=INK)
fig.savefig(out, dpi=110, bbox_inches="tight")
print(f"figure -> {out}")
