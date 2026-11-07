"""The rebuild's first unwrapped pairs against the lab's (2026-11-02; unw_check.py's npz): top, one pair's unwrapped phase
as LOS change (ours, the lab's, ours - lab); bottom, the closure of a triplet through that pair (ours, the lab's: a
whole cycle wrong somewhere in the triangle shows as +-1) and, for every sampled triplet, the share of pixels a cycle
off, ours against the lab's. Radar coordinates of the merged grid (every 4th cell).
usage (Mac, insar_sim env): python unw_check_fig.py CHECK.npz OUT.png"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, GRID = "#1f1f1f", "#e6e6e6"
OURS, LAB = "#1baf7a", "#eb6834"
MM = 55.465763 / (4 * np.pi)                    # mm of LOS change per radian (C-band)

z = np.load(sys.argv[1])
out = sys.argv[2]
pair = str(z["figpair"])
def ymd(d):                                     # GMTSAR's yyyyddd (day 0 = 1 January) -> YYYY-MM-DD
    from datetime import date, timedelta
    return (date(int(d[:4]), 1, 1) + timedelta(days=int(d[4:]))).isoformat()


pa, pb = (ymd(x) for x in pair.split("_"))
PAIR = f"{pa} x {pb[5:] if pb[:4] == pa[:4] else pb}"


def los(u):
    v = -u * MM
    return v - np.nanmedian(v)


fig, ax = plt.subplots(2, 3, figsize=(15, 12.5), layout="constrained")
HB = dict(orientation="horizontal", fraction=0.045, pad=0.02, aspect=40)
o, l = los(z["ours"]), los(z["lab"])
lim = np.nanpercentile(np.abs(np.r_[o[np.isfinite(o)], l[np.isfinite(l)]]), 99)
for a, v, name in ((ax[0, 0], o, "the rebuild"), (ax[0, 1], l, "the lab's")):
    im = a.imshow(v, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
    a.set_title(f"{name}: unwrapped {PAIR}", fontsize=10, color=INK)
fig.colorbar(im, ax=list(ax[0, :2]), **HB).set_label("LOS change [mm] (median removed)", fontsize=8.5)
d = (z["ours"] - z["lab"])
d = d - np.nanmedian(d)
im = ax[0, 2].imshow(d / (2 * np.pi), cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto", interpolation="nearest")
cyc = np.mean(np.abs(np.rint(d[np.isfinite(d)] / (2 * np.pi))) >= 1) * 100
ax[0, 2].set_title(f"rebuild - lab: {cyc:.2f} % of pixels a cycle apart", fontsize=10, color=INK)
fig.colorbar(im, ax=ax[0, 2], **HB).set_label("rebuild - lab [cycles]", fontsize=8.5)
if "clo_ours" in z.files:
    t = " + ".join(str(x) for x in z["trip"][:2])
    a_, b_, c_ = (str(x) for x in z["trip"])
    for a, key, name in ((ax[1, 0], "clo_ours", "the rebuild"), (ax[1, 1], "clo_lab", "the lab's")):
        if key not in z.files:
            a.axis("off"); continue
        c = z[key] / (2 * np.pi)
        im = a.imshow(c, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto", interpolation="nearest")
        frac = np.mean(np.abs(c[np.isfinite(c)]) > 0.5) * 100
        A, B, C = ymd(a_), ymd(b_)[5:], ymd(c_)[5:]
        a.set_title(f"{name}: closure of {A} / {B} / {C}\n{frac:.2f} % of pixels a cycle off", fontsize=9.5, color=INK)
    fig.colorbar(im, ax=list(ax[1, :2]), **HB).set_label("closure [cycles]", fontsize=8.5)
for a in ax.flat[:5]:
    a.set_xlabel("range cell", fontsize=8); a.set_ylabel("azimuth cell", fontsize=8); a.tick_params(labelsize=7)
s = ax[1, 2]
to, tl = z["trip_ours"], z["trip_lab"]
m = np.isfinite(to) & np.isfinite(tl)
s.scatter(tl[m], to[m], s=28, color=OURS, edgecolor="white", linewidth=0.6, zorder=3)
top = max(1.0, np.nanpercentile(tl[m], 95) * 1.3)
s.plot([0, top], [0, top], color="#9a9a9a", lw=1, ls="--", zorder=2)
off = m & (tl > top)
for x, y in zip(tl[off], to[off]):
    s.annotate(f"lab {x:.1f} %", xy=(top * 0.98, y), fontsize=7.5, color=INK, ha="right", va="bottom")
    s.scatter([top * 0.99], [y], marker=">", s=30, color=OURS, zorder=3)
s.set_xlim(0, top); s.set_ylim(0, top)
s.set_xlabel("the lab's closure: pixels a cycle off [%]", fontsize=8.5, color=INK)
s.set_ylabel("the rebuild's [%]", fontsize=8.5, color=INK)
s.set_title(f"{int(m.sum())} triplets: the rebuild lower on {int(np.sum(to[m] < tl[m]))}\n(median {np.median(to[m]):.2f} "
            f"vs {np.median(tl[m]):.2f} %)", fontsize=10, color=INK)
s.grid(color=GRID, lw=0.6); s.set_axisbelow(True); s.tick_params(labelsize=7.5)
fig.suptitle(f"T173's first unwrapped interferograms: the rebuild against the lab's (sampled from the 626 stored, "
             f"2014-2019)\n{len(z['cyc'])} pairs against the lab's same pair: a whole cycle apart at "
             f"{np.median(z['cyc']):.2f} % of pixels (median; max {np.max(z['cyc']):.2f} %)", fontsize=11.5, color=INK)
fig.savefig(out, dpi=110, bbox_inches="tight")
print(f"figure -> {out}")
