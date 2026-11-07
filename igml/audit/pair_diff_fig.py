"""Ours - the lab's for pairs of a replaced track (2026-11-06, pair_diff_data.py's npz): per pair, (a) unwrap.grd as stored on the same
nodes, (b) the filtered phase with the lab's content moved onto ours (the placement fix out), (c) (b) averaged along range per subswath over
one stretch of the track. Radar coordinates of the merged grid (descending: north up, west right).
usage (insar_sim env): python pair_diff_fig.py DATA.npz OUT.png [ROW0 ROW1]"""
import json
import sys
from datetime import date, timedelta

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

INK, MUTED, GRID, NODATA = "#1f1f1f", "#6b6a66", "#e6e6e6", "#c9c8c2"
SERIES = {"F1": "#2a78d6", "F2": "#eb6834", "F3": "#1baf7a"}          # the dataviz reference palette, categorical slots 1-3
DIV = LinearSegmentedColormap.from_list("div", ["#0d366b", "#2a78d6", "#9ec5f4", "#f0efec", "#f4b0ae", "#e34948", "#7a1c1c"], N=256)
DIV.set_bad(NODATA)
LIM_MM, LIM_DEG = 3.0, 10.0

z = np.load(sys.argv[1])
out = sys.argv[2]
R0, R1 = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (1800, 3000)
pairs = [str(p) for p in z["pairs"]]
K = int(z["K"])
OFFS = {"2016087_2016111": "68 / 174 / 120 m", "2022163_2022175": "69 / 174 / 115 m", "2024345_2024357": "71 / 173 / 106 m"}


def ymd(c):
    return date(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))


fig = plt.figure(figsize=(15, 15.5), layout="constrained")
gs = fig.add_gridspec(3, len(pairs), height_ratios=[6, 6, 1.7])
A = [fig.add_subplot(gs[0, i]) for i in range(len(pairs))]
B = [fig.add_subplot(gs[1, i]) for i in range(len(pairs))]
C = [fig.add_subplot(gs[2, i]) for i in range(len(pairs))]
lab_cm = ListedColormap([INK])
for i, p in enumerate(pairs):
    m = json.loads(str(z[f"{p}_meta"]))
    a, b = ymd(p[:7]), ymd(p[8:])
    u, ph = z[f"{p}_unw"], z[f"{p}_ph"]
    ext = (0, u.shape[1] * K, u.shape[0] * K, 0)
    im1 = A[i].imshow(u, cmap=DIV, vmin=-LIM_MM, vmax=LIM_MM, extent=ext, aspect="auto", interpolation="nearest")
    lo = np.ma.masked_where(~z[f"{p}_lab_only"], np.ones_like(u))
    A[i].imshow(lo, cmap=lab_cm, extent=ext, aspect="auto", interpolation="nearest", vmin=0, vmax=1)
    A[i].set_title(f"{a.isoformat()} × {b.strftime('%m-%d')}  ({p})\nas stored: sd {m['sd_mm']:.2f} mm; {m['cyc']:.3f} % of pixels a cycle apart",
                   fontsize=10, color=INK)
    im2 = B[i].imshow(ph, cmap=DIV, vmin=-LIM_DEG, vmax=LIM_DEG, extent=ext, aspect="auto", interpolation="nearest")
    v = ph[np.isfinite(ph)] * 55.465763 / 2 / 360
    q = np.percentile(v, [16, 84])
    B[i].set_title(f"lab moved onto ours (F1 / F2 / F3 {OFFS.get(p, '')}): sd {(q[1] - q[0]) / 2:.2f} mm", fontsize=10, color=INK)
    print(p, "as stored sd", round(m["sd_mm"], 2), "mm; placement out sd", round((q[1] - q[0]) / 2, 2), "mm", flush=True)
    s = [0] + m["seams"] + [u.shape[1] * K]
    for ax in (A[i], B[i]):
        for c in m["seams"]:
            ax.axvline(c, color=MUTED, lw=0.6, ls=(0, (4, 4)))
        ax.tick_params(labelsize=7.5, colors=MUTED, length=2)
        for sp in ax.spines.values():
            sp.set_color(GRID)
    for k, n in enumerate(("F1", "F2", "F3")):
        A[i].text((s[k] + s[k + 1]) / 2, 90, n, ha="center", va="top", fontsize=9, color=INK)
    B[i].axhspan(R0, R1, xmin=0, xmax=0.012, color=INK, lw=0)
    if i == 0:
        A[i].annotate("lab's zero-phase rows", xy=(330, 4757), xytext=(40, 4560), fontsize=8, color=INK,
                      arrowprops=dict(arrowstyle="-", color=INK, lw=0.7))
    B[i].set_xlabel("range cell (west →)", fontsize=8, color=MUTED)
    if i == 0:
        A[i].set_ylabel("(a) unwrapped phase\nazimuth cell (south →)", fontsize=9, color=INK)
        B[i].set_ylabel("(b) filtered phase, coherence > 0.3\nazimuth cell (south →)", fontsize=9, color=INK)
        C[i].set_ylabel("(c) along-range mean [deg]", fontsize=9, color=INK)
    y = np.arange(R0, R1)
    for n in ("F1", "F2", "F3"):
        C[i].plot(y, z[f"{p}_prof_{n}"][R0:R1], color=SERIES[n], lw=1.3, label=n)
    C[i].axhline(0, color=MUTED, lw=0.6)
    C[i].set_ylim(-12, 12); C[i].set_xlim(R0, R1)
    C[i].set_xlabel(f"azimuth cell (rows {R0}-{R1} of (b), about 100 km; the black tick there)", fontsize=8, color=MUTED)
    C[i].grid(color=GRID, lw=0.6); C[i].set_axisbelow(True)
    C[i].tick_params(labelsize=7.5, colors=MUTED, length=2)
    for sp in C[i].spines.values():
        sp.set_color(GRID)
    if i == 0:
        C[i].legend(fontsize=8, ncol=3, loc="upper right", frameon=True, framealpha=0.95, edgecolor=GRID)
cb = fig.colorbar(im1, ax=A, orientation="horizontal", fraction=0.03, pad=0.012, aspect=60, extend="both")
cb.set_label(f"ours − lab, LOS [mm] (median removed; whole cycles, ±27.7 mm, saturate). Black: only the lab unwrapped (outside mask_def). Gray: no data", fontsize=8.5, color=INK)
cb.ax.tick_params(labelsize=7.5, colors=MUTED)
cb = fig.colorbar(im2, ax=B, orientation="horizontal", fraction=0.03, pad=0.012, aspect=60, extend="both")
cb.set_label("ours − lab, phase [deg] (10 deg = 0.8 mm LOS). Gray: no data or coherence < 0.3", fontsize=8.5, color=INK)
cb.ax.tick_params(labelsize=7.5, colors=MUTED)
fig.suptitle("T173: our interferograms minus the lab's originals\n(a) as stored, the placement fix shows wherever fringes are steep, over the burst stripes; "
             "(b, c) with placement taken out, only the burst-step sawtooth is left", fontsize=12, color=INK)
fig.savefig(out, dpi=110, bbox_inches="tight")
print("figure ->", out)
