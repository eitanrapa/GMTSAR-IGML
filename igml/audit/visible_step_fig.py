"""One pair's interferogram from pass-1 (geometry-only) and pass-2 (network spectral diversity) images, side by side
(2026-11-02): the full subswath and zooms on three burst switches, where pass 1's misregistration shows as a jump in the
wrapped phase and pass 2's does not; under each zoom the phase averaged along range per line.
usage (Mac, insar_sim env): python visible_step_fig.py PASS1_phasefilt.grd PASS2_phasefilt.grd cuts OUT.png [LABEL [NAME1 NAME2 TITLE]]"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/Users/erapaportbruck/Documents/IGML")
from igml.real.grd import read_grd  # noqa: E402

p1, p2, cuts, out = sys.argv[1], sys.argv[2], [int(c) for c in sys.argv[3].split(",")], sys.argv[4]
label = sys.argv[5] if len(sys.argv) > 5 else ""
names = (sys.argv[6], sys.argv[7]) if len(sys.argv) > 7 else ("pass 1: geometry-only alignment", "pass 2: network spectral diversity")
title = sys.argv[8] if len(sys.argv) > 8 else "the same pair from pass-1 and pass-2 images (F1)"
G = [read_grd(p1), read_grd(p2)]
pick = [cuts[len(cuts) // 4], cuts[len(cuts) // 2], cuts[3 * len(cuts) // 4]]
fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(4, 4, height_ratios=[3, 1.1, 3, 1.1], width_ratios=[1.3, 1, 1, 1])
for i, (g, name) in enumerate(zip(G, names)):
    z = np.where(np.isfinite(g.z) & (g.z != 0), g.z, np.nan)
    a = fig.add_subplot(gs[2 * i:2 * i + 2, 0])
    a.imshow(z, cmap="hsv", vmin=-np.pi, vmax=np.pi, extent=[g.x[0], g.x[-1], g.y[-1], g.y[0]], aspect="auto",
             interpolation="nearest")
    for c in pick:
        a.axhline(c, color="k", lw=0.8, ls="--")
    a.set_title(name, fontsize=11)
    a.set_xlabel("range sample", fontsize=8); a.set_ylabel("azimuth line", fontsize=8); a.tick_params(labelsize=7)
    for k, c in enumerate(pick):
        r = (g.y > c - 300) & (g.y < c + 300)
        b = fig.add_subplot(gs[2 * i, k + 1])
        b.imshow(z[r], cmap="hsv", vmin=-np.pi, vmax=np.pi, extent=[g.x[0], g.x[-1], g.y[r][-1], g.y[r][0]],
                 aspect="auto", interpolation="nearest")
        b.axhline(c, color="k", lw=0.8, ls="--")
        b.set_title(f"burst switch at line {c}", fontsize=9); b.tick_params(labelsize=6)
        pr = fig.add_subplot(gs[2 * i + 1, k + 1])
        w = np.nansum(np.exp(1j * np.nan_to_num(z[r])) * np.isfinite(z[r]), axis=1)
        q = np.degrees(np.unwrap(np.angle(w)))
        pr.plot(g.y[r], q - np.median(q), lw=0.9, color="#2166ac")
        pr.axvline(c, color="k", lw=0.8, ls="--")
        pr.set_ylabel("phase [deg]", fontsize=7); pr.tick_params(labelsize=6)
fig.suptitle(f"T173 {label}: {title}", fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=110)
print(f"figure -> {out}")
