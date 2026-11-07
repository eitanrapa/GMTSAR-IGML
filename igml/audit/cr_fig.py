"""Piñon Flat corner reflectors in a geocoded product, the lab's against the rebuild (2026-11-02; user: "let's see if we
can plot the geolocation issue resolved" -> "Maybe the Pinion Flats reflector?"): per product the geocoded coherence
around the reflectors (cr_geocode_any.py's npz, 0.0002-deg cells), the surveyed apex (UCSD's survey) and the centroid of
the reflector's coherence blob, with the offset in metres; the product's cells drawn as tiles (nearest cell). The SLC's own geometry (cr_peak.py) goes in the subtitle.
usage (Mac, insar_sim env): python cr_fig.py LAB.npz OURS.npz OUT.png [SUBTITLE]"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import distance_transform_edt
from matplotlib.lines import Line2D

INK = "#1f1f1f"
SURVEY, MEASURED = "#2a78d6", "#eb6834"          # categorical slots 1 / 2 of the dataviz reference palette
GRAYS = LinearSegmentedColormap.from_list("grays", ["#f4f4f2", "#2b2b2b"])
HALO = [pe.withStroke(linewidth=3, foreground="white")]

files, out = sys.argv[1:3], sys.argv[3]
sub = sys.argv[4] if len(sys.argv) > 4 else ""
fig, axes = plt.subplots(1, 2, figsize=(11.5, 6.2), sharex=True, sharey=True)
for ax, f, name in zip(axes, files, ("the lab's product", "the rebuild")):
    z = np.load(f)
    g, E, N, crs, ids, res = z["g"], z["E"], z["N"], z["crs"], z["ids"], z["res"]
    lat0, lon0 = crs[:, 0].mean(), crs[:, 1].mean()
    kx, ky = 111320.0 * np.cos(np.radians(lat0)), 110950.0
    x, y = (E - lon0) * kx, (N - lat0) * ky                    # metres east / north of the reflectors' mean position
    dx, dy = abs(x[1] - x[0]) / 2, abs(y[1] - y[0]) / 2
    # the product's cells (~75 x 85 m) are coarser than the 0.0002-deg target grid, which holds a value only at each
    # cell's nearest node: every node takes its nearest cell's value (within 60 m), so the map shows the cells as tiles
    dist, (ii, jj) = distance_transform_edt(~np.isfinite(g), sampling=(2 * dy, 2 * dx), return_indices=True)
    g = np.where(dist < 60, g[ii, jj], np.nan)
    im = ax.imshow(g, cmap=GRAYS, vmin=0.3, vmax=1.0, interpolation="nearest",
                   extent=[x[0] - dx, x[-1] + dx, y[-1] - dy, y[0] + dy])
    off = []
    for i, cid in enumerate(ids[:len(res)]):
        sx, sy = (crs[i, 1] - lon0) * kx, (crs[i, 0] - lat0) * ky
        mx, my = sx + res[i, 0], sy + res[i, 1]
        d = float(np.hypot(res[i, 0], res[i, 1]))
        off.append(f"{cid} {d:.0f} m")
        ax.plot([sx, mx], [sy, my], color=INK, lw=1.0)
        ax.plot(sx, sy, "+", color=SURVEY, ms=14, mew=2.5)
        ax.plot(mx, my, "x", color=MEASURED, ms=10, mew=2.5)
        ax.text(sx + 18, sy + 22, f"{cid}: {d:.0f} m", color=INK, fontsize=10, path_effects=HALO)
    ax.set_title(f"{name}: {', '.join(off)} from the survey", fontsize=11, color=INK)
    ax.set_xlim(-480, 480); ax.set_ylim(-480, 480)
    ax.set_xlabel("east of the reflectors [m]", fontsize=9, color=INK)
    ax.tick_params(labelsize=8, colors=INK)
axes[0].set_ylabel("north [m]", fontsize=9, color=INK)
cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
cb.set_label("coherence", fontsize=9, color=INK)
cb.ax.tick_params(labelsize=8)
axes[1].legend(handles=[Line2D([], [], ls="", marker="+", color=SURVEY, ms=12, mew=2.5, label="surveyed apex (UCSD)"),
                        Line2D([], [], ls="", marker="x", color=MEASURED, ms=9, mew=2.5,
                               label="reflector in the product")],
               loc="lower right", fontsize=9, framealpha=0.9)
fig.suptitle("T173 Piñon Flat corner reflectors D1 / D2 in the geocoded super-master pair (2016-10-05 x 10-29)" +
             (f"\n{sub}" if sub else ""), fontsize=11.5, color=INK)
fig.savefig(out, dpi=120, bbox_inches="tight")
print(f"figure -> {out}")
