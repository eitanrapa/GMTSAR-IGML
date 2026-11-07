"""The Salton Sea's shoreline in the lab's merged product and in the rebuild (2026-11-02; user: "Could we do a shoreline
of Salton Sea figure?"): the geocoded coherence of one pair (shoreline_geocode.py's maps; the water decorrelates, the
desert stays coherent) under the Sentinel-2 water line of a day between the pair's dates (shoreline_clean.py: the lake
alone). Top left: the whole lake in the rebuild, with the zoom boxes; right: the west shore (Desert Shores) and the east
shore (north of Bombay Beach), the lab's product against the rebuild; bottom left: coherence across the long shores of
the desert mid-section (33.25-33.47 N: the ends are farmland and wetland, incoherent inland) against the distance from
the water line (trans.dat's DEM nodes, median per 25 m), with each curve's edge (half-way between water and land).
usage (insar_sim env): python shoreline_fig.py LAKE.npz COH.npz OUT.png [SUBTITLE]"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

INK, MUTED = "#1f1f1f", "#6b6b6b"
SHORE, LAB, OURS = "#2a78d6", "#eb6834", "#1baf7a"     # categorical slots 1-3 of the dataviz reference palette
GRAYS = LinearSegmentedColormap.from_list("grays", ["#f4f4f2", "#2b2b2b"])
ZOOMS = {"A  west-facing shore, Desert Shores": (-116.052, -116.028, 33.428, 33.448),
         "B  east shore, north of Bombay Beach": (-115.852, -115.828, 33.432, 33.452)}
# the long shores of the desert mid-section (33.25-33.47 N; the ends are farmland and wetland, incoherent inland)
FACES = {"west shore": (202.5, 292.5, "-"), "east shore": (22.5, 112.5, "--")}
MID = (33.25, 33.47)

LK, CO = np.load(sys.argv[1]), np.load(sys.argv[2])
out = sys.argv[3]
sub = sys.argv[4] if len(sys.argv) > 4 else ""
lon_m, lat_m, dm = CO["lon_map"], CO["lat_map"], CO["lon_map"][1] - CO["lon_map"][0]
sl, st, lake = LK["lon"], LK["lat"], LK["lake"].astype(np.float32)
dfs = sl[1] - sl[0]


def show(ax, z, box=None):
    w, e, s, n = box or (lon_m[0], lon_m[-1], lat_m[-1], lat_m[0])
    j = (lon_m >= w - dm) & (lon_m <= e + dm)
    i = (lat_m >= s - dm) & (lat_m <= n + dm)
    ax.imshow(z[np.ix_(i, j)], cmap=GRAYS, vmin=0.0, vmax=1.0, interpolation="nearest", aspect="auto",
              extent=[lon_m[j][0] - dm / 2, lon_m[j][-1] + dm / 2, lat_m[i][-1] - dm / 2, lat_m[i][0] + dm / 2])
    jj = (sl >= w - 0.003) & (sl <= e + 0.003)
    ii = (st >= s - 0.003) & (st <= n + 0.003)
    k = 1 if box else 4
    ax.contour(sl[jj][::k], st[ii][::k], lake[np.ix_(ii, jj)][::k, ::k], levels=[0.5], colors=SHORE,
               linewidths=2.2 if box else 1.0)
    ax.set_xlim(w, e); ax.set_ylim(s, n)
    ax.set_aspect(1 / np.cos(np.radians((s + n) / 2)))
    ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(4))
    fmt = FormatStrFormatter("%.3f" if box else "%.1f")
    ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
    ax.tick_params(labelsize=7, colors=INK)


def at_nodes(a):
    j = np.rint((CO["node_lon"] - sl[0]) / dfs).astype(np.int64)
    i = np.rint((st[0] - CO["node_lat"]) / dfs).astype(np.int64)
    ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1])
    v = np.full(j.shape, np.nan, np.float32)
    v[ok] = a[i[ok], j[ok]]
    return v


fig = plt.figure(figsize=(15.5, 10.5))
gs = fig.add_gridspec(2, 3, width_ratios=[1.3, 1, 1], wspace=0.22, hspace=0.28)
ax0 = fig.add_subplot(gs[0, 0])
show(ax0, CO["ours_map"])
for k, (name, (w, e, s, n)) in enumerate(ZOOMS.items()):
    ax0.add_patch(Rectangle((w, s), e - w, n - s, fill=False, ec=INK, lw=1.6))
    ax0.text(e + 0.008, n + 0.004, name[0], color=INK, fontsize=12, fontweight="bold", va="center")
ax0.set_title("the rebuild's coherence, 2016-10-05 x 10-29, T173 F2", fontsize=10, color=INK)
ax0.legend(handles=[Line2D([], [], color=SHORE, lw=2, label=f"water line, Sentinel-2 {LK['date']}")],
           loc="lower left", fontsize=8, framealpha=0.9)
for r, (name, box) in enumerate(ZOOMS.items()):
    for c, (key, label) in enumerate((("lab", "the lab's product"), ("ours", "the rebuild"))):
        ax = fig.add_subplot(gs[r, c + 1])
        show(ax, CO[f"{key}_map"], box)
        ax.set_title(f"{name}: {label}", fontsize=9.5, color=INK)
# coherence across the shore
ax = fig.add_subplot(gs[1, 0])
dist, bearing = at_nodes(LK["dist"]), at_nodes(LK["bearing"])
bins = np.arange(-800, 801, 25)
mid = (bins[:-1] + bins[1:]) / 2
for face, (lo, hi, ls) in FACES.items():
    sel = (bearing >= lo) & (bearing < hi) & np.isfinite(dist) & (CO["node_lat"] > MID[0]) & (CO["node_lat"] < MID[1])
    for key, color, label in (("lab", LAB, "lab"), ("ours", OURS, "rebuild")):
        coh = CO[f"{key}_coh"]
        k = np.digitize(dist[sel], bins) - 1
        med = np.array([np.nanmedian(coh[sel][k == b]) if np.sum(k == b) > 30 else np.nan for b in range(mid.size)])
        lo_lvl, hi_lvl = np.nanmedian(med[mid < -500]), np.nanmedian(med[mid > 500])
        half = (lo_lvl + hi_lvl) / 2
        up = np.where((med[:-1] < half) & (med[1:] >= half))[0]
        edge = np.nan
        if up.size:
            u = up[np.argmin(np.abs(mid[up]))]
            edge = mid[u] + (half - med[u]) / (med[u + 1] - med[u]) * (mid[u + 1] - mid[u])
        ax.plot(mid, med, ls=ls, color=color, lw=2, label=f"{label}, {face}: edge at {edge:+.0f} m")
ax.axvline(0, color=SHORE, lw=1.6)
ax.text(8, 0.08, "water line\n(Sentinel-2)", color=INK, fontsize=8)
ax.set_xlabel("distance from the water line [m]   (water < 0 < land)", fontsize=8.5, color=INK)
ax.set_ylabel("coherence (median at DEM nodes)", fontsize=8.5, color=INK)
ax.set_title("coherence across the shore (desert mid-section)", fontsize=10, color=INK)
ax.set_xlim(-800, 800); ax.set_ylim(0, 0.9)
ax.grid(color="#e6e6e6", lw=0.6); ax.set_axisbelow(True)
ax.tick_params(labelsize=7.5, colors=INK)
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.95)
cb = fig.colorbar(matplotlib.cm.ScalarMappable(cmap=GRAYS, norm=Normalize(0, 1)), ax=fig.axes[1:5],
                  fraction=0.02, pad=0.02)
cb.set_label("coherence (water decorrelates)", fontsize=9, color=INK); cb.ax.tick_params(labelsize=8)
fig.suptitle("The Salton Sea's shoreline in the geocoded T173 product: the lab's against the rebuild" +
             (f"\n{sub}" if sub else ""), fontsize=12, color=INK)
fig.savefig(out, dpi=120, bbox_inches="tight")
print(f"figure -> {out}")
