"""The Rosamond corner-reflector geolocation figure (2026-10-29; user: "Can I see the geolocation problem with the
reflectors?"), from fig_cr_geoloc_data.py's npz. Rows: the SLC at full resolution (the truth), the lab's merged
product, our conversion (M1 fixed), the proposed correction (preview). Left: the descending-facing reflectors 13-27;
right: a zoom on 4.8-m reflector 25. Rings = surveyed positions; arrows = surveyed -> imaged (the brightness centroid
of each 4.8-m reflector); the offset printed is the median over reflectors 23-27, measured on that panel's image.

usage (insar_sim env): python fig_cr_geoloc.py DATA.npz CR.csv OUT.png
"""
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402

INK, MUTED, RING, ARROW = "#1f2328", "#59636e", "#ffb000", "#ff4fd8"


def local(lon, lat, lon0, lat0):
    return (lon - lon0) * 111320.0 * np.cos(np.radians(lat0)), (lat - lat0) * 110950.0


if __name__ == "__main__":
    d = np.load(sys.argv[1])
    crs = {}
    for line in open(sys.argv[2]).read().splitlines()[1:]:
        f = line.split(",")
        if 13 <= int(f[0]) <= 27:
            crs[int(f[0])] = (float(f[2]), float(f[1]))
    lon0, lat0 = crs[25]
    gx, gy = local(d["lon"], d["lat"], lon0, lat0)
    X, Y = np.meshgrid(gx, gy)
    cx = {k: local(v[0], v[1], lon0, lat0) for k, v in crs.items()}
    rows = [("slc", "1  SLC at full resolution (the truth)"),
            ("lab", "2  Lab product (GMTSAR's merged grid)"),
            ("ours", "3  Our conversion today (F2 column fix)"),
            ("fix", "4  Proposed correction (preview only)")]
    fig, axes = plt.subplots(4, 2, figsize=(13.5, 10.2), gridspec_kw={"width_ratios": [4.6, 0.95], "hspace": 0.30,
                                                                         "wspace": 0.05})
    for (key, title), (axw, axz) in zip(rows, axes):
        img = d[key].astype(float)
        db = 20 * np.log10(np.maximum(img, 1e-12) / np.nanmedian(img))
        lo, hi = 0, np.nanpercentile(db, 99.97)                         # the brightest cell white, its neighbours grey
        # measured offset: brightness centroid (pixels within 6 dB of the peak, 250 m around) of each 4.8-m reflector
        vec = {}
        for k in (23, 24, 25, 26, 27):
            ex, ny = cx[k]
            m = np.hypot(X - ex, Y - ny) < 250
            pk = np.nanmax(db[m])
            h = m & (db >= pk - 6)
            w = 10 ** (db[h] / 10)
            vec[k] = (np.sum(w * X[h]) / w.sum() - ex, np.sum(w * Y[h]) / w.sum() - ny)
        v = np.median(np.array(list(vec.values())), axis=0)
        for ax, (x0, x1, y0, y1) in ((axw, (-1650, 3250, -330, 450)), (axz, (-230, 130, -180, 180))):
            ax.imshow(db, extent=(gx[0], gx[-1], gy[-1], gy[0]), cmap="gray", vmin=lo, vmax=hi, interpolation="nearest",
                      origin="upper")
            for k, (ex, ny) in cx.items():
                big = k >= 23
                ax.scatter(ex, ny, s=110 if big else 60, facecolors="none", edgecolors=RING, linewidths=1.4, zorder=3)
                if k in vec and (ax is axz or big):
                    dx, dy = vec[k]
                    if np.hypot(dx, dy) > 8:
                        ax.annotate("", xy=(ex + dx, ny + dy), xytext=(ex, ny),
                                    arrowprops=dict(arrowstyle="-|>", color=ARROW, lw=1.6, mutation_scale=11), zorder=4)
            ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.tick_params(colors=MUTED, labelsize=8)
            for s in ax.spines.values():
                s.set_color("#c9ced4")
        axz.set_title(f"reflector 25: {np.hypot(*vec[25]):.0f} m", fontsize=9, color=MUTED, loc="left")
        axw.set_title(f"{title}:  reflectors imaged {np.hypot(*v):.0f} m from their surveyed positions "
                      f"({abs(v[0]):.0f} m {'west' if v[0] < 0 else 'east'}, {abs(v[1]):.0f} m "
                      f"{'south' if v[1] < 0 else 'north'})", fontsize=10.5, color=INK, loc="left")
        axw.set_ylabel("north [m]", fontsize=8.5, color=MUTED)
        print(f"{key}: median offset of reflectors 23-27 = {v[0]:+.1f} m E, {v[1]:+.1f} m N ({np.hypot(*v):.1f} m); "
              + ", ".join(f"{k}: {np.hypot(*vv):.0f}" for k, vv in vec.items()))
    axes[-1][0].set_xlabel("east of reflector 25 [m]", fontsize=8.5, color=MUTED)
    axes[-1][1].set_xlabel("east [m]", fontsize=8.5, color=MUTED)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=RING, markersize=9, mew=1.4,
                      label="surveyed reflector position (JPL Rosamond array; large = 4.8 m, small = 2.4 m)"),
               Line2D([], [], color=ARROW, lw=1.6, marker=">", markersize=6, label="surveyed → imaged (4.8-m reflectors)")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.12, 0.035), ncol=2, frameon=False, fontsize=9,
               labelcolor=INK)
    fig.subplots_adjust(top=0.88)
    fig.text(0.125, 0.965, "Rosamond corner reflectors in T071 (descending), subswath F2", fontsize=13, color=INK,
             weight="bold")
    fig.text(0.125, 0.948,"Amplitude in dB above the lake bed. Every pixel takes its radar position from GMTSAR's "
             "trans.dat, the lookup GMTSAR and our conversion both geocode with.\nRows 2-4: pair 2018-10-13 / "
             "2019-01-05, decimated cells of 32 samples × 6 lines (~120 × 85 m on the ground); row 1: the "
             "super-master SLC, 2019-03-18.", fontsize=8.5, color=MUTED, va="top")
    fig.savefig(sys.argv[3], dpi=160, bbox_inches="tight", facecolor="white")
    print("wrote", sys.argv[3])
