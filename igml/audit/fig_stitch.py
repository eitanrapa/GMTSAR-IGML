"""The stitch_tops figure (2026-10-29; user: "can you show me the stitch_tops problem?"), from fig_stitch_data.py's npz.
a: 2026-01-29 stitched by GMTSAR against the same scenes stitched with the fix (identical through frame 1, then the TOPS
signature of a one-line azimuth shift); b: the 12-day pair's coherence along track with each stitch; c: the line shift
of 2026-01-29 that maximises that coherence, per block; d, e: zooms on the start of frame 2; the code at the foot.

usage (insar_sim env): python fig_stitch.py DATA.npz OUT.png
"""
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402

INK, MUTED, EDGE = "#1f2328", "#59636e", "#c9ced4"
GM, FX = "#D55E00", "#0072B2"                                          # GMTSAR's stitch_tops, ours
KM_LINE = 0.01397                                                      # along-track ground spacing of a line [km]
KM_LOOK = 16 * 0.0042                                                  # 16 range samples on the ground in F1 [km]


def style(ax):
    ax.tick_params(colors=MUTED, labelsize=8)
    for s in ax.spines.values():
        s.set_color(EDGE)


def smooth(y, n=50):
    k = np.ones(n)
    ok = np.isfinite(y)
    return np.convolve(np.where(ok, y, 0), k, "same") / np.maximum(np.convolve(ok, k, "same"), 1)


if __name__ == "__main__":
    d = np.load(sys.argv[1])
    s, ly = d["self"], int(d["ly"])
    ph, cs = np.angle(s), np.abs(s)
    cg, cf = d["coh_g"], d["coh_f"]
    nr, nc = s.shape
    xk = lambda rows: rows * ly * KM_LINE
    L, Y = xk(nr), nc * KM_LOOK
    med = np.nanmedian(cs, axis=1)
    b0 = int(np.argmax((med < 0.95) & (np.arange(nr) > 50)))           # the first look-row the two stitches differ
    xb = xk(b0)
    fig = plt.figure(figsize=(13.5, 12.4))
    gs = fig.add_gridspec(4, 3, height_ratios=[2.2, 1.25, 0.85, 2.5], width_ratios=[1.3, 1, 1], hspace=0.55,
                          wspace=0.12)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :], sharex=ax_a)
    ax_c = fig.add_subplot(gs[2, :], sharex=ax_a)
    ax_d, ax_e1, ax_e2 = fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1]), fig.add_subplot(gs[3, 2])
    cyc = plt.get_cmap("twilight_shifted").copy()                     # zero phase light: "identical" reads as empty
    cyc.set_bad("#e4e7eb")
    grey = plt.get_cmap("gray").copy()
    grey.set_bad("white")

    # a: the same scenes stitched two ways
    ax_a.imshow(ph.T, extent=(0, L, Y, 0), cmap=cyc, vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest")
    ax_a.set_title("a   2026-01-29 stitched by GMTSAR's stitch_tops, against the same two frames stitched with the fix "
                   "(phase of their product)", fontsize=9.5, color=INK, loc="left")
    box = dict(boxstyle="round,pad=0.3", fc="white", ec=EDGE, lw=0.6)
    ax_a.text(xb / 2, Y / 2, "frame 1: identical (zero phase)", ha="center", va="center", fontsize=9, color=INK, bbox=box)
    ax_a.text((xb + L) / 2, Y / 2, "frame 2: one line apart, a fringe every 1.9 km", ha="center", va="center",
              fontsize=9, color=INK, bbox=box)
    ax_a.set_ylabel("range [km]", fontsize=8.5, color=MUTED)

    # b: the 12-day pair's coherence along track
    xp = (np.arange(d["prof_g"].size) + 0.5) * 4 * KM_LINE
    for y, c, ls, lab in ((d["prof_f"], FX, "-", "with the fix"), (d["prof_g"], GM, "--", "GMTSAR's stitch_tops")):
        ax_b.plot(xp, smooth(y), color=c, lw=1.6, ls=ls, label=lab)
    for lo, hi, where in ((0, xb, "frame 1"), (xb, L, "frame 2")):
        m = (xp >= lo) & (xp < hi)
        g, f = np.nanmedian(d["prof_g"][m]), np.nanmedian(d["prof_f"][m])
        ax_b.text((lo + hi) / 2, 0.08, f"{where}: {g:.2f} with GMTSAR's, {f:.2f} with the fix", ha="center",
                  fontsize=8.5, color=INK)
    ax_b.set_ylim(0, 0.85)
    ax_b.set_ylabel("coherence", fontsize=8.5, color=MUTED)
    ax_b.legend(loc="center left", bbox_to_anchor=(0.0, 0.42), frameon=False, fontsize=8.5, labelcolor=INK, ncol=2)
    ax_b.set_title("b   Coherence of the 12-day pair 2026-01-17 / 2026-01-29, median across the subswath (4 x 4-pixel "
                   "windows, smoothed over 2.8 km)", fontsize=9.5, color=INK, loc="left")

    # c: the shift that maximises it
    sh, mid = d["shifts"], d["shift_mid"] * KM_LINE
    for key, c, mk, dy, lab in (("shift_f", FX, "s", -0.12, "with the fix"), ("shift_g", GM, "o", 0.12,
                                                                               "GMTSAR's stitch_tops")):
        best = sh[np.nanargmax(d[key], axis=1)]
        ax_c.scatter(mid, best + dy, s=34, marker=mk, facecolors="none" if mk == "o" else c, edgecolors=c,
                     linewidths=1.4, label=lab, zorder=3)
    ax_c.set_ylim(-2.5, 2.5)
    ax_c.set_yticks([-2, -1, 0, 1, 2])
    ax_c.grid(axis="y", color=EDGE, lw=0.6)
    ax_c.set_ylabel("lines", fontsize=8.5, color=MUTED)
    ax_c.set_xlabel("along track [km from the north end of the subswath]", fontsize=8.5, color=MUTED)
    ax_c.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=INK, ncol=2)
    ax_c.set_title("c   Line shift of 2026-01-29 that maximises its coherence with 2026-01-17, per block of 1,179 lines",
                   fontsize=9.5, color=INK, loc="left")
    for ax in (ax_a, ax_b, ax_c):
        ax.axvline(xb, color=INK, lw=0.9, ls=":")
        style(ax)
    ax_a.set_xlim(0, L)

    # d, e: zooms on the start of frame 2
    r0, r1 = max(b0 - 340, 0), min(b0 + 340, nr)
    ext = (xk(r0), xk(r1), Y, 0)
    ax_d.imshow(ph[r0:r1].T, extent=ext, cmap=cyc, vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest")
    ax_d.set_title("d   a, zoomed: ~10 fringes per 19-km burst, the burst's\nDoppler sweep times a one-line shift",
                   fontsize=9.5, color=INK, loc="left")
    ax_d.set_ylabel("range [km]", fontsize=8.5, color=MUTED)
    for ax, z, t in ((ax_e1, cg, "e   12-day coherence, GMTSAR's stitch_tops"), (ax_e2, cf, "    ... and with the fix")):
        ax.imshow(z[r0:r1].T, extent=ext, cmap=grey, vmin=0, vmax=0.8, aspect="auto", interpolation="nearest")
        ax.set_title(t, fontsize=9.5, color=INK, loc="left")
        ax.set_yticklabels([])
    for ax in (ax_d, ax_e1, ax_e2):
        ax.axvline(xb, color=GM if ax is not ax_e2 else FX, lw=1.0, ls=":")
        ax.set_xlabel("along track [km]", fontsize=8.5, color=MUTED)
        style(ax)

    code = ("stitch_tops.c, the last frame (GMTSAR's master, and the lab's build):\n"
            "    n_start = n_start - line_offset_total;   /* the overlap shifts so far, this one included: right   */\n"
            "    n_start = n_start - line_offset;         /* this frame's shift AGAIN; its printf is commented out */\n"
            "2026-01-29: ashift = -1 in the overlap, so frame 2 was read from its line 729 instead of 728. The\n"
            "alignment that follows is geometric (orbits + DEM), so nothing downstream measures or undoes it.")
    fig.text(0.125, ax_d.get_position().y0 - 0.05, code, family="monospace", fontsize=8.6, color=INK, va="top")
    fig.text(0.125, 0.955, "GMTSAR's stitch_tops shifts the last frame twice: T173 2026-01-29, subswath F1",
             fontsize=13, color=INK, weight="bold")
    fig.text(0.125, 0.94, "Both versions aligned to the super-master (2016-10-06) with the same recipe; only the stitch "
             "differs. One line is ~14 m along track, about 2/3 of the azimuth resolution.", fontsize=8.8,
             color=MUTED, va="top")
    fig.savefig(sys.argv[2], dpi=160, bbox_inches="tight", facecolor="white")
    print(f"frame 2 starts at look-row {b0} (line ~{b0 * ly}, {xb:.1f} km)")
    for key in ("shift_g", "shift_f"):
        print(key, "best:", sh[np.nanargmax(d[key], axis=1)].tolist())
        print("  coherence at 0 / best:", " ".join(f"{d[key][i, 2]:.2f}/{np.nanmax(d[key][i]):.2f}"
                                                 for i in range(d[key].shape[0])))
    print("wrote", sys.argv[2])
