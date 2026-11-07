"""Our interferogram of a pair against the lab's, all three subswaths (2026-11-01; user: "let's check on those first
interferograms ... and let me see them"): rows = our filtered phase, the lab's, their difference (ours - lab, degrees:
the misregistration the network spectral diversity removed shows as a sawtooth along azimuth that jumps at every
burst switch, dashed), and that difference averaged along range per line. Radar coordinates (x = range sample,
y = SLC line).
usage (Mac, insar_sim env): python intf_compare_fig.py DIR PAIR OUT.png cutsF1 cutsF2 cutsF3 [LABEL]"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/Users/erapaportbruck/Documents/IGML")
from igml.real.grd import read_grd  # noqa: E402

d, pair, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
cuts = [[int(c) for c in a.split(",")] for a in sys.argv[4:7]]
label = sys.argv[7] if len(sys.argv) > 7 else pair
fig, ax = plt.subplots(4, 3, figsize=(13, 13), gridspec_kw={"height_ratios": [3, 3, 3, 1.6]})
for j, n in enumerate((1, 2, 3)):
    po, pl = read_grd(d / f"ours_F{n}_phasefilt.grd"), read_grd(d / f"lab_F{n}_phasefilt.grd")
    ok = np.isfinite(po.z) & np.isfinite(pl.z) & (po.z != 0) & (pl.z != 0)
    ext = [po.x[0], po.x[-1], po.y[-1], po.y[0]]
    for i, (z, name) in enumerate(((po.z, "ours"), (pl.z, "lab"))):
        ax[i, j].imshow(np.where(ok, z, np.nan), cmap="hsv", vmin=-np.pi, vmax=np.pi, extent=ext, aspect="auto",
                        interpolation="nearest")
        ax[i, j].set_title(f"{name}: F{n} phase", fontsize=10)
    dif = np.degrees(np.angle(np.exp(1j * (po.z - pl.z))))
    dif = np.where(ok, dif, np.nan)
    v = max(5.0, float(np.nanpercentile(np.abs(dif - np.nanmedian(dif)), 95)))
    im = ax[2, j].imshow(dif, cmap="RdBu_r", vmin=-v, vmax=v, extent=ext, aspect="auto", interpolation="nearest")
    ax[2, j].set_title(f"ours - lab, F{n} [deg]", fontsize=10)
    fig.colorbar(im, ax=ax[2, j], fraction=0.04, pad=0.02)
    for c in cuts[j]:
        ax[2, j].axhline(c, color="k", lw=0.35, ls="--")
    w = np.where(ok, np.exp(1j * np.radians(np.nan_to_num(dif))), 0).sum(axis=1)
    prof = np.degrees(np.angle(w))
    ax[3, j].plot(po.y, prof, lw=0.8, color="#2166ac")
    for c in cuts[j]:
        ax[3, j].axvline(c, color="0.5", lw=0.4, ls="--")
    ax[3, j].set_xlim(po.y[0], po.y[-1])
    ax[3, j].set_ylim(np.nanpercentile(prof, 1) - 2, np.nanpercentile(prof, 99) + 2)
    ax[3, j].set_xlabel("azimuth line", fontsize=8)
    ax[3, j].set_ylabel("ours - lab [deg]", fontsize=8)
    for a in ax[:, j]:
        a.tick_params(labelsize=7)
fig.suptitle(f"T173 {label}: our interferogram against the lab's (dashed: burst switches)", fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=110)
print(f"figure -> {out}")
