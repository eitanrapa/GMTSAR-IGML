"""The first interferograms of a track's rebuild against the lab's of the same pair (2026-11-01): per subswath the
filtered phase and the coherence in radar coordinates (ours, the lab's), and a crop around one of the super-master's
burst-switch lines (where A1 decorrelated a band and A5 stepped the phase).

usage (Mac, insar_sim env): python first_intf_fig.py DIR PAIR OUT.png [cut_F1,cut_F2,cut_F3]
  DIR holds ours_F<n>_phasefilt.grd / ours_F<n>_corr.grd and, where the lab has the pair, lab_F<n>_*.grd
  cut_F<n>: a burst-switch line (SLC lines) per subswath for the crop
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/Users/erapaportbruck/Documents/IGML")
from igml.real.grd import read_grd  # noqa: E402

d, pair, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
cuts = [int(c) for c in sys.argv[4].split(",")] if len(sys.argv) > 4 else [None, None, None]
subs = [n for n in (1, 2, 3) if (d / f"ours_F{n}_phasefilt.grd").exists()]
has_lab = all((d / f"lab_F{n}_phasefilt.grd").exists() for n in subs)
srcs = ["ours", "lab"] if has_lab else ["ours"]
ncol = 2 * len(srcs)
fig, ax = plt.subplots(len(subs) * 2, ncol, figsize=(3.2 * ncol, 4.6 * len(subs)),
                       gridspec_kw={"height_ratios": [3, 1] * len(subs)})
ax = np.atleast_2d(ax)
for i, n in enumerate(subs):
    for j, s in enumerate(srcs):
        ph, co = read_grd(d / f"{s}_F{n}_phasefilt.grd"), read_grd(d / f"{s}_F{n}_corr.grd")
        y, z = ph.y, np.where(np.isfinite(ph.z) & (ph.z != 0), ph.z, np.nan)
        ext = [ph.x[0], ph.x[-1], ph.y[-1], ph.y[0]]
        a = ax[2 * i, 2 * j]
        a.imshow(z, cmap="hsv", vmin=-np.pi, vmax=np.pi, extent=ext, aspect="auto", interpolation="nearest")
        a.set_title(f"{s}: F{n} phase", fontsize=9)
        b = ax[2 * i, 2 * j + 1]
        b.imshow(co.z, cmap="gray", vmin=0, vmax=1, extent=[co.x[0], co.x[-1], co.y[-1], co.y[0]], aspect="auto",
                 interpolation="nearest")
        b.set_title(f"{s}: F{n} coherence (median {np.nanmedian(co.z):.2f})", fontsize=9)
        c = cuts[n - 1]
        for a_ in (a, b):
            a_.tick_params(labelsize=6)
            if c is not None:
                a_.axhline(c, color="w", lw=0.6, ls="--")
        # crop: +-200 SLC lines around the cut, the whole range
        ca, cb = ax[2 * i + 1, 2 * j], ax[2 * i + 1, 2 * j + 1]
        if c is not None:
            r = (y > c - 200) & (y < c + 200)
            cz, cc = z[r], co.z[(co.y > c - 200) & (co.y < c + 200)]
            ca.imshow(cz, cmap="hsv", vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest",
                      extent=[ph.x[0], ph.x[-1], y[r][-1], y[r][0]])
            cb.imshow(cc, cmap="gray", vmin=0, vmax=1, aspect="auto", interpolation="nearest",
                      extent=[co.x[0], co.x[-1], co.y[(co.y > c - 200) & (co.y < c + 200)][-1],
                              co.y[(co.y > c - 200) & (co.y < c + 200)][0]])
            for a_ in (ca, cb):
                a_.axhline(c, color="k", lw=0.6, ls="--")
                a_.tick_params(labelsize=6)
            ca.set_title(f"burst switch at line {c}: phase", fontsize=8)
            cb.set_title("coherence", fontsize=8)
        else:
            ca.axis("off"); cb.axis("off")
fig.suptitle(f"T173 {pair}: the rebuild's interferogram" + (" against the lab's" if has_lab else ""), fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=110)
print(f"figure -> {out}")
