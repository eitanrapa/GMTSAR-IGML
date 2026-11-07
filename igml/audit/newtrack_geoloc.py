"""Absolute geolocation check of a NEW track's merged grid (10-29): water in GSHHG's full shoreline (landmask_ra.grd,
GMTSAR's landmask.csh: the coastline projected to radar coordinates through the frame's trans.dat) against water seen
in the data (low coherence of a pair, corr.grd), inside a lon / lat box (a lake, a coast), for shifts of the GSHHG mask
by -SHIFT..SHIFT cells in azimuth and range: the best agreement should be at 0 / 0 (one cell = 32 samples x 6 lines,
~120 x 85 m on the ground). Agreement = correlation of the GSHHG water mask with 1 - coherence (sub-cell peak by a
quadratic), and the Jaccard index of the two water masks (coherence < COH_WATER). Read only.

usage (igml env): python newtrack_geoloc.py corr.grd landmask_ra.grd trans.dat LON0 LON1 LAT0 LAT1
"""
import sys

import numpy as np

from igml.real.grd import read_grd

SHIFT, COH_WATER = 4, 0.2

c = read_grd(sys.argv[1])
m = read_grd(sys.argv[2])
lon0, lon1, lat0, lat1 = map(float, sys.argv[4:8])
t = np.memmap(sys.argv[3], dtype="<f8", mode="r").reshape(-1, 5)
inbox = (t[:, 3] >= lon0) & (t[:, 3] <= lon1) & (t[:, 4] >= lat0) & (t[:, 4] <= lat1)
r, a = t[inbox, 0], t[inbox, 1]
x0, x1, y0, y1 = r.min(), r.max(), a.min(), a.max()
print(f"box {lon0}..{lon1} E, {lat0}..{lat1} N -> range {x0:.0f}..{x1:.0f}, azimuth {y0:.0f}..{y1:.0f} "
      f"({inbox.sum()} DEM nodes)")


def on(g, xs, ys):
    """Grid g's values at nodes (xs, ys) of another grid (nearest node)."""
    jx = np.clip(np.round((xs - g.x[0]) / (g.x[1] - g.x[0])).astype(int), 0, g.x.size - 1)
    iy = np.clip(np.round((ys - g.y[0]) / (g.y[1] - g.y[0])).astype(int), 0, g.y.size - 1)
    return g.z[np.ix_(iy, jx)].astype(float)


lm = on(m, c.x, c.y)
water_m = ~np.isfinite(lm) | (lm < 0.5)
cz = c.z.astype(float)
fin = np.isfinite(cz)
water_c = fin & (cz < COH_WATER)
win = (c.x[None, :] >= x0) & (c.x[None, :] <= x1) & (c.y[:, None] >= y0) & (c.y[:, None] <= y1) & fin
print(f"window: {win.sum()} cells with coherence; GSHHG water {water_m[win].mean():.3f}, coherence < {COH_WATER} "
      f"{water_c[win].mean():.3f}, median coherence on GSHHG water {np.median(cz[win & water_m]):.2f} / land "
      f"{np.median(cz[win & ~water_m]):.2f}")
C = np.full((2 * SHIFT + 1, 2 * SHIFT + 1), np.nan)
J = np.full_like(C, np.nan)
for di in range(-SHIFT, SHIFT + 1):
    for dj in range(-SHIFT, SHIFT + 1):
        s = np.roll(np.roll(water_m, di, axis=0), dj, axis=1)
        C[di + SHIFT, dj + SHIFT] = np.corrcoef(s[win].astype(float), 1 - cz[win])[0, 1]
        J[di + SHIFT, dj + SHIFT] = (s & water_c & win).sum() / max(((s | water_c) & win).sum(), 1)
print("correlation of GSHHG water with 1 - coherence, rows = GSHHG moved by -4..4 cells in azimuth, columns in range:")
for di in range(-SHIFT, SHIFT + 1):
    print(f"  {di:+d}: " + " ".join(f"{v:.3f}" for v in C[di + SHIFT]))
k = np.unravel_index(np.nanargmax(C), C.shape)
q = lambda f: 0.5 * (f[0] - f[2]) / (f[0] - 2 * f[1] + f[2])
fy = q(C[k[0] - 1:k[0] + 2, k[1]]) if 0 < k[0] < 2 * SHIFT else 0.0
fx = q(C[k[0], k[1] - 1:k[1] + 2]) if 0 < k[1] < 2 * SHIFT else 0.0
kj = np.unravel_index(np.nanargmax(J), J.shape)
dy, dx = c.y[1] - c.y[0], c.x[1] - c.x[0]
print(f"best: GSHHG moved {k[0] - SHIFT + fy:+.2f} cells in azimuth ({(k[0] - SHIFT + fy) * dy:+.1f} lines), "
      f"{k[1] - SHIFT + fx:+.2f} in range ({(k[1] - SHIFT + fx) * dx:+.1f} samples), correlation {C[k]:.3f} (at 0 / 0 "
      f"{C[SHIFT, SHIFT]:.3f}); Jaccard best at {kj[0] - SHIFT:+d} / {kj[1] - SHIFT:+d} ({J[kj]:.3f}, at 0 / 0 "
      f"{J[SHIFT, SHIFT]:.3f})")
