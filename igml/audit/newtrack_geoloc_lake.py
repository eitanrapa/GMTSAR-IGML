"""Absolute geolocation of a NEW track's grid at a lake (10-29): SRTM's water bodies are flattened, so the DEM nodes at
the lake's level (the most common height in a lon / lat box, +- TOL m) outline the lake on the ground; trans.dat puts
each node on the radar grid. Their water mask is cross-correlated with 1 - coherence of a pair (open water decorrelates)
for shifts of -L..L cells; the peak (quadratic, sub-cell) is where the data's lake sits relative to where trans.dat
says it is. GSHHG's mask (landmask_ra.grd, landmask.csh) is too coarse for this (blocks of ~20 x 10 cells). Read only.

usage (igml env): python newtrack_geoloc_lake.py corr.grd trans.dat dem.grd LON0 LON1 LAT0 LAT1 [PAIR_CORR ...]
"""
import sys

import numpy as np

from igml.real.grd import read_grd

TOL, L = 0.5, 6

t = np.memmap(sys.argv[2], dtype="<f8", mode="r").reshape(-1, 5)
dem = read_grd(sys.argv[3])
lon0, lon1, lat0, lat1 = map(float, sys.argv[4:8])
box = (t[:, 3] >= lon0) & (t[:, 3] <= lon1) & (t[:, 4] >= lat0) & (t[:, 4] <= lat1)
r, a, lon, lat = t[box, 0], t[box, 1], t[box, 3], t[box, 4]
jx = np.clip(np.round((lon - dem.x[0]) / (dem.x[1] - dem.x[0])).astype(int), 0, dem.x.size - 1)
iy = np.clip(np.round((lat - dem.y[0]) / (dem.y[1] - dem.y[0])).astype(int), 0, dem.y.size - 1)
h = dem.z[iy, jx].astype(np.float64)
vals, cnt = np.unique(np.round(h, 1), return_counts=True)
level = vals[np.argmax(cnt)]
lake = np.abs(h - level) <= TOL
print(f"box {lon0}..{lon1} E, {lat0}..{lat1} N: {box.sum()} DEM nodes, lake level {level:.1f} m ({cnt.max()} nodes at it), "
      f"{lake.sum()} lake nodes (+-{TOL} m), radar range {r.min():.0f}..{r.max():.0f}, azimuth {a.min():.0f}..{a.max():.0f}")
for f in [sys.argv[1]] + sys.argv[8:]:
    c = read_grd(f)
    dx, dy = c.x[1] - c.x[0], c.y[1] - c.y[0]
    col = np.floor((r - (c.x[0] - dx / 2)) / dx).astype(int)          # pixel-registered cells, centres declared
    row = np.floor((a - (c.y[0] - dy / 2)) / dy).astype(int)
    inn = (col >= L) & (col < c.x.size - L) & (row >= L) & (row < c.y.size - L)   # nodes L cells inside the grid
    col, row, lk = col[inn], row[inn], lake[inn]
    c0, c1, r0, r1 = col.min() - L, col.max() + L + 1, row.min() - L, row.max() + L + 1
    W = np.zeros((r1 - r0, c1 - c0))
    N = np.zeros_like(W)
    np.add.at(W, (row - r0, col - c0), lk.astype(float))
    np.add.at(N, (row - r0, col - c0), 1.0)
    water = np.where(N > 0, W / np.maximum(N, 1), np.nan)                 # share of the cell's DEM nodes on the lake
    z = c.z[max(r0, 0):r1, max(c0, 0):c1].astype(np.float64)
    if z.shape != water.shape:
        print(f"{f}: window outside the grid -- skipped")
        continue
    d = 1 - z
    C = np.full((2 * L + 1, 2 * L + 1), np.nan)
    core = (slice(L, -L), slice(L, -L))
    wv = water[core]
    for i in range(-L, L + 1):
        for j in range(-L, L + 1):
            dv = d[L + i:d.shape[0] - L + i, L + j:d.shape[1] - L + j]      # data at (cell + shift)
            ok = np.isfinite(wv) & np.isfinite(dv)
            C[i + L, j + L] = np.corrcoef(wv[ok], dv[ok])[0, 1]
    k = np.unravel_index(np.nanargmax(C), C.shape)
    q = lambda v: 0.5 * (v[0] - v[2]) / (v[0] - 2 * v[1] + v[2])
    fy = q(C[k[0] - 1:k[0] + 2, k[1]]) if 0 < k[0] < 2 * L else 0.0
    fx = q(C[k[0], k[1] - 1:k[1] + 2]) if 0 < k[1] < 2 * L else 0.0
    sy, sx = k[0] - L + fy, k[1] - L + fx
    zz = z[core]
    print(f"{f.split('/')[-2]}: the data's lake sits {sy:+.2f} cells in azimuth ({sy * dy:+.1f} lines), {sx:+.2f} in range "
          f"({sx * dx:+.1f} samples) from where trans.dat puts SRTM's; correlation {C[k]:.3f} (at 0 / 0 {C[L, L]:.3f}); "
          f"coherence on lake cells (> 90 % lake) median {np.nanmedian(zz[wv > 0.9]):.2f}, on dry cells "
          f"{np.nanmedian(zz[wv == 0]):.2f}")
