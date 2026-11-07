"""Where surveyed corner reflectors land in a GEOCODED merged product (2026-11-01; generalises cr_geocode.py beyond
Rosamond): a merged pair's corr.grd through igml.real.gmtsar.Geocoder and the product's merge/trans.dat (the
conversion's own path) on a 0.0002-deg grid around the reflectors. A trihedral dominates its decimated cell (coherence
~1 against the ground's), so the half-maximum centroid of that blob minus the surveyed position is the product's
geolocation error there, in metres east / north, and through trans.dat's local Jacobian in radar samples / lines.
Writes the geocoded window to an npz for the figure. Read only.

usage (igml env): python cr_geocode_any.py PRODUCT PAIR CRS OUT.npz [LABEL]
  CRS: "id,lat,lon;id,lat,lon;..."  e.g. Pinon Flat D1 / D2 (UCSD's survey, WGS-84)
"""
import sys
from pathlib import Path

import numpy as np

from igml.real.grd import read_coords, read_grd
from igml.real.gmtsar import Geocoder, TargetGrid

# 11-06: R_M was 150 m, which cut the blob of a product more than ~100 m off and pulled its centroid back to the survey: stock T173 read 81 / 62 m for a true 116 / 110 m.
# Prefer display_amp in radar coordinates (notebooks/google_earth/src/make_geolocation_kmz.py, Layer.blob); the coherence blob sits on a coherent background.
STEP, PAD, R_M = 0.0002, 0.004, 350.0


def jacobian(trans, lon0, lat0, rad=0.01):
    """d(r, a) / d(east m, north m) from trans.dat nodes near (lon0, lat0)."""
    t = np.memmap(trans, dtype="<f8", mode="r").reshape(-1, 5)
    out = []
    for s in range(0, t.shape[0], 5_000_000):
        b = np.asarray(t[s:s + 5_000_000])
        out.append(b[(np.abs(b[:, 3] - lon0) < rad) & (np.abs(b[:, 4] - lat0) < rad)])
    b = np.concatenate(out)
    e = (b[:, 3] - lon0) * 111320.0 * np.cos(np.radians(lat0))
    n = (b[:, 4] - lat0) * 110950.0
    A = np.c_[e, n, np.ones_like(e)]
    cr, ca = np.linalg.lstsq(A, b[:, 0], rcond=None)[0], np.linalg.lstsq(A, b[:, 1], rcond=None)[0]
    return np.array([[cr[0], cr[1]], [ca[0], ca[1]]])


prod, pair, out = Path(sys.argv[1]), sys.argv[2], sys.argv[4]
label = sys.argv[5] if len(sys.argv) > 5 else prod.name
crs = [(f[0], float(f[1]), float(f[2])) for f in (c.split(",") for c in sys.argv[3].split(";"))]
lat0, lat1 = min(c[1] for c in crs) - PAD, max(c[1] for c in crs) + PAD
lon0, lon1 = min(c[2] for c in crs) - PAD, max(c[2] for c in crs) + PAD
grid = TargetGrid(lon0, lat1, STEP, int(round((lat1 - lat0) / STEP)), int(round((lon1 - lon0) / STEP)))
pdir = prod / "merge" / pair
x, y = read_coords(pdir / "corr.grd")
trans = prod / "merge" / "trans.dat"
geo = Geocoder.build(trans, grid, x, y)
J = jacobian(trans, np.mean([c[2] for c in crs]), np.mean([c[1] for c in crs]))
z = read_grd(pdir / "corr.grd", rows=geo.rows).z.astype(np.float32)
g = geo.to_grid(z[:, geo.cols[0]:geo.cols[1]], min_frac=0.0)
E = np.arange(grid.width) * STEP + lon0
N = lat1 - np.arange(grid.length) * STEP
res = []
for cid, lat, lon in crs:
    de = np.broadcast_to((E[None, :] - lon) * 111320.0 * np.cos(np.radians(lat)), g.shape)
    dn = np.broadcast_to((N[:, None] - lat) * 110950.0, g.shape)
    m = (np.hypot(de, dn) < R_M) & np.isfinite(g)
    if m.sum() < 10:
        print(f"  {cid}: outside the product")
        continue
    pk, bg = g[m].max(), np.median(g[(np.hypot(de, dn) < 3 * R_M) & np.isfinite(g)])
    h = m & (g >= bg + 0.5 * (pk - bg))
    w = g[h] - bg
    ce, cn = np.sum(w * de[h]) / w.sum(), np.sum(w * dn[h]) / w.sum()
    dra = J @ np.array([ce, cn])
    res.append((cid, ce, cn, dra[0], dra[1], pk, bg))
    print(f"  {label} {pair} {cid}: product - surveyed = east {ce:+6.1f} m, north {cn:+6.1f} m "
          f"({np.hypot(ce, cn):.0f} m; radar {dra[0]:+.1f} samples, {dra[1]:+.2f} lines; peak coh {pk:.2f}, ground {bg:.2f})")
np.savez(out, g=g, E=E, N=N, crs=np.array([(c[1], c[2]) for c in crs]), ids=np.array([c[0] for c in crs]),
         res=np.array([r[1:] for r in res]) if res else np.zeros((0, 6)), label=label, pair=pair)
print(f"-> {out}")
