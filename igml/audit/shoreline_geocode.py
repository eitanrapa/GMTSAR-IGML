"""Coherence of merged products at trans.dat's DEM nodes around a box (2026-11-02, the Salton Sea shoreline test): per
product (DIR/merge/PAIR/corr.grd, DIR/merge/trans.dat) every DEM node in the box -> its radar cell (the products' own
x / y coordinates, pixel-registered as igml.real.gmtsar.Geocoder) -> coherence, plus the node's merged-grid column (for
the subswath) and a geocoded map at STEP degrees (Geocoder, nearest cell where a target cell holds no node). The nodes
and maps go to an npz for shoreline_fig.py. Read only.
usage (insar_sim env: scipy; PYTHONPATH=IGML): python shoreline_geocode.py W S E N STEP OUT.npz PAIR NAME=DIR [NAME=DIR ...]"""
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

from igml.real.grd import read_coords, read_grd
from igml.real.gmtsar import Geocoder, TargetGrid

(w, s, e, n), step, out, pair = map(float, sys.argv[1:5]), float(sys.argv[5]), sys.argv[6], sys.argv[7]
prods = [a.split("=", 1) for a in sys.argv[8:]]
grid = TargetGrid(w, n, step, int(round((n - s) / step)) + 1, int(round((e - w) / step)) + 1)
res = {"lon_map": w + np.arange(grid.width) * step, "lat_map": n - np.arange(grid.length) * step}
nodes = None
for name, d in prods:
    d = Path(d)
    trans, corr = d / "merge/trans.dat", d / f"merge/{pair}/corr.grd"
    if nodes is None:                                # the DEM nodes in the box (trans.dat is the lab's, shared)
        t = np.memmap(trans, dtype="<f8", mode="r").reshape(-1, 5)
        keep = []
        for k in range(0, t.shape[0], 5_000_000):
            b = np.asarray(t[k:k + 5_000_000])
            keep.append(b[(b[:, 3] >= w) & (b[:, 3] <= e) & (b[:, 4] >= s) & (b[:, 4] <= n)])
        nodes = np.concatenate(keep)
        res["node_lon"], res["node_lat"] = nodes[:, 3], nodes[:, 4]
        print(f"{nodes.shape[0]} DEM nodes in the box", flush=True)
    x, y = read_coords(corr)
    dx, dy = x[1] - x[0], y[1] - y[0]
    col = np.floor((nodes[:, 0] - (x[0] - dx / 2)) / dx).astype(np.int64)
    row = np.floor((nodes[:, 1] - (y[0] - dy / 2)) / dy).astype(np.int64)
    on = (col >= 0) & (col < x.size) & (row >= 0) & (row < y.size)
    z = read_grd(corr).z
    v = np.full(nodes.shape[0], np.nan, np.float32)
    v[on] = z[row[on], col[on]]
    res[f"{name}_coh"], res[f"{name}_col"] = v, np.where(on, col, -1)
    geo = Geocoder.build(trans, grid, x, y)
    g = geo.to_grid(read_grd(corr, rows=geo.rows, cols=geo.cols).z, min_frac=0.0)
    dist, (ii, jj) = distance_transform_edt(~np.isfinite(g), return_indices=True)
    res[f"{name}_map"] = np.where(dist * step * 111e3 < 150, g[ii, jj], np.nan).astype(np.float32)
    print(f"{name}: {int(on.sum())} nodes on the grid, coherence median {np.nanmedian(v):.2f}; map {g.shape}", flush=True)
np.savez_compressed(out, **res, pair=pair, names=np.array([p[0] for p in prods]))
print(f"-> {out}")
