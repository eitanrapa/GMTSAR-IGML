"""Where the Rosamond corner reflectors land in OUR geocoded product (2026-10-29, the end-to-end check of the
geolocation chain): a lab pair's merged corr.grd through igml.real.gmtsar.Geocoder and merge/trans.dat (the
conversion's own code path), with and without the F2 column fix, on a 0.0002-deg grid around the array. A trihedral
dominates its radar cell, so the cell's coherence is ~0.9-1 against ~0.4 for the playa; the half-maximum centroid of
that blob, minus the reflector's surveyed position, is the product's geolocation error there, in metres east / north.
The same error predicted from radar offsets (declared - true, samples / lines) through the local Jacobian of
trans.dat is printed beside it. Read only.

usage (igml env): python cr_geocode.py TRACK_PRODUCT PAIR CR_CSV [IDS]
"""
import sys
from pathlib import Path

import numpy as np

from igml.real.grd import read_coords, read_grd
from igml.real.gmtsar import Geocoder, TargetGrid, shift_f2

LON0, LON1, LAT0, LAT1, STEP = -118.115, -118.035, 34.785, 34.825, 0.0002
R_M = 150.0                                                    # search radius around each reflector [m]


def jacobian(trans, lon0, lat0, rad=0.01):
    """d(r, a) / d(east m, north m) from trans.dat nodes near (lon0, lat0)."""
    t = np.memmap(trans, dtype="<f8", mode="r").reshape(-1, 5)
    out = []
    for s in range(0, t.shape[0], 5_000_000):
        b = np.asarray(t[s:s + 5_000_000])
        m = (np.abs(b[:, 3] - lon0) < rad) & (np.abs(b[:, 4] - lat0) < rad)
        out.append(b[m])
    b = np.concatenate(out)
    e = (b[:, 3] - lon0) * 111320.0 * np.cos(np.radians(lat0))
    n = (b[:, 4] - lat0) * 110950.0
    A = np.c_[e, n, np.ones_like(e)]
    cr = np.linalg.lstsq(A, b[:, 0], rcond=None)[0]
    ca = np.linalg.lstsq(A, b[:, 1], rcond=None)[0]
    return np.array([[cr[0], cr[1]], [ca[0], ca[1]]]), len(b)


if __name__ == "__main__":
    prod, pair, csv = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    keep = set(sys.argv[4].split(",")) if len(sys.argv) > 4 else None
    crs = []
    for line in open(csv).read().splitlines()[1:]:
        f = line.split(",")
        if keep is None or f[0] in keep:
            crs.append((f[0], float(f[2]), float(f[1])))
    grid = TargetGrid(LON0, LAT1, STEP, int(round((LAT1 - LAT0) / STEP)), int(round((LON1 - LON0) / STEP)))
    pdir = prod / "merge" / pair
    x, y = read_coords(pdir / "corr.grd")
    geo = Geocoder.build(prod / "merge" / "trans.dat", grid, x, y)
    J, nn = jacobian(prod / "merge" / "trans.dat", np.mean([c[1] for c in crs]), np.mean([c[2] for c in crs]))
    print(f"Jacobian from {nn} trans.dat nodes: d range / d (E, N) = {J[0, 0]:+.4f}, {J[0, 1]:+.4f} samples/m; "
          f"d azimuth / d (E, N) = {J[1, 0]:+.4f}, {J[1, 1]:+.4f} lines/m", flush=True)
    Ji = np.linalg.inv(J)
    log = (pdir / "merge_log").read_text()
    ovl = [int(v) for v in log.split("ovl12,23:")[1].split("\n")[0].split(",")]
    n1 = int(log.split("n1 =")[1].split()[0]); n2 = int(log.split("n2 =")[1].split()[0])
    def width(n):
        try:
            return read_coords(prod / f"F{n}/intf_all" / pair / "corr.grd")[0].size
        except FileNotFoundError:                               # subswath grids deleted (T064): the PRM's width
            P = {k.strip(): v.strip() for k, s_, v in (l.partition("=") for l in open(prod / f"F{n}/topo/master.PRM")) if s_}
            return int(P["num_rng_bins"]) // 32
    w1, w2 = width(1), width(2)
    s12, s23 = w1 - ovl[0] + n1, w1 + w2 - ovl[0] - 1 - (ovl[1] - n2)
    full = read_grd(pdir / "corr.grd", rows=geo.rows).z.astype(np.float32)
    E = (np.arange(grid.width) * STEP + LON0)
    N = (LAT1 - np.arange(grid.length) * STEP)
    for tag in ("lab merge (no fix)", "F2 column fix"):
        z = full if tag.startswith("lab") else shift_f2(full.copy(), s12, s23)
        g = geo.to_grid(z[:, geo.cols[0]:geo.cols[1]], min_frac=0.0)
        res = []
        for cid, lon, lat in crs:
            de = (E[None, :] - lon) * 111320.0 * np.cos(np.radians(lat))
            dn = (N[:, None] - lat) * 110950.0 + 0 * de
            de = de + 0 * dn
            m = (np.hypot(de, dn) < R_M) & np.isfinite(g)
            if m.sum() < 10:
                continue
            v = g[m]
            pk, bg = v.max(), np.median(g[(np.hypot(de, dn) < 3 * R_M) & np.isfinite(g)])
            if pk - bg < 0.25:
                continue
            h = m & (g >= bg + 0.5 * (pk - bg))
            w = g[h] - bg
            ce, cn = np.sum(w * de[h]) / w.sum(), np.sum(w * dn[h]) / w.sum()
            res.append((cid, ce, cn, pk, bg))
        print(f"\n{tag}, s12 {s12}, s23 {s23}:", flush=True)
        for cid, ce, cn, pk, bg in res:
            print(f"  CR {cid}: product - surveyed = east {ce:+6.1f} m, north {cn:+6.1f} m (peak coh {pk:.2f}, background {bg:.2f})")
        if res:
            a = np.array([(r[1], r[2]) for r in res])
            med = np.median(a, axis=0)
            dra = J @ med
            print(f"  median of {len(a)}: east {med[0]:+.1f} m, north {med[1]:+.1f} m ({np.hypot(*med):.0f} m) = radar "
                  f"range {dra[0]:+.1f} samples, azimuth {dra[1]:+.2f} lines")
    for dr_, da_ in ((16, 3), (42, 6), (10, 6), (25, 6), (31, 3)):
        en = Ji @ np.array([dr_, da_])
        print(f"prediction: declared - true = {dr_:+d} samples, {da_:+d} lines -> product - truth = east {en[0]:+.1f} m, "
              f"north {en[1]:+.1f} m")
