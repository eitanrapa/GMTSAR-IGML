"""The geolocation test of geoloc_amp.py at FULL resolution (2026-10-28, audit I2): the super-master's own SLC (no
resampling, no conv decimation) against the range slope of topo_ra.grd, both averaged in OUR boxes of 32 samples x 6
lines ([32 m, 32 m + 32) x [6 k, 6 k + 6), centres known; 8 x 2 boxes leave too much speckle to correlate). What
geoloc_amp.py measures on the decimated grids = this + whatever conv's node convention adds; the half-cell of conv.c
(+16 samples, +3 lines) is absent here by construction, a timing / range calibration or a DEM offset is not. Sample i,
line l of the SLC sit at range sample i and line l (GMTSAR's PRM convention); a block's centre is its mean index.
Read only (the lab's stack over the Mac's NFS mount).

usage (igml env): python geoloc_full.py TRACK_PRODUCT SUPERMASTER_YYYYMMDD [LINE0 NLINES]
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from geoloc_amp import TILE, fshift, highpass, offsets            # noqa: E402
from igml.real.grd import read_grd                               # noqa: E402

prod, sm = Path(sys.argv[1]), sys.argv[2]
L0 = int(sys.argv[3]) if len(sys.argv) > 3 else 7998       # a multiple of 6 and of 2
NL = int(sys.argv[4]) if len(sys.argv) > 4 else 12000
BX, BY = 32, 6          # boxes of our own over the raw SLC: the speckle of GMTSAR's 32 x 6 grids, centres known


def refine(obs, sim, min_std, it=4):
    d = np.zeros(2)
    for _ in range(it):
        e = offsets(obs, highpass(fshift(sim, *d)), min_std)
        if len(e) < 10:
            return None, e
        d += np.median(e[:, :2], axis=0)
    return d, e


for n in (1, 2, 3):
    prm = {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(prod / f"F{n}/raw/S1_{sm}_ALL_F{n}.PRM"))
           if s}
    W = int(prm["bytes_per_line"]) // 4
    slc = np.memmap(prod / f"F{n}/raw/S1_{sm}_ALL_F{n}.SLC", dtype="<i2", mode="r")
    nl = slc.size // (2 * W)
    l1 = min(L0 + NL, nl) // BY * BY
    a = slc[L0 * 2 * W:l1 * 2 * W].reshape(-1, W, 2).astype(np.float32)
    p = a[..., 0] ** 2 + a[..., 1] ** 2
    del a
    nx = W // BX
    p = p[:, :nx * BX].reshape(-1, BY, nx, BX).mean(axis=(1, 3))       # block (i, j): lines L0+2i.., samples 8j..
    with np.errstate(divide="ignore"):
        obs = np.log(p)
    obs[~np.isfinite(obs)] = np.nan
    t = read_grd(prod / f"F{n}/topo/topo_ra.grd")
    tx, ty = t.x[1] - t.x[0], t.y[1] - t.y[0]                          # 2 x 2: node 2k+1 covers [2k, 2k+2)
    c0, r0 = int(round((0 - (t.x[0] - tx / 2)) / tx)), int(round((L0 - (t.y[0] - ty / 2)) / ty))
    kx, ky = int(BX / tx), int(BY / ty)
    h = t.z.astype(np.float64)
    slope = np.full_like(h, np.nan)
    slope[:, :-1] = np.diff(h, axis=1) / tx
    rows = obs.shape[0]
    S = slope[r0:r0 + ky * rows, c0:c0 + kx * nx]
    rows = min(rows, S.shape[0] // ky)
    obs, S = obs[:rows], S[:ky * rows]
    sim = np.nanmean(S.reshape(rows, ky, nx, kx), axis=(1, 3))
    sim[~np.isfinite(sim)] = 0.0
    sim_h = highpass(sim)
    min_std = np.percentile([sim_h[i:i + TILE, j:j + TILE].std() for i in range(0, rows - TILE, TILE)
                             for j in range(0, nx - TILE, TILE)], 70)
    d0, e0 = refine(highpass(fshift(sim, 0.25, -0.3)), sim, min_std)
    o = highpass(obs)
    o[~np.isfinite(obs)] = 0.0
    d, e = refine(o, sim, min_std)
    if d is None:
        print(f"F{n}: only {len(e)} tiles correlate -- no estimate", flush=True)
        continue
    print(f"F{n} full resolution, lines {L0}..{L0 + BY * rows}, blocks {BY} lines x {BX} samples: estimator check "
          f"(+0.25, -0.30) -> ({d0[0]:+.3f}, {d0[1]:+.3f}); observed - simulated = azimuth {d[0] * BY:+.2f} lines, range "
          f"{d[1] * BX:+.2f} samples ({len(e)} tiles, median peak corr {np.median(e[:, 2]):.2f})", flush=True)
