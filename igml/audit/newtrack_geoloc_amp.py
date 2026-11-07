"""geoloc_amp.py's test on a NEW track's subswaths (10-29): the observed amplitude of a pair (amp1 + amp2, decimated by
our conv like the phase) against the range slope of the DEM in radar coordinates (F<n>/topo/topo_ra.grd, from the
base's trans.dat), per 32 x 6 cell centred on its declared node. amp1 / amp2 / topo_ra.grd are stored REVERSED in
azimuth (AUDIT.md traps: the row at node y holds line y_min + y_max - y) and their extents differ, so both are flipped
to line order first (geoloc_amp.py compared them reversed: an azimuth term of the extents' difference). The estimator is
geoloc_amp.py's (tiles, quadratic peak, refined median), checked on the simulation against itself; its own bias is
measured on a grid the Rosamond reflectors place at ~0 (T071 F2 of the 10-29 fix test). Read only.

usage (igml env): python newtrack_geoloc_amp.py INTF_ROOT SUBSWATHS PAIR [PAIR ...]
       INTF_ROOT holds F<n>/intf_all/<pair>/amp{1,2}.grd and F<n>/topo/topo_ra.grd; SUBSWATHS e.g. 3 or 1,2
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from geoloc_amp import TILE, fshift, highpass, refine            # noqa: E402
from igml.real.grd import read_grd                               # noqa: E402

root, subs, pairs = Path(sys.argv[1]), [int(n) for n in sys.argv[2].split(",")], sys.argv[3:]
for n in subs:
    t = read_grd(root / f"F{n}/topo/topo_ra.grd")
    tx, ty = t.x[1] - t.x[0], t.y[1] - t.y[0]
    h = t.z.astype(np.float64)[::-1]                                         # line order (stored reversed)
    slope = np.full_like(h, np.nan)
    slope[:, :-1] = np.diff(h, axis=1) / tx                                  # d h / d range sample
    for pair in pairs:
        d = root / f"F{n}/intf_all/{pair}"
        g1, g2 = read_grd(d / "amp1.grd"), read_grd(d / "amp2.grd")
        x, y = g1.x, g1.y
        dx, dy = x[1] - x[0], y[1] - y[0]
        kx, ky = int(round(dx / tx)), int(round(dy / ty))
        c0 = int(round((x[0] - dx / 2 - (t.x[0] - tx / 2)) / tx))
        r0 = int(round((y[0] - dy / 2 - (t.y[0] - ty / 2)) / ty))
        S = slope[r0:r0 + ky * y.size, c0:c0 + kx * x.size]
        ny, nx = S.shape[0] // ky, S.shape[1] // kx                          # cells the topography fully covers
        S = np.nanmean(S[:ny * ky, :nx * kx].reshape(ny, ky, nx, kx), axis=(1, 3))
        S[~np.isfinite(S)] = 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            obs = (0.5 * (np.log(g1.z[::-1]) + np.log(g2.z[::-1])))[:ny, :nx]   # line order
        good = np.isfinite(obs)                                              # amp 0 (no data) -> log -inf
        obs[~good] = np.nan
        o = highpass(obs)
        o[~good] = 0.0
        base = S
        hb = highpass(base)
        min_std = np.percentile([hb[i:i + TILE, j:j + TILE].std() for i in range(0, ny - TILE, TILE)
                                 for j in range(0, nx - TILE, TILE)], 70)
        print(f"F{n} {pair}: grid x {x[0]:g}..{x[-1]:g} step {dx:g}, y {y[0]:g}..{y[-1]:g} step {dy:g}; topo_ra "
              f"{t.z.shape} at {tx:g} x {ty:g} -> {kx} x {ky} per cell from column {c0}, row {r0}", flush=True)
        for sy, sx in ((0.0, 0.0), (0.5, 0.5)):
            dd, e = refine(highpass(fshift(base, sy, sx)), base, min_std)
            print(f"  estimator check: shift ({sy:+.2f}, {sx:+.2f}) cells -> ({dd[0]:+.3f}, {dd[1]:+.3f}) from "
                  f"{len(e)} tiles", flush=True)
        dd, e = refine(o, base, min_std)
        if dd is None:
            print(f"  only {len(e)} tiles correlate -- no estimate", flush=True)
            continue
        sy = 1.4826 * np.median(np.abs(e[:, 0] - np.median(e[:, 0]))) / np.sqrt(len(e))
        sx = 1.4826 * np.median(np.abs(e[:, 1] - np.median(e[:, 1]))) / np.sqrt(len(e))
        print(f"  observed - simulated = azimuth {dd[0]:+.3f} +- {sy:.3f} cells ({dd[0] * dy:+.2f} lines), range "
              f"{dd[1]:+.3f} +- {sx:.3f} cells ({dd[1] * dx:+.1f} samples); {len(e)} tiles, median peak corr "
              f"{np.median(e[:, 2]):.2f}", flush=True)
