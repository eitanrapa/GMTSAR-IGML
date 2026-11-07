"""Where conv.c centres its decimated cells (2026-10-28, audit I2), with no terrain model and no calibration involved:
the SAME image twice -- GMTSAR's amp1.grd of a pair (the reference SLC through filter.csh / conv, nodes declared at
32 m + 16, 6 r + 3) against our own boxes of 32 samples x 6 lines over that reference's SLC ([32 m, 32 m + 32) x
[6 r, 6 r + 6): centres 32 m + 15.5, 6 r + 2.5, known). Same speckle, so the offset (amp1 - ours, in cells) is where
conv's filter sits relative to the declared node: the half-cell of the audit (filter centred on 32 m, 6 r) gives
+0.48 cell in range and +0.42 in azimuth; a filter centred on the declared node gives -0.02 / -0.08. Tiles of 64 x 64
cells everywhere there is data, refined as geoloc_amp.py; the estimator is checked against itself shifted. Read only.

usage (igml env): python geoloc_conv.py TRACK_PRODUCT PAIR [n ...]
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from geoloc_amp import fshift, highpass, offsets                  # noqa: E402
from igml.real.grd import read_grd                               # noqa: E402

prod, pair = Path(sys.argv[1]), sys.argv[2]
subs = [int(n) for n in sys.argv[3:]] or [1, 2, 3]


def refine(obs, sim, it=4):
    d = np.zeros(2)
    for _ in range(it):
        e = offsets(obs, highpass(fshift(sim, *d)), 0.0)
        if len(e) < 10:
            return None, e
        d += np.median(e[:, :2], axis=0)
    return d, e


for n in subs:
    d = prod / f"F{n}/intf_all/{pair}"
    a = min(p.name[3:11] for p in d.glob(f"S1_*_ALL_F{n}.PRM"))
    P = {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(prod / f"F{n}/raw/S1_{a}_ALL_F{n}.PRM")) if s}
    W = int(P["bytes_per_line"]) // 4
    g = read_grd(d / "amp1.grd")
    rows, nx = g.z.shape
    slc = np.memmap(prod / f"F{n}/raw/S1_{a}_ALL_F{n}.SLC", dtype="<i2", mode="r")
    rows = min(rows, slc.size // (2 * W) // 6)
    blk = slc[:rows * 6 * 2 * W].reshape(-1, W, 2)
    own = np.empty((rows, nx))
    for i in range(rows):
        b = blk[6 * i:6 * i + 6, :nx * 32].astype(np.float32)
        own[i] = (b[..., 0] ** 2 + b[..., 1] ** 2).reshape(6, nx, 32).mean(axis=(0, 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        o = np.log(own)
        m = np.log(g.z[:rows].astype(np.float64))
    ok = np.isfinite(o) & np.isfinite(m)
    o, m = np.where(ok, o, np.nan), np.where(ok, m, np.nan)
    oh, mh = highpass(o), highpass(m)
    oh[~ok], mh[~ok] = 0.0, 0.0
    chk, _ = refine(highpass(fshift(np.nan_to_num(o, nan=np.nanmean(o)), 0.25, -0.3)), np.nan_to_num(o, nan=np.nanmean(o)))
    dd, e = refine(mh, np.nan_to_num(o, nan=np.nanmean(o)))
    if dd is None:
        print(f"F{n}: only {len(e)} tiles -- no estimate", flush=True)
        continue
    print(f"F{n} {pair} (reference {a}): amp1.grd - our boxes = azimuth {dd[0]:+.3f} cells ({dd[0] * 6:+.2f} lines), range "
          f"{dd[1]:+.3f} cells ({dd[1] * 32:+.1f} samples); {len(e)} tiles, median peak corr {np.median(e[:, 2]):.2f}; "
          f"estimator check (+0.25, -0.30) -> ({chk[0]:+.3f}, {chk[1]:+.3f})", flush=True)
