"""Where conv.c centres its decimated cells, measured again with the azimuth ORIENTATION handled (2026-10-29; replaces
geoloc_conv.py, whose 10-28 "refutation" of I2 is withdrawn: the Rosamond corner reflectors, cr_peak.py, sit at +16.0
samples / +3.0 lines from their SLC position in T071 F2's amp1.grd / amp2.grd).

GMTSAR writes radar grids in two azimuth orientations: in the lab's T071 subswath grids (amp1, corr, phasefilt) and
topo_ra.grd (dem2topo_ra.csh ends with FLIPUD), the row at declared y holds line YMAX - y; elsewhere the row at y holds
line y. A comparison that ignores this compares an image with its own mirror. Here: GMTSAR's amp1.grd of a pair
against our OWN 32 x 6 boxes of the same SLC (lines [6k, 6k + 6), samples [32 m, 32 m + 32): centres 6k + 2.5,
32 m + 15.5), in BOTH orientations; the one that correlates at zero lag is the grid's, the other ~0. In that
orientation (rows put in natural line order), the offset amp1 - boxes in cells: the half-cell of I2 (cell holds the
filter centred on line 6k, sample 32 m; node declared at 6k + 3, 32 m + 16) gives +0.42 / +0.48, a filter centred on
the declared node -0.08 / -0.02. Read only.

usage (igml env): python geoloc_conv2.py TRACK_PRODUCT PAIR [n ...]
"""
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from geoloc_amp import fshift, highpass, offsets                  # noqa: E402
from igml.real.grd import read_grd                               # noqa: E402

K0, K1 = (int(os.environ.get(k, v)) for k, v in (("K0", 1200), ("K1", 3300)))  # decimated rows used (blocks of 6 lines;
                                                                  # 10-29: env K0 / K1 for shorter subswaths, T100 F3)


def refine(obs, sim, it=4):
    d = np.zeros(2)
    for _ in range(it):
        e = offsets(obs, highpass(fshift(sim, *d)), 0.0)
        if len(e) < 10:
            return None, e
        d += np.median(e[:, :2], axis=0)
    return d, e


if __name__ == "__main__":
    prod, pair = Path(sys.argv[1]), sys.argv[2]
    for n in [int(s) for s in sys.argv[3:]] or [1, 2, 3]:
        d = prod / f"F{n}/intf_all/{pair}"
        a = min(p.name[3:11] for p in d.glob(f"S1_*_ALL_F{n}.PRM"))
        P = {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(d / f"S1_{a}_ALL_F{n}.PRM")) if s}
        W = int(P["bytes_per_line"]) // 4
        g = read_grd(d / "amp1.grd")
        ny, nx = g.z.shape
        dy = g.y[1] - g.y[0]
        ymax = g.y[-1] + dy / 2
        slc = np.memmap(prod / f"F{n}/raw/S1_{a}_ALL_F{n}.SLC", dtype="<i2", mode="r")
        lines = slc.size // (2 * W)
        blk = slc[:lines * 2 * W].reshape(lines, W, 2)
        own = np.empty((K1 - K0, nx))
        for i, k in enumerate(range(K0, K1)):
            b = blk[6 * k:6 * k + 6, :nx * 32].astype(np.float32)
            own[i] = (b[..., 0] ** 2 + b[..., 1] ** 2).reshape(6, nx, 32).mean(axis=(0, 2))
        with np.errstate(divide="ignore", invalid="ignore"):
            o = np.log(own)
            lz = np.log(g.z.astype(np.float64))
        o[~np.isfinite(o)] = np.nan
        lz[~np.isfinite(lz)] = np.nan
        nat = lz[K0:K1]                                             # row k = line block k
        rev = lz[::-1][K0:K1] if ny * dy == ymax else None          # row k = line block k if stored reversed
        cc = {}
        for tag, m in (("natural", nat), ("reversed", rev)):
            if m is None or m.shape != o.shape:
                continue
            ok = np.isfinite(o) & np.isfinite(m)
            cc[tag] = np.corrcoef(highpass(np.where(ok, o, np.nan))[ok], highpass(np.where(ok, m, np.nan))[ok])[0, 1]
        tag = max(cc, key=cc.get)
        m = nat if tag == "natural" else rev
        ok = np.isfinite(o) & np.isfinite(m)
        mh = highpass(np.where(ok, m, np.nan))
        mh[~ok] = 0.0
        fill = np.nan_to_num(o, nan=np.nanmean(o))
        chk, _ = refine(highpass(fshift(fill, 0.25, -0.3)), fill)
        dd, e = refine(mh, fill)
        if dd is None:
            print(f"F{n} {pair}: orientation {tag} (zero-lag corr {cc}), only {len(e)} tiles -- no estimate", flush=True)
            continue
        print(f"F{n} {pair} (reference {a}): stored {tag} (zero-lag corr " + ", ".join(f"{k} {v:.2f}" for k, v in cc.items())
              + f"); amp1 - our boxes = azimuth {dd[0]:+.3f} cells ({dd[0] * 6:+.2f} lines), range {dd[1]:+.3f} cells "
              f"({dd[1] * 32:+.1f} samples); declared node - true centre = azimuth {dd[0] * 6 + 0.5:+.2f} lines, range "
              f"{dd[1] * 32 + 0.5:+.1f} samples; {len(e)} tiles, median peak corr {np.median(e[:, 2]):.2f}; estimator "
              f"check (+0.25, -0.30) -> ({chk[0]:+.3f}, {chk[1]:+.3f})", flush=True)
