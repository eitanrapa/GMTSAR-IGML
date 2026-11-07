"""Where GMTSAR's decimated grids really are (2026-10-28, audit I2; user: "geolocation error might not be right by the
way. It seems that all Sentinel processors have a 60 meter shift, it's just a calibration thing"). Per subswath: the
observed amplitude (amp1 + amp2 of a pair, decimated by conv like the phase) against a SIMULATED amplitude from the DEM
through GMTSAR's own geometry -- F<n>/topo/topo_ra.grd, the DEM in radar coordinates (2-sample cells, from trans.dat):
its slope along range, d h / d range, is what makes slopes facing the radar bright (foreshortening), averaged over each
32 x 6 decimated cell CENTRED ON ITS DECLARED NODE (32 m + 16, 6 r + 3), the position igml.real.gmtsar.Geocoder gives
the cell, so the measured offset (observed - simulated, in cells) is exactly the error our geocoding carries: the half-cell of conv.c predicts +0.5 cell in range (+16 samples)
AND +0.5 in azimuth (+3 lines); a range-delay or timing calibration moves one direction by an amount unrelated to the
32 x 6 decimation; a DEM georeferencing error moves both. Tiles of 64 x 64 cells over textured terrain, normalised
cross-correlation at -3..3 cells, a quadratic through the peak, the median over tiles refined by shifting the simulation
by it until the residual is ~0 (where the quadratic is unbiased); the estimator is checked on the simulation against
itself shifted by known fractions. Read only (the lab's stack over the Mac's NFS mount).

usage (igml env): python geoloc_amp.py TRACK_PRODUCT PAIR [PAIR ...]
"""
import sys
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

TILE, LAG = 64, 3


def smooth(a, s):
    """Gaussian blur (FFT, wrap-free by padding)."""
    from numpy.fft import irfft2, rfft2
    p = int(4 * s)
    b = np.pad(a, p, mode="reflect")
    ky = np.fft.fftfreq(b.shape[0])[:, None]
    kx = np.fft.rfftfreq(b.shape[1])[None, :]
    g = np.exp(-2 * (np.pi * s) ** 2 * (kx ** 2 + ky ** 2))
    return irfft2(rfft2(b) * g, s=b.shape)[p:-p, p:-p]


def highpass(a, s=8.0):
    a = np.where(np.isfinite(a), a, np.nanmean(a))
    return a - smooth(a, s)


def fshift(a, dy, dx):
    """a shifted by (dy, dx) cells (Fourier), for the estimator check."""
    ky = np.fft.fftfreq(a.shape[0])[:, None]
    kx = np.fft.fftfreq(a.shape[1])[None, :]
    return np.real(np.fft.ifft2(np.fft.fft2(a) * np.exp(-2j * np.pi * (ky * dy + kx * dx))))


def offsets(obs, sim, min_std):
    """Per tile, the (dy, dx) in cells by which obs is sim shifted (obs(i) = sim(i - d))."""
    out = []
    H, W = obs.shape
    for i in range(LAG, H - TILE - LAG, TILE):
        for j in range(LAG, W - TILE - LAG, TILE):
            s = sim[i:i + TILE, j:j + TILE]
            if s.std() < min_std:
                continue
            c = np.full((2 * LAG + 1, 2 * LAG + 1), np.nan)
            for u in range(-LAG, LAG + 1):
                for v in range(-LAG, LAG + 1):
                    o = obs[i + u:i + u + TILE, j + v:j + v + TILE]
                    c[u + LAG, v + LAG] = np.corrcoef(o.ravel(), s.ravel())[0, 1]
            k = np.unravel_index(np.nanargmax(c), c.shape)
            if not (0 < k[0] < 2 * LAG and 0 < k[1] < 2 * LAG) or c[k] < 0.2:
                continue
            fy = c[k[0] - 1:k[0] + 2, k[1]]
            fx = c[k[0], k[1] - 1:k[1] + 2]
            q = lambda f: 0.5 * (f[0] - f[2]) / (f[0] - 2 * f[1] + f[2])
            out.append((k[0] - LAG + q(fy), k[1] - LAG + q(fx), c[k]))
    return np.array(out)


def refine(obs, sim, min_std, it=4):
    """Median offset over tiles, refined: shift the simulation by the estimate and measure what is left."""
    d = np.zeros(2)
    for _ in range(it):
        e = offsets(obs, highpass(fshift(sim, *d)), min_std)
        if len(e) < 10:
            return None, e
        d += np.median(e[:, :2], axis=0)
    return d, e


if __name__ == "__main__":
    prod, pairs = Path(sys.argv[1]), sys.argv[2:]
    for n in (1, 2, 3):
        topo = None
        for pair in pairs:
            d = prod / f"F{n}/intf_all/{pair}"
            g1, g2 = read_grd(d / "amp1.grd"), read_grd(d / "amp2.grd")
            x, y = g1.x, g1.y
            dx, dy = x[1] - x[0], y[1] - y[0]
            with np.errstate(divide="ignore", invalid="ignore"):
                obs = 0.5 * (np.log(g1.z) + np.log(g2.z))
            obs[~np.isfinite(obs)] = np.nan
            if topo is None:
                t = read_grd(prod / f"F{n}/topo/topo_ra.grd")
                tx, ty = t.x[1] - t.x[0], t.y[1] - t.y[0]
                kx, ky = int(round(dx / tx)), int(round(dy / ty))          # topo cells per decimated cell
                # the decimated cell m spans samples [x_m - dx/2, x_m + dx/2): topo columns from there
                c0 = int(round((x[0] - dx / 2 - (t.x[0] - tx / 2)) / tx))
                r0 = int(round((y[0] - dy / 2 - (t.y[0] - ty / 2)) / ty))
                h = t.z.astype(np.float64)
                slope = np.full_like(h, np.nan)
                slope[:, :-1] = np.diff(h, axis=1) / tx                          # d h / d range sample
                S = slope[r0:r0 + ky * y.size, c0:c0 + kx * x.size]
                S = np.nanmean(S.reshape(y.size, ky, x.size, kx), axis=(1, 3))
                S[~np.isfinite(S)] = 0.0
                base = S
                topo = highpass(base)
                min_std = np.percentile([topo[i:i + TILE, j:j + TILE].std() for i in range(0, topo.shape[0] - TILE, TILE)
                                         for j in range(0, topo.shape[1] - TILE, TILE)], 70)
                print(f"F{n}: topo_ra {t.z.shape} at {tx:g} x {ty:g}, {kx} x {ky} per decimated cell from column {c0}, row "
                      f"{r0}", flush=True)
                for sy, sx in ((0.0, 0.0), (0.5, 0.5), (0.25, -0.3)):
                    dd, e = refine(highpass(fshift(base, sy, sx)), base, min_std)
                    print(f"F{n} estimator check: shift ({sy:+.2f}, {sx:+.2f}) cells -> ({dd[0]:+.3f}, {dd[1]:+.3f}) from "
                          f"{len(e)} tiles", flush=True)
            o = highpass(obs)
            o[~np.isfinite(obs)] = 0.0
            dd, e = refine(o, base, min_std)
            if dd is None:
                print(f"F{n} {pair}: only {len(e)} tiles correlate -- no estimate", flush=True)
                continue
            sy = 1.4826 * np.median(np.abs(e[:, 0] - np.median(e[:, 0]))) / np.sqrt(len(e))
            sx = 1.4826 * np.median(np.abs(e[:, 1] - np.median(e[:, 1]))) / np.sqrt(len(e))
            print(f"F{n} {pair}: observed - simulated = azimuth {dd[0]:+.3f} +- {sy:.3f} cells ({dd[0] * dy:+.2f} lines), "
                  f"range {dd[1]:+.3f} +- {sx:.3f} cells ({dd[1] * dx:+.1f} samples); {len(e)} tiles, median peak corr "
                  f"{np.median(e[:, 2]):.2f}", flush=True)
