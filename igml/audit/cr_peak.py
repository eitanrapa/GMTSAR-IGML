"""Where the Rosamond corner reflectors really are in a GMTSAR SLC (2026-10-29; user: "let's dig into the geolocation
problem and find out what is really happening before doing anything"). The 28 trihedrals of JPL's Rosamond array
(surveyed, ellipsoidal heights; uavsar.jpl.nasa.gov calibration page) through GMTSAR's own geometry (SAT_llt2rat with
the stack's master.PRM + LED, run beforehand into cr_rat.txt) give a predicted (range sample, azimuth line); the SLC's
own peak near it, found on the complex samples (spectrum centred per axis, then FFT-oversampled 32x), is where the
reflector is. observed - predicted is GMTSAR's end-to-end geolocation error at a point target, with no DEM, no
decimation and no correlation estimator in between. A path delay (troposphere ~2.2 m zenith at 660 m) makes the
observed range LARGER by ~1 sample; nothing physical makes it 15-20 samples. Read only.

usage (seasat venv): python cr_peak.py PRM SLC CR_RAT [CR_IDS]
  CR_RAT: SAT_llt2rat output (r a z lon lat per line), in the order of CR_IDS (one id per line)
"""
import sys

import numpy as np

HALF_R, HALF_A, CHIP, OVS = 64, 24, 16, 32


def prm(path):
    return {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(path)) if s}


def centred(chip):
    """The chip with each axis' spectral centroid moved to zero (TOPS azimuth spectra are not at baseband)."""
    out = chip
    for ax in (0, 1):
        n = chip.shape[ax]
        S = np.abs(np.fft.fft(out, axis=ax)) ** 2
        S = S.sum(axis=1 - ax)
        f = np.fft.fftfreq(n)
        fc = np.angle(np.sum(S * np.exp(2j * np.pi * f))) / (2 * np.pi)
        k = np.arange(n).reshape((-1, 1) if ax == 0 else (1, -1))
        out = out * np.exp(-2j * np.pi * fc * k)
    return out


def subpixel(chip):
    """Peak of |chip| after FFT oversampling, in chip pixel units."""
    n0, n1 = chip.shape
    F = np.fft.fftshift(np.fft.fft2(centred(chip)))
    P = np.zeros((n0 * OVS, n1 * OVS), complex)
    o0, o1 = (n0 * OVS - n0) // 2, (n1 * OVS - n1) // 2
    P[o0:o0 + n0, o1:o1 + n1] = F
    up = np.abs(np.fft.ifft2(np.fft.ifftshift(P)))
    i, j = np.unravel_index(np.argmax(up), up.shape)
    return i / OVS, j / OVS


if __name__ == "__main__":
    P = prm(sys.argv[1])
    W, L = int(P["bytes_per_line"]) // 4, int(P["num_lines"])
    slc = np.memmap(sys.argv[2], dtype="<i2", mode="r")
    L = min(L, slc.size // (2 * W))
    slc = slc[:L * W * 2].reshape(L, W, 2)
    rat = np.loadtxt(sys.argv[3], ndmin=2)
    ids = [l.split()[-1] for l in open(sys.argv[4])] if len(sys.argv) > 4 else [str(i) for i in range(len(rat))]
    rows = []
    for cid, (rp, ap, *_rest) in zip(ids, rat):
        r0, a0 = int(round(rp)) - HALF_R, int(round(ap)) - HALF_A
        if r0 < 0 or a0 < 0 or r0 + 2 * HALF_R >= W or a0 + 2 * HALF_A >= L:
            print(f"CR {cid}: predicted ({rp:.1f}, {ap:.1f}) too near the edge", flush=True)
            continue
        w = slc[a0:a0 + 2 * HALF_A + 1, r0:r0 + 2 * HALF_R + 1].astype(np.float32)
        z = w[..., 0] + 1j * w[..., 1]
        I = np.abs(z) ** 2
        med = np.median(I[I > 0]) if np.any(I > 0) else np.nan
        # local maxima >= 20 dB above the window median; the one nearest the prediction is this reflector
        pk = []
        for i in range(1, I.shape[0] - 1):
            for j in range(1, I.shape[1] - 1):
                v = I[i, j]
                if v > 100 * med and v == I[i - 1:i + 2, j - 1:j + 2].max():
                    pk.append((i, j, v))
        if not pk:
            print(f"CR {cid}: no peak 20 dB above the clutter near ({rp:.1f}, {ap:.1f})", flush=True)
            continue
        pk = np.array(pk)
        dist = np.hypot((pk[:, 1] + r0 - rp) * 2.33, (pk[:, 0] + a0 - ap) * 14.0)       # ~metres, to rank
        k = int(np.argmin(dist))
        i, j = int(pk[k, 0]), int(pk[k, 1])
        h = CHIP // 2
        if i - h < 0 or j - h < 0 or i + h > I.shape[0] or j + h > I.shape[1]:
            print(f"CR {cid}: peak at the window edge", flush=True)
            continue
        di, dj = subpixel(z[i - h:i + h, j - h:j + h])
        ao, ro = a0 + i - h + di, r0 + j - h + dj
        snr = 10 * np.log10(pk[k, 2] / med)
        nxt = np.sort(dist)[1] if len(dist) > 1 else np.inf
        rows.append((ro - rp, ao - ap, snr))
        print(f"CR {cid}: predicted r {rp:9.2f} a {ap:9.2f} | observed r {ro:9.2f} a {ao:9.2f} | observed - predicted "
              f"r {ro - rp:+7.2f} samples a {ao - ap:+6.2f} lines | peak {snr:4.1f} dB over the median, "
              f"{len(pk)} peaks in the window, next nearest {nxt:.0f} m", flush=True)
    if rows:
        d = np.array(rows)
        good = d[:, 2] >= 30
        for tag, m in (("all", np.ones(len(d), bool)), (">= 30 dB", good)):
            if m.sum():
                print(f"{tag}: {m.sum()} reflectors, median observed - predicted: range {np.median(d[m, 0]):+.2f} samples "
                      f"(sd {d[m, 0].std():.2f}), azimuth {np.median(d[m, 1]):+.2f} lines (sd {d[m, 1].std():.2f})", flush=True)
