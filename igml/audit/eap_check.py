"""Our EAP-corrected aligned SLC of an IPF 002.36 date against the lab's uncorrected one (2026-10-29): the ratio
ours * conj(lab) should be the elevation-antenna-pattern phase alone -- a smooth function of range per subswath, the
same in every burst -- and |ours| / |lab| = 1. Lines from the first quarter (the first frame: the lab's last frame
may sit a line off, audit A0). Prints the phase across range (deg) per subswath and the jump between subswaths at
their common slant range. Read only.

usage (seasat venv): python eap_check.py OURS_DIR LAB_RAW_PREFIX DATE   (LAB_RAW_PREFIX = .../product_T/F)
"""
import sys
from pathlib import Path

import numpy as np

ours, labp, date = Path(sys.argv[1]), sys.argv[2], sys.argv[3]


def prm(p):
    return {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(p)) if s}


prof = {}
for n in (1, 2, 3):
    po = prm(ours / f"F{n}/S1_{date}_ALL_F{n}.PRM")
    W, L = int(po["bytes_per_line"]) // 4, int(po["num_lines"])
    a = np.memmap(ours / f"F{n}/S1_{date}_ALL_F{n}.SLC", dtype="<i2", mode="r")
    b = np.memmap(f"{labp}{n}/raw/S1_{date}_ALL_F{n}.SLC", dtype="<i2", mode="r")
    r0, r1 = L // 8, L // 8 + 400
    za = a[r0 * W * 2:r1 * W * 2].astype(np.float32).view(np.complex64).reshape(-1, W)
    zb = b[r0 * W * 2:r1 * W * 2].astype(np.float32).view(np.complex64).reshape(-1, W)
    ok = (np.abs(za) > 0) & (np.abs(zb) > 0)
    rat = np.where(ok, za * np.conj(zb), 0)
    ph = np.degrees(np.angle(rat.sum(0)))                         # per range sample, over 400 lines
    amp = np.median(np.abs(za[ok]) / np.abs(zb[ok]))
    good = ok.sum(0) > 100
    x = np.flatnonzero(good)
    phu = np.degrees(np.unwrap(np.radians(ph[good])))
    k = np.linspace(0, x.size - 1, 9).astype(int)
    prof[n] = (float(po["near_range"]), float(po["rng_samp_rate"]), x, phu)
    print(f"F{n}: |ours| / |lab| median {amp:.4f}; EAP phase removed across range (deg, 9 points): "
          + " ".join(f"{phu[i]:+.1f}" for i in k) + f"; spread {phu.max() - phu.min():.1f} deg", flush=True)
c = 299792458.0
for n in (1, 2):
    nr1, fs, x1, p1 = prof[n]; nr2, _, x2, p2 = prof[n + 1]
    r_end = nr1 + x1[-40:] * c / 2 / fs                           # the far end of F_n
    r_beg = nr2 + x2[:40] * c / 2 / fs                            # the near end of F_n+1
    print(f"F{n}/F{n + 1}: removed phase at F{n}'s far end {np.median(p1[-40:]):+.1f} deg "
          f"({r_end.mean() / 1e3:.2f} km), at F{n + 1}'s near end {np.median(p2[:40]):+.1f} deg ({r_beg.mean() / 1e3:.2f} km) "
          f"-> step {np.median(p2[:40]) - np.median(p1[-40:]):+.1f} deg", flush=True)
