"""A5: the interferometric phase step across each burst boundary, from a1_extract.py's windows (npz): each date's
switch line from its Doppler jump, then the mean phase per line on either side of the switch band (+-30 lines, the
band and a line of margin left out), one linear gradient fitted across both sides plus a step. 2026-10-30.
usage (insar_sim env): python burst_jump.py windows.npz"""
import sys
import numpy as np

z = np.load(sys.argv[1]); H = int(z["H"])
for k, c in enumerate(z["centres"]):
    r0 = c - H
    da, db = z[f"da{k}"], z[f"db{k}"]
    ja = np.abs(np.angle(np.exp(1j * np.diff(da)))); jb = np.abs(np.angle(np.exp(1j * np.diff(db))))
    sa, sb = r0 + 1 + int(np.argmax(ja)), r0 + 1 + int(np.argmax(jb))
    lo, hi = min(sa, sb) - 1, max(sa, sb) + 1
    ph = z[f"ph{k}"]
    lines = r0 + np.arange(ph.shape[0])
    p = np.exp(1j * ph).mean(axis=1)
    sel = ((lines >= lo - 30) & (lines < lo)) | ((lines > hi) & (lines <= hi + 30))
    l, q = lines[sel], np.unwrap(np.angle(p[sel]))
    A = np.c_[np.ones(l.size), l - c, (l > hi).astype(float)]
    coef, res, *_ = np.linalg.lstsq(A, q, rcond=None)
    err = np.sqrt(res[0] / (l.size - 3)) * np.sqrt(np.linalg.inv(A.T @ A)[2, 2])
    print(f"boundary near line {c}: switches {sa}/{sb}; step {np.degrees(coef[2]):+6.2f} deg (+- {np.degrees(err):.2f})")
