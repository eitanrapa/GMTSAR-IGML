"""A5 (2026-10-31): the phase step of a pair across each burst boundary, measured in range blocks across the WHOLE swath
(burst_jump.py measures one 2048-sample cut at mid-swath; spectral diversity is a swath-wide median): per block, each
date's switch line from its Doppler jump, the dominant fringe removed, the mean phase per line on either side of the
switch band (+-30 lines, the band and a line of margin left out), one linear gradient + a step; then the median over
blocks and their spread.
usage (seasat, the venv): python steps_wide.py A.SLC B.SLC c1,c2,... [H]"""
import sys
import numpy as np


def prm(p):
    return {a.strip(): b.strip() for a, s, b in (l.partition("=") for l in open(p)) if s}


A, B = sys.argv[1], sys.argv[2]
centres = [int(c) for c in sys.argv[3].split(",")]
H = int(sys.argv[4]) if len(sys.argv) > 4 else 80
P = prm(A[:-4] + ".PRM"); w, nl = int(P["num_rng_bins"]), int(P["num_lines"])
ma = np.memmap(A, dtype="<i2", mode="r", shape=(nl, w, 2)); mb = np.memmap(B, dtype="<i2", mode="r", shape=(nl, w, 2))
allsteps = []
for c in centres:
    a = ma[c - H:c + H].astype(np.float32); a = a[..., 0] + 1j * a[..., 1]
    b = mb[c - H:c + H].astype(np.float32); b = b[..., 0] + 1j * b[..., 1]
    ok = (np.abs(a) > 0) & (np.abs(b) > 0)
    cols = np.where(ok.mean(axis=0) > 0.99)[0]
    steps = []
    for x0 in range(cols.min(), cols.max() - 2048, 2048):
        aa, bb = a[:, x0:x0 + 2048], b[:, x0:x0 + 2048]
        ja = np.abs(np.angle(np.exp(1j * np.angle((aa[1:] * np.conj(aa[:-1])).sum(axis=1)))))
        jb = np.abs(np.angle(np.exp(1j * np.angle((bb[1:] * np.conj(bb[:-1])).sum(axis=1)))))
        jd = np.abs(np.diff(np.angle((aa[1:] * np.conj(aa[:-1])).sum(axis=1))))
        sa = c - H + 1 + int(np.argmax(np.abs(np.angle(np.exp(1j * np.diff(np.angle((aa[1:] * np.conj(aa[:-1])).sum(axis=1)))))))) + 1
        sb = c - H + 1 + int(np.argmax(np.abs(np.angle(np.exp(1j * np.diff(np.angle((bb[1:] * np.conj(bb[:-1])).sum(axis=1)))))))) + 1
        lo, hi = min(sa, sb) - 1, max(sa, sb) + 1
        ifg = aa * np.conj(bb)
        F = np.fft.fft2(ifg); ky, kx = np.unravel_index(np.argmax(np.abs(F)), F.shape)
        yy, xx = np.mgrid[0:ifg.shape[0], 0:ifg.shape[1]]
        ifg = ifg * np.exp(-2j * np.pi * (ky * yy / ifg.shape[0] + kx * xx / ifg.shape[1]))
        z = ifg.reshape(ifg.shape[0], -1, 8).sum(axis=2)
        coh = np.abs(z.sum()) / np.abs(z).sum()
        p = np.exp(1j * np.angle(z)).mean(axis=1)
        lines = c - H + np.arange(ifg.shape[0])
        sel = ((lines >= lo - 30) & (lines < lo)) | ((lines > hi) & (lines <= hi + 30))
        l, q = lines[sel], np.unwrap(np.angle(p[sel]))
        M = np.c_[np.ones(l.size), l - c, (l > hi).astype(float)]
        coef = np.linalg.lstsq(M, q, rcond=None)[0]
        if np.abs(p[sel]).mean() > 0.2:                 # a block coherent enough to read
            steps.append(np.degrees(coef[2]))
    s = np.array(steps)
    allsteps += list(s)
    print(f"boundary near line {c}: {len(s)} blocks, step median {np.median(s):+6.1f} deg, "
          f"spread (IQR) {np.percentile(s, 25):+6.1f} .. {np.percentile(s, 75):+6.1f}")
print(f"all boundaries: median {np.median(allsteps):+.1f} deg over {len(allsteps)} blocks")
