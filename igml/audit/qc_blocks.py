# The qc's 12 blocks per subswath for three dates, every pair: complex coherence p90 (5 x 20 windows) and amplitude
# correlation per block (|a| vs |b|, normalised; robust to lost coherence, ~0 when the content differs)
import sys
import numpy as np
from pathlib import Path

def prm(p):
    return {a.strip(): b.strip() for a, s, b in (l.partition("=") for l in open(p)) if s}

def blocks(slc, nblk=12, rows=50):
    p = prm(slc.with_suffix(".PRM")); w, nl = int(p["bytes_per_line"]) // 4, int(p["num_lines"])
    out = []
    with open(slc, "rb") as fh:
        for r in np.linspace(0.02 * nl, 0.98 * nl - rows, nblk).astype(int):
            fh.seek(int(r) * w * 4)
            a = np.fromfile(fh, dtype="<i2", count=rows * w * 2).astype(np.float32).reshape(-1, w, 2)
            out.append((int(r), a[..., 0] + 1j * a[..., 1]))
    return out

def coh_p90(a, b, ly=5, lx=20):
    n, m = min(len(a), len(b)) // ly * ly, min(a.shape[1], b.shape[1]) // lx * lx
    A = a[:n, :m].reshape(n // ly, ly, m // lx, lx); B = b[:n, :m].reshape(n // ly, ly, m // lx, lx)
    pa, pb = (abs(A) ** 2).sum(axis=(1, 3)), (abs(B) ** 2).sum(axis=(1, 3)); ok = (pa > 0) & (pb > 0)
    c = abs((A * B.conj()).sum(axis=(1, 3)))[ok] / np.sqrt(pa[ok] * pb[ok])
    return float(np.percentile(c, 90)) if c.size else float("nan"), ok.mean()

def amp_corr(a, b):
    m = min(a.shape[1], b.shape[1]); A, B = abs(a[:, :m]), abs(b[:, :m]); ok = (A > 0) & (B > 0)
    if ok.sum() < 1000:
        return float("nan")
    # 4 x 16 looks first: the amplitude pattern, not the speckle
    k = lambda x: x[: x.shape[0] // 5 * 5, : x.shape[1] // 20 * 20].reshape(x.shape[0] // 5, 5, -1, 20).mean(axis=(1, 3))
    A, B, ok = k(A), k(B), k(ok.astype(float)) > 0.99
    return float(np.corrcoef(A[ok], B[ok])[0, 1])

root = Path(sys.argv[1]); lab = Path(sys.argv[2]); dates = sys.argv[3:]
import os
for f in os.environ.get("SUBS", "F1 F2 F3").split():
    paths = {d: (root if d == "20191113" else lab) / f / "raw" / f"S1_{d}_ALL_{f}.SLC" for d in dates}
    B = {d: blocks(p) for d, p in paths.items()}
    print(f"== {f}")
    pairs = [(dates[i], dates[j]) for i in range(len(dates)) for j in range(i + 1, len(dates))]
    print("block  line   " + "   ".join(f"{a[4:]}x{b[4:]} coh/amp" for a, b in pairs))
    for k in range(12):
        row = []
        for a, b in pairs:
            ra, xa = B[a][k]; rb, xb = B[b][k]
            c, _ = coh_p90(xa, xb); r = amp_corr(xa, xb)
            row.append(f"{c:5.2f}/{r:5.2f}      ")
        print(f"{k:5d} {B[dates[0]][k][0]:6d}   " + "".join(row))
