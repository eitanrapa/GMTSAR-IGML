# Per block along an aligned SLC pair: coherence p90 (5 x 20 windows) with the second image shifted -2..+2 lines;
# plus an interferogram cut (6 x 24 looks) over rows [c0, c1). npz -> stdout.
import io, sys
import numpy as np
A, B = sys.argv[1], sys.argv[2]; nl, w = int(sys.argv[3]), int(sys.argv[4]); c0, c1 = int(sys.argv[5]), int(sys.argv[6])
ma = np.memmap(A, dtype="<i2", mode="r", shape=(nl, w, 2)); mb = np.memmap(B, dtype="<i2", mode="r", shape=(nl, w, 2))
cx = lambda m, r0, r1: (lambda x: x[..., 0] + 1j * x[..., 1])(m[r0:r1].astype(np.float32))
def coh_p90(a, b, ly=5, lx=20):
    n, m = min(len(a), len(b)) // ly * ly, min(a.shape[1], b.shape[1]) // lx * lx
    A_ = a[:n, :m].reshape(n // ly, ly, m // lx, lx); B_ = b[:n, :m].reshape(n // ly, ly, m // lx, lx)
    pa, pb = (abs(A_) ** 2).sum(axis=(1, 3)), (abs(B_) ** 2).sum(axis=(1, 3)); ok = (pa > 0) & (pb > 0)
    c = abs((A_ * B_.conj()).sum(axis=(1, 3)))[ok] / np.sqrt(pa[ok] * pb[ok])
    return float(np.percentile(c, 90)) if c.size > 50 else np.nan
starts = np.linspace(0.01 * nl, 0.99 * nl - 60, 48).astype(int)
C = np.full((len(starts), 5), np.nan)
for i, r in enumerate(starts):
    a = cx(ma, r, r + 54); b = cx(mb, r, r + 54)
    for j, s in enumerate((-2, -1, 0, 1, 2)):
        aa, bb = (a[2 + s:52 + s], b[2:52])          # a's line k + s against b's line k
        C[i, j] = coh_p90(aa, bb)
a = cx(ma, c0, c1); b = cx(mb, c0, c1)
xm = w // 2 - 1024
ifg = a[:, xm:xm + 2048] * np.conj(b[:, xm:xm + 2048])
lk = lambda z: z[: z.shape[0] // 6 * 6, : z.shape[1] // 24 * 24].reshape(z.shape[0] // 6, 6, -1, 24).sum(axis=(1, 3))
I = lk(ifg); pa = lk(np.abs(a[:, xm:xm + 2048]) ** 2); pb = lk(np.abs(b[:, xm:xm + 2048]) ** 2)
buf = io.BytesIO()
np.savez_compressed(buf, starts=starts, C=C, ph=np.angle(I).astype(np.float32),
                    coh=(np.abs(I) / np.sqrt(pa * pb + 1e-30)).astype(np.float32), c0=c0, c1=c1, xm=xm)
sys.stdout.buffer.write(buf.getvalue())
