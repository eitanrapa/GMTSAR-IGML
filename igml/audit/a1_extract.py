# Windows around burst boundaries of two aligned GMTSAR SLCs (int16 I/Q): each date's line-to-line phase (its Doppler
# centroid, which jumps where that date switches bursts), the pair's coherence per line, and an image cut. npz -> stdout.
import io, sys
import numpy as np

def prm(p):
    return {a.strip(): b.strip() for a, s, b in (l.partition("=") for l in open(p)) if s}

A, B = sys.argv[1], sys.argv[2]
centres = [int(c) for c in sys.argv[3].split(",")]
H = int(sys.argv[4]) if len(sys.argv) > 4 else 80
P = prm(A[:-4] + ".PRM"); w, nl = int(P["num_rng_bins"]), int(P["num_lines"])
prf = float(P["PRF"])
ma = np.memmap(A, dtype="<i2", mode="r", shape=(nl, w, 2)); mb = np.memmap(B, dtype="<i2", mode="r", shape=(nl, w, 2))
out = {"prf": prf, "width": w, "centres": np.array(centres), "H": H}
for k, c in enumerate(centres):
    r0, r1 = c - H, c + H
    a = ma[r0:r1].astype(np.float32); a = a[..., 0] + 1j * a[..., 1]
    b = mb[r0:r1].astype(np.float32); b = b[..., 0] + 1j * b[..., 1]
    ok = (np.abs(a) > 0) & (np.abs(b) > 0)
    cols = np.where(ok.mean(axis=0) > 0.99)[0]
    x0, x1 = cols.min(), cols.max() + 1                     # the range span both cover in the whole window
    a, b, ok = a[:, x0:x1], b[:, x0:x1], ok[:, x0:x1]
    # each date's line-to-line phase, summed over range
    da = np.angle((a[1:] * np.conj(a[:-1])).sum(axis=1)); db = np.angle((b[1:] * np.conj(b[:-1])).sum(axis=1))
    # the pair's coherence per line: 3 lines x 32 samples, median over range
    nx = (x1 - x0) // 32 * 32
    ifg = a[:, :nx] * np.conj(b[:, :nx]); pa = np.abs(a[:, :nx]) ** 2; pb = np.abs(b[:, :nx]) ** 2
    s = lambda z: z.reshape(z.shape[0], -1, 32).sum(axis=2)
    I, PA, PB = s(ifg), s(pa), s(pb)
    I3 = I[:-2] + I[1:-1] + I[2:]; PA3 = PA[:-2] + PA[1:-1] + PA[2:]; PB3 = PB[:-2] + PB[1:-1] + PB[2:]
    coh_line = np.median(np.abs(I3) / np.sqrt(PA3 * PB3 + 1e-30), axis=1)
    # image cut: 2048 samples at mid-swath, dominant fringe removed, 2 x 8 looks for the phase, 3 x 16 coherence
    xm = (x1 - x0) // 2 - 1024
    cut = ifg[:, xm:xm + 2048]; ca = pa[:, xm:xm + 2048]; cb = pb[:, xm:xm + 2048]
    F = np.fft.fft2(cut); ky, kx = np.unravel_index(np.argmax(np.abs(F)), F.shape)
    yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
    cut = cut * np.exp(-2j * np.pi * (ky * yy / cut.shape[0] + kx * xx / cut.shape[1]))
    lk = lambda z, ly, lx: z[: z.shape[0] // ly * ly, : z.shape[1] // lx * lx].reshape(z.shape[0] // ly, ly, -1, lx).sum(axis=(1, 3))
    ph = np.angle(lk(cut, 1, 8))
    cm = np.abs(lk(cut, 3, 16)) / np.sqrt(lk(ca, 3, 16) * lk(cb, 3, 16) + 1e-30)
    out.update({f"da{k}": da, f"db{k}": db, f"coh{k}": coh_line, f"ph{k}": ph.astype(np.float32),
                f"cm{k}": cm.astype(np.float32), f"x{k}": np.array([x0, x1, xm])})
    print(f"window {k}: lines {r0}-{r1}, range {x0}-{x1}", file=sys.stderr)
buf = io.BytesIO(); np.savez_compressed(buf, **out); sys.stdout.buffer.write(buf.getvalue())
