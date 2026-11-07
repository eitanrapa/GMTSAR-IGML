"""How far apart two merged products of the SAME pair place the same content (2026-11-02): per subswath (column ranges
of the merged grid) the cross-correlation of the two corr.grd arrays (FFT, sub-cell peak by parabola), turned into the
offset of A's labels against B's in range samples / lines (e = x_A[0] - x_B[0] + k * dx), and into metres east / north
at a point through trans.dat's local Jacobian (as cr_geocode_any.py). With B the rebuild (validated at corner
reflectors), this is A's geolocation error -- the same pixels, so no reference shoreline and no edge bias in it.
usage (insar_sim env, PYTHONPATH=IGML): python relshift.py A_corr.grd B_corr.grd TRANS LON LAT c1,c2 [c1,c2 ...]"""
import sys

import numpy as np

from igml.real.grd import read_coords, read_grd

A, B, trans = read_grd(sys.argv[1]), read_grd(sys.argv[2]), sys.argv[3]
lon0, lat0 = float(sys.argv[4]), float(sys.argv[5])
ranges = [tuple(int(v) for v in r.split(",")) for r in sys.argv[6:]]
dx, dy = A.x[1] - A.x[0], A.y[1] - A.y[0]
assert np.isclose(dx, B.x[1] - B.x[0]) and np.isclose(dy, B.y[1] - B.y[0]), "cell sizes differ"


def jacobian(rad=0.02):
    t = np.memmap(trans, dtype="<f8", mode="r").reshape(-1, 5)
    out = []
    for s in range(0, t.shape[0], 5_000_000):
        b = np.asarray(t[s:s + 5_000_000])
        out.append(b[(np.abs(b[:, 3] - lon0) < rad) & (np.abs(b[:, 4] - lat0) < rad)])
    b = np.concatenate(out)
    e = (b[:, 3] - lon0) * 111320.0 * np.cos(np.radians(lat0))
    n = (b[:, 4] - lat0) * 110950.0
    M = np.c_[e, n, np.ones_like(e)]
    cr, ca = np.linalg.lstsq(M, b[:, 0], rcond=None)[0], np.linalg.lstsq(M, b[:, 1], rcond=None)[0]
    return np.array([[cr[0], cr[1]], [ca[0], ca[1]]])


J = jacobian()
Ji = np.linalg.inv(J)
print(f"cells {dx:g} samples x {dy:g} lines; A origin ({A.x[0]:g}, {A.y[0]:g}), B origin ({B.x[0]:g}, {B.y[0]:g}); "
      f"A {A.z.shape}, B {B.z.shape}")
for c1, c2 in ranges:
    rows = min(A.z.shape[0], B.z.shape[0])
    a = A.z[:rows, c1:c2].astype(np.float64)
    b = B.z[:rows, c1:c2].astype(np.float64)
    m = np.isfinite(a) & np.isfinite(b)
    a = np.where(m, a - a[m].mean(), 0.0)
    b = np.where(m, b - b[m].mean(), 0.0)
    F = np.fft.rfft2(a) * np.conj(np.fft.rfft2(b))
    cc = np.fft.irfft2(F, s=a.shape)
    cc = np.fft.fftshift(cc)
    i, j = np.unravel_index(np.argmax(cc), cc.shape)

    def par(m1, c0, p1):
        d = m1 - 2 * c0 + p1
        return 0.0 if d == 0 else 0.5 * (m1 - p1) / d
    fi = i + par(cc[i - 1, j], cc[i, j], cc[i + 1, j]) - a.shape[0] // 2
    fj = j + par(cc[i, j - 1], cc[i, j], cc[i, j + 1]) - a.shape[1] // 2
    # cc peaks at (fi, fj) where a[r, c] ~ b[r - fi, c - fj]: A's content sits fj columns / fi rows later than B's
    er = A.x[0] - B.x[0] + fj * dx
    ea = A.y[0] - B.y[0] + fi * dy
    en = Ji @ np.array([er, ea])
    print(f"columns {c1}-{c2}: A - B = {fj:+.3f} cells / {fi:+.3f} rows -> A's labels {er:+.1f} samples, {ea:+.2f} lines off "
          f"B's = east {en[0]:+.0f} m, north {en[1]:+.0f} m ({np.hypot(*en):.0f} m) at ({lon0}, {lat0}); peak {cc[i, j] / np.sqrt((a * a).sum() * (b * b).sum()):.3f}")
