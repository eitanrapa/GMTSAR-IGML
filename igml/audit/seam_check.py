"""Where merge_swath put each subswath, and where the subswaths really overlap (T173 super-master pair).
Grids streamed with grd2xyz -ZTLf (top row = last azimuth line first); flipped so row 0 = first line."""
import sys
import numpy as np

D = sys.argv[1]
shape = {"F1": (4715, 666), "F2": (4714, 789), "F3": (4938, 762), "M": (4939, 2129)}


def load(name, g):
    a = np.fromfile(f"{D}/{name}_{g}.bin", dtype=np.float32).reshape(shape[name])
    return a[::-1]          # row 0 = first azimuth line (y = 0)


# merge_swath's own geometry (from the PRMs printed in the run)
dr = 299792458.0 / 2 / 64345238.125714
near = {"F1": 799876.347427, "F2": 846316.167874, "F3": 901899.519560}
clock = {"F1": 279.572009173409, "F2": 279.572020212499, "F3": 279.571999419028}
prf = 486.486310
d12 = (near["F2"] - near["F1"]) / dr / 32
d23 = (near["F3"] - near["F2"]) / dr / 32
print(f"geometric column offsets: d12 {d12:.3f}  d23 {d23:.3f}  d13 {d12 + d23:.3f}")
h = {k: (clock[k] - clock["F1"]) * 86400 * prf / 6 for k in clock}
print("geometric row offsets vs F1:", {k: round(v, 3) for k, v in h.items()})

for g in ("phasefilt", "corr"):
    M = load("M", g)
    for name in ("F1", "F2", "F3"):
        S = load(name, g)
        ns, nc = S.shape
        # find the (row, col) placement of S inside M by exact equality of a central block
        r0, c0 = ns // 2, nc // 2
        blk = S[r0:r0 + 20, c0:c0 + 20]
        hits = []
        for dy in range(-200, 200):
            for dx in range(-50, 1450):
                rr, cc = r0 + dy, c0 + dx
                if 0 <= rr and rr + 20 <= M.shape[0] and 0 <= cc and cc + 20 <= M.shape[1]:
                    if np.array_equal(M[rr:rr + 20, cc:cc + 20], blk, equal_nan=True):
                        hits.append((dy, dx))
        print(f"{g} {name}: placed in merged grid at row offset / col offset {hits}")


def best_shift(A, B, rows, cols_a, cols_b, shifts, rshifts, use_phase):
    """corr of A[r, ca] with B[r + dy, cb + s] over the overlap, for column shifts s and row shifts dy."""
    out = {}
    for dy in rshifts:
        for s in shifts:
            a = A[rows][:, cols_a]
            b = B[rows + dy][:, cols_b + s]
            ok = np.isfinite(a) & np.isfinite(b)
            if use_phase:
                z = np.exp(1j * a[ok]) * np.exp(-1j * b[ok])
                out[(dy, s)] = abs(z.mean())
            else:
                x, y = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
                out[(dy, s)] = float((x * y).mean() / x.std() / y.std())
    return out


# TRUE overlap: F1 col (623 + c) vs F2 col c ; F2 col (746 + c) vs F3 col c.  Rows: F2 row = F1 row - 77 (head2 77.33),
# F3 row = F1 row + 68 (head3 -68.33); F3 row = F2 row + 145 (145.67)
for g, ph in (("phasefilt", True), ("corr", False)):
    F1, F2, F3 = load("F1", g), load("F2", g), load("F3", g)
    rows = np.arange(400, 4200)
    cb = np.arange(14, 30)                    # inside F2's valid near-range part, inside F1's far-range valid part
    res = best_shift(F1, F2, rows, cb + 623, cb, shifts=range(-3, 4), rshifts=range(-79, -74), use_phase=ph)
    best = max(res, key=res.get)
    print(f"{g} F1 vs F2 (F2 row = F1 row + dy, F2 col = F1 col - 623 + s): best (dy, s) {best} {res[best]:.4f}")
    for s in range(-3, 4):
        print("   s", s, " ".join(f"{res[(dy, s)]:.4f}" for dy in range(-79, -74)))
    rows2 = np.arange(400, 4200)
    res = best_shift(F2, F3, rows2, cb + 746, cb, shifts=range(-3, 4), rshifts=range(143, 149), use_phase=ph)
    best = max(res, key=res.get)
    print(f"{g} F2 vs F3 (F3 row = F2 row + dy, F3 col = F2 col - 746 + s): best (dy, s) {best} {res[best]:.4f}")
    for s in range(-3, 4):
        print("   s", s, " ".join(f"{res[(dy, s)]:.4f}" for dy in range(143, 149)))
