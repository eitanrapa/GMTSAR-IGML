"""In the F1/F2 and F2/F3 overlaps of the T173 super-master pair: phase of the two subswaths at the same ground (true
registration) and at the merged registration (F2 one column off), as circular mean / spread, in rad and mm LOS."""
import sys
import numpy as np

D = sys.argv[1]
shape = {"F1": (4715, 666), "F2": (4714, 789), "F3": (4938, 762)}
ld = lambda n: np.fromfile(f"{D}/{n}_phasefilt.bin", dtype=np.float32).reshape(shape[n])[::-1]
F1, F2, F3 = ld("F1"), ld("F2"), ld("F3")
mm = 0.0554658 / (4 * np.pi) * 1e3


def stats(a, b, tag):
    ok = np.isfinite(a) & np.isfinite(b)
    z = np.exp(1j * (a[ok] - b[ok]))
    m = np.angle(z.mean())
    d = np.angle(z * np.exp(-1j * m))
    print(f"{tag}: n {ok.sum()}  circ mean {m:+.4f} rad ({m * mm:+.2f} mm)  rms about mean {np.sqrt((d ** 2).mean()):.4f} rad"
          f" ({np.sqrt((d ** 2).mean()) * mm:.2f} mm)  |coh| {abs(z.mean()):.4f}")


rows = np.arange(400, 4200)
cs = np.arange(14, 30)
# F1 row R <-> F2 row R - 77 ; F2 col c <-> F1 col c + 623 (measured)
a = F1[rows][:, cs + 623]
stats(a, F2[rows - 77][:, cs], "F1 vs F2, true registration       ")
stats(a, F2[rows - 77][:, cs - 1], "F1 vs F2, merged registration (+1)")
# F2 row R <-> F3 row R + 146 ; F3 col c <-> F2 col c + 746
a = F2[rows][:, cs + 746]
stats(a, F3[rows + 146][:, cs], "F2 vs F3, true registration       ")
stats(a, F3[rows + 146][:, cs + 1], "F2 vs F3, merged registration     ")
# per-azimuth-block mean offsets (burst structure)
a = F1[rows][:, cs + 623]
b = F2[rows - 77][:, cs]
blk = []
for k in range(0, len(rows), 200):
    x, y = a[k:k + 200], b[k:k + 200]
    ok = np.isfinite(x) & np.isfinite(y)
    blk.append(np.angle(np.exp(1j * (x[ok] - y[ok])).mean()))
print("F1-F2 true-reg block means (rad):", " ".join(f"{v:+.3f}" for v in blk))
