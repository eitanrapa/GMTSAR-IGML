"""Where along azimuth a lab date is misregistered (2026-10-28, after lab_shift_survey flagged T144's block 4 -- not the
last frame that stitch_tops' double offset shifts). Per NB blocks of 50 lines over the aligned SLC of date d and of
neighbour x (subswath n): p90 coherence unshifted and at -2..+2 lines, the best shift, and the line. Read only.

usage (seasat): python shift_profile.py T d x n [NB]
"""
import sys

sys.path.insert(0, "/zdata/datasets/IGML/code/scratchpad/s1_gapfill")
T, d, x, n = sys.argv[1:5]
NB = int(sys.argv[5]) if len(sys.argv) > 5 else 96
sys.argv = sys.argv[:1]
import numpy as np                                               # noqa: E402
import process_track as P                                        # noqa: E402

raw = P.Track(T).labp / f"F{n}/raw"


def blocks(date):
    f = raw / f"S1_{date}_ALL_F{n}.SLC"
    p = P._prm(f.with_suffix(".PRM"))
    w, nl = int(p["bytes_per_line"]) // 4, int(p["num_lines"])
    rows = np.linspace(0.005 * nl, 0.995 * nl - 50, NB).astype(int)
    out = []
    with open(f, "rb") as fh:
        for r in rows:
            fh.seek(int(r) * w * 4)
            a = np.fromfile(fh, dtype="<i2", count=50 * w * 2).astype(np.float32).reshape(-1, w, 2)
            out.append(a[..., 0] + 1j * a[..., 1])
    return rows, nl, out


p90 = lambda v: float(np.percentile(v[np.isfinite(v)], 90)) if np.isfinite(v).any() else np.nan
rows, nl, A = blocks(d)
_, _, B = blocks(x)
print(f"{T} F{n} {d} vs {x}: {nl} lines; line, p90 at shift -2 -1 0 +1 +2, best", flush=True)
for r, a, b in zip(rows, A, B):
    c = {s: p90(P._windows(*P._shifted(a, b, s))[2]) if s else p90(P._windows(a, b)[2]) for s in (-2, -1, 0, 1, 2)}
    best = max(c, key=lambda s: c[s] if np.isfinite(c[s]) else -1)
    flag = " <-- SHIFT" if best != 0 and c[best] > c[0] + P.QC_SHIFT_GAIN else ""
    print(f"{r:6d} {r / nl:5.3f} " + " ".join(f"{c[s]:.3f}" for s in (-2, -1, 0, 1, 2)) + f"  {best:+d}{flag}", flush=True)
