"""Synthetic test of stitch_tops (2026-10-29): frames cut from one random scene, each frame's clock_start annotated
e_k lines late; a correct stitch reproduces the scene line for line. For each case and each binary, prints the
offset (output line - true scene line) on each stitched piece: 0 everywhere = correct.

usage: python stitch_synth.py BIN [BIN ...]
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

PRF, W, N, STEP = 486.486310, 64, 1500, 1000                # frames of 1,500 lines every 1,000: 500-line overlaps
rng = np.random.default_rng(0)
SCENE = (rng.normal(0, 800, (STEP * 4 + N, W, 2))).astype("<i2")
KEY = {SCENE[i, :8].tobytes(): i for i in range(SCENE.shape[0])}


def prm(path, a, e):
    t = 100.0 + (a + e) / PRF / 86400.0
    path.write_text(f"SLC_file = {path.with_suffix('.SLC').name}\nPRF = {PRF}\nrng_samp_rate = 64345238.125714\n"
                    f"num_rng_bins = {W}\nbytes_per_line = {4 * W}\ngood_bytes_per_line = {4 * W}\nnum_lines = {N}\n"
                    f"nrows = {N}\nnum_valid_az = {N}\nnear_range = 800000.0\nclock_start = {t:.12f}\n"
                    f"clock_stop = {t + N / PRF / 86400:.12f}\nSC_clock_start = {2026000 + t:.12f}\n"
                    f"SC_clock_stop = {2026000 + t + N / PRF / 86400:.12f}\nashift = 0\nsub_int_a = 0\n")


def run(binary, errors):
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        stems = []
        for k, e in enumerate(errors):
            s = d / f"frame{k}"
            SCENE[k * STEP:k * STEP + N].tofile(s.with_suffix(".SLC"))
            prm(s.with_suffix(".PRM"), k * STEP, e)
            stems.append(s.name)
        (d / "list").write_text("\n".join(stems) + "\n")
        r = subprocess.run([binary, "list", "out"], cwd=d, capture_output=True, text=True)   # main() returns 1 on success
        if not (d / "out.SLC").exists():
            raise RuntimeError(r.stdout + r.stderr)
        out = np.fromfile(d / "out.SLC", dtype="<i2").reshape(-1, W, 2)
    off = np.array([i - KEY.get(out[i, :8].tobytes(), -10**6) for i in range(out.shape[0])])
    pieces, start = [], 0                                         # runs of constant offset
    for i in range(1, off.size + 1):
        if i == off.size or off[i] != off[start]:
            pieces.append(f"lines {start}-{i - 1}: {off[start]:+d}")
            start = i
    return pieces


CASES = {"2 frames, frame 1 annotated 1 line late": (0, 1),
         "3 frames, frame 2 (the last) 1 line late": (0, 0, 1),
         "3 frames, frame 1 (the middle) 1 line late": (0, 1, 0),
         "3 frames, frames 1 and 2 both 1 line late": (0, 1, 1)}
for name, errors in CASES.items():
    print(f"== {name} {errors}")
    for b in sys.argv[1:]:
        print(f"   {Path(b).parent.name:6s}", "; ".join(run(b, errors)))
