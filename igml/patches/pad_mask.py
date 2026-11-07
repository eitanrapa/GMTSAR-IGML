"""Extend a merged-frame mask grid (the lab's landmask_ra.grd / mask_def.grd) to our merged grid's size (2026-10-29).
With every subswath placed exactly (patched merge_swath) a merged grid gains up to a column at far range and a row at
the end; its nodes are the lab's lattice (origin 0, 32 x 6), so the mask is padded by copying its last column / row
(the new cells only; every existing value, NaN included, stays). Checked: OUT cut back to SRC's region == SRC.

usage (seasat venv, GMT on PATH): python pad_mask.py SRC.grd REF.grd OUT.grd
"""
import subprocess
import sys
from pathlib import Path

import numpy as np


def info(g):
    f = subprocess.run(["gmt", "grdinfo", "-C", g], capture_output=True, text=True, check=True).stdout.split()
    x0, x1, y0, y1 = map(float, f[1:5]); dx, dy = float(f[7]), float(f[8]); nx, ny = int(f[9]), int(f[10])
    return x0, x1, y0, y1, dx, dy, nx, ny


def read_bf(g, tmp):
    subprocess.run(["gmt", "grdconvert", g, f"{tmp}=bf"], check=True)
    raw = np.fromfile(tmp, dtype=np.float32, offset=892)
    return raw


if __name__ == "__main__":
    src, ref, out = sys.argv[1:4]
    s, r = info(src), info(ref)
    if s[4:6] != r[4:6]:                       # landmask_ra.grd is 4 x 8: onto the pairs' 32 x 6, as snaphu_interp.csh
        rs = str(Path(out).with_suffix(".rs.grd"))   # does per pair (grdsample to phase_patch's increment, own region)
        subprocess.run(["gmt", "grdsample", src, f"-I{r[4]:g}/{r[5]:g}", f"-G{rs}"], check=True)
        src, s = rs, info(rs)
    assert s[0] == r[0] == 0 and s[2] == r[2] == 0 and s[4:6] == r[4:6], (s, r)
    k, m = r[6] - s[6], r[7] - s[7]
    assert 0 <= k <= 2 and 0 <= m <= 2, (k, m)
    t = Path(out).with_suffix(".tmp")
    a = read_bf(src, str(t)).reshape(s[7], s[6])                 # top row (y max) first
    b = np.pad(a, ((m, 0), (0, k)), mode="edge")                 # new rows at the top (y end), new columns at the right
    subprocess.run(["gmt", "grdmath", f"-R0/{r[1]:g}/0/{r[3]:g}", f"-I{r[4]:g}/{r[5]:g}", "-r", "0", "=", f"{t}=bf"],
                   check=True)
    hdr = open(t, "rb").read(892)
    assert Path(t).stat().st_size == 892 + b.size * 4, (Path(t).stat().st_size, b.shape)
    with open(t, "wb") as f:
        f.write(hdr); f.write(b.astype(np.float32).tobytes())
    subprocess.run(["gmt", "grdconvert", f"{t}=bf", out], check=True)
    t.unlink()
    if src.endswith(".rs.grd"):
        Path(src).unlink()
    # check: OUT on SRC's region equals SRC, NaN for NaN
    c = Path(out).with_suffix(".chk")
    subprocess.run(["gmt", "grdcut", out, f"-R0/{s[1]:g}/0/{s[3]:g}", f"-G{c}=bf"], check=True)
    back = np.fromfile(c, dtype=np.float32, offset=892); c.unlink()
    same = np.array_equal(back, a.ravel(), equal_nan=True)
    print(f"{out}: {s[6]} x {s[7]} -> {r[6]} x {r[7]} (+{k} columns, +{m} rows); original cells unchanged: {same}")
    assert same
