"""The unwrapped phase of given pairs, three ways (2026-11-06): the lab's unwrap.grd, the lab's unwrap_notropo.grd (SNAPHU's output,
per notes/gapfill_plan.md) and ours (SNAPHU's output on the fixed interferogram). Read only, over the Mac's mount.
For each: the difference's spread after its median, in rad and mm LOS, and the share of pixels beyond half a cycle; is the lab's
unwrap - notropo what ours - lab is?
usage (insar_sim env, PYTHONPATH=IGML): python pair_unwrap_tropo.py T PAIR [PAIR ...]"""
import sys
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

LABS = {"T173": "seasatHD3"}
T, pairs = sys.argv[1], sys.argv[2:]
ARCH = Path(f"/seasat/IGML/lab_originals_20261029/{T}/product_{T}/merge")
OURS = Path(f"/seasat/{LABS[T]}/InSAR_Processing/{T}/product_{T}/merge")
MM = 55.465763 / (4 * np.pi)


def load(f):
    g = read_grd(f)
    z = g.z.astype(np.float64)
    z[z == 0] = np.nan
    return z, g


def stats(d, name):
    m = np.isfinite(d)
    v = d[m] - np.median(d[m])
    q = np.percentile(v, [1, 16, 84, 99])
    print(f"   {name}: {m.sum()} px; robust sd {(q[2] - q[1]) / 2:.3f} rad = {(q[2] - q[1]) / 2 * MM:.1f} mm; p1..p99 {q[0]:+.2f}..{q[3]:+.2f} rad = "
          f"{q[0] * MM:+.1f}..{q[3] * MM:+.1f} mm; beyond pi {np.mean(np.abs(v) > np.pi) * 100:.3f} %", flush=True)
    return v


for pair in pairs:
    print(f"=== {pair}", flush=True)
    L, gl = load(ARCH / pair / "unwrap.grd")
    O, go = load(OURS / pair / "unwrap.grd")
    f = ARCH / pair / "unwrap_notropo.grd"
    print(f"   lab unwrap {L.shape} x {gl.x[0]:g}..{gl.x[-1]:g} y {gl.y[0]:g}..{gl.y[-1]:g}; ours {O.shape}", flush=True)
    r, c = min(L.shape[0], O.shape[0]), min(L.shape[1], O.shape[1])
    stats(O[:r, :c] - L[:r, :c], "ours - lab unwrap.grd")
    if f.exists():
        N, gn = load(f)
        print(f"   lab notropo {N.shape} x {gn.x[0]:g}..{gn.x[-1]:g} y {gn.y[0]:g}..{gn.y[-1]:g}", flush=True)
        if N.shape == L.shape:
            a = stats(L - N, "lab unwrap.grd - unwrap_notropo.grd")
            b = stats(O[:r, :c] - N[:r, :c], "ours - lab unwrap_notropo.grd")
            d = np.isfinite(L - N) & np.isfinite(O[:r, :c] - L[:r, :c])
            print(f"   corr(lab unwrap - notropo, ours - lab unwrap) = {np.corrcoef((L - N)[d], (O[:r, :c] - L[:r, :c])[d])[0, 1]:+.3f}", flush=True)
        else:
            ys = np.searchsorted(gn.y, gl.y).clip(0, N.shape[0] - 1)
            xs = np.searchsorted(gn.x, gl.x).clip(0, N.shape[1] - 1)
            NN = N[np.ix_(ys, xs)]
            print(f"   (notropo on another grid: nearest node, {N.shape} -> {L.shape})", flush=True)
            a = stats(L - NN, "lab unwrap.grd - unwrap_notropo.grd")
            b = stats(O[:r, :c] - NN[:r, :c], "ours - lab unwrap_notropo.grd")
            d = np.isfinite(L - NN) & np.isfinite(O[:r, :c] - L[:r, :c])
            d[r:, :] = False
            A = (L - NN)[:r, :c]; B = (O[:r, :c] - L[:r, :c])
            dd = np.isfinite(A) & np.isfinite(B)
            print(f"   corr(lab unwrap - notropo, ours - lab unwrap) = {np.corrcoef(A[dd], B[dd])[0, 1]:+.3f}", flush=True)
    else:
        print("   no unwrap_notropo.grd", flush=True)
