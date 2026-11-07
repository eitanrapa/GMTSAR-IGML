"""A whole-cycle unwrapping error in the lab's OWN pairs, for a figure (2026-10-28, user: "Can I see one of the whole
cycle offsets in one of the lab's processed interferograms"). One lab triangle (i, k, j): the closure
u_ik + u_kj - u_ij in whole cycles; and per leg the shortest other triangle that contains that leg and neither of the
other two (the leg whose other triangle repeats a patch is the leg that carries it). Unwrapped = unwrap_notropo.grd
(SNAPHU's output; the lab's unwrap.grd has its troposphere correction subtracted). Saves the int8 cycle maps
(-128 = no data), the three legs' unwrapped / wrapped phase (float16) and coherence (x 250, uint8), and a coarse
radar -> lon / lat table from merge/trans.dat (range, azimuth, height, lon, lat). Read only.

usage (seasat): python cycle_example.py T i k j       (GMTSAR dates YYYYDDD, i < k < j)
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

T, i, k, j = sys.argv[1:5]
VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
lab = Path(f"/zdata/datasets/seasatHD{VOL[T]}/InSAR_Processing/{T}/product_{T}/merge")
PAT = re.compile(r"^\d{7}_\d{7}$")
net = {d.name for d in lab.iterdir() if PAT.match(d.name) and (d / "unwrap.grd").exists()}
dates = sorted({x for p in net for x in p.split("_")})
legs = {"A": f"{i}_{k}", "B": f"{k}_{j}", "C": f"{i}_{j}"}
SIGN = {"A": 1, "B": 1, "C": -1}                   # each leg's sign in the main closure
out = Path(f"/zdata/datasets/IGML/work/chk/cycle_example/{T}_{i}_{k}_{j}.npz")
out.parent.mkdir(parents=True, exist_ok=True)


def load(path):
    f = subprocess.run(["gmt", "grdinfo", "-C", str(path)], capture_output=True, text=True, check=True).stdout.split()
    nx, ny = int(f[9]), int(f[10])
    raw = subprocess.run(["gmt", "grd2xyz", str(path), "-ZTLf"], capture_output=True, check=True).stdout
    return np.frombuffer(raw, dtype=np.float32).reshape(ny, nx), [float(v) for v in f[1:7]]


cache = {}


def unw(p):
    if p not in cache:
        f = lab / p / "unwrap_notropo.grd"
        cache[p] = load(f if f.exists() else lab / p / "unwrap.grd")[0]
    return cache[p]


def cycles(x, y, z):
    c = unw(f"{x}_{y}") + unw(f"{y}_{z}") - unw(f"{x}_{z}")
    ok = np.isfinite(c)
    n = np.full(c.shape, -128, np.int8)
    n[ok] = np.clip(np.rint((c[ok] - np.median(c[ok])) / (2 * np.pi)), -127, 127)
    return n


def day(d):
    return int(d[:4]) * 365.25 + int(d[4:])


def other_triangle(p):
    a, b = p.split("_")
    best = None
    for m in dates:
        if m in (a, b):
            continue
        x, y, z = sorted((a, b, m))
        tri = {f"{x}_{y}", f"{y}_{z}", f"{x}_{z}"}
        if tri <= net and len(tri & set(legs.values())) == 1 and (best is None or day(z) - day(x) < best[0]):
            best = (day(z) - day(x), (x, y, z))
    return best[1] if best else None


res = {"n_main": cycles(i, k, j), "dates": np.array([i, k, j])}
for L, p in legs.items():
    tri = other_triangle(p)
    x, y, z = tri
    s = -1 if p == f"{x}_{z}" else 1
    res[f"n_{L}"] = cycles(*tri)
    res[f"tri_{L}"] = np.array(tri)
    res[f"ratio_{L}"] = np.int8(s * SIGN[L])        # in L's patch: n_L = ratio_L * n_main
    print(f"leg {L} {p}: other triangle {x}_{y} + {y}_{z} - {x}_{z}", flush=True)
for L, p in legs.items():
    res[f"unw_{L}"] = unw(p).astype(np.float16)
    res[f"phase_{L}"] = load(lab / p / "phasefilt.grd")[0].astype(np.float16)
    c = load(lab / p / "corr.grd")[0]
    res[f"corr_{L}"] = np.where(np.isfinite(c), np.clip(np.rint(c * 250), 0, 250), 255).astype(np.uint8)
res["header"] = np.array(load(lab / legs["C"] / "phasefilt.grd")[1])    # x0 x1 y0 y1 dx dy
tr = np.memmap(lab / "trans.dat", dtype=np.float64, mode="r").reshape(-1, 5)
res["trans"] = np.asarray(tr[::40][:, [0, 1, 3, 4]], dtype=np.float32)   # range, azimuth, lon, lat
np.savez_compressed(out, **res)
n = res["n_main"]
print(f"{out}: closure off by >= 1 cycle on {100 * np.mean(n[n != -128] != 0):.2f} % of pixels", flush=True)
