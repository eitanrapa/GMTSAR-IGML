"""The rebuild's first unwrapped pairs, checked (2026-11-02; user: "When will first unwrapped interferograms be made, to
check they are good"), on seasat, read only:
  1. each sampled stored pair against the lab's same pair: the share of pixels where the two unwrapped phases differ by
     a whole cycle (|round((ours - lab - median) / 2 pi)| >= 1, on pixels both have; same nominal grid, the lab's
     content sits 0.5-1.5 cells off -- AUDIT 11-02 -- which moves no pixel a cycle except on steep fringes);
  2. closure of sampled triplets (a, b) + (b, c) - (a, c) among the stored pairs, ours and, where it has all three, the
     lab's: the share of pixels with |closure - median| > pi (a whole cycle wrong somewhere in the triangle).
Also writes, for FIGPAIR and a triplet through it, decimated arrays for the figure (unw_check_fig.py).
usage (Mac, igml env, PYTHONPATH=IGML; seasat's venv has no h5py): python unw_check.py T NPAIRS NTRIP FIGPAIR OUT.npz
  ROOT=/seasat (default: the Mac's automount of seasat:/zdata/datasets)"""
import os
import re
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

T, npairs, ntrip, figpair, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
LABS = {"T173": "seasatHD3", "T071": "seasatHD1", "T064": "seasatHD2", "T166": "seasatHD4"}
ROOT = os.environ.get("ROOT", "/seasat")
OURS = Path(f"{ROOT}/IGML/InSAR_Processing/{T}/product_{T}/merge")
LAB = Path(f"{ROOT}/{LABS[T]}/InSAR_Processing/{T}/product_{T}/merge")
TWO_PI = 2 * np.pi


@lru_cache(maxsize=8)
def unw(root, pair):
    f = Path(root) / pair / "unwrap.grd"
    if not f.exists():
        return None
    z = read_grd(f).z.astype(np.float32)
    z[z == 0] = np.nan                        # GMTSAR leaves 0 where it did not unwrap
    return z


def common(*zs):
    r = min(z.shape[0] for z in zs); c = min(z.shape[1] for z in zs)
    return [z[:r, :c] for z in zs]


ours = sorted(p.name for p in OURS.iterdir() if re.match(r"^\d{7}_\d{7}$", p.name) and (p / "unwrap.grd").exists())
print(f"{T}: {len(ours)} stored pairs", flush=True)
# 1. against the lab
both = [p for p in ours if (LAB / p / "unwrap.grd").exists()]
pick = [both[int(i)] for i in np.linspace(0, len(both) - 1, min(npairs, len(both)))]
rows = []
for p in pick:
    a, b = common(unw(str(OURS), p), unw(str(LAB), p))
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 1000:
        continue
    d = (a - b)[m]
    k = np.rint((d - np.median(d)) / TWO_PI)
    rows.append((p, float(np.mean(k != 0) * 100), float(m.mean() * 100), float(np.mean(np.isfinite(a)) * 100),
                 float(np.mean(np.isfinite(b)) * 100)))
    print(f"  {p}: {rows[-1][1]:6.3f} % of common pixels a cycle apart (common {rows[-1][2]:.1f} %, ours "
          f"{rows[-1][3]:.1f} %, lab {rows[-1][4]:.1f} % unwrapped)", flush=True)
v = np.array([r[1] for r in rows])
print(f"ours vs lab, {len(rows)} pairs: pixels a whole cycle apart median {np.median(v):.3f} %, p90 "
      f"{np.percentile(v, 90):.3f} %, max {v.max():.3f} % ({rows[int(np.argmax(v))][0]})", flush=True)
# 2. triplets
S = set(ours)
by_a = {}
for p in ours:
    a, b = p.split("_")
    by_a.setdefault(a, []).append(b)
trip = [(a, b, c) for a in by_a for b in by_a[a] for c in by_a.get(b, []) if f"{a}_{c}" in S]
trip.sort()
sel = [trip[int(i)] for i in np.linspace(0, len(trip) - 1, min(ntrip, len(trip)))]
print(f"{len(trip)} triplets among the stored pairs; {len(sel)} sampled", flush=True)


def closure(root, a, b, c):
    zs = [unw(root, f"{a}_{b}"), unw(root, f"{b}_{c}"), unw(root, f"{a}_{c}")]
    if any(z is None for z in zs):
        return None, None
    ab, bc, ac = common(*zs)
    cl = ab + bc - ac
    m = np.isfinite(cl)
    if m.sum() < 1000:
        return None, None
    cl = cl - np.median(cl[m])
    return float(np.mean(np.abs(cl[m]) > np.pi) * 100), cl


res = []
for a, b, c in sel:
    fo, _ = closure(str(OURS), a, b, c)
    fl, _ = closure(str(LAB), a, b, c)
    res.append((f"{a}_{b}_{c}", fo, fl))
    print(f"  {a} {b} {c}: closure a cycle off, ours {fo if fo is None else round(fo, 3)} %, lab "
          f"{fl if fl is None else round(fl, 3)} %", flush=True)
fo = np.array([r[1] for r in res if r[1] is not None])
pl = [(r[1], r[2]) for r in res if r[1] is not None and r[2] is not None]
print(f"closure, ours ({fo.size} triplets): median {np.median(fo):.3f} %, p90 {np.percentile(fo, 90):.3f} %, "
      f"max {fo.max():.3f} %", flush=True)
if pl:
    pl = np.array(pl)
    print(f"closure, the {len(pl)} triplets the lab also has: ours median {np.median(pl[:, 0]):.3f} % (p90 "
          f"{np.percentile(pl[:, 0], 90):.3f}), lab median {np.median(pl[:, 1]):.3f} % (p90 {np.percentile(pl[:, 1], 90):.3f}); "
          f"ours lower on {int(np.sum(pl[:, 0] < pl[:, 1]))}, higher on {int(np.sum(pl[:, 0] > pl[:, 1]))}", flush=True)
# 3. the figure's arrays (every 4th row / column)
a, b = figpair.split("_")
fig = {}
o, l = common(unw(str(OURS), figpair), unw(str(LAB), figpair))
fig["ours"], fig["lab"] = o[::4, ::4], l[::4, ::4]
cs = [c for c in by_a.get(b, []) if f"{a}_{c}" in S]
if cs:
    c = cs[0]
    _, clo = closure(str(OURS), a, b, c)
    _, cll = closure(str(LAB), a, b, c)
    fig["trip"] = np.array([a, b, c])
    if clo is not None:
        fig["clo_ours"] = clo[::4, ::4]
    if cll is not None:
        fig["clo_lab"] = cll[::4, ::4]
np.savez_compressed(out, pairs=np.array([r[0] for r in rows]), cyc=v, trip_names=np.array([r[0] for r in res]),
                    trip_ours=np.array([np.nan if r[1] is None else r[1] for r in res]),
                    trip_lab=np.array([np.nan if r[2] is None else r[2] for r in res]), figpair=figpair, **fig)
print(f"-> {out}")
