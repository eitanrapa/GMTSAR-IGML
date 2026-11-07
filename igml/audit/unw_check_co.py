"""The rebuild's unwrapped pairs across an earthquake, and one batch's others, checked (2026-11-02; user: "Okay, let's
check" -- batch 3 and the pairs spanning Ridgecrest, unwrapped in SNAPHU's DEFO mode, defomax 60). Read only, over the
Mac's mount of seasat (seasat's venv has no h5py). As unw_check.py:
  1. every stored pair spanning the track's earthquake (process_track.unwrap_cmd's rule: a <= the event's last day and
     its first day <= b) against the lab's same pair -- the share of pixels a whole cycle apart, over the grid and in its
     north-west quarter (rows < 1/4, columns > 1/2 of the merged grid: descending, so north / west), nearest Ridgecrest;
  2. closure of triplets with at least one leg across the earthquake, ours and the lab's;
  3. the same two checks on a sample of one batch's other pairs (second dates B0..B1).
Writes an npz (+ decimated arrays of the shortest coseismic pair and one coseismic triplet's closure).
usage (Mac, igml env, PYTHONPATH=IGML): python unw_check_co.py T B0 B1 NPAIRS NTRIP OUT.npz   (B0 / B1: YYYYMMDD)"""
import json
import os
import re
import sys
from datetime import date, timedelta
from functools import lru_cache
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

T, B0, B1, npairs, ntrip, out = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
LABS = {"T173": "seasatHD3", "T071": "seasatHD1", "T064": "seasatHD2", "T166": "seasatHD4"}
ROOT = os.environ.get("ROOT", "/seasat")
OURS = Path(f"{ROOT}/IGML/InSAR_Processing/{T}/product_{T}/merge")
LAB = Path(f"{ROOT}/{LABS[T]}/InSAR_Processing/{T}/product_{T}/merge")
EV = json.load(open(Path(__file__).resolve().parents[2] / "coseismic_events.json"))["tracks"].get(T, [])
TWO_PI = 2 * np.pi


def ymd(d):                                     # GMTSAR's yyyyddd (day 0 = 1 January) -> YYYYMMDD
    return (date(int(d[:4]), 1, 1) + timedelta(days=int(d[4:]))).strftime("%Y%m%d")


def coseismic(p):
    a, b = (ymd(x) for x in p.split("_"))
    return any(a <= e.get("end_date", e["date"]) and e["date"] <= b for e in EV)


@lru_cache(maxsize=8)
def unw(root, pair):
    f = Path(root) / pair / "unwrap.grd"
    if not f.exists():
        return None
    z = read_grd(f).z.astype(np.float32)
    z[z == 0] = np.nan
    return z


def common(*zs):
    r = min(z.shape[0] for z in zs); c = min(z.shape[1] for z in zs)
    return [z[:r, :c] for z in zs]


def vs_lab(p):
    o, l = unw(str(OURS), p), unw(str(LAB), p)
    if o is None or l is None:
        return None
    a, b = common(o, l)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 1000:
        return None
    d = np.where(m, a - b, np.nan)
    k = np.rint((d - np.nanmedian(d)) / TWO_PI)
    nw = np.zeros_like(m)
    nw[: m.shape[0] // 4, m.shape[1] // 2:] = True
    return float(np.mean(k[m] != 0) * 100), float(np.mean(k[m & nw] != 0) * 100) if (m & nw).sum() else np.nan


def closure(root, a, b, c, keep=False):
    zs = [unw(root, f"{a}_{b}"), unw(root, f"{b}_{c}"), unw(root, f"{a}_{c}")]
    if any(z is None for z in zs):
        return None, None
    ab, bc, ac = common(*zs)
    cl = ab + bc - ac
    m = np.isfinite(cl)
    if m.sum() < 1000:
        return None, None
    cl = cl - np.median(cl[m])
    return float(np.mean(np.abs(cl[m]) > np.pi) * 100), (cl if keep else None)


ours = sorted(p.name for p in OURS.iterdir() if re.match(r"^\d{7}_\d{7}$", p.name) and (p / "unwrap.grd").exists())
S = set(ours)
co = [p for p in ours if coseismic(p)]
print(f"{T}: {len(ours)} stored pairs, {len(co)} spanning {[e['date'] + '..' + e.get('end_date', e['date']) for e in EV]}",
      flush=True)
res = {}
# 1. the coseismic pairs against the lab
rows = []
for p in co:
    r = vs_lab(p)
    if r is None:
        print(f"  {p} ({ymd(p[:7])} x {ymd(p[8:])}): no lab pair", flush=True)
        continue
    rows.append((p, *r))
    print(f"  {p} ({ymd(p[:7])} x {ymd(p[8:])}): a cycle apart {r[0]:.3f} % of pixels, NW quarter {r[1]:.3f} %", flush=True)
v = np.array([r[1] for r in rows]); w = np.array([r[2] for r in rows])
print(f"coseismic pairs vs lab ({len(rows)}): median {np.median(v):.3f} % (max {v.max():.3f}), NW quarter median "
      f"{np.nanmedian(w):.3f} % (max {np.nanmax(w):.3f}, {rows[int(np.nanargmax(w))][0]})", flush=True)
res.update(co_pairs=np.array([r[0] for r in rows]), co_cyc=v, co_cyc_nw=w)
# 2. triplets with a coseismic leg
by_a = {}
for p in ours:
    a, b = p.split("_")
    by_a.setdefault(a, []).append(b)
trip = [(a, b, c) for a in by_a for b in by_a[a] for c in by_a.get(b, []) if f"{a}_{c}" in S]
tco = [t for t in trip if any(coseismic(f"{x}_{y}") for x, y in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])))]
sel = [tco[int(i)] for i in np.linspace(0, len(tco) - 1, min(ntrip, len(tco)))] if tco else []
print(f"{len(tco)} triplets with a coseismic leg; {len(sel)} sampled", flush=True)
tr = []
for a, b, c in sel:
    fo, _ = closure(str(OURS), a, b, c); fl, _ = closure(str(LAB), a, b, c)
    tr.append((f"{a}_{b}_{c}", fo, fl))
    print(f"  {a} {b} {c}: closure a cycle off, ours {fo if fo is None else round(fo, 3)} %, lab "
          f"{fl if fl is None else round(fl, 3)} %", flush=True)
res.update(cot_names=np.array([t[0] for t in tr]), cot_ours=np.array([np.nan if t[1] is None else t[1] for t in tr]),
           cot_lab=np.array([np.nan if t[2] is None else t[2] for t in tr]))
# 3. one batch's other pairs
bat = [p for p in ours if B0 <= ymd(p[8:]) <= B1 and not coseismic(p)]
pick = [bat[int(i)] for i in np.linspace(0, len(bat) - 1, min(npairs, len(bat)))]
br = []
for p in pick:
    r = vs_lab(p)
    if r is not None:
        br.append((p, *r))
bv = np.array([r[1] for r in br])
print(f"batch {B0}..{B1}: {len(bat)} non-coseismic pairs; vs lab ({len(br)} sampled with a lab pair): median "
      f"{np.median(bv):.3f} %, p90 {np.percentile(bv, 90):.3f} %, max {bv.max():.3f} % ({br[int(np.argmax(bv))][0]})",
      flush=True)
bt = [t for t in trip if B0 <= ymd(t[2]) <= B1 and t not in tco]
bsel = [bt[int(i)] for i in np.linspace(0, len(bt) - 1, min(ntrip, len(bt)))]
btr = []
for a, b, c in bsel:
    fo, _ = closure(str(OURS), a, b, c); fl, _ = closure(str(LAB), a, b, c)
    btr.append((f"{a}_{b}_{c}", fo, fl))
fo = np.array([t[1] for t in btr if t[1] is not None])
pl = np.array([(t[1], t[2]) for t in btr if t[1] is not None and t[2] is not None])
print(f"batch closure ({fo.size} triplets): ours median {np.median(fo):.3f} % (p90 {np.percentile(fo, 90):.3f}, max "
      f"{fo.max():.3f})" + (f"; with the lab ({len(pl)}): ours {np.median(pl[:, 0]):.3f} vs lab {np.median(pl[:, 1]):.3f} %, "
      f"ours lower on {int(np.sum(pl[:, 0] < pl[:, 1]))}" if len(pl) else ""), flush=True)
res.update(b_pairs=np.array([r[0] for r in br]), b_cyc=bv, bt_names=np.array([t[0] for t in btr]),
           bt_ours=np.array([np.nan if t[1] is None else t[1] for t in btr]),
           bt_lab=np.array([np.nan if t[2] is None else t[2] for t in btr]))
# figure arrays: the shortest coseismic pair with a lab pair, and the coseismic triplet with the worst lab closure
if rows:
    fp = min((r[0] for r in rows), key=lambda p: int(ymd(p[8:])) - int(ymd(p[:7])))
    o, l = common(unw(str(OURS), fp), unw(str(LAB), fp))
    res.update(fig_pair=fp, fig_ours=o[::4, ::4], fig_lab=l[::4, ::4])
ok = [t for t in tr if t[1] is not None and t[2] is not None]
if ok:
    tw = max(ok, key=lambda t: t[2])[0]
    a, b, c = tw.split("_")
    _, co_ = closure(str(OURS), a, b, c, keep=True); _, cl_ = closure(str(LAB), a, b, c, keep=True)
    res.update(fig_trip=tw, fig_clo_ours=co_[::4, ::4], fig_clo_lab=cl_[::4, ::4])
np.savez_compressed(out, **res)
print(f"-> {out}")
