"""Why T071's pass 1 set 29 of 427 dates aside (2026-11-06). 25 of them fail only _content_blocks (the amplitude pattern of
one 50-line block matches neither neighbour while the neighbours match each other), T173 none. Read only, on the pass-1
images (OUT/T/superseded_nesd1/raw/F<n>/). For every flagged (date, subswath, block), with the two neighbours qc used:
  amplitude correlation with each neighbour and between them, as qc computes it;
  coherence against the nearer neighbour: p90 / p99 and the share of windows above 0.5, at shift 0 and the noise floor
  (the same with the neighbour QC_SHIFT lines off): coherent windows prove the block holds this ground;
  the same numbers for the date's other blocks in that subswath (median), as a control.
For one block per (subswath, block) also every date's amplitude correlation against the stack's median date, and the
amplitude strips of the worst cases for a figure (npz).
usage (seasat venv): python qc_content_check.py T OUT.npz"""
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, out = sys.argv[1], sys.argv[2]
tr = pt.Track(T)
RAW = tr.out / "superseded_nesd1/raw"
s = tr.state
q, bad = s["qc_nesd1"], s["bad_align_nesd1"]
recut = sorted(s.get("recut", {}))
order = sorted(q, key=lambda d: (d not in recut, d))   # qc's own order: re-cut dates first, then time


def files(d, n):
    if d == tr.master:
        p = tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
        p = p if p.exists() else tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
        return p, p.with_suffix(".PRM")
    p = RAW / f"F{n}/S1_{d}_ALL_F{n}.SLC"
    return p, p.with_suffix(".PRM")


def neighbours(d):
    """As qc: the nearest dates before / after d that had already passed when d was judged."""
    i = order.index(d)
    ref = sorted(x for x in order[:i] if q[x].get("ok"))
    before, after = [x for x in ref if x < d], [x for x in ref if x > d]
    pool = sorted(set(q) | {tr.master})
    if not before:
        before = [x for x in pool if x < d]
    if not after:
        after = [x for x in pool if x > d]
    return ([before[-1]] if before else []) + ([after[0]] if after else [])


def blocks(d, n, cache={}):
    k = (d, n)
    if k not in cache:
        if len(cache) > 40:
            cache.clear()
        cache[k] = pt._blocks(*files(d, n))
    return cache[k]


def strip(a, ly=5, lx=20):
    n, m = len(a) // ly * ly, a.shape[1] // lx * lx
    return np.abs(a[:n, :m]).reshape(n // ly, ly, m // lx, lx).mean(axis=(1, 3))


def coh_stats(a, b):
    c = pt._windows(a, b)[2]; c = c[np.isfinite(c)]
    f = pt._windows(a[:-pt.QC_SHIFT], b[pt.QC_SHIFT:])[2]; f = f[np.isfinite(f)]
    pc = lambda x, p: float(np.percentile(x, p)) if x.size else np.nan
    return {"p90": pc(c, 90), "p99": pc(c, 99), "gt05": float(np.mean(c > 0.5)) if c.size else np.nan,
            "f90": pc(f, 90), "f99": pc(f, 99), "fgt05": float(np.mean(f > 0.5)) if f.size else np.nan}


rows, figs = [], {}
print(f"{T}: {len(q)} dates judged, {len(bad)} set aside; re-cut (judged first) {recut}", flush=True)
print("date       sw blk neighbours            amp a-b a-c b-c | coh p90 p99 >0.5  floor p90 p99 >0.5 | other blocks p90 p99 >0.5")
flagged = [(d, n, k) for d in sorted(bad) for n in tr.subs for k in (q[d].get(f"F{n}") or {}).get("content", [])]
for d, n, k in flagged:
    nb = neighbours(d)
    A = blocks(d, n)
    B = [blocks(x, n) for x in nb]
    ab, ac, bc = pt._amp_corr(A[k], B[0][k]), pt._amp_corr(A[k], B[1][k]), pt._amp_corr(B[0][k], B[1][k])
    c = coh_stats(A[k], B[0][k])
    oth = [coh_stats(A[j], B[0][j]) for j in range(len(A)) if j != k]
    om = {x: float(np.nanmedian([o[x] for o in oth])) for x in ("p90", "p99", "gt05")}
    rows.append({"d": d, "n": n, "k": k, "nb": nb, "ab": ab, "ac": ac, "bc": bc, **c, "other": om})
    print(f"{d} F{n} {k:3d} {nb[0]} {nb[1]}  {ab:+.2f} {ac:+.2f} {bc:+.2f} | {c['p90']:.3f} {c['p99']:.3f} {c['gt05']:.3f}  "
          f"{c['f90']:.3f} {c['f99']:.3f} {c['fgt05']:.3f} | {om['p90']:.3f} {om['p99']:.3f} {om['gt05']:.3f}", flush=True)
    key = f"{d}_F{n}_b{k}"
    figs[key + "_a"], figs[key + "_b"], figs[key + "_c"] = strip(A[k]), strip(B[0][k]), strip(B[1][k])
    figs[key + "_coh"] = pt._windows(A[k], B[0][k])[2]
# every date against the stack's median date, in each flagged (subswath, block)
dates = sorted(d for d in q)
ref = dates[len(dates) // 2]
for n, k in sorted({(n, k) for _, n, k in flagged}):
    R = pt._blocks(*files(ref, n))[k]
    v = []
    for d in dates:
        a = pt._blocks(*files(d, n))[k]
        v.append(pt._amp_corr(a, R))
    v = np.array(v)
    fl = np.array([d in bad and k in (q[d].get(f"F{n}") or {}).get("content", []) for d in dates])
    print(f"F{n} block {k}, amplitude correlation with {ref}: all dates percentiles 5 / 25 / 50 / 75 / 95 "
          f"{np.round(np.nanpercentile(v, [5, 25, 50, 75, 95]), 2)}; flagged dates median {np.nanmedian(v[fl]):.2f}, "
          f"unflagged {np.nanmedian(v[~fl]):.2f}", flush=True)
    figs[f"stack_F{n}_b{k}"] = v
np.savez_compressed(out, rows=json.dumps(rows), dates=np.array(dates), **figs)
print("->", out)
