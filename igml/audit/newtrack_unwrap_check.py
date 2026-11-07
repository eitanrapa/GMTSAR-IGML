"""A NEW track's stored pairs checked without a lab stack to compare with (10-29): per pair the coherence (corr.grd,
median of the finite cells), the share of cells unwrapped, the unwrapped phase's spread (p1..p99 about its median),
congruence with the wrapped phase it was unwrapped from (unwrap - phasefilt a whole number of cycles where both are
finite), the connected components; per triplet of consecutive dates (a, b, c) with all three pairs stored, the closure
u_ab + u_bc - u_ac less its median: its robust sd and the share of cells more than half a cycle off (whole-cycle
unwrapping errors; T173's proc1 batch: 0.03-0.28 % a cycle off against the lab's pairs). Read only.

usage (igml env): python newtrack_unwrap_check.py PRODUCT_MERGE_DIR [MAX_PAIRS]
"""
import re
import sys
from datetime import datetime, timedelta
from itertools import combinations
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

root = Path(sys.argv[1])
ymd = lambda c: (datetime(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d")
pairs = sorted((d.name for d in root.iterdir() if re.match(r"^\d{7}_\d{7}$", d.name) and (d / "unwrap.grd").exists()),
               key=lambda n: (n[8:], n[:7]))
if len(sys.argv) > 2:
    pairs = pairs[:int(sys.argv[2])]
U, rows = {}, []
for p in pairs:
    d = root / p
    u = read_grd(d / "unwrap.grd").z.astype(np.float64)
    c = read_grd(d / "corr.grd").z.astype(np.float64)
    f = read_grd(d / "phasefilt.grd").z.astype(np.float64)
    cc = read_grd(d / "conncomp.grd").z if (d / "conncomp.grd").exists() else None
    fin = np.isfinite(u)
    U[p] = u
    both = fin & np.isfinite(f)
    k = (u - f)[both] / (2 * np.pi)
    off = np.abs(k - np.round(k)) > 0.02                                    # not a whole number of cycles
    q = np.percentile(u[fin] - np.median(u[fin]), [1, 99]) if fin.any() else (np.nan, np.nan)
    comp = ""
    if cc is not None:
        lab, n = np.unique(cc[fin & (cc > 0)], return_counts=True) if (fin & (cc > 0)).any() else ([], [])
        comp = f", {sum(x > 0.01 * fin.sum() for x in n)} components > 1 % (largest {max(n) / fin.sum():.0%})" if len(n) else ""
    dt = (datetime.strptime(ymd(p[8:]), "%Y%m%d") - datetime.strptime(ymd(p[:7]), "%Y%m%d")).days
    rows.append((p, dt, np.nanmedian(c), fin.mean(), q, off.mean() if both.any() else np.nan))
    print(f"{p} ({ymd(p[:7])}..{ymd(p[8:])}, {dt:3d} d): coherence median {np.nanmedian(c):.2f}, unwrapped "
          f"{fin.mean():.0%} of the grid, p1..p99 {q[0]:+.1f}..{q[1]:+.1f} rad, not congruent with phasefilt "
          f"{off.mean() if both.any() else np.nan:.2%}{comp}", flush=True)
dates = sorted({x for p in pairs for x in (p[:7], p[8:])})
trip = []
for a, b, c in zip(dates, dates[1:], dates[2:]):
    ab, bc, ac = f"{a}_{b}", f"{b}_{c}", f"{a}_{c}"
    if all(x in U for x in (ab, bc, ac)):
        cl = U[ab] + U[bc] - U[ac]
        ok = np.isfinite(cl)
        r = cl[ok] - np.median(cl[ok])
        sd = 1.4826 * np.median(np.abs(r))
        trip.append((a, b, c, sd, np.mean(np.abs(r) > np.pi), ok.mean()))
        print(f"closure {ymd(a)} {ymd(b)} {ymd(c)}: robust sd {sd:.3f} rad, {np.mean(np.abs(r) > np.pi):.2%} of "
              f"{ok.sum()} cells > half a cycle off", flush=True)
if trip:
    print(f"{len(trip)} triplets: closure robust sd median {np.median([t[3] for t in trip]):.3f} rad, cells a cycle off "
          f"median {np.median([t[4] for t in trip]):.2%} (max {max(t[4] for t in trip):.2%})")
by = {}
for p, dt, c, *_ in rows:
    by.setdefault(dt, []).append(c)
print("coherence median by pair length: " + ", ".join(f"{dt} d {np.median(v):.2f} ({len(v)})" for dt, v in sorted(by.items())))
