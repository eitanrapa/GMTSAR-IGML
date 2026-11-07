"""Burst-boundary steps before / after network spectral diversity (2026-11-01), from qc's records in a track's
state.json: pass 1 (geometry-only, kept as qc_nesd1) against pass 2 (the network's shifts, qc). Each record holds the
step of a date against its qc neighbour at up to 4 of the super-master's burst-switch lines, per subswath.
usage: python steps_before_after.py state.json"""
import json, sys
from pathlib import Path
import numpy as np

s = json.load(open(sys.argv[1]))
sat = {k: Path(v[0]).name[:3] for k, v in s.get("stitched", {}).items() if v}
sat[s.get("master", "")] = sat.get(s.get("master", ""), "S1A")


def rows(q):
    out = []
    for d, rec in (q or {}).items():
        for n in ("F1", "F2", "F3"):
            v = rec.get(n) or {}
            if v.get("steps"):
                out.append((d, v["vs"], n, np.array(v["steps"], float)))
    return out


for name, q in (("pass 1 (geometry only)", s.get("qc_nesd1")), ("pass 2 (network shifts)", s.get("qc"))):
    r = rows(q)
    if not r:
        print(f"{name}: no records")
        continue
    a = np.concatenate([x[3] for x in r])
    print(f"{name}: {len(r)} date x subswath records, {a.size} boundary steps; |step| median {np.median(np.abs(a)):.1f} deg, "
          f"p90 {np.percentile(np.abs(a), 90):.1f}, p99 {np.percentile(np.abs(a), 99):.1f}, max {np.abs(a).max():.1f}")
    for n in ("F1", "F2", "F3"):
        b = np.concatenate([x[3] for x in r if x[2] == n])
        print(f"   {n}: |step| median {np.median(np.abs(b)):.1f}, p90 {np.percentile(np.abs(b), 90):.1f}")
    g = {}
    for d, x, n, st in r:
        g.setdefault("-".join(sorted({sat.get(d, "?"), sat.get(x, "?")})), []).append(st)
    for k, v in sorted(g.items()):
        v = np.concatenate(v)
        print(f"   pairs {k:8s}: {v.size:5d} steps, |step| median {np.median(np.abs(v)):5.1f}, p90 {np.percentile(np.abs(v), 90):5.1f}")
    years = {}
    for d, x, n, st in r:
        years.setdefault(d[:4], []).append(st)
    print("   |step| p90 by year: " + "  ".join(f"{y} {np.percentile(np.abs(np.concatenate(v)), 90):.0f}" for y, v in sorted(years.items())))
