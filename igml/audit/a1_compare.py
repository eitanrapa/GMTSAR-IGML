"""A1 / A5 validation (2026-10-30), after a1_validate.py's runs: which lines of each aligned image differ between two
runs (or a run and the pipeline's own image), as line ranges. B (bin_a1, IGML_A1 unset) must equal the production
image byte for byte; C (IGML_A1) may differ only between each date's own and the super-master's cut lines; the
super-master never.

usage (seasat, the venv): python a1_compare.py IMG1 IMG2 [IMG1b IMG2b ...]   (GMTSAR .SLC paths, pairs)
"""
import sys
from pathlib import Path

import numpy as np


def prm(p):
    return {a.strip(): b.strip() for a, s, b in (l.partition("=") for l in open(p)) if s}


def ranges(idx):
    out, start, prev = [], None, None
    for i in idx:
        if start is None:
            start = prev = i
        elif i == prev + 1:
            prev = i
        else:
            out.append((start, prev))
            start = prev = i
    if start is not None:
        out.append((start, prev))
    return out


a = sys.argv[1:]
for p1, p2 in zip(a[0::2], a[1::2]):
    p1, p2 = Path(p1), Path(p2)
    s1, s2 = p1.stat().st_size, p2.stat().st_size
    w = int(prm(p1.with_suffix(".PRM"))["bytes_per_line"]) // 4
    if s1 != s2:
        print(f"{p1.name}: sizes differ {s1} vs {s2}")
        continue
    n = s1 // (4 * w)
    m1 = np.memmap(p1, dtype="<i4", mode="r", shape=(n, w))
    m2 = np.memmap(p2, dtype="<i4", mode="r", shape=(n, w))
    diff = []
    for r in range(0, n, 2000):
        d = np.any(m1[r:r + 2000] != m2[r:r + 2000], axis=1)
        diff += [r + i for i in np.nonzero(d)[0]]
    rg = ranges(diff)
    print(f"{p1.parent.parent.name}/{p1.name} vs {p2.parent.parent.name}/{p2.name}: {len(diff)} of {n} lines differ"
          + (f" in {len(rg)} ranges: " + ", ".join(f"{x}-{y}" for x, y in rg[:30]) + (" ..." if len(rg) > 30 else "")
             if rg else " (identical)"), flush=True)
