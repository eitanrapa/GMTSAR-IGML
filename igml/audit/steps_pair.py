"""Swath-wide burst-switch steps of ONE pair from its aligned SLCs, in three versions: pass 1 (OUT/<T>/superseded_nesd1),
pass 2 (the product's raw/) and the lab's own images (product_<T>/F<n>/raw), 2026-11-01. As steps_sample.py.
usage (seasat, the venv): python steps_pair.py T DATE_A DATE_B"""
import importlib.util
import sys
from pathlib import Path
import numpy as np

sys.argv, _a = sys.argv[:2] + ["0", "0"], sys.argv
spec = importlib.util.spec_from_file_location("ss", "/zdata/datasets/IGML/work/a1_test/steps_sample.py")
ss = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(ss)
except SystemExit:
    pass
T, a, b = _a[1], _a[2], _a[3]
tr, pt = ss.tr, ss.pt
for n in tr.subs:
    cs = pt._master_cut_lines(tr, n)
    cs = [cs[int(i)] for i in np.linspace(0, len(cs) - 1, min(8, len(cs)))]
    lab = lambda d: tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
    res = {}
    for name, A, B in (("pass 1", ss.path(a, n, 1), ss.path(b, n, 1)), ("pass 2", ss.path(a, n, 2), ss.path(b, n, 2)),
                       ("lab", lab(a), lab(b))):
        try:
            s = ss.steps(A, B, cs)
            res[name] = f"median {np.median(s):+6.1f} (|max| {np.max(np.abs(s)):5.1f}, n {len(s)})" if s else "-"
        except Exception as e:
            res[name] = f"error {str(e)[:60]}"
    print(f"{T} {a} x {b} F{n}: " + " | ".join(f"{k} {v}" for k, v in res.items()), flush=True)
