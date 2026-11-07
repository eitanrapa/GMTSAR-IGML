"""Per-boundary burst-switch steps of early (2014-15) pairs after pass 2, against pass 1 (2026-11-02; the open "early dates
~10 deg" of AUDIT 11-01): every burst switch of the super-master grid, each subswath, from the aligned SLCs (steps_sample.py's
steps()). A constant residual along the swath = a misregistration the network left (its shift for a date is off); one that
changes along the swath = within-frame variation a constant per frame cannot remove.
usage (seasat, the venv): python steps_early.py T DATE_A:DATE_B [DATE_A:DATE_B ...]"""
import importlib.util
import sys
import numpy as np

sys.argv, _a = sys.argv[:2] + ["0", "0"], sys.argv
spec = importlib.util.spec_from_file_location("ss", "/zdata/datasets/IGML/work/a1_test/steps_sample.py")
ss = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(ss)
except SystemExit:
    pass
T, pairs = _a[1], [p.split(":") for p in _a[2:]]
tr, pt = ss.tr, ss.pt
for a, b in pairs:
    for n in tr.subs:
        cs = pt._master_cut_lines(tr, n)
        out = []
        for ps in (1, 2):
            try:
                s = ss.steps(ss.path(a, n, ps), ss.path(b, n, ps), cs)
                out.append(f"pass {ps}: " + " ".join(f"{v:+.0f}" for v in s))
            except Exception as e:
                out.append(f"pass {ps}: error {str(e)[:60]}")
        print(f"{T} {a} x {b} F{n} ({len(cs)} switches)\n   " + "\n   ".join(out), flush=True)
