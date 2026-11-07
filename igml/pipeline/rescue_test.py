"""Test of the rescue stage on seasat (2026-10-27): stitch partial-swath dates with clipped pins, align them to the
lab's super-master and check them against their nearest lab-aligned neighbours, all in
/zdata/datasets/IGML/work/<T>_rescue_test. The track's state is only read.

usage (seasat): /zdata/datasets/IGML/venv/bin/python rescue_test.py T064 20190117 [20160812 ...]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_track as pt  # noqa: E402

T, dates = sys.argv[1], sys.argv[2:]
tr = pt.Track(T)
tr.save = lambda: None
plan = tr.state["plan"]
tr.out = tr.work = pt.WORK_ROOT / f"{T}_rescue_test"
tr.outp, tr.outd = tr.work / f"product_{T}", tr.work / f"data_{T}"
tr.work.mkdir(parents=True, exist_ok=True)
tr.state = {"plan": {d: plan[d] for d in dates}, "stitched": {}, "aligned": []}
pt.rescue(tr, dates)
ok = [d for d in dates if tr.state["rescued"][d]["ok"]]
if ok:
    pt.align(tr, ok, batch=len(ok))
lab = sorted(d for d in tr.lab_dates() if tr.lab_aligned(d))
for d in ok:
    nb = [x for x in lab if x < d][-1:] + [x for x in lab if x > d][:1]
    for n in (1, 2, 3):
        mine = pt._blocks(tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC")
        for x in nb:
            cov, coh, flo = pt._qc_pair(mine, pt._blocks(tr.labp / f"F{n}/raw/S1_{x}_ALL_F{n}.SLC"))
            good = cov >= pt.QC_COVER and coh >= flo + pt.QC_MARGIN
            pt.log(f"rescue test {T} {d} F{n} vs {x}: cover {cov:.3f}, coh p90 {coh:.3f}, floor {flo:.3f} -> "
                   f"{'ok' if good else 'FAIL'}")
print(json.dumps(tr.state["rescued"], indent=1))
