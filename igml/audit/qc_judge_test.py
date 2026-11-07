"""The 11-06 qc (_qc_judge: nearest neighbours, water blocks, votes from neighbours that see the ground) on dates whose
answer is known, on the LAB's images (read only; nothing saved):
  must FAIL: T144 2016-06-30, 2016-11-21, 2020-10-25, 2024-04-01 (misregistered, lab survey 10-28; re-cuts cured them);
             T071 2015-11-22 (the last third of every subswath holds other content, 11-06);
  must PASS: T071's 2017-02-20, 2020-03-06, 2025-02-02 (the old qc's open-water flags), 2022-09-22's neighbours;
             controls next to each failing date.
usage (seasat venv): python qc_judge_test.py"""
import importlib.util
import os
import re

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
P.slc_src = lambda tr, d, n: tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"


def images(tr):
    n = tr.subs[0]
    raw = tr.labp / f"F{n}/raw"
    out = []
    for f in os.listdir(raw):
        m = re.match(rf"S1_(\d{{8}})_ALL_F{n}\.SLC$", f)
        if m and (raw / f).stat().st_size > 0:
            out.append(m.group(1))
    return sorted(out)


def run(T, dates):
    tr = P.Track(T)
    L = images(tr)
    print(f"{T}: {len(L)} lab images, subswaths {tr.subs}", flush=True)
    for d, want in dates:
        if d not in L:
            print(f" {d}: no lab image", flush=True)
            continue
        i = L.index(d)
        before, after = L[:i], L[i + 1:]
        near = sorted(before[-3:] + after[:3], key=lambda x: abs((P.day(x) - P.day(d)).days))
        nb = before[-1:] + after[:1]
        rec, ok = P._qc_judge(tr, d, nb, P.QC_COVER)
        first = ok
        spare = [x for x in near if x not in nb][:P.QC_NEIGHBOURS - len(nb)]
        if not ok and spare:
            rec, ok = P._qc_judge(tr, d, nb + spare, P.QC_COVER)
        s = " ".join(f"{k}:{v['coh']:.2f}/{v['floor']:.2f}" + (f" shifted {v['shifted']}" if "shifted" in v else "") +
                     (f" content {v['content']}" if "content" in v else "") for k, v in rec.items() if k.startswith("F") and v)
        verdict = "as expected" if ok == want else "*** NOT AS EXPECTED ***"
        print(f" {d}: {'ok' if ok else 'FAILED'} (first judge {'ok' if first else 'failed'}) want {'ok' if want else 'fail'} "
              f"{verdict} | {s}", flush=True)


run("T144", [("20160630", False), ("20161121", False), ("20201025", False), ("20240401", False),
             ("20160606", True), ("20161103", True), ("20201013", True), ("20240320", True)])
run("T071", [("20151122", False), ("20151029", True), ("20170220", True), ("20200306", True), ("20250202", True),
             ("20220910", True)])
