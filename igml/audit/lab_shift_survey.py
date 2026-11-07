"""The lab's OWN aligned dates through qc v2's line-shift test (2026-10-28, GMTSAR audit: stitch_tops has applied the
last frame's overlap line offset twice since 9d55e66, 2025-01-21, and the lab aligned all its multi-frame stacks after
that date -- T137 2025-02, the rest 2026-08..10). Per track, every lab date whose aligned image holds data, against the
nearest such dates before and after, per subswath: a date is MISREGISTERED when some block's coherence peaks at a 1-2
line shift (beats 0 by QC_SHIFT_GAIN) against both neighbours -- process_track.qc's rule. Read only (the lab's SLCs);
results appended to work/chk/lab_shift_<T>.json as they come.

TRAP (10-28): the first / last date has ONE neighbour, so a misregistered neighbour flags it too (T173 2014-11-10 via
2014-12-04: clean against 12-28). Re-check a flagged end date with shift_profile.py against a further date.
TRAP 2 (10-28): a block over OPEN WATER can flag a good date (T035 2026-03-15, F2 block 1 at line 4,038: p90 0.21-0.23
unshifted, +0.03-0.04 at -2 against both neighbours, evenly across the whole range width, gone 72 lines later; every
coherent stretch peaks at 0 in shift_profile, and the lab's pairs with it are as coherent as the controls in
corr_profile.py). A real misregistration peaks well above the noise floor (T173 0.71, T144 0.63).

usage (seasat): python lab_shift_survey.py T [T ...]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, "/zdata/datasets/IGML/code/scratchpad/s1_gapfill")
sys.argv, tracks = sys.argv[:1], sys.argv[1:]
import process_track as P                                       # noqa: E402

for T in tracks:
    tr = P.Track(T)
    raw = {n: tr.labp / f"F{n}/raw" for n in (1, 2, 3)}
    dates = sorted(f.name[3:11] for f in raw[1].glob("S1_*_ALL_F1.SLC"))
    have = [d for d in dates if all(P._holds_data(raw[n] / f"S1_{d}_ALL_F{n}.SLC") for n in (1, 2, 3))]
    out_path = Path(f"/zdata/datasets/IGML/work/chk/lab_shift_{T}.json")
    res = json.load(open(out_path)) if out_path.exists() else {}
    P.log(f"{T} lab shift survey: {len(have)} of {len(dates)} aligned dates hold data, {len(res)} done before")
    cache = {}

    def blocks(d, n):
        if (d, n) not in cache:
            cache[(d, n)] = P._blocks(raw[n] / f"S1_{d}_ALL_F{n}.SLC")
        return cache[(d, n)]

    for i, d in enumerate(have):
        if d in res:
            continue
        nb = ([have[i - 1]] if i else []) + ([have[i + 1]] if i + 1 < len(have) else [])
        rec = {"vs": nb}
        for n in (1, 2, 3):
            off, c0s = None, []
            for x in nb:
                blk = P._qc_pair(blocks(d, n), blocks(x, n), shifts=True)[3]
                mis = {k: s for k, (c0, cs, s) in enumerate(blk) if cs > c0 + P.QC_SHIFT_GAIN}
                off = mis if off is None else {k: s for k, s in off.items() if k in mis}
                c0s.append([round(c0, 3) for c0, _, _ in blk])
            rec[f"F{n}"] = {"shifted": {str(k): s for k, s in sorted((off or {}).items())}, "coh": c0s}
        rec["misregistered"] = any(rec[f"F{n}"]["shifted"] for n in (1, 2, 3))
        res[d] = rec
        json.dump(res, open(out_path, "w"))
        if rec["misregistered"]:
            P.log(f"{T} lab shift survey {d}: MISREGISTERED " +
                  " ".join(f"F{n} {rec[f'F{n}']['shifted']}" for n in (1, 2, 3) if rec[f"F{n}"]["shifted"]))
        for k in [k for k in cache if k[0] < d]:          # the next date needs d's blocks and its own only
            del cache[k]
    bad = sorted(d for d, v in res.items() if v["misregistered"])
    P.log(f"{T} lab shift survey done: {len(res)} dates, {len(bad)} misregistered: {bad}")
