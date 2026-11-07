"""merge_swath's seam columns for every merged pair (2026-10-28, GMTSAR audit M1). The 3-subswath branch of
merge_swath.c writes G2[jj - N1 + ovl12 - 1] at output column jj (the 2-subswath branch and F3 have no -1), so F2's block
[s12, s23) sits one decimated column too far in range; the conversion moves it back and needs each pair's s12 / s23:

  s12 = N1 - ovl12 + n1                 s23 = N1 + N2 - ovl12 - 1 - (ovl23 - n2)          width = N1 + N2 + N3 - ovl12 - ovl23 - 2
  ovl12 = N1 - round((near2 - near1) / (c / fs / 2) / incx)   ovl23 = N2 - round((near3 - near2) / (c / fs / 2) / incx)
  n1 = max(10, ceil((-rshift + first_sample + 150) / incx)) with F2's values, n2 the same with F3's

near / fs from the SUPER-MASTER's PRMs (merge_batch.csh copies the first line's reference PRM as every pair's
supermaster.PRM), N1 N2 N3 the subswath grids' widths, incx their x increment (32 samples); rshift the SECONDARY's (the
last rshift line of its aligned PRM, merge_unwrap_geocode_tops.csh L61), first_sample the largest of the super-master's,
the reference's and the secondary's (merge_batch.csh L68-71, merge_unwrap_geocode_tops.csh L62-66).

Every lab pair carries merge_log; ours stored before 10-28 do not (store() copies it since). check: the rule against the
lab's merge_logs, with the PRMs the merge used (F<n>/intf_all/<pair>/) and with the dates' aligned PRMs in F<n>/raw/
(how ours are reconstructed), + our work pairs that still hold a merge_log. backfill: merge_seams.json into OUR stored
pairs without merge_log (our tree only; a dry run unless --apply). Read only otherwise.

usage (seasat): python merge_seams.py check T [N]   |   python merge_seams.py backfill T [--apply]
"""
import json
import math
import random
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, "/zdata/datasets/IGML/code/scratchpad/s1_gapfill")
MODE, T, REST = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.argv = sys.argv[:1]
import process_track as P                                       # noqa: E402

C = 299792458.0
PAIR = re.compile(r"^\d{7}_\d{7}$")
FIELDS = ("ovl12", "ovl23", "n1", "n2", "width")


def prm(path):
    """A PRM as {name: value}; the LAST occurrence wins, as GMTSAR's reader (aligned PRMs carry rshift twice)."""
    out = {}
    for line in open(path):
        k, s, v = line.partition("=")
        if s:
            out[k.strip()] = v.strip()
    return out


def cround(x):
    """C's round(): halves away from zero."""
    return math.floor(x + 0.5) if x >= 0 else -math.floor(-x + 0.5)


def grid(path):
    """(n_columns, x increment) of a GMT grid."""
    f = subprocess.run(["gmt", "grdinfo", "-C", str(path), "--GMT_HISTORY=false"], capture_output=True, text=True,
                       check=True).stdout.split()
    return int(f[9]), float(f[7])


def from_log(d):
    t = (d / "merge_log").read_text()
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    return {"ovl12": int(o[1]), "ovl23": int(o[2]), "width": int(re.search(r"Size\((\d+)x\d+\)", t)[1]),
            "n1": int(re.search(r"Stitching location n1 = (-?\d+)", t)[1]),
            "n2": int(re.search(r"Stitching location n2 = (-?\d+)", t)[1])}


def seams(sm, ref, sec, N, incx):
    """merge_swath's numbers from {n: PRM} of the super-master / reference / secondary and {n: grid width}."""
    dr = C / float(sm[1]["rng_samp_rate"]) / 2
    near = {n: float(sm[n]["near_range"]) for n in (1, 2, 3)}
    ovl12 = N[1] - cround((near[2] - near[1]) / dr / incx)
    ovl23 = N[2] - cround((near[3] - near[2]) / dr / incx)

    def ncol(n):
        first = max(int(sm[n]["first_sample"]), int(ref[n]["first_sample"]), int(sec[n]["first_sample"]))
        rshift = int(float(sec[n]["rshift"]))               # an int in GMTSAR's PRM struct
        return max(10, math.ceil((-rshift + first + 150.0) / incx))

    out = {"ovl12": ovl12, "ovl23": ovl23, "n1": ncol(2), "n2": ncol(3),
           "width": N[1] + N[2] + N[3] - ovl12 - ovl23 - 2}
    return out


def spans(s, N):
    """F2's output columns [s12, s23) in the merged grid."""
    return {"s12": N[1] - s["ovl12"] + s["n1"], "s23": N[1] + N[2] - s["ovl12"] - 1 - (s["ovl23"] - s["n2"])}


tr = P.Track(T)
sm_date = tr.master_pair[0]          # merge_batch.csh: the FIRST line's reference PRM is every pair's supermaster.PRM
raw = lambda n, d: next((p for p in (tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.PRM", tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.PRM")
                         if p.exists()), None)
SM = {n: prm(raw(n, sm_date)) for n in (1, 2, 3)}
smdir = f"{P.doy0(tr.master_pair[0])}_{P.doy0(tr.master_pair[1])}"
N0, INC = {}, None
for cand in [smdir] + sorted(d.name for d in (tr.labp / "F1/intf_all").iterdir() if PAIR.match(d.name)):
    try:                                 # the super-master pair's grids, or the first readable pair's (T064's are
        got = [grid(tr.labp / f"F{n}/intf_all/{cand}/phasefilt.grd") for n in (1, 2, 3)]      # emptied)
    except (subprocess.CalledProcessError, FileNotFoundError):
        continue
    N0, INC = {n: g[0] for n, g in zip((1, 2, 3), got)}, got[0][1]
    break
if not N0:              # every subswath grid deleted (T064): the super-master's num_rng_bins // incx, which equals the
    INC = grid(next(d for d in sorted((tr.labp / "merge").iterdir()) if PAIR.match(d.name)    # grids' widths on the
                    and (d / "unwrap.grd").exists()) / "unwrap.grd")[1]                         # eight other tracks
    N0 = {n: int(SM[n]["num_rng_bins"]) // int(round(INC)) for n in (1, 2, 3)}
print(f"{T}: super-master {sm_date} (first merge line {smdir}), subswath widths {N0}, incx {INC}", flush=True)

if MODE == "check":
    n_max = int(REST[0]) if REST else 40
    labs = sorted(d for d in (tr.labp / "merge").iterdir() if PAIR.match(d.name) and (d / "merge_log").exists()
                  and (d / "unwrap.grd").exists())
    random.seed(1)
    sample = random.sample(labs, min(n_max, len(labs)))
    work = [d for d in (tr.work / "intf/merge").iterdir() if PAIR.match(d.name) and (d / "merge_log").exists()] \
        if (tr.work / "intf/merge").is_dir() else []
    bad = {"used": 0, "raw": 0}
    widths = set()
    for d in sample + work:
        a, b = (P.ymd0(x) for x in d.name.split("_"))
        lab = d.parent == tr.labp / "merge"
        base = (lambda n: tr.labp / f"F{n}/intf_all/{d.name}") if lab else (lambda n: tr.work / f"intf/F{n}/intf_all/{d.name}")
        try:
            N = {n: grid(base(n) / "phasefilt.grd")[0] for n in (1, 2, 3)}
        except subprocess.CalledProcessError:                # the lab deleted this pair's subswath grids
            continue
        widths.add(tuple(N.values()))
        want = from_log(d)
        got = {}
        if lab:
            got["used"] = seams(SM, {n: prm(base(n) / f"S1_{a}_ALL_F{n}.PRM") for n in (1, 2, 3)},
                                {n: prm(base(n) / f"S1_{b}_ALL_F{n}.PRM") for n in (1, 2, 3)}, N, INC)
        try:
            got["raw"] = seams(SM, {n: prm(raw(n, a)) for n in (1, 2, 3)}, {n: prm(raw(n, b)) for n in (1, 2, 3)},
                               N0, INC)
        except TypeError:                                   # a date with no aligned PRM left in raw/
            got["raw"] = None
        line = [f"{d.name} [{'lab' if lab else 'ours'}] log {[want[k] for k in FIELDS]}"]
        for k, g in got.items():
            if g is None:
                line.append(f"{k}: no PRM")
                continue
            diff = {f: (g[f], want[f]) for f in FIELDS if g[f] != want[f]}
            bad[k] += bool(diff)
            line.append(f"{k}: {'ok' if not diff else diff}")
        print("  ".join(line), flush=True)
    print(f"{T} check: {len(sample)} lab pairs + {len(work)} of ours; mismatches {bad}; subswath widths seen {widths}")

elif MODE == "backfill":
    apply = "--apply" in REST
    todo = sorted(d for d in (tr.outp / "merge").iterdir() if PAIR.match(d.name) and not (d / "merge_log").exists()
                  and not (d / "merge_seams.json").exists())
    n_ok = 0
    for d in todo:
        a, b = (P.ymd0(x) for x in d.name.split("_"))
        s = seams(SM, {n: prm(raw(n, a)) for n in (1, 2, 3)}, {n: prm(raw(n, b)) for n in (1, 2, 3)}, N0, INC)
        w = grid(d / "unwrap.grd")[0]
        if w != s["width"]:
            print(f"{d.name}: width {s['width']} from the PRMs, {w} in unwrap.grd -- NOT written", flush=True)
            continue
        rec = dict(s, **spans(s, N0), N=[N0[1], N0[2], N0[3]], incx=INC, source="PRMs (merge_seams.py)",
                   prms={"supermaster": str(raw(1, sm_date)).replace("F1/", "F<n>/").replace("_F1.", "_F<n>."),
                         "reference": str(raw(1, a)).replace("F1/", "F<n>/").replace("_F1.", "_F<n>."),
                         "secondary": str(raw(1, b)).replace("F1/", "F<n>/").replace("_F1.", "_F<n>.")})
        if apply:
            (d / "merge_seams.json").write_text(json.dumps(rec, indent=1))
        n_ok += 1
    print(f"{T} backfill: {len(todo)} of our pairs without merge_log / merge_seams.json, {n_ok} "
          f"{'written' if apply else 'would be written (dry run)'}")
else:
    sys.exit(__doc__)
