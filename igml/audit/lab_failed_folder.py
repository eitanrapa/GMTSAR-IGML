"""What is in the lab's merge/failed folder of a track, and did we remake those pairs? (2026-11-06; user: "is there an interferogram data on the failed subfolder from the original labs
T173 that we never remade?") Read only. For every entry of product_<T>/merge/failed: its files, whether it holds real grids, its dates (satellite, days apart), and whether the pair is in our
stored list, in the lab's merged pairs now, in the archive's, whether both dates passed our qc, and whether the lab's subswath folders hold it.
usage (seasat venv): python lab_failed_folder.py T"""
import importlib.util
import os
import re
import sys
from collections import Counter
from datetime import date, timedelta
from pathlib import Path

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T = sys.argv[1]
tr = pt.Track(T)
s = tr.state
F = tr.labp / "merge/failed"
ARCH = pt.LAB_ARCHIVE / T / f"product_{T}"
print(f"{F}: {len(os.listdir(F))} entries; modified {__import__('time').strftime('%F', __import__('time').localtime(F.stat().st_mtime))}")


def ymd(c):
    return (date(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d")


entries = sorted(os.listdir(F))
pairs = [e for e in entries if re.match(r"^\d{7}_\d{7}$", e)]
print("not pair folders:", [e for e in entries if e not in pairs][:20])
stored = set(s.get("stored", []))
q = s.get("qc", {})
sat = {d: Path(fr[0]).name[:3] for d, fr in s.get("stitched", {}).items() if fr}
lab_now = {p.name for p in (tr.labp / "merge").iterdir() if re.match(r"^\d{7}_\d{7}$", p.name)}
arch = {p.name for p in (ARCH / "merge").iterdir() if re.match(r"^\d{7}_\d{7}$", p.name)}
sw = {n: set(os.listdir(tr.labp / f"F{n}/intf_all")) for n in tr.subs}
asw = {n: set(os.listdir(ARCH / f"F{n}/intf_all")) for n in tr.subs}
oursm = {p.name for p in (tr.outp / "merge").iterdir() if re.match(r"^\d{7}_\d{7}$", p.name)}
files = Counter()
rows = []
grids = Counter()
for p in pairs:
    d = F / p
    names = sorted(os.listdir(d))
    for n in names:
        files[n] += 1
    a, b = ymd(p[:7]), ymd(p[8:])
    size = {n: (d / n).stat().st_size if (d / n).exists() else -1 for n in ("phasefilt.grd", "corr.grd", "unwrap.grd", "merge_log")}
    for n, v in size.items():
        grids[n] += v > 1000 if n != "merge_log" else v > 0
    days = (date(int(b[:4]), int(b[4:6]), int(b[6:])) - date(int(a[:4]), int(a[4:6]), int(a[6:]))).days
    rows.append({"p": p, "a": a, "b": b, "days": days, "sa": sat.get(a), "sb": sat.get(b), "stored": p in stored, "lab_now": p in lab_now, "arch": p in arch, "ours_merge": p in oursm,
                 "qa": q.get(a, {}).get("ok"), "qb": q.get(b, {}).get("ok"), "sw": [p in sw[n] for n in tr.subs], "asw": [p in asw[n] for n in tr.subs], "n": len(names), "size": size})
print(f"{len(pairs)} pair folders; with phasefilt / corr / unwrap / merge_log over 1 kB: {grids['phasefilt.grd']} / {grids['corr.grd']} / {grids['unwrap.grd']} / {grids['merge_log']}")
print("files per folder (top):", files.most_common(14))
print("in our stored list:", sum(r["stored"] for r in rows), "| same pair in lab merge/ now (ours):", sum(r["lab_now"] for r in rows), "| in the archive's merge/:", sum(r["arch"] for r in rows),
      "| in OUT merge/:", sum(r["ours_merge"] for r in rows))
print("pairs NOT remade:", sum(not r["stored"] and not r["lab_now"] and not r["ours_merge"] for r in rows))
print("both dates passed our qc:", sum(bool(r["qa"]) and bool(r["qb"]) for r in rows), "| a date with no qc record:", sum(r["qa"] is None or r["qb"] is None for r in rows))
print("lab tree F<n>/intf_all holds it (F1, F2, F3 all):", sum(all(r["sw"]) for r in rows), "| archive's:", sum(all(r["asw"]) for r in rows))
print("days apart:", sorted(Counter(r["days"] for r in rows).items())[:24])
print("years of the first date:", sorted(Counter(r["a"][:4] for r in rows).items()))
print("satellites:", Counter(f"{r['sa']} x {r['sb']}" for r in rows).most_common(8))
for r in rows[:60]:
    print(f"  {r['p']} {r['a']} x {r['b']} {r['days']:3d} d {r['sa']}x{r['sb']} files {r['n']:2d} phasefilt {r['size']['phasefilt.grd']:>9d} unwrap {r['size']['unwrap.grd']:>9d} "
          f"stored {int(r['stored'])} lab-now {int(r['lab_now'])} arch {int(r['arch'])} OUT {int(r['ours_merge'])} qc {r['qa']} {r['qb']} sw {r['sw']}")
all_pairs = {f"{pt.doy0(a)}_{pt.doy0(b)}" for a, b in pt.all_pairs(tr)}
print("in today's planned pair list:", sum(r["p"] in all_pairs for r in rows))
