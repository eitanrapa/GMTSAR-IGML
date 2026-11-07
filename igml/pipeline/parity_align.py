"""Parity test of the align stage on seasat (2026-10-27): re-align dates the LAB has aligned too, from our stitched
frames, with seasat's GMTSAR (rebuilt 10-15, after the lab's 09-24 run), and compare with the lab's SLC / PRM byte by
byte. The Mac pilot against the same lab files: PRMs identical, ~25,000 of 2.4e9 SLC bytes differ (by 1). Writes only
to /zdata/datasets/IGML/work/<T>_parity; the track's state is not touched.

usage (seasat): /zdata/datasets/IGML/venv/bin/python parity_align.py T173 20150310 [20250615 ...]
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_track as pt  # noqa: E402

T, dates = sys.argv[1], sys.argv[2:]
tr = pt.Track(T)
tr.save = lambda: None
tr.work = pt.WORK_ROOT / f"{T}_parity"
tr.outp = tr.work / "out"
tr.state = {"stitched": tr.state.get("stitched", {}), "plan": {}, "aligned": []}
for d in dates:
    first = Path(tr.state["stitched"][d][0]).name
    orb = pt.orbit_for(tr, d, first[:3], first[17:32])
    if orb is None:
        sys.exit(f"{d}: no orbit")
    tr.state["plan"][d] = {"orbit": str(orb)}
pt.align(tr, dates, batch=len(dates))
for d in dates:
    for n in (1, 2, 3):
        a, b = tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}", tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}"
        prm = subprocess.run(["cmp", "-s", f"{a}.PRM", f"{b}.PRM"]).returncode == 0
        nb = int(subprocess.run(f"cmp -l {a}.SLC {b}.SLC | wc -l", shell=True, capture_output=True,
                                text=True).stdout.strip() or 0)
        pt.log(f"parity {T} {d} F{n}: PRM {'identical' if prm else 'DIFFERS'}, "
               f"SLC {nb} of {Path(f'{a}.SLC').stat().st_size} bytes differ")
