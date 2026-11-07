"""Fetch the precise orbits (POEORB) of every date of the extra tracks (T013 / T100 / T108; 2026-10-29: no lab stack,
so no state.json plan to read yet), into /palsar/IGML/orbits, as prefetch_orbits.py does for the nine tracks: the dates
are the scenes in the archive (/palsar/S1<A-D>/CA/<year>/<T>/, listed over ssh on seasat), the Earthdata login stays in
this Mac's ~/.netrc. Dates younger than ~3 weeks have no POEORB yet.

usage (Mac, repo root): python3 -I scratchpad/s1_gapfill/prefetch_orbits_extra.py [T100 T013 T108] [--dates d1,d2]
"""
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_track as pt  # noqa: E402

args = [a for a in sys.argv[1:] if not a.startswith("--")]
only = None
if "--dates" in sys.argv:
    only = set(sys.argv[sys.argv.index("--dates") + 1].split(","))
    args = [a for a in args if a != sys.argv[sys.argv.index("--dates") + 1]]
TRACKS = args or ["T100", "T013", "T108"]
nobody = SimpleNamespace(labd=Path("/nonexistent"))
got, none = 0, []
for T in TRACKS:
    names = subprocess.run(["ssh", "-o", "BatchMode=yes", "erapaportbruck@seasat.ucsd.edu",
                            f"ls -d /palsar/S1?/CA/20*/{T}/*.SAFE"], capture_output=True, text=True, check=True).stdout.split()
    first = {}
    for n in sorted(Path(x).name for x in names):
        first.setdefault(n[17:25], n)               # the date's earliest scene: satellite and start time
    todo = sorted(d for d in first if only is None or d in only)
    for d in todo:
        f = pt.orbit_for(nobody, d, first[d][:3], first[d][17:32])
        if f:
            got += 1
        else:
            none.append(f"{T} {d} {first[d][:3]}")
    pt.log(f"{T}: {len(todo)} dates")
for p in Path("/palsar/IGML/orbits").glob("._*"):
    p.unlink()
pt.log(f"orbits found or fetched for {got} dates; none yet for {len(none)}: {none}")
