"""Fetch the precise orbits (POEORB) seasat's plan is waiting for, into /palsar/IGML/orbits (2026-10-27: processing runs
on seasat, but ASF s1qc needs the Earthdata login, which stays in this Mac's ~/.netrc). Reads each track's state.json
from seasat over ssh; a date waits for its orbit when all its scenes are on disk and no orbit was found. Then re-run the
plan stage on seasat. Dates younger than ~3 weeks have no POEORB yet (restituted orbits are off: USE_RESORB).

usage (Mac, repo root): python3 -I scratchpad/s1_gapfill/prefetch_orbits.py [T173 T071 ...]   (default: the 9 tracks)
"""
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_track as pt  # noqa: E402

TRACKS = sys.argv[1:] or "T173 T071 T064 T166 T137 T144 T035 T042 T115".split()
nobody = SimpleNamespace(labd=Path("/nonexistent"))          # the plan on seasat already looked in the lab's data_<T>
got, none = 0, []
for T in TRACKS:
    s = json.loads(subprocess.run(["ssh", "-o", "BatchMode=yes", "erapaportbruck@seasat.ucsd.edu",
                                   f"cat /zdata/datasets/IGML/InSAR_Processing/{T}/state.json"],
                                  capture_output=True, text=True, check=True).stdout)
    wait = {d: v for d, v in s.get("plan", {}).items() if not v["ready"] and v["safes"] and all(v["safes"])}
    for d, v in sorted(wait.items()):
        first = Path(v["safes"][0]).name
        f = pt.orbit_for(nobody, d, first[:3], first[17:32])
        if f:
            got += 1
        else:
            none.append(f"{T} {d} {first[:3]}")
    pt.log(f"{T}: {len(wait)} dates waiting for an orbit")
for p in Path("/palsar/IGML/orbits").glob("._*"):        # macOS AppleDouble companions of the new files
    p.unlink()
pt.log(f"orbits found or fetched for {got} dates; none yet for {len(none)}: {none}")
