"""Network spectral diversity test (2026-10-31) on a collecting run of a1_validate.py (IGML_ESD_COLLECT): every pair of
the run's frames aligned to one super-master frame through GMTSAR's spectral_diversity (process_track._sd_pair), the
triangle closures, the network solution (process_track._nesd_solve, super-master = 0) against each date's estimate
against the super-master, and a shifts file for the pass-2 run (IGML_ESD_SHIFTS).

usage (seasat, the venv): python nesd_test.py RUN_DIR PROCESS_TRACK.py SHIFTS_OUT [K]
  RUN_DIR  e.g. /zdata/datasets/IGML/work/a1_test/T166/F/F1 (nesd_frames.txt, masterlist, the overlap files)
  K        pair each node with its next K nodes in time (default: every pair)
"""
import importlib.util
import itertools
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

run, ptf, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
K = int(sys.argv[4]) if len(sys.argv) > 4 else 99
spec = importlib.util.spec_from_file_location("pt", ptf)
pt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pt)


class T:                                       # the bits of a Track _sd_pair uses
    work = run.parent.parent / "nesd_work"


masters = [s.strip() for s in (run / "masterlist").read_text().split("\n") if s.strip()]
frames = {}
for line in (run / "nesd_frames.txt").read_text().split("\n"):
    p = line.split()
    if len(p) == 3:
        frames[p[0]] = {"master": p[1], "tmp_da": int(p[2]),
                        **pt._sd_parse((run / f"{p[0]}.sd").read_text())}
share = "/zdata/datasets/IGML/gmtsar_igml/install/share/gmtsar"
filt = f"{share}/filters/{pt.SPECDIV_FILTER}"
rows = []
for m in masters:
    nodes = [m] + sorted(s for s, v in frames.items() if v["master"] == m)      # time order (stems hold the date)
    prs = [(a, b) for i, a in enumerate(nodes) for b in nodes[i + 1:i + 1 + K]]
    with ThreadPoolExecutor(max_workers=10) as ex:
        obs = dict(zip(prs, ex.map(lambda ab: pt._sd_pair(T, 0, ab[0], run, ab[1], run, filt), prs)))
    for (a, b), r in obs.items():
        print(f"{a} -> {b}: shift {r.get('shift', float('nan')):+.5f} +- {r.get('sigma', float('nan')):.5f} lines, "
              f"coherent sum {r.get('residual_shift', float('nan')):+.5f}, rows {r['n']}, corr {r.get('corr')}",
              flush=True)
    sep = np.median([r["spectral_spectrationXdta"] for r in obs.values()])
    deg = 360.0 * sep
    for a, b, c in itertools.combinations(nodes, 3):
        if (a, b) in obs and (b, c) in obs and (a, c) in obs:
            cl = obs[(a, b)]["shift"] + obs[(b, c)]["shift"] - obs[(a, c)]["shift"]
            print(f"closure {a[3:11]} {b[3:11]} {c[3:11]}: {cl:+.5f} lines ({cl * deg:+.1f} deg)")
    sol, res, scale = pt._nesd_solve(m, [(a, b, r["shift"], r["sigma"]) for (a, b), r in obs.items()
                                        if r.get("n", 0) >= pt.NESD_MIN_ROWS and "shift" in r])
    print(f"frame {m}: robust z scale {scale:.2f}; 1 line = {deg:.0f} deg of burst-boundary step")
    for s in nodes[1:]:
        x, sx, npr = sol.get(s, (float("nan"), float("nan"), 0))
        e = frames[s].get("esd1_shift", float("nan"))
        print(f"  {s}: network {x:+.5f} +- {sx:.5f} ({npr} pairs) | against the super-master {e:+.5f} | "
              f"difference {(x - e) * deg:+.1f} deg")
        rows.append((s, x))
out.write_text("".join(f"{s} {x:.9f}\n" for s, x in rows))
print(f"shifts -> {out}")
