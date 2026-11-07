"""The lab's own interferograms along azimuth (2026-10-28, after lab_shift_survey flagged T035 2026-03-15 in one F2
block): corr.grd per azimuth bin for each pair, side by side, so pairs that contain a misregistered date can be
read against pairs that skip it (a misregistration halves coherence over its stretch). p90 per bin, as qc. Read only.

usage (seasat): python corr_profile.py T n pair [pair ...]      (pairs as the lab's intf_all names, e.g. 2026061_2026073)
"""
import sys

sys.path.insert(0, "/zdata/datasets/IGML/code/scratchpad/s1_gapfill")
T, n, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.argv = sys.argv[:1]
import subprocess                                                # noqa: E402

import numpy as np                                               # noqa: E402
import process_track as P                                        # noqa: E402

NBIN = 48
intf = P.Track(T).labp / f"F{n}/intf_all"


def load(pair):
    """y (line) per row and the grid, through GMT (seasat's venv has no netCDF reader)."""
    g = str(intf / pair / "corr.grd")
    gmt = lambda *a: subprocess.run(["gmt", *a, "--GMT_HISTORY=false"], capture_output=True, check=True).stdout
    w, e, s, n, _, _, dx, dy, nx, ny = gmt("grdinfo", "-C", g).split()[1:11]
    z = np.frombuffer(gmt("grd2xyz", g, "-ZTLf"), dtype="<f4").reshape(int(float(ny)), int(float(nx)))
    return np.linspace(float(n), float(s), z.shape[0])[:, None] * np.ones((1, z.shape[1])), z


cols = {}
for pr in pairs:
    y, z = load(pr)
    edges = np.linspace(y.min(), y.max(), NBIN + 1)
    k = np.clip(np.searchsorted(edges, y, side="right") - 1, 0, NBIN - 1)
    med = []
    for b in range(NBIN):
        v = z[k == b]
        v = v[np.isfinite(v) & (v > 0)]
        med.append(np.percentile(v, 90) if v.size else np.nan)
    cols[pr] = (edges, np.array(med))
edges = cols[pairs[0]][0]
print(f"{T} F{n} corr.grd p90 per azimuth bin (line = bin centre; p90 as qc: water rows pin a median)")
print("  line  " + " ".join(f"{p:>16}" for p in pairs))
for b in range(NBIN):
    print(f"{0.5 * (edges[b] + edges[b + 1]):6.0f}  " + " ".join(f"{cols[p][1][b]:16.3f}" for p in pairs))
