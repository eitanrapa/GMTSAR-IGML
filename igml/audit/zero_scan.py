"""Where merged phasefilt.grd is EXACTLY zero (2026-10-28, GMTSAR audit I1: phasefilt.c leaves its patch-grid edges at
0 phase, mask 1, corr valid). Per pair: the rows and the columns that are >= 50 % exact zeros among their finite
pixels -- the azimuth-end strip is expected; a zero COLUMN inside the data (a subswath edge landing inside the merged
grid) would be a zero stripe along the whole track. Read only; GMT reads, numpy counts.

usage (seasat): python zero_scan.py [N_PAIRS_PER_TRACK [T ...]]  (the lab's latest pairs + ours; every track by default)
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
PAT = re.compile(r"^\d{7}_\d{7}$")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2
TRACKS = sys.argv[2:]


def load(path):
    f = subprocess.run(["gmt", "grdinfo", "-C", str(path)], capture_output=True, text=True, check=True).stdout.split()
    nx, ny = int(f[9]), int(f[10])
    raw = subprocess.run(["gmt", "grd2xyz", str(path), "-ZTLf"], capture_output=True, check=True).stdout
    return np.frombuffer(raw, dtype=np.float32).reshape(ny, nx)          # row 0 = the top = the azimuth END


def runs(idx):
    """[a-b] runs of consecutive indices."""
    out, s = [], None
    for k, i in enumerate(idx):
        if s is None:
            s = i
        if k + 1 == len(idx) or idx[k + 1] != i + 1:
            out.append(f"{s}-{i}" if i > s else f"{s}")
            s = None
    return ",".join(out)


for T, v in VOL.items():
    if TRACKS and T not in TRACKS:
        continue
    roots = [("lab", Path(f"/zdata/datasets/seasatHD{v}/InSAR_Processing/{T}/product_{T}/merge")),
             ("ours", Path(f"/zdata/datasets/IGML/InSAR_Processing/{T}/product_{T}/merge"))]
    for who, root in roots:
        if not root.is_dir():
            continue
        # complete pairs only (unwrap.grd, as lab_dates): 10-28, the lab's T042 2025-26 pairs (252) are all-NaN, no unwrap
        prs = sorted(d for d in root.iterdir() if PAT.match(d.name) and (d / "phasefilt.grd").exists()
                     and (d / "unwrap.grd").exists())[-N:]
        for d in prs:
            p = load(d / "phasefilt.grd")
            fin = np.isfinite(p)
            if not fin.any():                       # an all-NaN pair (10-28: crashed the first run on T042)
                print(f"{T} {who:4s} {d.name} {p.shape[0]}x{p.shape[1]}  no finite pixels", flush=True)
                continue
            z = fin & (p == 0)
            rf = z.sum(1) / np.maximum(fin.sum(1), 1)
            cf = z.sum(0) / np.maximum(fin.sum(0), 1)
            zr = np.flatnonzero((rf >= 0.5) & (fin.sum(1) > 0.05 * p.shape[1]))
            zc = np.flatnonzero((cf >= 0.5) & (fin.sum(0) > 0.05 * p.shape[0]))
            c0, c1 = np.flatnonzero(fin.any(0))[[0, -1]]
            inside = [c for c in zc if c0 + 3 < c < c1 - 3]
            print(f"{T} {who:4s} {d.name} {p.shape[0]}x{p.shape[1]}  zero rows [{runs(list(zr))}]  zero cols "
                  f"[{runs(list(zc))}]  inside the data: {len(inside)}  all-zero share {z.sum() / max(fin.sum(), 1):.4f}",
                  flush=True)
