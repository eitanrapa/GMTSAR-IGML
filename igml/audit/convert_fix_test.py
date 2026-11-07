"""The conversion with and without the GMTSAR fixes on the lab's real pairs (2026-10-28, audit I1 / M2 / M1): the same
pairs through igml.real.gmtsar.build twice, fixes off and on, compared cell by cell -- cells the fixes blank (the zero
strip, mask_def) and cells whose value moves (the F2 column), by position across the swath. Reads the lab's stack over
the Mac's read-only NFS mount; writes only OUT.

usage (igml env): python convert_fix_test.py TRACK_ROOT START END MAX_DAYS OUT [STEP]
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

from igml.real.gmtsar import TargetGrid, build

if __name__ == "__main__":          # macOS spawns the geocoding workers by re-importing this file
    root, start, end, maxd, out = sys.argv[1:6]
    step = float(sys.argv[6]) if len(sys.argv) > 6 else 0.002
    lon0, lon1, lat0, lat1 = -119.4167, -113.6087, 31.7878, 36.2142          # the SoCal product box
    grid = TargetGrid(lon0, lat1, step, int(round((lat1 - lat0) / step)), int(round((lon1 - lon0) / step)))
    res = {}
    for tag, fixes in (("off", ()), ("on", None)):
        o = Path(out) / tag
        if not (o / "inputs/ifgramStack.h5").exists():
            build(root, o, grid, start, end, max_days=int(maxd), gnss_tables=None, weather=False, phase_sign=1,
                  workers=4, require_connected=False, missing_bperp="nan", stream=False, fixes=fixes)
        with h5py.File(o / "inputs/ifgramStack.h5", "r") as f:
            res[tag] = f["unwrapPhase"][:]
        res[tag + "_m"] = json.load(open(o / "manifest.json"))
    A, B = res["off"], res["on"]
    pairs = res["on_m"]["pairs"]
    print("fixes:", json.dumps(res["on_m"]["gmtsar_fixes"]))
    for p, a, b in zip(pairs, A, B):
        fa, fb = np.isfinite(a), np.isfinite(b)
        blank = fa & ~fb
        moved = fa & fb & (a != b)
        rows = np.flatnonzero(blank.any(1))
        d = np.abs(a - b)[moved]
        print(f"{p}: finite {fa.sum():,} -> {fb.sum():,}; blanked {blank.sum():,} cells in rows "
              f"{rows.min() if rows.size else '-'}..{rows.max() if rows.size else '-'}; newly finite {(fb & ~fa).sum():,}; "
              f"value moved {moved.sum():,} cells ({100 * moved.sum() / max(fb.sum(), 1):.1f} %), |change| median "
              f"{np.median(d) if d.size else 0:.3f} rad, p99 {np.percentile(d, 99) if d.size else 0:.3f} rad", flush=True)
    np.savez_compressed(Path(out) / "diff_first_pair.npz", off=A[0], on=B[0])
