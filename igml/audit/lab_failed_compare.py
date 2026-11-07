"""The lab's merge/failed pairs of T173 against our remade ones (2026-11-06), over the Mac's mount, read only: median and p90 coherence, the share of the grid with phase and unwrapped phase.
usage (insar_sim env, PYTHONPATH=IGML): python lab_failed_compare.py"""
import os
import re
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

F = Path("/seasat/seasatHD3/InSAR_Processing/T173/product_T173/merge/failed")
O = Path("/seasat/IGML/InSAR_Processing/T173/product_T173/merge")
rows = []
for p in sorted(x for x in os.listdir(F) if re.match(r"^\d{7}_\d{7}$", x)):
    def info(d):
        c = read_grd(d / "corr.grd").z.astype(np.float32)
        ph = read_grd(d / "phasefilt.grd").z
        u = d / "unwrap.grd"
        un = np.nan
        if u.exists():
            z = read_grd(u).z
            un = float(np.mean(np.isfinite(z) & (z != 0)) * 100)
        return float(np.nanmedian(c)), float(np.nanpercentile(c, 90)), float(np.mean(np.isfinite(ph) & (ph != 0)) * 100), un
    a, b = info(F / p), info(O / p)
    rows.append((p, a, b))
    print(f"{p} lab corr median {a[0]:.3f} p90 {a[1]:.3f} phase {a[2]:5.1f} % unwrap {a[3]:5.1f} % | ours corr median {b[0]:.3f} p90 {b[1]:.3f} phase {b[2]:5.1f} % unwrap {b[3]:5.1f} %", flush=True)
A = np.array([r[1] for r in rows]); B = np.array([r[2] for r in rows])
print(f"{len(rows)} pairs. Median over pairs: corr p90 lab {np.median(A[:, 1]):.3f}, ours {np.median(B[:, 1]):.3f}; phase lab {np.median(A[:, 2]):.1f} %, ours {np.median(B[:, 2]):.1f} %; "
      f"lab has unwrap on {np.sum(np.isfinite(A[:, 3]))}, ours on {np.sum(np.isfinite(B[:, 3]))}; ours' lowest unwrap share {np.nanmin(B[:, 3]):.1f} %; corr p90 ours - lab: min {np.min(B[:, 1] - A[:, 1]):+.3f}, median {np.median(B[:, 1] - A[:, 1]):+.3f}, max {np.max(B[:, 1] - A[:, 1]):+.3f}")
