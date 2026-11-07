"""WHERE triangle closure fails by whole cycles (2026-10-28, GMTSAR audit): in phasefilt's zero strip (I1), at burst
boundaries (A1: each date keeps its own burst-switch lines, ~1342 lines = ~223.6 merged rows apart), or elsewhere.
Triangles as closure_check.py (our stored pairs with the lab's legs; or the lab's alone with LAB_ONLY=1). Per triangle:
the share of pixels a cycle or more off the closure's median, in the strip rows (exactly zero phasefilt in any leg)
vs outside; and, outside the strip, the off pixels' rows folded at the burst period per subswath third, as the peak
bin's excess over uniform. Read only; GMT reads, numpy counts.

usage (seasat): python closure_rows.py T [N]
"""
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

T, N = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 12
LAB_ONLY = os.environ.get("LAB_ONLY") == "1"
VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
ours = Path(f"/zdata/datasets/IGML/InSAR_Processing/{T}/product_{T}/merge")
lab = Path(f"/zdata/datasets/seasatHD{VOL[T]}/InSAR_Processing/{T}/product_{T}/merge")
PAT = re.compile(r"^\d{7}_\d{7}$")
PERIOD, BINS = 1342 / 6, 32


def load(path):
    f = subprocess.run(["gmt", "grdinfo", "-C", str(path)], capture_output=True, text=True, check=True).stdout.split()
    nx, ny = int(f[9]), int(f[10])
    raw = subprocess.run(["gmt", "grd2xyz", str(path), "-ZTLf"], capture_output=True, check=True).stdout
    return np.frombuffer(raw, dtype=np.float32).reshape(ny, nx)


have = {}
for root in ((lab,) if LAB_ONLY else (lab, ours)):
    if root.is_dir():
        for d in root.iterdir():
            if PAT.match(d.name) and (d / "unwrap.grd").exists() and (d / "phasefilt.grd").exists():
                have.setdefault(d.name, d)
test = sorted(have) if LAB_ONLY else sorted(d.name for d in ours.iterdir() if PAT.match(d.name))
random.seed(1)
random.shuffle(test)
done, agg = 0, []
for p in test:
    i, j = p.split("_")
    legs = None
    for k in sorted({x for q in have for x in q.split("_")}):
        if i < k < j and f"{i}_{k}" in have and f"{k}_{j}" in have:
            legs = (have[f"{i}_{k}"], have[f"{k}_{j}"], have[p])            # (i,k) + (k,j) - (i,j)
            break
    if not legs:
        continue
    u = [load(x / "unwrap.grd") for x in legs]
    z = [load(x / "phasefilt.grd") == 0 for x in legs]
    if len({a.shape for a in u}) != 1:
        continue
    c = u[0] + u[1] - u[2]
    ok = np.isfinite(c)
    off = ok & (np.abs(c - np.median(c[ok])) > np.pi)
    strip = z[0] | z[1] | z[2]          # exact zeros in some leg, per PIXEL: each subswath ends at its own row, so a
    S = ok & strip                      # row test over the whole width misses the strip (found 10-28 on T042: F1 only)
    O = ok & ~strip
    nx = c.shape[1]
    exc = []
    for a, b in ((0, nx // 3), (nx // 3, 2 * nx // 3), (2 * nx // 3, nx)):   # ~F1 / F2 / F3
        rows = np.nonzero(off[:, a:b] & O[:, a:b])[0]
        if rows.size < 50:
            exc.append(float("nan"))
            continue
        h = np.histogram((rows % PERIOD) / PERIOD, bins=BINS, range=(0, 1))[0]
        exc.append(h.max() / h.mean())                                  # 1 = uniform; a burst boundary -> >> 1
    who = "".join("L" if str(x).startswith(str(lab)) else "o" for x in legs)
    r = (100 * off[S].mean() if S.any() else float("nan"), 100 * off[O].mean(), S.sum() / ok.sum())
    agg.append(r + tuple(exc))
    print(f"{p} [{who}] off >= 1 cycle: in the strip {r[0]:.1f} %, outside {r[1]:.2f} % (strip = {100 * r[2]:.2f} % of "
          f"pixels; of all off pixels {100 * off[S].sum() / max(off.sum(), 1):.0f} % in it); burst-period peak / mean "
          f"per third {' '.join('%.2f' % e for e in exc)}", flush=True)
    done += 1
    if done >= N:
        break
a = np.array(agg)
print(f"{T} median of {len(a)}: in the strip {np.nanmedian(a[:, 0]):.1f} %, outside {np.nanmedian(a[:, 1]):.2f} %; "
      f"burst-period peak / mean {np.nanmedian(a[:, 3:], axis=0).round(2).tolist()} (uniform ~{1 + 2 / np.sqrt(BINS):.1f}"
      f" at most by chance for ~1000 pixels)")
