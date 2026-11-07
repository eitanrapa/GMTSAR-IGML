"""The F2 fix on the lab's real grids (2026-10-28, audit M1): in a merged phasefilt.grd every column is one subswath
grid's column, value for value (merge_swath copies; the unwrap scripts only landmask it in place), so the seams from
merge_log + the subswath widths and the one-column move of igml.real.gmtsar.shift_f2 can be checked exactly: per
block, the share of the merged column's finite values found in the subswath column it should hold -- F1 col jj, F2 col
jj - (N1 - ovl12) [after the move; jj - (N1 - ovl12) - 1 before], F3 col jj - (N1 - ovl12) - (N2 - ovl23). Read only.

usage (seasat): python seam_identity.py T [N_PAIRS]
"""
import random
import re
import subprocess
import sys

import numpy as np

T, NP = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2
VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
P = f"/zdata/datasets/seasatHD{VOL[T]}/InSAR_Processing/{T}/product_{T}"


def load(path):
    f = subprocess.run(["gmt", "grdinfo", "-C", path, "--GMT_HISTORY=false"], capture_output=True, text=True,
                       check=True).stdout.split()
    raw = subprocess.run(["gmt", "grd2xyz", path, "-ZTLf", "--GMT_HISTORY=false"], capture_output=True,
                         check=True).stdout
    return np.frombuffer(raw, dtype=np.float32).reshape(int(f[10]), int(f[9]))


def share(col, ref):
    v = col[np.isfinite(col) & (col != 0)]
    return np.isin(v, ref[np.isfinite(ref)]).mean() if v.size else np.nan


ls = subprocess.run(["ls", f"{P}/merge"], capture_output=True, text=True).stdout.split()
pairs = [p for p in ls if re.match(r"^\d{7}_\d{7}$", p)]
random.seed(2)
done = 0
for p in random.sample(pairs, len(pairs)):
    try:
        t = open(f"{P}/merge/{p}/merge_log").read()
        G = [load(f"{P}/F{n}/intf_all/{p}/phasefilt.grd") for n in (1, 2, 3)]
        M = load(f"{P}/merge/{p}/phasefilt.grd")
    except (OSError, subprocess.CalledProcessError):
        continue
    if not np.isfinite(M).any():
        continue
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    ovl12, ovl23 = int(o[1]), int(o[2])
    n1 = int(re.search(r"n1 = (-?\d+)", t)[1])
    n2 = int(re.search(r"n2 = (-?\d+)", t)[1])
    N1, N2, N3 = (g.shape[1] for g in G)
    off12, off23 = N1 - ovl12, N2 - ovl23
    s12, s23 = off12 + n1, off12 + off23 + n2 - 1
    Mf = M.copy()
    Mf[:, s12 - 1:s23 - 1] = M[:, s12:s23].copy()                  # igml.real.gmtsar.shift_f2
    rng = np.random.default_rng(0)
    pick = lambda a, b, k=12: rng.choice(np.arange(a, b), size=min(k, b - a), replace=False)
    f1 = np.nanmean([share(M[:, j], G[0][:, j]) for j in pick(5, s12 - 1)])
    f2_before = np.nanmean([share(M[:, j], G[1][:, j - off12 - 1]) for j in pick(s12, s23)])
    f2_true = np.nanmean([share(M[:, j], G[1][:, j - off12]) for j in pick(s12, s23)])
    f2_after = np.nanmean([share(Mf[:, j], G[1][:, j - off12]) for j in pick(s12 - 1, s23 - 1)])
    f3 = np.nanmean([share(M[:, j], G[2][:, j - off12 - off23]) for j in pick(s23, M.shape[1] - 2)])
    print(f"{T} {p}: widths {N1} {N2} {N3}, seams s12 {s12} s23 {s23}; merged values found in the subswath column "
          f"F1 {f1:.3f} | F2 at the code's index {f2_before:.3f}, at the true offset {f2_true:.3f}, after shift_f2 "
          f"{f2_after:.3f} | F3 {f3:.3f}", flush=True)
    done += 1
    if done >= NP:
        break
