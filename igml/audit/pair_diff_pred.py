"""Is the as-stored ours - lab (pair_diff_data.py) the placement fix acting on the fringes? Per subswath: the difference predicted from our own
unwrapped phase and the measured offset, fi dU/drow + fj dU/dcol, against the actual difference (median removed; |d| < 5 mm, so whole cycles
stay out). Read only. usage (insar_sim env, PYTHONPATH=IGML): python pair_diff_pred.py T DATA.npz"""
import json
import sys
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

T, data = sys.argv[1], sys.argv[2]
z = np.load(data)
ARCH = Path(f"/seasat/IGML/lab_originals_20261029/{T}/product_{T}/merge")
OURS = Path(f"/seasat/seasatHD3/InSAR_Processing/{T}/product_{T}/merge")
MM = 55.465763 / (4 * np.pi)
for p in (str(x) for x in z["pairs"]):
    m = json.loads(str(z[f"{p}_meta"]))
    U = read_grd(OURS / p / "unwrap.grd").z.astype(np.float64)
    L = read_grd(ARCH / p / "unwrap.grd").z.astype(np.float64)
    U[U == 0] = np.nan; L[L == 0] = np.nan
    r, c = min(U.shape[0], L.shape[0]), min(U.shape[1], L.shape[1])
    U, L = U[:r, :c], L[:r, :c]
    d = (U - L); d = (d - np.nanmedian(d)) * MM
    gy, gx = np.gradient(U)
    s = [0] + m["seams"] + [c]
    for k, n in enumerate(("F1", "F2", "F3")):
        a, b = s[k] + 4, s[k + 1] - 4
        fi, fj = m["offsets"][n]
        pred = (fi * gy[:, a:b] + fj * gx[:, a:b]) * MM
        act = d[:, a:b]
        ok = np.isfinite(pred) & np.isfinite(act) & (np.abs(act) < 5)
        rr = np.corrcoef(pred[ok], act[ok])[0, 1]
        left = act[ok] - pred[ok]
        print(f"{p} {n}: corr(predicted, actual) {rr:+.3f}; sd actual {act[ok].std():.2f} mm, sd after taking the prediction out {left.std():.2f} mm", flush=True)
