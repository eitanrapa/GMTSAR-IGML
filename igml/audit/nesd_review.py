"""Review of a track's network spectral diversity solution (2026-11-01), before pass 2 uses it: per subswath the
residuals (deg of burst-boundary step), the worst pairs, the fallbacks; per date its shift x in each subswath and frame
against the others (orbit timing is common to a date's frames and subswaths, the ionosphere varies along track); and a
figure of x against time by satellite.

usage (Mac, insar_sim env): python nesd_review.py DIR OUT.png
  DIR  holds nesd_F<n>.json and pairs_F<n>.jsonl (OUT/<T>/esd/F<n>/) and state.json (the track's)
"""
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

d, out = Path(sys.argv[1]), sys.argv[2]
st = json.loads((d / "state.json").read_text())
sat = {k: Path(v[0]).name[:3] for k, v in st.get("stitched", {}).items() if v}
X = {}                                           # (n, frame index, date) -> x [lines]
deg = {}
for n in (1, 2, 3):
    f = d / f"nesd_F{n}.json"
    if not f.exists():
        continue
    J = json.loads(f.read_text())
    S, fac = J["summary"], J.get("factors", {})
    sol = J["solution"]
    print(f"F{n}: {S}")
    res = J["residuals"]
    z = np.array([r["z"] for r in res]); r_ = np.array([r["r"] for r in res])
    pj = d / f"pairs_F{n}.jsonl"                 # the subswath's spec_sep x dta (spectral_diversity's), else ~11
    sxd = np.median([json.loads(l)["spectral_spectrationXdta"] for l in pj.read_text().split("\n")
                     if l.strip() and "spectral_spectrationXdta" in l]) if pj.exists() else 11.0
    k = 360.0 * sxd / np.median(list(fac.values()) or [1.0])     # deg of step per line (A6-corrected)
    deg[n] = k
    print(f"  residuals: |r| median {np.median(np.abs(r_)) * k:.2f} deg, p95 {np.percentile(np.abs(r_), 95) * k:.2f}, "
          f"max {np.abs(r_).max() * k:.1f}; |z| > 5: {int((np.abs(z) > 5).sum())} of {len(z)}")
    for r in sorted(res, key=lambda r: -abs(r["z"]))[:6]:
        print(f"    worst {r['a'][3:11]} -> {r['b'][3:11]} ({r['a'][-2:]}): r {r['r'] * k:+.1f} deg, z {r['z']:+.1f}")
    if J["fallback"]:
        print(f"  fallback (estimate against the super-master): {sorted(J['fallback'])}")
    # frames of one date in time order -> index 0, 1
    byd = defaultdict(list)
    for s in sol:
        byd[s[3:11]].append(s)
    for date, ss in byd.items():
        for i, s in enumerate(sorted(ss)):
            X[(n, i, date)] = sol[s]["x"]
dates = sorted({k[2] for k in X})
print("\nconsistency of a date's shift (deg of step):")
for (a, b) in (((1, 0), (1, 1)), ((2, 0), (2, 1)), ((3, 0), (3, 1)), ((1, 0), (2, 0)), ((2, 0), (3, 0)), ((1, 0), (3, 0))):
    pr = [(X[(*a, t)], X[(*b, t)]) for t in dates if (*a, t) in X and (*b, t) in X]
    if pr:
        p = np.array(pr)
        dd = (p[:, 0] - p[:, 1]) * deg.get(a[0], 3960)
        print(f"  F{a[0]} frame {a[1] + 1} vs F{b[0]} frame {b[1] + 1}: {len(p)} dates, corr {np.corrcoef(p.T)[0, 1]:+.3f}, "
              f"difference median {np.median(dd):+.1f}, rms {np.sqrt(np.mean(dd ** 2)):.1f}")
fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
col = {"S1A": "#1f77b4", "S1B": "#d62728", "S1C": "#2ca02c", "S1D": "#9467bd"}
for n in (1, 2, 3):
    a = ax[n - 1]
    for i, mk in ((0, "o"), (1, "x")):
        for sname, c in col.items():
            pts = [(datetime.strptime(t, "%Y%m%d"), X[(n, i, t)] * deg.get(n, 3960)) for t in dates
                   if (n, i, t) in X and sat.get(t) == sname]
            if pts:
                a.plot(*zip(*pts), mk, ms=3.5, color=c, alpha=0.8, label=f"{sname} frame {i + 1}" if n == 1 else None)
    a.axhline(0, color="0.5", lw=0.6)
    a.set_ylabel(f"F{n}: step corrected [deg]")
ax[0].legend(ncol=4, fontsize=8, loc="upper left")
ax[0].set_title("Network spectral diversity: each date's azimuth shift against the super-master (as burst-boundary step)")
fig.tight_layout()
fig.savefig(out, dpi=130)
print(f"figure -> {out}")
