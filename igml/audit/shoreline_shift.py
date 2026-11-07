"""The geolocation offset of each product from the Salton Sea's shoreline (2026-11-02): at trans.dat's DEM nodes within
~1 km of the Sentinel-2 water line (shoreline_s2.py), the shift (east, north) of the optical land fraction -- smoothed
like the coherence estimate (Gaussian, SIGMA m) -- that best correlates with each product's coherence
(shoreline_geocode.py), searched over +-500 m (20-m steps, then 5-m), per subswath (merged-grid column against the
seams). A product that shows a feature at its true position + d has its best shift at d (Piñon Flat convention:
product - surveyed).
usage (insar_sim env): python shoreline_shift.py S2.npz COH.npz F2_START F3_START [STITCH] [OUT.json]"""
import json
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

S2, CO = np.load(sys.argv[1]), np.load(sys.argv[2])
f2, f3 = int(sys.argv[3]), int(sys.argv[4])
stitch = int(sys.argv[5]) if len(sys.argv) > 5 else 13
out = sys.argv[6] if len(sys.argv) > 6 else None
SIGMA = 100.0
lon, lat, water = S2["lon"], S2["lat"], S2["water"]
step = lon[1] - lon[0]
lat0 = lat.mean()
kx, ky = 111320.0 * np.cos(np.radians(lat0)), 110950.0
land = (~water).astype(np.float32)
sm = gaussian_filter(land, (SIGMA / (step * ky), SIGMA / (step * kx)))
wide = uniform_filter(land, (int(1000 / (step * ky)), int(1000 / (step * kx))))
nl, nt = CO["node_lon"], CO["node_lat"]


def sample(a, dlon, dlat):
    j = np.rint((nl - dlon - lon[0]) / step).astype(np.int64)
    i = np.rint((lat[0] - (nt - dlat)) / step).astype(np.int64)
    ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1])
    v = np.full(nl.shape, np.nan, np.float32)
    v[ok] = a[i[ok], j[ok]]
    return v


band = sample(wide, 0, 0)
band = (band > 0.05) & (band < 0.95)
res = {}
for name in CO["names"]:
    coh, col = CO[f"{name}_coh"], CO[f"{name}_col"]
    sub = np.where(col < 0, 0, np.where(col < f2 + stitch, 1, np.where(col < f3 + stitch, 2, 3)))
    for sw in (1, 2, 3):
        m = band & np.isfinite(coh) & (sub == sw)
        if m.sum() < 500:
            continue
        c = coh[m] - coh[m].mean()

        def score(de, dn):
            v = sample(sm, de / kx, dn / ky)[m]
            ok = np.isfinite(v)
            a, b = c[ok] - c[ok].mean(), v[ok] - v[ok].mean()
            return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))

        best = max(((score(de, dn), de, dn) for de in range(-500, 501, 20) for dn in range(-500, 501, 20)))
        best = max(((score(de, dn), de, dn) for de in range(best[1] - 20, best[1] + 21, 5)
                    for dn in range(best[2] - 20, best[2] + 21, 5)))
        r0 = score(0, 0)
        res[f"{name}_F{sw}"] = dict(east=best[1], north=best[2], dist=float(np.hypot(best[1], best[2])),
                                    r_best=best[0], r_zero=r0, nodes=int(m.sum()))
        print(f"{name} F{sw}: product - true = east {best[1]:+4d} m, north {best[2]:+4d} m ({np.hypot(best[1], best[2]):.0f} m)"
              f"; r {best[0]:.3f} (at zero shift {r0:.3f}), {int(m.sum())} nodes near the shore", flush=True)
if out:
    json.dump(res, open(out, "w"), indent=1)
