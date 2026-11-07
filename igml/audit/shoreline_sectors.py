"""The shoreline offset per stretch of shore, as a check of shoreline_shift.py's 2-D result (2026-11-02): DEM nodes near
the Sentinel-2 water line binned by the direction of the shore normal (8 sectors of 45 deg; the land fraction's
gradient smoothed at 300 m), per sector the 1-D shift along its mean normal that best correlates the optical land
fraction (Gaussian 100 m) with the coherence (+-400 m, 5-m steps), then the 2-D offset d that fits every sector's shift
(s_k = d . n_k + e, weighted by nodes; e = a common edge bias along every normal) and the sectors' rms misfit --
the method's precision.
usage (insar_sim env): python shoreline_sectors.py S2.npz COH.npz"""
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

S2, CO = np.load(sys.argv[1]), np.load(sys.argv[2])
lon, lat, water = S2["lon"], S2["lat"], S2["water"]
step = lon[1] - lon[0]
kx, ky = 111320.0 * np.cos(np.radians(lat.mean())), 110950.0
land = (~water).astype(np.float32)
sm = gaussian_filter(land, (100 / (step * ky), 100 / (step * kx)))
big = gaussian_filter(land, (300 / (step * ky), 300 / (step * kx)))
gn, ge = np.gradient(big)                       # rows go south, so north = -d/drow
gn, ge = -gn / (step * ky), ge / (step * kx)
wide = uniform_filter(land, (int(1000 / (step * ky)), int(1000 / (step * kx))))
nl, nt = CO["node_lon"], CO["node_lat"]


def sample(a, dlon=0.0, dlat=0.0):
    j = np.rint((nl - dlon - lon[0]) / step).astype(np.int64)
    i = np.rint((lat[0] - (nt - dlat)) / step).astype(np.int64)
    ok = (i >= 0) & (i < a.shape[0]) & (j >= 0) & (j < a.shape[1])
    v = np.full(nl.shape, np.nan, np.float32)
    v[ok] = a[i[ok], j[ok]]
    return v


w = sample(wide)
band = (w > 0.05) & (w < 0.95)
ang = np.degrees(np.arctan2(sample(ge), sample(gn))) % 360      # bearing of the normal, pointing to land
sector = ((ang + 22.5) // 45).astype(int) % 8
names = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
for name in CO["names"]:
    coh = CO[f"{name}_coh"]
    rows = []
    for k in range(8):
        m = band & np.isfinite(coh) & (sector == k)
        if m.sum() < 1500:
            continue
        a0 = np.radians(np.degrees(np.arctan2(np.nanmean(sample(ge)[m]), np.nanmean(sample(gn)[m]))))
        ne, nn = np.sin(a0), np.cos(a0)
        c = coh[m]

        def score(sh):
            v = sample(sm, sh * ne / kx, sh * nn / ky)[m]
            ok = np.isfinite(v)
            p, q = c[ok] - c[ok].mean(), v[ok] - v[ok].mean()
            return float((p * q).sum() / np.sqrt((p * p).sum() * (q * q).sum()))

        sc = [(score(sh), sh) for sh in range(-400, 401, 5)]
        r, sh = max(sc)
        rows.append((names[k], int(m.sum()), np.degrees(a0) % 360, sh, r, ne, nn))
    # s_k = d . n_k + e: e is a common EDGE bias along every normal (the coherence edge seaward / landward of the optical
    # water line, e.g. wet mudflats: NDWI water, coherent to the radar) -- isotropic, so a full circle of sectors
    # separates it from the geolocation offset d (a dipole)
    A = np.array([[x[5], x[6], 1.0] for x in rows]); b = np.array([x[3] for x in rows], float)
    W = np.sqrt(np.array([x[1] for x in rows], float))
    d = np.linalg.lstsq(A * W[:, None], b * W, rcond=None)[0]
    mis = b - A @ d
    rms = np.sqrt(np.average(mis ** 2, weights=W ** 2))
    cov = np.linalg.inv((A * W[:, None]).T @ (A * W[:, None])) * rms ** 2 * W.size / max(W.size - 3, 1) * np.mean(W ** 2)
    se = np.sqrt(np.diag(cov))
    print(f"{name}: offset east {d[0]:+.0f} +- {se[0]:.0f} m, north {d[1]:+.0f} +- {se[1]:.0f} m ({np.hypot(d[0], d[1]):.0f} m); "
          f"edge bias {d[2]:+.0f} m along the land normal; sectors' rms misfit {rms:.0f} m")
    for (nm, cnt, az, sh, r, *_), mm in zip(rows, mis):
        print(f"    shore facing {nm:>2} (normal {az:5.1f} deg, {cnt:6d} nodes): shift along the normal {sh:+4d} m, r {r:.3f}, "
              f"misfit {mm:+.0f}")
