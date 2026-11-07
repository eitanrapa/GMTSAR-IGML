"""The lake alone from shoreline_s2.py's NDWI mask (2026-11-02): the largest connected water body, holes filled, specks
under ~3 cells opened / closed away -- one clean water line -- plus the signed distance to it (m; + on land) and the
direction of the shore normal (bearing toward land, from the distance's gradient). Prints where shores facing east /
west (a north-south water line, which shows an east-west offset) are densest. Writes an npz next to the input.
usage (insar_sim env): python shoreline_clean.py S2.npz OUT.npz"""
import sys

import numpy as np
from scipy import ndimage as ndi

S2 = np.load(sys.argv[1])
lon, lat, w = S2["lon"], S2["lat"], S2["water"]
step = lon[1] - lon[0]
kx, ky = 111320.0 * np.cos(np.radians(lat.mean())), 110950.0
lab, n = ndi.label(w)
lake = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)
lake = ndi.binary_fill_holes(ndi.binary_closing(ndi.binary_opening(lake, iterations=2), iterations=3))
d_land = ndi.distance_transform_edt(~lake, sampling=(step * ky, step * kx))
d_water = ndi.distance_transform_edt(lake, sampling=(step * ky, step * kx))
dist = np.where(lake, -d_water, d_land).astype(np.float32)
sm = ndi.gaussian_filter(dist, 15)
gn, ge = np.gradient(sm)
bearing = (np.degrees(np.arctan2(ge / kx, -gn / ky)) % 360).astype(np.float32)
print(f"lake {lake.mean():.3f} of the box ({lake.sum() * step * step * kx * ky / 1e6:.0f} km2), from {n} water bodies")
near = np.abs(dist) < 300
for face, lo, hi in (("west-facing (land to the W)", 247.5, 292.5), ("east-facing (land to the E)", 67.5, 112.5)):
    m = near & (bearing >= lo) & (bearing < hi)
    H, ye, xe = np.histogram2d(np.broadcast_to(lat[:, None], m.shape)[m], np.broadcast_to(lon[None, :], m.shape)[m],
                               bins=[np.arange(lat.min(), lat.max() + 0.02, 0.02), np.arange(lon.min(), lon.max() + 0.02, 0.02)])
    top = np.argsort(H.ravel())[::-1][:4]
    print(face + ":", ", ".join(f"({xe[j % H.shape[1]] + 0.01:.3f}, {ye[j // H.shape[1]] + 0.01:.3f}) {int(H.ravel()[j])}" for j in top))
np.savez_compressed(sys.argv[2], lake=lake, dist=dist, bearing=bearing, lon=lon, lat=lat, date=S2["date"])
print(f"-> {sys.argv[2]}")
