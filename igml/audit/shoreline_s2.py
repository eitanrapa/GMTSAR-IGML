"""The Salton Sea's water line on a Sentinel-2 L2A scene (2026-11-02; user: "Could we do a shoreline of Salton Sea
figure?"): Planetary Computer's STAC (anonymous SAS token), every tile of one acquisition that meets the box, B03 / B08
(10 m) read as windows from the COGs, NDWI = (G - NIR) / (G + NIR) > 0 = water, resampled by nearest onto a lon / lat
grid of STEP degrees. SRTM cannot be the reference: the lake has receded since 2000. Writes an npz (water mask, NDWI,
lon, lat, scene ids).
usage (remote_sensing env, python -I): python -I shoreline_s2.py DATE W S E N STEP OUT.npz   e.g. 2016-10-13"""
import json
import sys
import urllib.request

import numpy as np
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

date, (w, s, e, n), step, out = sys.argv[1], map(float, sys.argv[2:6]), float(sys.argv[6]), sys.argv[7]
api = "https://planetarycomputer.microsoft.com/api"
body = json.dumps({"collections": ["sentinel-2-l2a"], "bbox": [w, s, e, n], "limit": 50,
                   "datetime": f"{date}T00:00:00Z/{date}T23:59:59Z"}).encode()
req = urllib.request.Request(f"{api}/stac/v1/search", body, {"Content-Type": "application/json"})
items = json.load(urllib.request.urlopen(req, timeout=120))["features"]
token = json.load(urllib.request.urlopen(f"{api}/sas/v1/token/sentinel-2-l2a", timeout=60))["token"]
lon = np.arange(w + step / 2, e, step)
lat = np.arange(n - step / 2, s, -step)
LON, LAT = np.meshgrid(lon, lat)
ndwi = np.full(LON.shape, np.nan, np.float32)
used = []
for it in sorted(items, key=lambda f: f["properties"].get("eo:cloud_cover", 100)):
    with rasterio.open(it["assets"]["B03"]["href"] + "?" + token) as g, \
            rasterio.open(it["assets"]["B08"]["href"] + "?" + token) as r:
        tr = Transformer.from_crs("EPSG:4326", g.crs, always_xy=True)
        X, Y = tr.transform(LON, LAT)
        inside = (X >= g.bounds.left) & (X < g.bounds.right) & (Y > g.bounds.bottom) & (Y <= g.bounds.top) & np.isnan(ndwi)
        if not inside.any():
            continue
        win = from_bounds(X[inside].min() - 20, Y[inside].min() - 20, X[inside].max() + 20, Y[inside].max() + 20,
                          g.transform).round_offsets().round_lengths()
        G = g.read(1, window=win, boundless=True, fill_value=0).astype(np.float32)
        R = r.read(1, window=win, boundless=True, fill_value=0).astype(np.float32)
        t = g.window_transform(win)
        col = np.floor((X - t.c) / t.a).astype(int)
        row = np.floor((Y - t.f) / t.e).astype(int)
        ok = inside & (row >= 0) & (row < G.shape[0]) & (col >= 0) & (col < G.shape[1])
        gg, rr = G[row[ok], col[ok]], R[row[ok], col[ok]]
        good = (gg > 0) & (rr > 0)
        v = np.full(gg.shape, np.nan, np.float32)
        v[good] = (gg[good] - rr[good]) / (gg[good] + rr[good])
        sub = ndwi[ok]
        sub[np.isnan(sub)] = v[np.isnan(sub)]
        ndwi[ok] = sub
        used.append(it["id"])
        print(f"{it['id']}: cloud {it['properties'].get('eo:cloud_cover'):.2f} %, {int(ok.sum())} cells", flush=True)
water = ndwi > 0
print(f"grid {ndwi.shape}, {np.isfinite(ndwi).mean():.3f} covered, water {np.nanmean(water[np.isfinite(ndwi)]):.3f}")
np.savez_compressed(out, water=water, ndwi=ndwi, lon=lon, lat=lat, scenes=np.array(used), date=date)
print(f"-> {out}")
