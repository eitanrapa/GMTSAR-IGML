"""A check and a preview of the Google Earth files, drawn FROM the files (doc.kml's LatLonBox + the PNG inside each KMZ), so what it shows is what Google
Earth gets (2026-11-06). Rows: conv stage, merged stage. Columns: the Salton Sea's west shore, stock / fixed coherence under the Sentinel-2 water line of 2016-10-13; the
Pinon Flat reflectors, stock / fixed amplitude with the surveyed points. Offsets from measured_offsets.json.
usage (insar_sim env): python preview_geolocation_kmz.py KMZ_DIR OUT.png"""
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from PIL import Image

D, out = Path(sys.argv[1]), sys.argv[2]
NS = {"k": "http://www.opengis.net/kml/2.2"}
M = json.loads((D / "measured_offsets.json").read_text())
SRC = Path(__file__).resolve().parent
WL = np.load(SRC / "s2_waterline_20161013.npy")
REFL = {"D1": (33.612253, -116.457893), "D2": (33.607373, -116.451836)}
INK, MUTED, SURVEY, WATER = "#1f1f1f", "#6b6a66", "#2a78d6", "#eb6834"
FILES = {"conv": "gmtsar_geolocation_t173_conv_stage", "merged": "gmtsar_geolocation_t173_merged_stage"}
ROWS = (("conv", "conv stage (subswath grid)"), ("merged", "merged stage"))
SITE = {"salton": ((-116.075, -115.995, 33.405, 33.470), "corr"), "pinon": ((-116.475, -116.435, 33.595, 33.625), "amp")}


def layer(zf, tag, build, prod):
    root = ET.fromstring(zf.read("doc.kml"))
    for g in root.iter("{http://www.opengis.net/kml/2.2}GroundOverlay"):
        href = g.find("k:Icon/k:href", NS).text
        if href == f"files/{tag}_{build}_{prod}.png":
            b = g.find("k:LatLonBox", NS)
            box = [float(b.find(f"k:{x}", NS).text) for x in ("west", "east", "south", "north")]
            im = Image.open(BytesIO(zf.read(href))).convert("RGBA")
            return np.asarray(im), box
    raise KeyError((tag, build, prod))


fig, ax = plt.subplots(2, 4, figsize=(15, 8.6), layout="constrained")
for r, (st, stn) in enumerate(ROWS):
    with zipfile.ZipFile(D / f"{FILES[st]}.kmz") as zf:
        for c, (tag, b, bn) in enumerate((("salton", "stock", "Stock GMTSAR"), ("salton", "fixed", "Fixed"), ("pinon", "stock", "Stock GMTSAR"), ("pinon", "fixed", "Fixed"))):
            win, prod = SITE[tag]
            im, box = layer(zf, tag, b, prod)
            a = ax[r, c]
            a.set_facecolor("#c9c8c2")
            a.imshow(im, extent=box, interpolation="nearest", aspect=1 / np.cos(np.radians(0.5 * (win[2] + win[3]))))
            a.set_xlim(win[0], win[1]); a.set_ylim(win[2], win[3])
            a.xaxis.set_major_locator(MaxNLocator(4)); a.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            a.tick_params(labelsize=7, colors=MUTED, length=2)
            for sp in a.spines.values():
                sp.set_color("#e6e6e6")
            if tag == "salton":
                a.plot(WL[:, 0], WL[:, 1], color=WATER, lw=1.2)
                name = "Salton Sea"
                k = [v for key, v in M["sites"].items() if key.startswith(FILES[st]) and name in key][0]
                ttl = f"{bn}: west shore, coherence"
                if b == "stock":
                    ttl += f"\n{k['total_m']} m from fixed (east {k['east_m']:+d}, north {k['north_m']:+d})"
            else:
                for cid, (la, lo) in REFL.items():
                    a.plot(lo, la, "+", color=SURVEY, ms=11, mew=1.6)
                    a.annotate(cid, (lo, la), xytext=(6, 6), textcoords="offset points", fontsize=8, color=SURVEY)
                s = "F3" if st == "conv" else "merged"
                d1, d2 = M["reflectors"][f"{s} {b} D1"], M["reflectors"][f"{s} {b} D2"]
                ttl = f"{bn}: Piñon Flat reflectors, amplitude\nD1 {d1['total_m']} m, D2 {d2['total_m']} m from the survey"
            a.set_title(ttl, fontsize=9.5, color=INK)
            if c == 0:
                a.set_ylabel(stn, fontsize=10.5, color=INK)
fig.suptitle("T173, 2016-10-06 × 10-30, read back from the Google Earth files. Orange line: Sentinel-2 water line, 2016-10-13. Blue crosses: surveyed reflectors", fontsize=11.5, color=INK)
fig.savefig(out, dpi=120, bbox_inches="tight")
print("->", out)
# Pinon Flat alone: replaces pinon_flat_geolocation_t173.png, whose 81 / 62 m had the search radius cutting the stock product's blob
if len(sys.argv) > 3:
    fig, ax = plt.subplots(2, 2, figsize=(9.5, 8.2), layout="constrained")
    win = SITE["pinon"][0]
    for r, (st, stn) in enumerate(ROWS):
        with zipfile.ZipFile(D / f"{FILES[st]}.kmz") as zf:
            for c, (b, bn) in enumerate((("stock", "Stock GMTSAR"), ("fixed", "Fixed"))):
                im, box = layer(zf, "pinon", b, "amp")
                a = ax[r, c]
                a.imshow(im, extent=box, interpolation="nearest", aspect=1 / np.cos(np.radians(0.5 * (win[2] + win[3]))))
                a.set_xlim(win[0], win[1]); a.set_ylim(win[2], win[3])
                a.xaxis.set_major_locator(MaxNLocator(5)); a.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
                a.tick_params(labelsize=7.5, colors=MUTED, length=2)
                for sp in a.spines.values():
                    sp.set_color("#e6e6e6")
                for cid, (la, lo) in REFL.items():
                    a.plot(lo, la, "+", color=SURVEY, ms=13, mew=1.8)
                    a.annotate(cid, (lo, la), xytext=(7, 7), textcoords="offset points", fontsize=9, color=SURVEY)
                s = "F3" if st == "conv" else "merged"
                d1, d2 = M["reflectors"][f"{s} {b} D1"], M["reflectors"][f"{s} {b} D2"]
                a.set_title(f"{bn}, {stn}\nD1 {d1['total_m']} m, D2 {d2['total_m']} m from the survey", fontsize=10, color=INK)
    fig.suptitle("Piñon Flat corner reflectors, T173 F3, 2016-10-06 × 10-30: amplitude through GMTSAR's geocoding. Blue crosses: UCSD's survey", fontsize=11, color=INK)
    fig.savefig(sys.argv[3], dpi=120, bbox_inches="tight")
    print("->", sys.argv[3])
