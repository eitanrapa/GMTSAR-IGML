"""Earthquakes that make a pair COSEISMIC for unwrapping, per track (2026-10-26; user: "let's do smooth except for
coseismic pairs"): igml.catalog's rule -- merged ComCat events whose point-source static displacement reaches COSEIS_MM
somewhere in the track's swath (its outline in s1_acquisitions.json; the extra tracks': extra_tracks/subswaths.json). process_track.unwrap_cmd() runs SNAPHU in DEFO mode
for a pair spanning one of these (by date, both ends included) and in SMOOTH mode otherwise. Writes coseismic_events.json.
10-29: each merged event also gets the date of its LAST member that qualifies on its own ("end_date"; Ridgecrest: the
M7.1 of 2019-07-06, the merged time being the M6.4's of 07-04): T166 has a date between the two (2019-07-05), so a pair
from it spans the M7.1 -- a pair is coseismic if it overlaps [date, end_date] (user 10-29: "smooth unless there is an
earthquake").

usage (insar_sim env, repo root): python -I scratchpad/s1_gapfill/coseismic_events.py <comcat csv>
"""
import json
import sys
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import Polygon

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from igml import catalog  # noqa: E402

COSEIS_MM = 50.0                                              # ~ two C-band LOS fringes
events = catalog.read_comcat(sys.argv[1], min_mag=4.0)
S = json.load(open(REPO / "scratchpad/s1_inventory_20261015/s1_acquisitions.json"))["summary"]
out = {"threshold_mm": COSEIS_MM, "comcat": sys.argv[1], "tracks": {}}
byid = {e.id: e for e in events}


def last_member(ev):
    """The time of the merged event's last member that qualifies on its own (catalog.merge's own test for its first)."""
    own = [byid[i].time for i in ev.id.split("+") if i in byid and catalog.zone_radius_km(byid[i], COSEIS_MM) > 0]
    return max(own) if own else ev.time


# 10-29 (user: "queue the other tracks that were added"): the extra tracks, with no lab stack, by the outline we process
# (extra_tracks/subswaths.json: the California part of their subswaths); the nine tracks keep their hulls
X = json.load(open(REPO / "scratchpad/s1_inventory_20261015/extra_tracks/subswaths.json"))
HULLS = {T: s["hull"] for T, s in S.items()} | {T: s["outline"] for T, s in X.items()}
for T, hull in sorted(HULLS.items()):
    poly = Polygon(hull)
    x0, y0, x1, y1 = poly.bounds
    lon, lat = np.meshgrid(np.arange(x0, x1, 0.02), np.arange(y0, y1, 0.02))
    inside = shapely.contains_xy(poly, lon, lat)
    sel = catalog.select(events, lon, lat, inside, "20141001", "20261101", COSEIS_MM)
    out["tracks"][T] = [{"time": ev.time.isoformat(), "date": ev.time.strftime("%Y%m%d"),
                         "end_time": last_member(ev).isoformat(), "end_date": last_member(ev).strftime("%Y%m%d"),
                         "mw": round(float(ev.mw), 2), "lon": ev.lon, "lat": ev.lat, "depth_km": ev.depth_km,
                         "radius_km": round(r, 1)}
                        for ev, r in sel]
    print(T, [(e["date"], e["mw"], e["radius_km"]) for e in out["tracks"][T]])
json.dump(out, open(REPO / "scratchpad/s1_gapfill/coseismic_events.json", "w"), indent=1)
