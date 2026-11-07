"""Data for the Rosamond corner-reflector geolocation figure (2026-10-29; user: "Can I see the geolocation problem with
the reflectors?"). T071 (descending) F2, the lab's stack: on a 2-m lon / lat grid over the array, every pixel's radar
position comes from merge/trans.dat (the lookup our Geocoder and GMTSAR's proj_ra2ll use), interpolated bilinearly
between its DEM nodes, and the pixel shows
  slc   the super-master SLC's amplitude at that position (full resolution: where the reflectors really are)
  lab   F2's decimated amplitude (display_amp.grd of pair 2018285_2019004) in the cell the lab's merged grid puts
        there: merged column floor(r / 32) - 624, row floor(a / 6) - 77 (merge_swath's placement, value-for-value)
  ours  the same with M1 corrected (F2 one column back: floor(r / 32) - 623)
  fix   the proposed correction (a preview, applied nowhere): the cell whose CONTENT is centred nearest, in F2's own
        frame, round(r_F2 / 32), round(a_F2 / 6)
Read only (the lab's files over the Mac's NFS mount); writes one npz.

usage (igml env): python fig_cr_geoloc_data.py OUT.npz
"""
import sys

import numpy as np

from igml.real.grd import read_grd

T = "/Volumes/seasat2/seasatHD1/InSAR_Processing/T071/product_T071"
PAIR = "2018285_2019004"
LON0, LON1, LAT0, LAT1, D = -118.0920, -118.0300, 34.7995, 34.8090, 0.00002


def prm(p):
    return {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(p)) if s}


if __name__ == "__main__":
    # trans.dat nodes around the array -> r, a on the DEM's own lon / lat grid
    t = np.memmap(f"{T}/merge/trans.dat", dtype="<f8", mode="r").reshape(-1, 5)
    keep = []
    for s in range(0, t.shape[0], 4_000_000):
        b = np.asarray(t[s:s + 4_000_000])
        m = (b[:, 3] > LON0 - 0.01) & (b[:, 3] < LON1 + 0.01) & (b[:, 4] > LAT0 - 0.01) & (b[:, 4] < LAT1 + 0.01)
        keep.append(b[m])
    b = np.concatenate(keep)
    lons, lats = np.unique(np.round(b[:, 3], 7)), np.unique(np.round(b[:, 4], 7))
    R = np.full((lats.size, lons.size), np.nan)
    A = np.full_like(R, np.nan)
    i = np.searchsorted(lats, np.round(b[:, 4], 7)); j = np.searchsorted(lons, np.round(b[:, 3], 7))
    R[i, j], A[i, j] = b[:, 0], b[:, 1]
    print(f"{len(b)} trans.dat nodes, DEM step {np.median(np.diff(lons)) * 3600:.2f} x {np.median(np.diff(lats)) * 3600:.2f} arcsec")
    # bilinear onto the fine grid
    glon = np.arange(LON0, LON1 + D / 2, D)
    glat = np.arange(LAT1, LAT0 - D / 2, -D)
    fx = np.interp(glon, lons, np.arange(lons.size)); fy = np.interp(glat, lats, np.arange(lats.size))
    x0, y0 = np.floor(fx).astype(int), np.floor(fy).astype(int)
    wx, wy = (fx - x0)[None, :], (fy - y0)[:, None]
    X0, Y0 = x0[None, :], y0[:, None]

    def bil(Z):
        return ((1 - wy) * (1 - wx) * Z[Y0, X0] + (1 - wy) * wx * Z[Y0, X0 + 1]
                + wy * (1 - wx) * Z[Y0 + 1, X0] + wy * wx * Z[Y0 + 1, X0 + 1])

    Rm, Am = bil(R), bil(A)
    # merged -> F2 frame (continuous): near-range and start-time differences from the super-master PRMs
    P1, P2 = prm(f"{T}/F1/topo/master.PRM"), prm(f"{T}/F2/topo/master.PRM")
    dr = 299792458.0 / (2 * float(P1["rng_samp_rate"]))
    C2 = (float(P2["near_range"]) - float(P1["near_range"])) / dr
    D2 = (float(P2["clock_start"]) - float(prm(f"{T}/merge/{PAIR}/supermaster.PRM")["clock_start"])) * 86400 * float(P1["PRF"])
    r2, a2 = Rm - C2, Am - D2
    print(f"F2 offsets in the merged frame: {C2:.3f} samples, {D2:.3f} lines")
    # SLC amplitude (nearest sample)
    W = int(prm(f"{T}/F2/raw/S1_20190318_ALL_F2.PRM")["bytes_per_line"]) // 4
    slc = np.memmap(f"{T}/F2/raw/S1_20190318_ALL_F2.SLC", dtype="<i2", mode="r")
    slc = slc[:(slc.size // (2 * W)) * 2 * W].reshape(-1, W, 2)
    ra, rr = np.rint(a2).astype(int), np.rint(r2).astype(int)
    a_lo, a_hi, r_lo, r_hi = ra.min(), ra.max() + 1, rr.min(), rr.max() + 1
    win = slc[a_lo:a_hi, r_lo:r_hi].astype(np.float32)
    amp = np.hypot(win[..., 0], win[..., 1])
    slc_img = amp[ra - a_lo, rr - r_lo]
    # decimated amplitude, three placements
    G = read_grd(f"{T}/F2/intf_all/{PAIR}/display_amp.grd")
    z = G.z.astype(np.float64)
    cells = {"lab": (np.floor(Rm / 32).astype(int) - 624, np.floor(Am / 6).astype(int) - 77),
             "ours": (np.floor(Rm / 32).astype(int) - 623, np.floor(Am / 6).astype(int) - 77),
             "fix": (np.floor((r2 + 16) / 32).astype(int), np.floor((a2 + 3) / 6).astype(int))}
    out = {"lon": glon, "lat": glat, "slc": slc_img, "C2": C2, "D2": D2}
    for k, (cc, rw) in cells.items():
        out[k] = z[rw, cc]
    np.savez_compressed(sys.argv[1], **out)
    print("wrote", sys.argv[1], {k: v.shape for k, v in out.items() if hasattr(v, "shape")})
