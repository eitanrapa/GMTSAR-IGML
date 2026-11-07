"""Google Earth files of GMTSAR's geolocation error at the conv and merge_swath stages (2026-11-06; user: "Google Earth files where I can see
the geolocation problem, in the merged and conv stages ... I am going to show this at a presentation later so GMTSAR can fix upstream").

T173's pair 2016-10-06 x 10-30, stock GMTSAR (3dd43d3: conv, merge_swath as installed) against the fixed conv + merge_swath, through the SAME
trans.dat and GMTSAR's own lookup grids (proj_ra2ll.csh's raln.grd / ralt.grd, made by /zdata/datasets/IGML/work/geoloc_kml_seasat.sh):
every output pixel goes back through those lookups to fractional radar coordinates and takes the grid's value there (bilinear between the DECLARED
nodes), so a layer is exactly where a GMTSAR user's geocoding puts the cells, without the 100-m blockmedian that hides a sub-cell shift. Offsets are measured
in radar coordinates (cross-correlation of the two grids, as audit/relshift.py) and turned into metres through the lookups' local gradient.
A TRAP avoided: cross-correlating nearest-cell tiles of two grids on one lattice locks the peak to whole cells (F1 read 5-11 m for a true 69 m).
Reflector positions come from display_amp (peak ~7 x the ground), not coherence, whose blob sits on a coherent background.

Writes <OUT>/gmtsar_geolocation_t173_conv_stage.kmz and ..._merged_stage.kmz: per site a radio folder (one layer at a time: click between
Stock and Fixed to flicker), the Sentinel-2 water line of 2016-10-13 at the Salton Sea, the surveyed Pinon Flat reflectors with where each layer
puts them. Prints every measured number. usage (insar_sim env, PYTHONPATH=IGML): python make_geolocation_kmz.py OUT_DIR [--check PNG_DIR]"""
import html
import json
import sys
import zipfile
from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import label, map_coordinates

from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

OUT = Path(sys.argv[1])
CHECK = Path(sys.argv[sys.argv.index("--check") + 1]) if "--check" in sys.argv else None
P = "2016279_2016303"
W = Path("/seasat/IGML/work/geoloc_kml")
STOCK = Path("/seasat/IGML/lab_originals_20261029/T173/product_T173")
FIXED = Path("/seasat/seasatHD3/InSAR_Processing/T173/product_T173")
SRC = Path(__file__).resolve().parent
KX, KY = 111320.0, 110950.0
PAIR = "2016-10-06 x 2016-10-30"
REFL = {"D1": (33.612253, -116.457893), "D2": (33.607373, -116.451836)}     # UCSD's survey, WGS-84, apex (topex.ucsd.edu/geodynamics/pinon/pinon_corner.pdf)
# name, west, east, south, north, pixel [deg], conv-stage subswath, views: (label, lon, lat, range m)
SITES = [
    ("Salton Sea", -116.16, -115.52, 33.05, 33.56, 0.00025, 2, [("West shore, Desert Shores", -116.040, 33.438, 5000), ("East shore, north of Bombay Beach", -115.840, 33.442, 5000)]),
    ("Pinon Flat corner reflectors", -116.50, -116.41, 33.585, 33.635, 0.0001, 3, [("Reflectors D1 and D2", -116.4549, 33.6098, 1800)]),
    ("Lake Havasu", -114.47, -114.09, 34.27, 34.62, 0.00025, 1, [("Parker Dam end", -114.15, 34.31, 9000), ("Lake Havasu City", -114.36, 34.47, 9000)]),
    ("San Diego coast", -117.32, -117.05, 32.60, 32.92, 0.00025, 3, [("Mission Bay", -117.235, 32.775, 7000), ("Coronado, San Diego Bay", -117.17, 32.67, 9000)]),
    # across the F1 / F2 seam: merged stage only (no conv-stage subswath holds both sides)
    ("Mexicali seam F1-F2", -115.63, -115.46, 32.60, 32.74, 0.0002, None, [("Where the seam crosses the canal", -115.553, 32.664, 5000)]),
    ("Imperial Valley seam F1-F2", -115.55, -115.40, 32.95, 33.08, 0.0002, None, [("Where the seam crosses the river", -115.475, 33.020, 5000)]),
]
GRAYS = np.array([[0, 0, 0]] + [[int(v)] * 3 for v in np.linspace(10, 255, 255)], dtype=np.uint8)


def lookup(d):
    return read_grd(d / "raln.grd"), read_grd(d / "ralt.grd")


def at(g, xq, yq):
    """A lookup grid at fractional radar coordinates, bilinear (grdtrack -nl)."""
    dx, dy = g.x[1] - g.x[0], g.y[1] - g.y[0]
    return map_coordinates(g.z, [(yq - g.y[0]) / dy, (xq - g.x[0]) / dx], order=1, mode="nearest")


class Layer:
    """One product: its radar grid and where GMTSAR's lookups put each cell's declared node."""

    def __init__(self, grd, d, log=False):
        g = read_grd(grd)
        self.x, self.y, self.z = g.x, g.y, g.z.astype(np.float32)
        self.v = np.log10(np.where(self.z > 0, self.z, np.nan)) if log else self.z
        self.ln, self.lt = lookup(d)
        X, Y = np.meshgrid(self.x, self.y)
        self.lon, self.lat = at(self.ln, X, Y), at(self.lt, X, Y)

    def tile(self, w, e, s, n, step, thin=1):
        """North-up lon / lat image: each pixel through the lookups back to fractional radar coordinates (a linear interpolation of the cells'
        positions, every thin-th cell), then the grid's value there, bilinear; no value within a NaN cell's reach."""
        lon = np.arange(w + step / 2, e, step)
        lat = np.arange(n - step / 2, s, -step)
        m = (self.lon > w - 0.01) & (self.lon < e + 0.01) & (self.lat > s - 0.01) & (self.lat < n + 0.01)
        out = np.full((lat.size, lon.size), np.nan, dtype=np.float32)
        if m.sum() < 10:
            return lon, lat, out
        r, c = np.nonzero(m)
        if thin > 1:
            k = (r % thin == 0) & (c % thin == 0)
            r, c = r[k], c[k]
        f = LinearNDInterpolator(np.c_[self.lon[r, c], self.lat[r, c]], np.c_[r, c].astype(np.float64))
        LO, LA = np.meshgrid(lon, lat)
        rc = f(LO, LA)
        ok = np.isfinite(rc[..., 0])
        good = np.isfinite(self.v)
        val = map_coordinates(np.nan_to_num(self.v), [rc[ok, 0], rc[ok, 1]], order=1)
        wgt = map_coordinates(good.astype(np.float32), [rc[ok, 0], rc[ok, 1]], order=1)
        out[ok] = np.where(wgt > 0.75, val / np.maximum(wgt, 1e-6), np.nan)
        return lon, lat, out

    def blob(self, lat0, lon0, R=7, near=4):
        """Where this layer puts a corner reflector: the half-maximum centroid of its blob near the surveyed point, in radar coordinates, through the
        lookups. Returns (lon, lat)."""
        c = np.cos(np.radians(lat0))
        dd = ((self.lon - lon0) * KX * c) ** 2 + ((self.lat - lat0) * KY) ** 2
        i0, j0 = np.unravel_index(np.nanargmin(dd), dd.shape)
        win = self.z[i0 - R:i0 + R + 1, j0 - R:j0 + R + 1].astype(np.float64)
        nr = np.zeros(win.shape, bool); nr[R - near:R + near + 1, R - near:R + near + 1] = True
        ip, jp = np.unravel_index(np.nanargmax(np.where(nr, win, np.nan)), win.shape)
        bg = np.nanmedian(win)
        lab, _ = label(win >= bg + 0.5 * (win[ip, jp] - bg))
        sel = lab == lab[ip, jp]
        wt = win - bg
        ii, jj = np.nonzero(sel)
        fi = (ii * wt[sel]).sum() / wt[sel].sum() + i0 - R
        fj = (jj * wt[sel]).sum() / wt[sel].sum() + j0 - R
        xq, yq = self.x[0] + fj * (self.x[1] - self.x[0]), self.y[0] + fi * (self.y[1] - self.y[0])
        return float(at(self.ln, np.array([xq]), np.array([yq]))[0]), float(at(self.lt, np.array([xq]), np.array([yq]))[0])


def shift_cells(a, b, zoom=0.02, half=1.5):
    """(rows, columns, peak): a's content sits that much later than b's in the grids. FFT cross-correlation, then its true maximum on a fine lattice around
    the whole-cell peak (the matrix DFT of Guizar-Sicairos et al. 2008): a parabola through three lags reads +5..+17 samples for the same +16."""
    m = np.isfinite(a) & np.isfinite(b)
    A = np.where(m, a - a[m].mean(), 0.0); B = np.where(m, b - b[m].mean(), 0.0)
    S = np.fft.fft2(A) * np.conj(np.fft.fft2(B))
    cc = np.fft.ifft2(S).real
    i, j = np.unravel_index(np.argmax(cc), cc.shape)
    n, c = A.shape
    i0 = i - n if i > n // 2 else i
    j0 = j - c if j > c // 2 else j
    u = np.arange(-half, half + zoom / 2, zoom)
    kr, kc = np.fft.fftfreq(n), np.fft.fftfreq(c)
    R = np.exp(2j * np.pi * np.outer(i0 + u, kr))
    C = np.exp(2j * np.pi * np.outer(kc, j0 + u))
    fine = (R @ S @ C).real / (n * c)
    p, q = np.unravel_index(np.argmax(fine), fine.shape)
    return i0 + u[p], j0 + u[q], fine[p, q] / np.sqrt((A * A).sum() * (B * B).sum())


# How far stock GMTSAR's declared nodes sit from their content, T173, in samples and lines. conv: +16, +3 by construction (the cell holds a filter centred on sample 32 m,
# line 6 r; its node is declared at 32 m + 16, 6 r + 3). merge_swath adds its whole-cell placement: F2 one column (M1, +32) plus the residual of the fixed grid's origin (1, 4
# lines = -2 mod 6), F3 the residual of (13, 2). Measured over each whole subswath they come back within 1-2 samples (printed below).
EXACT = {"conv": {1: (16, 3), 2: (16, 3), 3: (16, 3)}, "merged": {1: (16, 3), 2: (49, 1), 3: (29, 5)}}


def to_metres(B, er, ea, r0, r1, c0, c1):
    """(east, north) metres of er samples, ea lines at the middle of B's rows r0:r1, columns c0:c1, through the lookups."""
    xm, ym = np.array([0.5 * (B.x[c0] + B.x[c1 - 1])]), np.array([0.5 * (B.y[r0] + B.y[r1 - 1])])
    hx, hy = 256.0, 128.0
    dx_lon = (at(B.ln, xm + hx, ym) - at(B.ln, xm - hx, ym))[0] / (2 * hx)
    dy_lon = (at(B.ln, xm, ym + hy) - at(B.ln, xm, ym - hy))[0] / (2 * hy)
    dx_lat = (at(B.lt, xm + hx, ym) - at(B.lt, xm - hx, ym))[0] / (2 * hx)
    dy_lat = (at(B.lt, xm, ym + hy) - at(B.lt, xm, ym - hy))[0] / (2 * hy)
    lat0 = at(B.lt, xm, ym)[0]
    return (dx_lon * er + dy_lon * ea) * KX * np.cos(np.radians(lat0)), (dx_lat * er + dy_lat * ea) * KY


def window_cells(B, w, e, s, n, cols=None):
    m = (B.lon > w) & (B.lon < e) & (B.lat > s) & (B.lat < n)
    if cols:
        m[:, :cols[0]] = False
        m[:, cols[1]:] = False
    r, c = np.nonzero(m)
    return r.min(), r.max() + 1, c.min(), c.max() + 1


def exact_m(stage, sw, B, w, e, s, n, cols=None):
    """The exact offset of subswath sw at this stage, in metres in this window: (east, north, samples, lines)."""
    er, ea = EXACT["merged" if stage == "merged" else "conv"][sw]
    return (*to_metres(B, er, ea, *window_cells(B, w, e, s, n, cols)), er, ea)


def radar_offset(A, B, w, e, s, n, cols=None):
    """How far layer A puts the ground from where layer B puts it, in the window (and the grid's columns cols): (east m, north m, samples, lines, peak)."""
    r0, r1, c0, c1 = window_cells(B, w, e, s, n, cols)
    fi, fj, peak = shift_cells(A.v[r0:r1, c0:c1].astype(np.float64), B.v[r0:r1, c0:c1].astype(np.float64))
    er = A.x[0] - B.x[0] + fj * (B.x[1] - B.x[0])
    ea = A.y[0] - B.y[0] + fi * (B.y[1] - B.y[0])
    return (*to_metres(B, er, ea, r0, r1, c0, c1), er, ea, peak)


def png(z, lo, hi):
    """North-up tile -> a 64-level gray paletted PNG, no data transparent."""
    q = np.clip((z - lo) / (hi - lo), 0, 1)
    idx = np.where(np.isfinite(z), 1 + np.rint(np.nan_to_num(q) * 62).astype(np.uint8), 0).astype(np.uint8)
    im = Image.fromarray(idx, mode="P")
    pal = [0, 0, 0] + [int(v) for g in np.linspace(12, 255, 63) for v in (g, g, g)]
    im.putpalette(pal + [0] * (768 - len(pal)))
    b = BytesIO()
    im.save(b, format="PNG", optimize=True, transparency=0)
    return b.getvalue()


def lookat(lon, lat, rng):
    return f"<LookAt><longitude>{lon}</longitude><latitude>{lat}</latitude><altitude>0</altitude><heading>0</heading><tilt>0</tilt><range>{rng}</range></LookAt>"


def overlay(name, href, lon, lat, step, vis, desc, view):
    w, e = lon[0] - step / 2, lon[-1] + step / 2
    n, s = lat[0] + step / 2, lat[-1] - step / 2
    return (f"<GroundOverlay><name>{html.escape(name)}</name><visibility>{int(vis)}</visibility><description><![CDATA[{desc}]]></description>{view}"
            f"<color>e6ffffff</color><drawOrder>1</drawOrder><Icon><href>{href}</href></Icon>"
            f"<LatLonBox><north>{n:.7f}</north><south>{s:.7f}</south><east>{e:.7f}</east><west>{w:.7f}</west></LatLonBox></GroundOverlay>\n")


def line(name, pts, color, width, vis=1, desc=""):
    co = " ".join(f"{a:.6f},{b:.6f},0" for a, b in pts)
    return (f"<Placemark><name>{html.escape(name)}</name><visibility>{int(vis)}</visibility><description><![CDATA[{desc}]]></description>"
            f"<Style><LineStyle><color>{color}</color><width>{width}</width></LineStyle></Style>"
            f"<LineString><tessellate>1</tessellate><altitudeMode>clampToGround</altitudeMode><coordinates>{co}</coordinates></LineString></Placemark>\n")


def point(name, lon, lat, color, desc="", scale=0.7, view="", vis=1):
    return (f"<Placemark><name>{html.escape(name)}</name><visibility>{int(vis)}</visibility><description><![CDATA[{desc}]]></description>{view}"
            f"<Style><IconStyle><color>{color}</color><scale>{scale}</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href></Icon></IconStyle>"
            f"<LabelStyle><scale>0.7</scale></LabelStyle></Style><Point><coordinates>{lon:.7f},{lat:.7f},0</coordinates></Point></Placemark>\n")


def water_line():
    """The Salton Sea's Sentinel-2 water line, 2016-10-13 (shoreline_s2.py + shoreline_clean.py): the lake mask's longest contour."""
    f = SRC / "s2_waterline_20161013.npy"
    if not f.exists():
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        L = np.load(SRC / "s2_lake_20161013.npz")
        fig, ax = plt.subplots()
        cs = ax.contour(L["lon"], L["lat"], L["lake"].astype(np.float32), levels=[0.5])
        segs = cs.allsegs[0]
        plt.close(fig)
        seg = max(segs, key=len)
        np.save(f, seg[::3])
    return np.load(f)


# ---- the layers
layers = {}
for b, R in (("stock", STOCK), ("fixed", FIXED)):
    layers["merged", b, "corr"] = Layer(W / f"merged_{b}/corr_ra.grd", W / f"merged_{b}")
    layers["merged", b, "amp"] = Layer(W / f"merged_{b}/display_amp_ra.grd", W / f"merged_{b}", log=True)
    for n in (1, 2, 3):
        layers[f"F{n}", b, "corr"] = Layer(R / f"F{n}/intf_all/{P}/corr.grd", W / f"conv_{b}_F{n}")
        layers[f"F{n}", b, "amp"] = Layer(R / f"F{n}/intf_all/{P}/display_amp.grd", W / f"conv_{b}_F{n}", log=True)
z = layers["merged", "stock", "amp"].v
amp_lo, amp_hi = np.nanpercentile(z[np.isfinite(z)], [2, 99.5])
print(f"amplitude stretch log10 {amp_lo:.2f}..{amp_hi:.2f}", flush=True)
measured = {}


def tiles(stage, w, e, s, n, step, thin=1):
    out = {}
    for b in ("stock", "fixed"):
        lon, lat, c = layers[stage, b, "corr"].tile(w, e, s, n, step, thin)
        _, _, a = layers[stage, b, "amp"].tile(w, e, s, n, step, thin)
        out[b] = (c, a)
    return lon, lat, out


# ---- reflectors: where each layer puts them
refl = {}
for stage in ("F3", "merged"):
    for b in ("stock", "fixed"):
        for k, (la, lo) in REFL.items():
            ln, lt = layers[stage, b, "amp"].blob(la, lo)
            e, n_ = (ln - lo) * KX * np.cos(np.radians(la)), (lt - la) * KY
            refl[stage, b, k] = (ln, lt, e, n_, float(np.hypot(e, n_)))
            print(f"reflector {k} {stage} {b}: east {e:+.0f} m, north {n_:+.0f} m, {np.hypot(e, n_):.0f} m from the survey", flush=True)

# ---- the merged grid's seams, and each subswath's offset over the whole track
SEAMS = [int(c) for c in merge_seams(STOCK / "merge" / P, tuple(read_grd(STOCK / f"F{n}/intf_all/{P}/phasefilt.grd").z.shape[1] for n in (1, 2, 3)))]
print("merged seams at columns", SEAMS, flush=True)
BIG = (-180.0, 180.0, -90.0, 90.0)
whole_offsets = {}
for stage, cols, nm in (("F1", None, "conv F1"), ("F2", None, "conv F2"), ("F3", None, "conv F3"), ("merged", (20, SEAMS[0] - 20), "merged F1"),
                        ("merged", (SEAMS[0] + 20, SEAMS[1] - 20), "merged F2"), ("merged", (SEAMS[1] + 20, layers["merged", "stock", "amp"].z.shape[1] - 20), "merged F3")):
    A, B = layers[stage, "stock", "amp"], layers[stage, "fixed", "amp"]
    ae, an, asm, al, ap = radar_offset(A, B, *BIG, cols=cols)
    ce, cn, cs, cl, cp = radar_offset(layers[stage, "stock", "corr"], layers[stage, "fixed", "corr"], *BIG, cols=cols)
    whole_offsets[nm] = {"amp": {"east_m": round(ae), "north_m": round(an), "total_m": round(float(np.hypot(ae, an))), "samples": round(float(asm), 1), "lines": round(float(al), 1), "peak": round(float(ap), 3)},
                         "corr": {"east_m": round(ce), "north_m": round(cn), "total_m": round(float(np.hypot(ce, cn))), "samples": round(float(cs), 1), "lines": round(float(cl), 1), "peak": round(float(cp), 3)}}
    sw_ = int(nm[-1])
    xe, xn, xs, xl = exact_m("merged" if stage == "merged" else "conv", sw_, B, *BIG, cols=cols)
    whole_offsets[nm]["exact"] = {"samples": xs, "lines": xl, "east_m": round(xe), "north_m": round(xn), "total_m": round(float(np.hypot(xe, xn)))}
    print(f"whole {nm}: exact {xs:+d} samples, {xl:+d} lines = east {xe:+.0f}, north {xn:+.0f} = {np.hypot(xe, xn):.0f} m; measured amplitude {asm:+.1f}, {al:+.1f} (peak {ap:.3f}), coherence {cs:+.1f}, {cl:+.1f} (peak {cp:.3f})", flush=True)


def kmz(path, title, intro, stages_for_site, whole=None, extra=""):
    files, body = {}, []
    for name, w, e, s, n, step, sw, views in SITES:
        stage = stages_for_site(sw)
        if stage is None:
            continue
        lon, lat, T = tiles(stage, w, e, s, n, step)
        lat0 = 0.5 * (s + n)
        if sw is None:                  # a seam site: each side has its own offset
            AS, BS = layers[stage, "stock", "amp"], layers[stage, "fixed", "amp"]
            m1, m2 = radar_offset(AS, BS, w, e, s, n, cols=(0, SEAMS[0] - 4)), radar_offset(AS, BS, w, e, s, n, cols=(SEAMS[0] + 4, SEAMS[1] - 4))
            e1, n1, s1, l1 = exact_m("merged", 1, BS, w, e, s, n, cols=(0, SEAMS[0] - 4))
            e2, n2, s2_, l2 = exact_m("merged", 2, BS, w, e, s, n, cols=(SEAMS[0] + 4, SEAMS[1] - 4))
            print(f"{path.stem} {name}: F1 side exact {np.hypot(e1, n1):.0f} m ({s1:+d}, {l1:+d}), measured {np.hypot(m1[0], m1[1]):.0f} m ({m1[2]:+.1f}, {m1[3]:+.1f}, peak {m1[4]:.2f}); "
                  f"F2 side exact {np.hypot(e2, n2):.0f} m ({s2_:+d}, {l2:+d}), measured {np.hypot(m2[0], m2[1]):.0f} m ({m2[2]:+.1f}, {m2[3]:+.1f}, peak {m2[4]:.2f})", flush=True)
            measured[path.stem, name] = {"stage": stage, "F1_side_m": round(float(np.hypot(e1, n1))), "F2_side_m": round(float(np.hypot(e2, n2))), "step_m": round(float(np.hypot(e2 - e1, n2 - n1)))}
            tag = name.split()[0].lower()
            desc = (f"<b>{name}</b>, merged stage.<br/>Stock GMTSAR moves the east side of the seam (F1) {np.hypot(e1, n1):.0f} m and the west side (F2) {np.hypot(e2, n2):.0f} m, so features step by about "
                    f"<b>{np.hypot(e2 - e1, n2 - n1):.0f} m</b> where they cross it. Look at the white seam line. Click between the layers to flicker.")
            c0 = lookat(0.5 * (w + e), lat0, max(e - w, n - s) * 111000 * 1.15)
            items = []
            for b, bname in (("stock", "Stock GMTSAR"), ("fixed", "Fixed")):
                for prod, pname, lo, hi, q in (("amp", "amplitude", amp_lo, amp_hi, 1), ("corr", "coherence", 0.0, 0.8, 0)):
                    href = f"files/{tag}_{b}_{prod}.png"
                    files[href] = png(T[b][q], lo, hi)
                    items.append(overlay(f"{bname}: {pname}", href, lon, lat, step, False, desc, c0))
            items.append("<Folder><name>Off (imagery only)</name></Folder>\n")
            vw = "".join(point(lb, lo_, la_, "00ffffff", "", 0.01, lookat(lo_, la_, rg), vis=0) for lb, lo_, la_, rg in views)
            body.append(f"<Folder><name>{html.escape(name)}: flicker these</name><open>1</open><description><![CDATA[{desc}]]></description>{c0}"
                        f"<Style><ListStyle><listItemType>radioFolder</listItemType></ListStyle></Style>\n{''.join(items)}</Folder>\n"
                        f"<Folder><name>{html.escape(name)}: go to</name><open>0</open>{vw}</Folder>\n")
            continue
        ce, cn, cs, cl, cp = radar_offset(layers[stage, "stock", "corr"], layers[stage, "fixed", "corr"], w, e, s, n)
        ae, an, asm, al, ap = radar_offset(layers[stage, "stock", "amp"], layers[stage, "fixed", "amp"], w, e, s, n)
        cov = np.mean(np.isfinite(T["fixed"][0])) * 100
        print(f"{path.stem} {name} ({stage}): stock from fixed, coherence east {ce:+.0f}, north {cn:+.0f} = {np.hypot(ce, cn):.0f} m ({cs:+.1f} samples, {cl:+.1f} lines, peak {cp:.2f}); "
              f"amplitude east {ae:+.0f}, north {an:+.0f} = {np.hypot(ae, an):.0f} m ({asm:+.1f}, {al:+.1f}, peak {ap:.2f}); {cov:.0f} % of the window has data", flush=True)
        xe, xn, xs, xl = exact_m(stage, sw, layers[stage, "fixed", "amp"], w, e, s, n)
        print(f"   exact there: {xs:+d} samples, {xl:+d} lines = east {xe:+.0f}, north {xn:+.0f} = {np.hypot(xe, xn):.0f} m", flush=True)
        measured[path.stem, name] = {"stage": stage, "east_m": round(xe), "north_m": round(xn), "total_m": round(float(np.hypot(xe, xn))), "samples": xs, "lines": xl,
                                     "measured_amplitude_m": round(float(np.hypot(ae, an))), "measured_coherence_m": round(float(np.hypot(ce, cn)))}
        tag = name.split()[0].lower()
        desc = (f"<b>{name}</b>, {stage if stage == 'merged' else 'subswath ' + stage + ', straight from conv'}.<br/>Stock GMTSAR's labels are {xs:+d} samples, {xl:+d} lines from their content here, "
                f"so it puts the ground <b>{np.hypot(xe, xn):.0f} m</b> from the fixed layer (east {xe:+.0f} m, north {xn:+.0f} m).<br/>"
                f"Click between the layers to flicker. Same trans.dat, same GMTSAR lookup tables for both.")
        c0 = lookat(0.5 * (w + e), lat0, max(e - w, n - s) * 111000 * 1.15)
        items = []
        for k, (b, bname) in enumerate((("stock", "Stock GMTSAR"), ("fixed", "Fixed"))):
            for prod, pname, lo, hi, q in (("corr", "coherence", 0.0, 0.8, 0), ("amp", "amplitude", amp_lo, amp_hi, 1)):
                href = f"files/{tag}_{b}_{prod}.png"
                files[href] = png(T[b][q], lo, hi)
                items.append(overlay(f"{bname}: {pname}", href, lon, lat, step, b == "stock" and prod == "corr", desc, c0))
        items.append("<Folder><name>Off (imagery only)</name></Folder>\n")
        vw = "".join(point(lb, lo_, la_, "00ffffff", "", 0.01, lookat(lo_, la_, rg), vis=0) for lb, lo_, la_, rg in views)
        body.append(f"<Folder><name>{html.escape(name)}: flicker these</name><open>1</open><description><![CDATA[{desc}]]></description>{c0}"
                    f"<Style><ListStyle><listItemType>radioFolder</listItemType></ListStyle></Style>\n{''.join(items)}</Folder>\n"
                    f"<Folder><name>{html.escape(name)}: go to</name><open>0</open>{vw}</Folder>\n")
        if CHECK:
            CHECK.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(CHECK / f"{path.stem}_{tag}.npz", lon=lon, lat=lat, sc=T["stock"][0], fc=T["fixed"][0], sa=T["stock"][1], fa=T["fixed"][1])
    if whole:
        w, e, s, n, step = whole
        lon, lat, T = tiles("merged", w, e, s, n, step, thin=6)
        items = []
        dsc = "The whole track. Each subswath moves by its own amount: zoom in to the seam lines."
        for b, bname in (("stock", "Stock GMTSAR"), ("fixed", "Fixed")):
            for prod, pname, lo, hi, q in (("corr", "coherence", 0.0, 0.8, 0), ("amp", "amplitude", amp_lo, amp_hi, 1)):
                href = f"files/track_{b}_{prod}.png"
                files[href] = png(T[b][q], lo, hi)
                items.append(overlay(f"{bname}: {pname}", href, lon, lat, step, False, dsc, lookat(0.5 * (w + e), 0.5 * (s + n), 560000)))
        items.append("<Folder><name>Off (imagery only)</name></Folder>\n")
        body.append(f"<Folder><name>Whole track: flicker these</name><open>0</open><description><![CDATA[{dsc}]]></description>"
                    f"<Style><ListStyle><listItemType>radioFolder</listItemType></ListStyle></Style>\n{''.join(items)}</Folder>\n")
    doc = (f'<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2"><Document>'
           f"<name>{html.escape(title)}</name><open>1</open><description><![CDATA[{intro}]]></description>\n" + extra + "".join(body) + "</Document></kml>\n")
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
        zf.writestr("doc.kml", doc, compress_type=zipfile.ZIP_DEFLATED)
        for k, v in files.items():
            zf.writestr(k, v)
    print(f"-> {path} {path.stat().st_size / 1e6:.0f} MB, {len(files)} images", flush=True)


wl = water_line()
print(f"Sentinel-2 water line: {len(wl)} vertices", flush=True)
s2 = line("Salton Sea water line, Sentinel-2 2016-10-13", [(a, b) for a, b in wl], "ff00ffff", 2, 1,
          "NDWI &gt; 0 on Sentinel-2 L2A, 2016-10-13, between the pair's two dates. Use this rather than today's imagery: the lake has receded since 2016.")


def reflectors(stage, slabel):
    out = ""
    for k, (la, lo) in REFL.items():
        out += point(f"{k} surveyed", lo, la, "ff00ffff", "UCSD's survey, WGS-84, apex.")
        for b, bn, col in (("stock", "stock", "ff3467eb"), ("fixed", "fixed", "ff7aaf1b")):
            ln, lt, e, n_, d = refl[stage, b, k]
            out += point(f"{k} {bn}: {d:.0f} m", ln, lt, col, f"Where the {bn} {slabel} puts reflector {k}: {d:.0f} m from the survey (east {e:+.0f} m, north {n_:+.0f} m). "
                         "Half-maximum centroid of its coherence blob, through GMTSAR's lookup tables.")
    return f"<Folder><name>Pinon Flat reflectors: survey vs each layer</name><open>1</open>{out}</Folder>\n"


how = ("T173 (descending), " + PAIR + ". <b>Stock GMTSAR</b> = conv and merge_swath as distributed. <b>Fixed</b> = centred conv cells and exact subswath placement. Both go through the same "
       "trans.dat and GMTSAR's own proj_ra2ll.csh lookup tables, and each radar cell is drawn at its declared position with no smoothing.<br/><br/>"
       "Each site is a radio list: click <b>Stock</b> then <b>Fixed</b> to flicker. Coherence shows water edges; amplitude shows towns and roads.")
kmz(OUT / "gmtsar_geolocation_t173_conv_stage.kmz", "GMTSAR geolocation: conv stage (T173)",
    "<b>Subswath grids straight from conv, before any merging.</b> conv's cell holds a filter centred on sample 32 m, line 6 r, but its node is declared at 32 m + 16, 6 r + 3.<br/><br/>" + how,
    lambda sw: f"F{sw}" if sw else None, extra=f"<Folder><name>Reference lines</name><open>1</open>{s2}</Folder>\n" + reflectors("F3", "F3 subswath grid"))
# the subswath seams of the merged grid, along its columns, through the lookups
M = layers["merged", "fixed", "corr"]
seam_lines = ""
for c, nm in zip(SEAMS, ("F1 / F2 seam", "F2 / F3 seam")):
    rows = np.arange(0, M.z.shape[0], 20)
    seam_lines += line(nm, list(zip(M.lon[rows, c], M.lat[rows, c])), "ffffffff", 1, 1, "Where the merged grid changes subswath. In stock GMTSAR the two sides are displaced differently, so features step here.")
kmz(OUT / "gmtsar_geolocation_t173_merged_stage.kmz", "GMTSAR geolocation: merged stage (T173)",
    "<b>Grids after merge_swath.</b> On top of conv's half cell, merge_swath rounds each subswath to whole cells and places F2 one column too far, so F1, F2 and F3 each sit off by a different amount.<br/><br/>" + how,
    lambda sw: "merged", whole=(-117.46, -114.0, 31.82, 35.89, 0.0008),
    extra=f"<Folder><name>Reference lines</name><open>1</open>{s2}{seam_lines}</Folder>\n" + reflectors("merged", "merged grid"))
(OUT / "measured_offsets.json").write_text(json.dumps({"subswaths": whole_offsets, "sites": {f"{a} | {b}": v for (a, b), v in measured.items()},
                                                       "reflectors": {f"{a} {b} {c}": {"east_m": round(v[2]), "north_m": round(v[3]), "total_m": round(v[4])} for (a, b, c), v in refl.items()}}, indent=1))
print("done", flush=True)
