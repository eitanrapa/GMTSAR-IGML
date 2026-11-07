"""Which of the rebuild's changes reached given pairs, at the level of the merged products (2026-11-06): the lab's original
pair (lab_originals_20261029) against ours (the lab's tree now), over the Mac's mount, read only:
  placement (I2, M1, M6): per subswath, how far the lab's content sits from ours (corr.grd cross-correlation, as relshift.py);
  zero strips (I1, I3): rows at the azimuth ends where phasefilt is exactly 0, rows / columns of data;
  the seams (M1, M6): phase step across the F1/F2 and F2/F3 seams, as eap_seam.py (A2 would show here);
  the phase itself: ours - lab, each subswath with the lab's moved onto ours by the measured offset: rms after the along-track
  profile, the along-track profile's spread (A5's staircase), the across-swath ramp and the offsets between subswaths (I4, orbits);
  unwrapping: pixels a whole cycle apart; pixels each unwrapped, and those where mask_def is NaN (M2).
usage (igml env, PYTHONPATH=IGML): python pair_changes_grd.py T PAIR [PAIR ...]"""
import re
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

LABS = {"T173": "seasatHD3", "T071": "seasatHD1", "T064": "seasatHD2", "T166": "seasatHD4"}
T, pairs = sys.argv[1], sys.argv[2:]
ROOT = "/seasat"
ARCH = Path(f"{ROOT}/IGML/lab_originals_20261029/{T}/product_{T}/merge")
TREE = Path(f"{ROOT}/{LABS[T]}/InSAR_Processing/{T}/product_{T}")
OURS = TREE / "merge"
TWO_PI = 2 * np.pi


def widths_from_log(d):
    t = (d / "merge_log").read_text()
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    f2 = re.search(r"F2 offset (-?[\d.]+) columns", t)
    f3 = re.search(r"F3 offset (-?[\d.]+) columns", t)
    w = int(re.search(r"Size\((\d+)x\d+\)", t)[1])
    ovl12, ovl23, c2, c3 = int(o[1]), int(o[2]), round(float(f2[1])), round(float(f3[1]))
    n1, n2 = ovl12 + c2, ovl23 + c3 - c2
    return n1, n2, w + ovl12 + ovl23 + 2 - n1 - n2


def cmean(z):
    z = z[np.isfinite(z)]
    m = z.mean()
    return np.angle(m), np.abs(m)


def step(ph, c_old, c_new, spacing):
    e = np.exp(1j * ph)
    g = 0.0
    if spacing:
        side = list(range(c_old - 6, c_old)) + list(range(c_new + 1, c_new + 7))
        g, _ = cmean(np.concatenate([e[:, k + 1] * np.conj(e[:, k]) for k in side]))
    a, coh = cmean(e[:, c_new] * np.conj(e[:, c_old]) * np.exp(-1j * spacing * g))
    return np.degrees(a), coh


def shift_cells(a, b):
    """(rows, columns) a's content sits later than b's: cross-correlation, sub-cell by parabola."""
    m = np.isfinite(a) & np.isfinite(b)
    A = np.where(m, a - a[m].mean(), 0.0)
    B = np.where(m, b - b[m].mean(), 0.0)
    cc = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(A) * np.conj(np.fft.rfft2(B)), s=A.shape))
    i, j = np.unravel_index(np.argmax(cc), cc.shape)

    def par(m1, c0, p1):
        d = m1 - 2 * c0 + p1
        return 0.0 if d == 0 else 0.5 * (m1 - p1) / d
    return i + par(cc[i - 1, j], cc[i, j], cc[i + 1, j]) - A.shape[0] // 2, j + par(cc[i, j - 1], cc[i, j], cc[i, j + 1]) - A.shape[1] // 2, \
        cc[i, j] / np.sqrt((A * A).sum() * (B * B).sum())


T0 = np.memmap(OURS / "trans.dat", dtype="<f8", mode="r").reshape(-1, 5)
def metres(er, ea):
    """(east, north) of an offset er samples / ea lines, from trans.dat's local Jacobian mid-track (the track's own)."""
    b = np.asarray(T0[T0.shape[0] // 2 - 300000:T0.shape[0] // 2 + 300000])
    lon0, lat0 = b[:, 3].mean(), b[:, 4].mean()
    e = (b[:, 3] - lon0) * 111320.0 * np.cos(np.radians(lat0))
    n = (b[:, 4] - lat0) * 110950.0
    M = np.c_[e, n, np.ones_like(e)]
    J = np.array([np.linalg.lstsq(M, b[:, 0], rcond=None)[0][:2], np.linalg.lstsq(M, b[:, 1], rcond=None)[0][:2]])
    return np.linalg.inv(J) @ np.array([er, ea])


maskL = read_grd(ARCH / "mask_def.grd")
maskO = read_grd(TREE / "merge/mask_def.grd")
for pair in pairs:
    dl, dO = ARCH / pair, OURS / pair
    Pl, Po = read_grd(dl / "phasefilt.grd"), read_grd(dO / "phasefilt.grd")
    Cl, Co = read_grd(dl / "corr.grd"), read_grd(dO / "corr.grd")
    Ul, Uo = read_grd(dl / "unwrap.grd"), read_grd(dO / "unwrap.grd")
    dx, dy = Po.x[1] - Po.x[0], Po.y[1] - Po.y[0]
    print(f"\n=== {pair}: lab {Pl.z.shape} origin ({Pl.x[0]:g}, {Pl.y[0]:g}); ours {Po.z.shape} origin ({Po.x[0]:g}, {Po.y[0]:g}); cell {dx:g} x {dy:g}", flush=True)
    wl = tuple(read_grd(ARCH.parent / f"F{n}/intf_all/{pair}/phasefilt.grd").z.shape[1] for n in (1, 2, 3))   # the lab's, archived
    s12l, s23l = merge_seams(dl, wl)
    s12o, s23o = merge_seams(dO, widths_from_log(dO))
    print(f" seams: lab columns {s12l}, {s23l}; ours {s12o}, {s23o}", flush=True)
    # zero strips
    for name, P in (("lab", Pl), ("ours", Po)):
        z = P.z
        zr = [int(np.sum(np.nansum(np.abs(z[r]), axis=0) == 0 and np.isfinite(z[r]).sum() > 100)) for r in range(z.shape[0])]
        ends = np.where(np.array(zr) > 0)[0]
        print(f" {name} phasefilt: rows exactly 0 at {list(ends[:3])}..{list(ends[-3:])} ({ends.size} rows); pixels == 0: {int(np.sum(z == 0))}", flush=True)
    # seams' phase steps
    pl = Pl.z.astype(np.float64); pl[pl == 0] = np.nan
    po = Po.z.astype(np.float64)
    a12, c12 = step(po, s12o - 1, s12o, 1); a23, c23 = step(po, s23o - 1, s23o, 1)
    b12, d12 = step(pl, s12l - 1, s12l, 0); b23, d23 = step(pl, s23l - 1, s23l, 2)
    print(f" seam steps [deg]: lab F1/F2 {b12:+.1f} (R {d12:.2f}), F2/F3 {b23:+.1f} (R {d23:.2f}); ours {a12:+.1f} ({c12:.2f}), {a23:+.1f} ({c23:.2f})", flush=True)
    # placement and phase, per subswath
    rows = min(Pl.z.shape[0], Po.z.shape[0])
    offs = []
    prof_all, lev = {}, {}
    for name, (c1l, c2l), (c1o, c2o) in (("F1", (0, s12l), (0, s12o)), ("F2", (s12l, s23l), (s12o, s23o)), ("F3", (s23l, Pl.z.shape[1]), (s23o, Po.z.shape[1]))):
        c1, c2 = max(c1l, c1o) + 12, min(c2l, c2o) - 12
        a, b = Cl.z[:rows, c1:c2].astype(np.float64), Co.z[:rows, c1:c2].astype(np.float64)
        fi, fj, peak = shift_cells(a, b)
        er = Pl.x[0] - Po.x[0] + fj * dx
        ea = Pl.y[0] - Po.y[0] + fi * dy
        en = metres(er, ea)
        offs.append((fi, fj))
        print(f" {name}: lab content {fj:+.2f} columns, {fi:+.2f} rows from ours in the grids; labels {er:+.1f} samples, {ea:+.1f} lines off -> east {en[0]:+.0f} m, north {en[1]:+.0f} m ({np.hypot(*en):.0f} m); corr peak {peak:.3f}", flush=True)
        # ours - lab on ours's cells, the lab's complex phase moved by the measured offset
        rr, cc = np.mgrid[0:rows, c1:c2]
        el = np.exp(1j * np.nan_to_num(pl[:rows]))
        el[~np.isfinite(pl[:rows])] = np.nan
        re_, im_ = (map_coordinates(np.nan_to_num(f(el)), [rr + fi, cc + fj], order=1) for f in (np.real, np.imag))
        okl = map_coordinates(np.isfinite(pl[:rows]).astype(float), [rr + fi, cc + fj], order=1) > 0.99
        d = np.exp(1j * po[:rows, c1:c2]) * np.conj(re_ + 1j * im_)
        d[~(okl & np.isfinite(po[:rows, c1:c2]))] = np.nan
        C = Co.z[:rows, c1:c2]
        good = np.isfinite(d) & (C > 0.3)
        dg = np.where(good, d, 0)
        n_row = good.sum(axis=1)
        prof = np.where(n_row > 30, dg.sum(axis=1), np.nan + 0j)
        lvl = np.nansum(prof)
        prof = prof * np.exp(-1j * np.angle(lvl))
        ang = np.degrees(np.angle(prof))
        prof_all[name] = ang
        lev[name] = np.degrees(np.angle(lvl))
        res = d * np.exp(-1j * np.radians(np.nan_to_num(ang)))[:, None]
        n_col = good.sum(axis=0)
        cprof = np.where(n_col > 30, (np.where(good, res * np.exp(-1j * np.angle(np.nansum(res[good]))), 0)).sum(axis=0), np.nan + 0j)
        cang = np.degrees(np.angle(cprof))
        x = np.arange(cang.size); m = np.isfinite(cang)
        ramp = np.polyfit(x[m], cang[m], 1)[0] * cang.size if m.sum() > 20 else np.nan
        rms = np.degrees(np.sqrt(np.mean(np.angle(res[good] * np.exp(-1j * np.angle(np.nansum(res[good])))) ** 2)))
        print(f"    ours - lab where corr > 0.3 ({good.mean() * 100:.0f} % of its pixels): mean {lev[name]:+.1f} deg; along-track profile sd {np.nanstd(ang):.1f}, range {np.nanmin(ang):+.1f}..{np.nanmax(ang):+.1f} deg; across-swath ramp {ramp:+.1f} deg; rms left {rms:.1f} deg", flush=True)
    # unwrapping
    a, b = Uo.z[:rows, :min(Ul.z.shape[1], Uo.z.shape[1])].astype(np.float64), Ul.z[:rows, :min(Ul.z.shape[1], Uo.z.shape[1])].astype(np.float64)
    a[a == 0] = np.nan; b[b == 0] = np.nan
    m = np.isfinite(a) & np.isfinite(b)
    k = np.rint((a - b - np.nanmedian((a - b)[m])) / TWO_PI)
    print(f" unwrap: ours finite {np.isfinite(Uo.z).sum()}, lab finite {np.isfinite(Ul.z).sum()}; both {m.sum()}; a whole cycle apart {np.mean(k[m] != 0) * 100:.3f} %; "
          f"max |cycles| {np.nanmax(np.abs(k[m])):.0f}", flush=True)
    ml, mo = maskL.z, maskO.z
    r, c = min(Ul.z.shape[0], ml.shape[0]), min(Ul.z.shape[1], ml.shape[1])
    lab_out = np.isfinite(Ul.z[:r, :c]) & ~np.isfinite(ml[:r, :c])
    r2, c2 = min(Uo.z.shape[0], mo.shape[0]), min(Uo.z.shape[1], mo.shape[1])
    our_out = np.isfinite(Uo.z[:r2, :c2]) & ~np.isfinite(mo[:r2, :c2])
    print(f" where mask_def is NaN: lab unwrapped {int(lab_out.sum())} px ({lab_out.sum() / np.isfinite(Ul.z).sum() * 100:.2f} %), ours {int(our_out.sum())} px", flush=True)
