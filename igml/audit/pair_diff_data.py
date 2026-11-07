"""Ours - the lab's, for pairs of a replaced track, as arrays for pair_diff_fig.py (2026-11-06; user: "Can you plot the difference
between the two?"). Read only, over the Mac's mount; the lab's originals from lab_originals_20261029, ours from the lab's tree.
  unw: unwrap.grd ours - lab as stored, on the same nodes, median removed, in mm LOS (every 4th block median); lab_only / ours_only pixels;
  ph : phasefilt ours - lab with the lab's content moved onto ours, per subswath, by the cross-correlated offset of the two corr.grd
       (the placement fix out), complex mean per 4 x 4 cells where coherence > 0.3, deg;
  prof_F<n>: ph's along-track profile at full resolution per subswath (deg), and the placement offsets.
usage (insar_sim env, PYTHONPATH=IGML): python pair_diff_data.py T OUT.npz PAIR [PAIR ...]"""
import json
import re
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

LABS = {"T173": "seasatHD3"}
T, out, pairs = sys.argv[1], sys.argv[2], sys.argv[3:]
ARCH = Path(f"/seasat/IGML/lab_originals_20261029/{T}/product_{T}/merge")
TREE = Path(f"/seasat/{LABS[T]}/InSAR_Processing/{T}/product_{T}")
OURS = TREE / "merge"
MM = 55.465763 / (4 * np.pi)
K = 4


def widths_from_log(d):
    t = (d / "merge_log").read_text()
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    f2 = re.search(r"F2 offset (-?[\d.]+) columns", t)
    f3 = re.search(r"F3 offset (-?[\d.]+) columns", t)
    w = int(re.search(r"Size\((\d+)x\d+\)", t)[1])
    ovl12, ovl23, c2, c3 = int(o[1]), int(o[2]), round(float(f2[1])), round(float(f3[1]))
    n1, n2 = ovl12 + c2, ovl23 + c3 - c2
    return n1, n2, w + ovl12 + ovl23 + 2 - n1 - n2


def shift_cells(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    A = np.where(m, a - a[m].mean(), 0.0)
    B = np.where(m, b - b[m].mean(), 0.0)
    cc = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(A) * np.conj(np.fft.rfft2(B)), s=A.shape))
    i, j = np.unravel_index(np.argmax(cc), cc.shape)

    def par(m1, c0, p1):
        d = m1 - 2 * c0 + p1
        return 0.0 if d == 0 else 0.5 * (m1 - p1) / d
    return i + par(cc[i - 1, j], cc[i, j], cc[i + 1, j]) - A.shape[0] // 2, j + par(cc[i, j - 1], cc[i, j], cc[i, j + 1]) - A.shape[1] // 2


def blocks(z, k=K):
    n, m = z.shape[0] // k * k, z.shape[1] // k * k
    return z[:n, :m].reshape(n // k, k, m // k, k)


res = {}
for pair in pairs:
    dl, dO = ARCH / pair, OURS / pair
    Ul, Uo = (read_grd(x / "unwrap.grd").z.astype(np.float64) for x in (dl, dO))
    Pl, Po = (read_grd(x / "phasefilt.grd").z.astype(np.float64) for x in (dl, dO))
    Cl, Co = (read_grd(x / "corr.grd").z.astype(np.float64) for x in (dl, dO))
    r, c = min(Ul.shape[0], Uo.shape[0]), min(Ul.shape[1], Uo.shape[1])
    Ul, Uo, Pl, Po, Cl, Co = (x[:r, :c] for x in (Ul, Uo, Pl, Po, Cl, Co))
    Ul[Ul == 0] = np.nan; Uo[Uo == 0] = np.nan; Pl[Pl == 0] = np.nan
    d = Uo - Ul
    med = np.nanmedian(d)
    d = (d - med) * MM
    k = np.rint((Uo - Ul - med) / (2 * np.pi))
    cyc = float(np.mean(k[np.isfinite(k)] != 0) * 100)
    v = d[np.isfinite(d)]
    q = np.percentile(v, [16, 84])
    res[f"{pair}_unw"] = np.nanmedian(blocks(d), axis=(1, 3)).astype(np.float32)
    lab_only, ours_only = np.isfinite(Ul) & ~np.isfinite(Uo), ~np.isfinite(Ul) & np.isfinite(Uo)
    res[f"{pair}_lab_only"] = blocks(lab_only).mean(axis=(1, 3)) > 0.5
    res[f"{pair}_ours_only"] = blocks(ours_only).mean(axis=(1, 3)) > 0.5
    wl = tuple(read_grd(ARCH.parent / f"F{n}/intf_all/{pair}/phasefilt.grd").z.shape[1] for n in (1, 2, 3))   # the lab's, archived
    sl, so = merge_seams(dl, wl), merge_seams(dO, widths_from_log(dO))
    ph = np.full((r, c), np.nan, dtype=np.complex64)
    meta = {"cyc": cyc, "sd_mm": float((q[1] - q[0]) / 2), "lab_only": int(lab_only.sum()), "ours_only": int(ours_only.sum()),
            "seams": [int(x) for x in so], "offsets": {}}
    el = np.exp(1j * np.nan_to_num(Pl))
    okl = np.isfinite(Pl).astype(np.float32)
    for n, (c1l, c2l), (c1o, c2o) in ((1, (0, sl[0]), (0, so[0])), (2, sl, so), (3, (sl[1], c), (so[1], c))):
        c1, c2 = max(c1l, c1o) + 12, min(c2l, c2o) - 12
        fi, fj = shift_cells(Cl[:, c1:c2], Co[:, c1:c2])
        meta["offsets"][f"F{n}"] = [float(fi), float(fj)]
        a1, a2 = max(c1o + 3, 0), min(c2o - 3, c)
        rr, cc = np.mgrid[0:r, a1:a2]
        co = [rr + fi, cc + fj]
        mv = map_coordinates(el.real, co, order=1) + 1j * map_coordinates(el.imag, co, order=1)
        ok = (map_coordinates(okl, co, order=1) > 0.99) & np.isfinite(Po[:, a1:a2]) & (Co[:, a1:a2] > 0.3)
        z = np.exp(1j * np.nan_to_num(Po[:, a1:a2])) * np.conj(mv)
        z = z / np.maximum(np.abs(z), 1e-6)
        z[~ok] = np.nan
        ph[:, a1:a2] = z
        s = np.where(ok, z, 0).sum(axis=1)
        n_row = ok.sum(axis=1)
        p = np.full(r, np.nan)
        p[n_row > 30] = np.degrees(np.angle(s[n_row > 30]))
        res[f"{pair}_prof_F{n}"] = p.astype(np.float32)
        print(f"{pair} F{n}: offset {fi:+.2f} rows, {fj:+.2f} columns; profile sd {np.nanstd(p):.1f} deg, {np.nanmin(p):+.1f}..{np.nanmax(p):+.1f}", flush=True)
    b = blocks(ph)
    f = np.isfinite(b.real)
    s = np.where(f, b, 0).sum(axis=(1, 3))
    a = np.full(s.shape, np.nan, dtype=np.float32)
    a[f.sum(axis=(1, 3)) >= 4] = np.degrees(np.angle(s[f.sum(axis=(1, 3)) >= 4]))
    res[f"{pair}_ph"] = a
    res[f"{pair}_meta"] = json.dumps(meta)
    print(f"{pair}: a cycle apart {cyc:.3f} %, robust sd {meta['sd_mm']:.2f} mm, lab only {meta['lab_only']}, ours only {meta['ours_only']}; seams {meta['seams']}", flush=True)
np.savez_compressed(out, K=K, pairs=np.array(pairs), **res)
print("->", out)
