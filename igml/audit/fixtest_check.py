"""Checks of the patched GMTSAR on T071's test pairs (work/fixtest/T071, run_fixtest.sh; 2026-10-29):
  1. F2's subswath grids: header origin = the subswath's alignment offset; the Rosamond reflectors in display_amp at
     declared - true ~ 0 (was +16.0 samples / +3.02 lines)
  2. the merged grids equal each subswath's at whole-cell offsets, and those offsets equal the continuous geometry
     exactly (merged declared - true at the reflectors ~ 0; was +42.1 / +0.02 for F2)
  3. phasefilt: no exact zeros left (I1), NaN only at the outermost rows / columns
  4. the phase against the lab's same pair: the lab's grid moved by half a cell (bilinear on the phasor) should match
Read only (over the Mac's NFS mount).

usage (igml env): python fixtest_check.py
"""
import numpy as np

from igml.real.grd import read_grd

W = "/Volumes/seasat2/IGML/work/fixtest/T071"
LAB = "/Volumes/seasat2/seasatHD1/InSAR_Processing/T071/product_T071"
PAIR = "2018285_2019004"
SLC = {13: (23385.72, 9967.62), 20: (22413.91, 9919.94), 21: (22260.25, 9912.72), 22: (22126.28, 9905.75),
       23: (23406.62, 9990.28), 24: (23173.72, 9976.84), 25: (22973.34, 9966.09), 26: (22763.69, 9956.34),
       27: (22552.75, 9944.75)}                                     # F2 frame, from cr_peak.py (the SLC's own peaks)
q = lambda m, c, p: 0.5 * (np.log(m) - np.log(p)) / (np.log(m) - 2 * np.log(c) + np.log(p))

if __name__ == "__main__":
    # 1. F2 subswath
    G = read_grd(f"{W}/F2/intf_all/{PAIR}/display_amp.grd")
    x, y, z = G.x, G.y, G.z.astype(float)
    print(f"F2 display_amp: x {x[0]:.0f}..{x[-1]:.0f} ({x.size}), y {y[0]:.0f}..{y[-1]:.0f} ({y.size}) "
          f"-> origin x {x[0] - 16:.0f}, y {y[0] - 3:.0f}")
    d = []
    for cid, (r, a) in SLC.items():
        j = int(np.argmin(abs(x - r))); i = int(np.argmin(abs(y - a)))
        w = z[i - 1:i + 2, j - 1:j + 2]
        ii, jj = np.unravel_index(np.argmax(w), w.shape); i, j = i - 1 + ii, j - 1 + jj
        xr = x[j] + 32 * q(z[i, j - 1], z[i, j], z[i, j + 1]); ya = y[i] + 6 * q(z[i - 1, j], z[i, j], z[i + 1, j])
        d.append((xr - r, ya - a))
    d = np.array(d)
    print(f"  reflectors (9): declared - true = range {np.median(d[:, 0]):+.2f} samples (sd {d[:, 0].std():.2f}), "
          f"azimuth {np.median(d[:, 1]):+.2f} lines (sd {d[:, 1].std():.2f})   [was +16.1 / +3.02]")
    # 2. merged placement
    M = read_grd(f"{W}/merge/{PAIR}/corr.grd")
    for n in (1, 2, 3):
        S = read_grd(f"{W}/F{n}/intf_all/{PAIR}/corr.grd")
        best = None
        for dr_ in range(-1, 160):
            for dc in ((0,) if n == 1 else range(620, 628) if n == 2 else range(1366, 1374)):
                a = S.z[200:400, 60:500]; b = M.z[200 + dr_:400 + dr_, 60 + dc:500 + dc]
                if b.shape != a.shape:
                    continue
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() > 1000:
                    f = np.mean(a[ok] == b[ok])
                    if best is None or f > best[0]:
                        best = (f, dr_, dc)
        f, dr_, dc = best
        # continuous: subswath cell (row, col) centre in merged-frame samples / lines vs the merged node there
        xs, ys = S.x[60], S.y[200]; xm, ym = M.x[60 + dc], M.y[200 + dr_]
        print(f"F{n}: in the merged grid at row +{dr_}, column +{dc} (identical {f:.3f}); node coordinates there: "
              f"subswath ({xs:.0f}, {ys:.0f}) -> merged ({xm:.0f}, {ym:.0f})")
    # merged-frame reflector check through F2's placement: merged declared - true
    dr = 299792458.0 / (2 * 64345238.125714)
    off_r = (846064.907962 - 799608.780580) / dr
    off_a = (76.577765563471 - 76.577754500585) * 86400 * 486.486310
    S2 = read_grd(f"{W}/F2/intf_all/{PAIR}/corr.grd")
    best = None
    for dr_ in range(60, 100):
        a = S2.z[200:400, 60:500]; b = M.z[200 + dr_:400 + dr_, 60 + 624:500 + 624]
        ok = np.isfinite(a) & np.isfinite(b)
        f = np.mean(a[ok] == b[ok]) if ok.sum() > 1000 else 0
        if best is None or f > best[0]:
            best = (f, dr_)
    mx = M.x[60 + 624] - S2.x[60]; my = M.y[200 + best[1]] - S2.y[200]
    print(f"merged-frame F2 node shift ({mx:.1f} samples, {my:.1f} lines) vs the continuous geometry ({off_r:.1f}, "
          f"{off_a:.1f}): merged declared - true at the reflectors = range {np.median(d[:, 0]) + mx - off_r:+.2f}, "
          f"azimuth {np.median(d[:, 1]) + my - off_a:+.2f}   [lab: +42.1 / +0.02]")
    # 3. phasefilt zeros / NaN
    Pm = read_grd(f"{W}/merge/{PAIR}/phasefilt.grd").z
    L = read_grd(f"{LAB}/merge/{PAIR}/phasefilt.grd").z
    print(f"merged phasefilt: exact zeros ours {np.sum(Pm == 0):,} vs the lab's {np.sum(L == 0):,}; NaN ours "
          f"{np.isnan(Pm).mean():.3%}, lab {np.isnan(L).mean():.3%}")
    # 4. phase vs the lab's: ours (F1 part) = the lab's moved by +0.5 cell in both directions
    zl = np.exp(1j * np.where(np.isfinite(L), L, 0)) * np.isfinite(L)
    zl_s = 0.25 * (zl[:-1, :-1] + zl[1:, :-1] + zl[:-1, 1:] + zl[1:, 1:])            # lab content at its cells' corners
    zo = np.exp(1j * np.where(np.isfinite(Pm), Pm, 0)) * np.isfinite(Pm)
    rows, cols = slice(500, 4300), slice(40, 600)                                     # F1, away from the edges
    for tag, ref in (("no shift", zl[rows, cols]), ("lab + half a cell", zl_s[rows, cols])):
        o = zo[rows, cols]
        c = np.abs(np.sum(o * np.conj(ref))) / np.sqrt(np.sum(np.abs(o) ** 2) * np.sum(np.abs(ref) ** 2))
        print(f"F1 phase, ours vs the lab's ({tag}): complex correlation {c:.3f}, mean phase difference "
              f"{np.degrees(np.angle(np.sum(o * np.conj(ref)))):+.1f} deg")
