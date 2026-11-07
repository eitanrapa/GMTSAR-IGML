"""Checks of a track's first batch made with the patched GMTSAR (PROC proc1; 2026-10-29), against the lab's same pairs:
  1. the F1/F2 and F2/F3 phase steps of the wrapped phasefilt (eap_seam.py's measure), ours and the lab's: a pair whose
     dates sit on either side of IPF 002.43 had ~-76 / +4 deg in the lab's (the elevation antenna pattern, audit A2);
     ours should be ~0 on every pair. Seams: merge_swath's regions (igml.real.gmtsar.merge_seams) with each grid's own
     subswath widths; the patched merge_swath puts every subswath on its true columns, so the last column before a seam
     and the first after it are adjacent ground (spacing 1), where the lab's had 0 (F1/F2: M1) and 2 (F2/F3)
  2. F1's phase, ours against the lab's moved by half a cell (I2: the lab's nodes sit +16 samples / +3 lines from
     their content): complex correlation and mean difference
  3. the unwrapped phase, ours less the lab's (both referenced to their common median): pixels a whole cycle apart
Read only (over the Mac's NFS mount).

usage (igml env): python proc1_check.py T173 [PAIR ...]
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

LABS = {"T173": "seasatHD3", "T071": "seasatHD1", "T064": "seasatHD2", "T166": "seasatHD4"}
T = sys.argv[1]
OURS = Path(f"/Volumes/seasat2/IGML/InSAR_Processing/{T}/product_{T}")
LAB = Path(f"/Volumes/seasat2/{LABS[T]}/InSAR_Processing/{T}/product_{T}")


def widths_from_log(d):
    """(N1, N2, N3) of the subswath grids a merged pair was made from, from its merge_log: our merge_swath prints the
    exact whole-cell offsets of F2 / F3, and ovl = N - offset."""
    t = (d / "merge_log").read_text()
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    f2 = re.search(r"F2 offset (-?[\d.]+) columns", t)
    f3 = re.search(r"F3 offset (-?[\d.]+) columns", t)
    w = int(re.search(r"Size\((\d+)x\d+\)", t)[1])
    ovl12, ovl23, c2, c3 = int(o[1]), int(o[2]), round(float(f2[1])), round(float(f3[1]))
    n1, n2 = ovl12 + c2, ovl23 + c3 - c2
    return n1, n2, w + ovl12 + ovl23 + 2 - n1 - n2


def lab_widths(pair):
    return tuple(read_grd(LAB / f"F{n}/intf_all/{pair}/phasefilt.grd").z.shape[1] for n in (1, 2, 3))


def cmean(z):
    z = z[np.isfinite(z)]
    m = z.mean()
    return np.angle(m), np.abs(m)


def step(ph, c_old, c_new, spacing):
    """eap_seam.step: the phase change from column c_old to c_new, less the local gradient times their spacing."""
    e = np.exp(1j * ph)
    g = 0.0
    if spacing:
        side = list(range(c_old - 6, c_old)) + list(range(c_new + 1, c_new + 7))
        g, _ = cmean(np.concatenate([e[:, k + 1] * np.conj(e[:, k]) for k in side]))
    a, coh = cmean(e[:, c_new] * np.conj(e[:, c_old]) * np.exp(-1j * spacing * g))
    return np.degrees(a), coh


def seams(d, ours):
    if ours:
        s12, s23 = merge_seams(d, widths_from_log(d))
        return (s12, 1), (s23, 1)
    s12, s23 = merge_seams(d, lab_widths(d.name))
    return (s12, 0), (s23, 2)


def half_cell(z):
    """The lab's grid sampled half a cell back (-16 samples / -3 lines = half a 32 x 6 cell), bilinear."""
    return 0.25 * (z[:-1, :-1] + z[1:, :-1] + z[:-1, 1:] + z[1:, 1:])


pairs = sys.argv[2:] or sorted(p.name for p in (OURS / "merge").iterdir() if re.match(r"^\d{7}_\d{7}$", p.name))
rows = []
for pair in pairs:
    d = OURS / "merge" / pair
    meta = json.loads((d / "igml_proc.json").read_text())
    ph = read_grd(d / "phasefilt.grd").z.astype(np.float64)
    (s12, sp12), (s23, sp23) = seams(d, True)
    a12, c12 = step(ph, s12 - 1, s12, sp12)
    a23, c23 = step(ph, s23 - 1, s23, sp23)
    out = f"{pair} EAP removed {','.join(meta['eap_removed']) or '-':19s} ours F1/F2 {a12:+6.1f} ({c12:.2f}) F2/F3 {a23:+6.1f} ({c23:.2f})"
    rec = {"pair": pair, "eap": meta["eap_removed"], "ours": (a12, a23)}
    ld = LAB / "merge" / pair
    if (ld / "phasefilt.grd").exists():
        lph = read_grd(ld / "phasefilt.grd").z.astype(np.float64)
        lph[lph == 0] = np.nan                                 # I1's strip and the lab's masked pixels
        (l12, lsp12), (l23, lsp23) = seams(ld, False)
        b12, _ = step(lph, l12 - 1, l12, lsp12)
        b23, _ = step(lph, l23 - 1, l23, lsp23)
        out += f" | lab {b12:+6.1f} {b23:+6.1f}"
        rec["lab"] = (b12, b23)
        # F1, ours vs the lab's half a cell back: F1 is at the merged grid's origin in both
        n = min(s12, l12) - 2
        zo = np.exp(1j * ph[:, :n])
        zl = np.exp(1j * lph[:, :n + 1])
        zl = np.exp(1j * np.angle(half_cell(np.pad(zl, ((0, 1), (0, 0)), constant_values=np.nan))))[:, :n]
        ok = np.isfinite(zo) & np.isfinite(zl)
        c = zo[ok] * np.conj(zl[ok])
        out += f" | F1 corr {np.abs(c.mean()):.3f} mean diff {np.degrees(np.angle(c.mean())):+5.1f}"
        rec["f1"] = (float(np.abs(c.mean())), float(np.degrees(np.angle(c.mean()))))
        if (d / "unwrap.grd").exists() and (ld / "unwrap.grd").exists():
            uo = read_grd(d / "unwrap.grd").z.astype(np.float64)
            ul = read_grd(ld / "unwrap.grd").z.astype(np.float64)
            # whole cells of F1 only (the merged grids differ by the subswath placement elsewhere)
            uo, ul = uo[:-1, :n], 0.25 * (ul[:-1, :n] + ul[1:, :n] + ul[:-1, 1:n + 1] + ul[1:, 1:n + 1])
            dd = uo - ul
            ok = np.isfinite(dd)
            dd = dd[ok] - np.median(dd[ok])
            k = np.round(dd / (2 * np.pi))
            out += f" | unwrap F1: {100 * np.mean(k != 0):.2f} % a cycle off ({ok.sum()} px)"
            rec["cyc"] = float(100 * np.mean(k != 0))
    print(out, flush=True)
    rows.append(rec)

for name, sel in (("across IPF 002.43", lambda r: len(r["eap"]) == 1),
                  ("both dates EAP-corrected", lambda r: len(r["eap"]) == 2),
                  ("neither (IPF >= 002.43)", lambda r: len(r["eap"]) == 0)):
    r = [x for x in rows if sel(x)]
    if not r:
        continue
    s = f"{name}: {len(r)} pairs, ours F1/F2 median {np.median([x['ours'][0] for x in r]):+.1f}, F2/F3 {np.median([x['ours'][1] for x in r]):+.1f}"
    lr = [x for x in r if "lab" in x]
    if lr:
        s += (f"; the lab's ({len(lr)}) {np.median([x['lab'][0] for x in lr]):+.1f} / {np.median([x['lab'][1] for x in lr]):+.1f}"
              f"; F1 corr median {np.median([x['f1'][0] for x in lr]):.3f}")
    cr = [x["cyc"] for x in r if "cyc" in x]
    if cr:
        s += f"; unwrap a cycle off median {np.median(cr):.2f} % (max {max(cr):.2f})"
    print(s)
