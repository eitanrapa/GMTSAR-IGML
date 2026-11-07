"""Is the uncorrected elevation antenna pattern (audit A2) in the lab's interferograms? (2026-10-28.) make_s1a_tops never
applies the EAP correction, which IPF < 002.43 data need (S1A to 2015-03-10 on our tracks): GMTSAR's own model
(eap_est.py) predicts ~146 deg across IW1, ~38 IW2, ~2 IW3, i.e. a ~70 deg phase STEP at the F1/F2 seam in pairs whose
two dates sit on either side of IPF 002.43 (each date's version from its manifests), and none where both dates are on one side (the patterns cancel). Per
merged pair: the step across each seam from the WRAPPED phasefilt (circular means over every row: the phase change from
the last F1 column to the first F2 column, less the local gradient times the columns between them), so unwrapping
cannot make or hide it (F1/F2: the same ground in both subswaths, merge_swath's duplicated column). Seams from merge_log (igml.real.gmtsar.merge_seams). Read only.

usage (igml env): python eap_seam.py TRACK_ROOT START END [MAX_DAYS]
                  python eap_seam.py TRACK_ROOT --dirs PAIR_DIR [PAIR_DIR ...]   (e.g. our pairs, with
                  merge_seams.json; TRACK_ROOT the lab's tree, for the subswath widths)
"""
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

from igml.real.gmtsar import Track, doy0_to_date, merge_seams
from igml.real.grd import read_grd

root = sys.argv[1]
dirs = [Path(x) for x in sys.argv[sys.argv.index("--dirs") + 1:]] if "--dirs" in sys.argv else None
if dirs is None:
    start, end = sys.argv[2:4]
    maxd = int(sys.argv[4]) if len(sys.argv) > 4 else 48
tr = Track.open(root)
widths = tr.subswath_widths()
day = lambda s: datetime.strptime(s, "%Y%m%d")
_IPF = {}


def ipf(ymd):
    """The date's IPF version(s) from its scenes' manifests in the archive (ESA reprocessed some scenes: T115's
    2015-02-22 is 002.43 between 002.36 dates, so the version is read per scene, never inferred from the date)."""
    if ymd not in _IPF:
        vs = set()
        for m in Path("/palsar").glob(f"S1?/CA/{ymd[:4]}/{tr.name}/S1?_IW_SLC__1S?V_{ymd}T*.SAFE/manifest.safe"):   # 2014-15: 1SSV
            hit = re.search(r'name="Sentinel-1 IPF" version="([0-9.]+)"', m.read_text(errors="replace"))
            if hit:
                vs.add(hit.group(1))
        _IPF[ymd] = sorted(vs)
    return _IPF[ymd]


def uncorrected(ymd):
    v = ipf(ymd)
    return None if not v else (all(x < "002.43" for x in v) if len(set(x < "002.43" for x in v)) == 1 else "mixed")


def cmean(z):
    z = z[np.isfinite(z)]
    m = z.mean()
    return np.angle(m), np.abs(m), z.size


def step(ph, c_old, c_new, spacing):
    """Phase step [deg] from column c_old (the old subswath) to c_new (the next), whose grounds are ``spacing`` columns
    apart, less the local gradient (circular mean of adjacent-column differences beside the seam) times the spacing.
    F1/F2: s12 - 1 (F1) and s12 (F2, merge_swath's extra column) image the SAME ground -> spacing 0, no gradient;
    F2/F3: s23 - 1 (F2: ground s23 - 2) and s23 (F3) -> spacing 2."""
    e = np.exp(1j * ph)
    g = 0.0
    if spacing:
        side = list(range(c_old - 6, c_old)) + list(range(c_new + 1, c_new + 7))
        g, _, _ = cmean(np.concatenate([e[:, k + 1] * np.conj(e[:, k]) for k in side]))
    a, coh, n = cmean(e[:, c_new] * np.conj(e[:, c_old]) * np.exp(-1j * spacing * g))
    return np.degrees(a), coh, n


rows = []
todo = [(doy0_to_date(d.name[:7]), doy0_to_date(d.name[8:]), d) for d in dirs] if dirs else \
    [(a, b, tr.pair_dir(a, b)) for a, b in tr.pairs(files=("phasefilt.grd", "merge_log"))
     if start <= a and b <= end and (day(b) - day(a)).days <= maxd]
for ref, sec, d in todo:
    s12, s23 = merge_seams(d, widths)
    ph = read_grd(d / "phasefilt.grd").z.astype(np.float64)
    ph[ph == 0] = np.nan                                  # the I1 strip and masked pixels
    a12, c12, n12 = step(ph, s12 - 1, s12, 0)
    a23, c23, n23 = step(ph, s23 - 1, s23, 2)
    ua, ub = uncorrected(ref), uncorrected(sec)
    span = None if None in (ua, ub) or "mixed" in (ua, ub) else ua != ub
    rows.append((ref, sec, span, a12, c12, a23, c23))
    print(f"{ref}_{sec} IPF {'/'.join(ipf(ref)) or '?'} -> {'/'.join(ipf(sec)) or '?'} "
          f"{'ACROSS' if span else ('one side' if span is False else 'unknown')}: F1/F2 step {a12:+6.1f} deg (|mean| "
          f"{c12:.2f}, {n12} rows)  F2/F3 step {a23:+6.1f} deg (|mean| {c23:.2f})", flush=True)
for span in (True, False):
    r = [x for x in rows if x[2] == span]
    if r:
        print(f"{'across' if span else 'one side'}: {len(r)} pairs, F1/F2 median {np.median([x[3] for x in r]):+.1f} deg, "
              f"F2/F3 median {np.median([x[5] for x in r]):+.1f} deg")
