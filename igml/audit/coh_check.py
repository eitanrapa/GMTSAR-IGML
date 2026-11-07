"""Does the patched GMTSAR (PROC proc1: the geolocation fixes I2 / M6 / M1, the edges I1, the EAP A2) change the
coherence? (2026-10-29; user: "Did we find this affects coherence? Or is that still to be determined".) A track's
first proc1 batch against the lab's same pairs, corr.grd:
  1. per subswath, the median / mean coherence of each grid over its own columns (merge_swath's seams, as
     proc1_check.py), every pixel that is finite and > 0: no pixel matching needed
  2. F1 pixel for pixel: ours against the lab's moved by half a cell (I2), as proc1_check.py does for the phase
  3. exact zeros in each grid (unfilled edges)
Read only (over the Mac's NFS mount).

usage (igml env): python coh_check.py T173 [EVERY]       (every EVERY-th stored pair the lab also has; default 6)
"""
import json
import sys
from pathlib import Path

import numpy as np

from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

sys.argv += [] if len(sys.argv) > 2 else ["6"]
sys.path.insert(0, str(Path(__file__).parent))
LABS = {"T173": "seasatHD3", "T071": "seasatHD1", "T064": "seasatHD2", "T166": "seasatHD4"}
T, EVERY = sys.argv[1], int(sys.argv[2])
OURS = Path(f"/Volumes/seasat2/IGML/InSAR_Processing/{T}/product_{T}/merge")
LAB = Path(f"/Volumes/seasat2/{LABS[T]}/InSAR_Processing/{T}/product_{T}")


def widths_from_log(d):                                                # = proc1_check.widths_from_log
    import re
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


def stats(z):
    v = z[np.isfinite(z) & (z > 0)]
    return np.median(v), v.mean(), np.mean(v > 0.3), v.size


pairs = sorted(p.name for p in OURS.iterdir() if (p / "igml_proc.json").exists() and (LAB / "merge" / p.name / "corr.grd").exists())
pairs = pairs[::EVERY]
rows = []
for pair in pairs:
    d, ld = OURS / pair, LAB / "merge" / pair
    eap = json.loads((d / "igml_proc.json").read_text())["eap_removed"]
    co = read_grd(d / "corr.grd").z.astype(np.float64)
    cl = read_grd(ld / "corr.grd").z.astype(np.float64)
    so, sl = merge_seams(d, widths_from_log(d)), merge_seams(ld, lab_widths(pair))
    rec = {"pair": pair, "eap": len(eap), "zeros": (int(np.sum(co == 0)), int(np.sum(cl == 0))),
           "shape": (co.shape, cl.shape)}
    out = f"{pair} eap {len(eap)}"
    for k, (ao, bo, al, bl) in enumerate(((0, so[0], 0, sl[0]), (so[0], so[1], sl[0], sl[1]),
                                          (so[1], co.shape[1], sl[1], cl.shape[1])), 1):
        mo, mno, fo, no = stats(co[:, ao:bo])
        ml, mnl, fl, nl = stats(cl[:, al:bl])
        rec[f"F{k}"] = (mo, ml, mno, mnl, fo, fl, no, nl)
        out += f" | F{k} median {mo:.3f} vs lab {ml:.3f} ({mo - ml:+.4f})"
    # F1 pixel for pixel: the lab's half a cell back (bilinear), proc1_check's comparison
    n = min(so[0], sl[0]) - 2
    zo = co[:-1, :n]
    zl = cl[:, :n + 1]
    zl = 0.25 * (zl[:-1, :-1] + zl[1:, :-1] + zl[:-1, 1:] + zl[1:, 1:])
    ok = np.isfinite(zo) & np.isfinite(zl) & (zo > 0) & (zl > 0)
    r = np.corrcoef(zo[ok], zl[ok])[0, 1]
    dd = zo[ok] - zl[ok]
    rec["f1px"] = (float(r), float(dd.mean()), float(np.median(np.abs(dd))))
    out += f" | F1 px r {r:.3f}, mean diff {dd.mean():+.4f} | zeros {rec['zeros'][0]} vs {rec['zeros'][1]}"
    print(out, flush=True)
    rows.append(rec)

print(f"\n{T}: {len(rows)} pairs (every {EVERY}th of the stored pairs the lab also has)")
for k in (1, 2, 3):
    a = np.array([r[f"F{k}"] for r in rows])
    dm = a[:, 0] - a[:, 1]
    print(f"F{k}: median coherence ours {np.median(a[:, 0]):.3f} vs lab {np.median(a[:, 1]):.3f}; per-pair difference "
          f"median {np.median(dm):+.4f} (range {dm.min():+.4f} .. {dm.max():+.4f}); mean coh diff "
          f"{np.median(a[:, 2] - a[:, 3]):+.4f}; fraction > 0.3: {np.median(a[:, 4]):.3f} vs {np.median(a[:, 5]):.3f}; "
          f"valid pixels ours / lab {np.median(a[:, 6] / a[:, 7]):.4f}")
p = np.array([r["f1px"] for r in rows])
print(f"F1 pixel for pixel (lab half a cell back): r median {np.median(p[:, 0]):.3f} ({p[:, 0].min():.3f}..{p[:, 0].max():.3f}), "
      f"mean diff median {np.median(p[:, 1]):+.4f}, median |diff| {np.median(p[:, 2]):.4f}")
for e in (0, 1, 2):
    s = [r for r in rows if r["eap"] == e]
    if s:
        print(f"  pairs with {e} EAP-corrected date(s): {len(s)}, F1/F2/F3 median differences "
              + " / ".join(f"{np.median([r[f'F{k}'][0] - r[f'F{k}'][1] for r in s]):+.4f}" for k in (1, 2, 3)))
z = np.array([r["zeros"] for r in rows])
print(f"exact zeros per grid: ours median {np.median(z[:, 0]):.0f}, lab {np.median(z[:, 1]):.0f}")
