# A2 (EAP) figure: the F1/F2 seam of T173 pairs across IPF 002.43, the lab's merged phase against our rebuild.
# Seams and the step measure as gmtsar_fix/audit/proc1_check.py; read only over the Mac's NFS mount (/seasat).
import json, re, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/erapaportbruck/Documents/IGML")
from igml.real.gmtsar import merge_seams
from igml.real.grd import read_grd

OURS = Path("/seasat/IGML/InSAR_Processing/T173/product_T173")
LAB = Path("/seasat/seasatHD3/InSAR_Processing/T173/product_T173")

def widths_from_log(d):
    t = (d / "merge_log").read_text()
    o = re.search(r"ovl12,23:\s*(-?\d+),\s*(-?\d+)", t)
    f2 = re.search(r"F2 offset (-?[\d.]+) columns", t); f3 = re.search(r"F3 offset (-?[\d.]+) columns", t)
    w = int(re.search(r"Size\((\d+)x\d+\)", t)[1])
    ovl12, ovl23, c2, c3 = int(o[1]), int(o[2]), round(float(f2[1])), round(float(f3[1]))
    n1, n2 = ovl12 + c2, ovl23 + c3 - c2
    return n1, n2, w + ovl12 + ovl23 + 2 - n1 - n2

def seam12(d, ours):
    if ours:
        return merge_seams(d, widths_from_log(d))[0], 1
    w = tuple(read_grd(LAB / f"F{n}/intf_all/{d.name}/phasefilt.grd").z.shape[1] for n in (1, 2, 3))
    return merge_seams(d, w)[0], 0

def step_rows(ph, c_old, c_new, spacing, block=None):
    """eap_seam.step per block of rows (block None: one value): phase change across the seam less the local gradient."""
    e = np.exp(1j * ph)
    side = list(range(c_old - 6, c_old)) + list(range(c_new + 1, c_new + 7))
    def one(rs):
        g = 0.0
        if spacing:
            z = np.concatenate([(e[rs, k + 1] * np.conj(e[rs, k])).ravel() for k in side]); z = z[np.isfinite(z)]
            g = np.angle(z.mean()) if z.size else np.nan
        z = (e[rs, c_new] * np.conj(e[rs, c_old]) * np.exp(-1j * spacing * g)); z = z[np.isfinite(z)]
        return np.degrees(np.angle(z.mean())) if z.size > 20 else np.nan
    if block is None:
        return one(slice(None))
    return np.array([one(slice(r, r + block)) for r in range(0, ph.shape[0] - block + 1, block)])

def load(pair, ours):
    d = (OURS if ours else LAB) / "merge" / pair
    ph = read_grd(d / "phasefilt.grd").z.astype(np.float64)
    if not ours:
        ph[ph == 0] = np.nan
    s, sp = seam12(d, ours)
    return ph, s, sp

INK, MUTED, CL, CO = "#2d3748", "#718096", "#c05621", "#2b6cb0"   # text; the lab (orange), ours (blue)
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False})
P = "2015020_2015092"                                   # 2015-01-21 (IPF 002.36) x 2015-04-03 (002.43)
lph, ls, lsp = load(P, False); oph, os_, osp = load(P, True)
H = 90                                                   # columns either side of the seam
fig = plt.figure(figsize=(12.5, 8.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.15], hspace=0.55, wspace=0.18)
for i, (ph, s, name, col) in enumerate(((lph, ls, "the lab's (released GMTSAR)", CL), (oph, os_, "ours (EAP corrected)", CO))):
    ax = fig.add_subplot(gs[i, :])
    cut = ph[:, s - H:s + H].T                           # range down, azimuth across
    ax.imshow(cut, cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest",
              extent=[0, ph.shape[0], H, -H])
    ax.axhline(-0.5, color=col, lw=0.8, ls=(0, (4, 3)))
    ax.set_ylabel("columns from\nthe F1/F2 seam")
    ax.set_title(f"{'ab'[i]}  {name}: wrapped phase of 2015-01-21 x 2015-04-03 around the F1/F2 seam (F1 above, F2 below)",
                 loc="left")
    if i == 1:
        ax.set_xlabel("merged row (azimuth)")
# c: the step along the seam, per 150 rows
ax = fig.add_subplot(gs[2, 0])
blk = 150
for ph, s, sp, col, name in ((lph, ls, lsp, CL, "the lab's"), (oph, os_, osp, CO, "ours")):
    st = step_rows(ph, s - 1, s, sp, blk)
    ax.plot((np.arange(len(st)) + 0.5) * blk, st, color=col, lw=1.3, marker="o", ms=3, label=f"{name}: {np.nanmedian(st):+.1f} deg")
ax.axhline(0, color=MUTED, lw=0.6)
ax.set_ylim(-120, 40); ax.set_xlabel("merged row (azimuth)"); ax.set_ylabel("phase step at the seam [deg]")
ax.set_title("c  The step along the seam (150-row blocks)", loc="left"); ax.legend(frameon=False, loc="center right")
# d: every early-2015 pair in both stacks
ax = fig.add_subplot(gs[2, 1])
pairs = sorted(p.name for p in (OURS / "merge").iterdir() if re.match(r"^2015(020|068|092|116|140)_2015(068|092|116|140|164)$", p.name))
xs, rows = [], []
for k, pr in enumerate(pairs):
    meta = json.loads((OURS / "merge" / pr / "igml_proc.json").read_text())
    across = len(meta.get("eap_removed", [])) == 1
    ph, s, sp = load(pr, True); so = step_rows(ph, s - 1, s, sp)
    lab = np.nan
    if (LAB / "merge" / pr / "phasefilt.grd").exists():
        ph, s, sp = load(pr, False); lab = step_rows(ph, s - 1, s, sp)
    rows.append((pr, across, lab, so))
    print(f"{pr}: {'across' if across else 'one side'}  lab {lab:+6.1f}  ours {so:+6.1f}", flush=True)
for j, (pr, across, lab, so) in enumerate(rows):
    ax.plot([j, j], [lab, so], color="#cbd5e0", lw=1, zorder=1)
    ax.scatter([j], [lab], color=CL, s=28, zorder=2, label="the lab's" if j == 0 else None)
    ax.scatter([j], [so], color=CO, s=28, zorder=3, label="ours" if j == 0 else None)
ax.set_xticks(range(len(rows)))
import datetime as _dt
md = lambda y: (_dt.date(int(y[:4]), 1, 1) + _dt.timedelta(days=int(y[4:]))).strftime("%m-%d")
ax.set_xticklabels([f"{md(p[:7])}\n{md(p[8:])}" + ("\nacross" if a else "") for p, a, _, _ in rows], fontsize=7)
ax.axhline(0, color=MUTED, lw=0.6); ax.set_ylim(-120, 40)
ax.set_ylabel("F1/F2 step [deg]"); ax.legend(frameon=False, loc="lower left")
ax.set_title("d  Every early-2015 pair (the lab made only one across the change)", loc="left")
fig.suptitle("GMTSAR A2: IPF 002.36 data never get the antenna-pattern correction, so pairs across the change step at the "
             "F1/F2 seam   (T173)", x=0.125, ha="left", color=INK, fontsize=11)
out = "/Users/erapaportbruck/Documents/IGML/notebooks/figures/a2_eap_seam_t173.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white"); print(out)
