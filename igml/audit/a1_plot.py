import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

z = np.load("t166_f1.npz")
H = int(z["H"]); C = z["centres"]
CA, CB, INK, MUTED = "#2b6cb0", "#dd6b20", "#2d3748", "#718096"   # date A, date B (CVD-safe pair), text
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False})

def info(k):
    r0 = C[k] - H
    da, db, coh = z[f"da{k}"], z[f"db{k}"], z[f"coh{k}"]
    ja = np.abs(np.angle(np.exp(1j * np.diff(da)))); jb = np.abs(np.angle(np.exp(1j * np.diff(db))))
    return r0, ja, jb, r0 + 1 + int(np.argmax(ja)), r0 + 1 + int(np.argmax(jb)), coh, r0 + 1 + np.arange(len(coh))

fig = plt.figure(figsize=(12.5, 8.2))
gs = fig.add_gridspec(2, 5, height_ratios=[1.25, 1], hspace=0.42, wspace=0.32)
k = 1
r0, ja, jb, sa, sb, coh, lines = info(k)
lo, hi = 6725, 6805                                        # zoom around the boundary
# (a) interferogram phase
ax = fig.add_subplot(gs[0, 0:2])
ph = z[f"ph{k}"]; x0, x1, xm = z[f"x{k}"]
ax.imshow(ph, cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest",
          extent=[x0 + xm, x0 + xm + ph.shape[1] * 8, r0 + ph.shape[0], r0])
ax.set_ylim(hi, lo); ax.set_xlabel("range sample (mid-swath; dominant fringe removed)"); ax.set_ylabel("line (super-master grid)")
ax.set_title("a  Interferogram phase", loc="left")
ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
for s, col in ((sa, CA), (sb, CB)):
    ax.axhline(s - 0.5, color=col, lw=1.2)
# (b) coherence per line
ax2 = fig.add_subplot(gs[0, 2], sharey=ax)
ax2.axhspan(min(sa, sb) - 0.5, max(sa, sb) - 0.5, color="#e2e8f0", lw=0)
ax2.plot(coh, lines, color=INK, lw=1.2)
ax2.set_xlim(0, 1); ax2.set_xlabel("coherence (3 lines x 32 samples,\nmedian over range)")
ax2.set_title("b  Pair coherence", loc="left")
plt.setp(ax2.get_yticklabels(), visible=False)
# (c) each date's Doppler jump
ax3 = fig.add_subplot(gs[0, 3:5], sharey=ax)
ax3.plot(ja, r0 + 1 + np.arange(len(ja)), color=CA, lw=1.2, label=f"2026-05-29: switches at line {sa}")
ax3.plot(jb, r0 + 1 + np.arange(len(jb)), color=CB, lw=1.2, label=f"2026-06-11: switches at line {sb}")
ax3.set_xlabel("line-to-line change of each date's Doppler phase (rad)")
ax3.set_title("c  Where each date switches burst", loc="left")
ax3.legend(frameon=False, loc="lower right", labelcolor=INK)
plt.setp(ax3.get_yticklabels(), visible=False)
# bottom: the other boundaries
others = [0, 2, 3, 4, 5]
for j, kk in enumerate(others):
    r0k, jak, jbk, sak, sbk, cohk, linesk = info(kk)
    axk = fig.add_subplot(gs[1, j])
    axk.axhspan(min(sak, sbk) - 0.5, max(sak, sbk) - 0.5, color="#e2e8f0", lw=0)
    axk.plot(cohk, linesk, color=INK, lw=1.1)
    axk.axhline(sak - 0.5, color=CA, lw=1); axk.axhline(sbk - 0.5, color=CB, lw=1, ls=(0, (3, 2)))
    mid = (sak + sbk) / 2
    axk.set_ylim(mid + 25, mid - 25); axk.set_xlim(0, 1)
    gap = abs(sbk - sak)
    axk.set_title(f"lines {sak} / {sbk}: {gap} apart" if gap else f"line {sak}: same line", loc="left", fontsize=8.5)
    axk.set_xlabel("coherence")
    if j == 0:
        axk.set_ylabel("line")
fig.text(0.125, 0.43, "Every other boundary of the same pair: the dip is as wide as the gap between the two dates' switch lines "
         "(blue 2026-05-29, orange 2026-06-11)", color=INK, fontsize=9)
fig.suptitle("GMTSAR A1: each date switches bursts at its own line, and the pair decorrelates in between   "
             "(T166 F1, 2026-05-29 x 2026-06-11, full resolution)", x=0.125, ha="left", color=INK, fontsize=11)
out = "/Users/erapaportbruck/Documents/IGML/notebooks/figures/a1_burst_switch_t166.png"
fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
print(out)
