import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
INK, MUTED, CL, CO = "#2d3748", "#718096", "#c05621", "#2b6cb0"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False})
L, O = np.load("a4_lab_f1.npz"), np.load("a4_ours_f1.npz")
S = np.array([-2, -1, 0, 1, 2])
mid = L["starts"] + 27
fig = plt.figure(figsize=(12.5, 10.2))
gs = fig.add_gridspec(5, 1, height_ratios=[0.8, 0.8, 1, 1, 1], hspace=0.66)
ax = fig.add_subplot(gs[0])
for z, col, mk, name, dy in ((L, CL, "o", "the lab's image (released GMTSAR)", 0.08), (O, CO, "s", "our re-cut (assemble_tops fixed)", -0.08)):
    best = S[np.nanargmax(z["C"], axis=1)]
    ax.scatter(mid, best + dy, s=22, marker=mk, facecolor="none" if mk == "o" else col, edgecolor=col, label=name)
ax.set_ylim(-2.5, 2.5); ax.set_yticks([-2, -1, 0, 1, 2]); ax.grid(axis="y", color="#e2e8f0", lw=0.8)
ax.set_xlim(0, 28292); ax.set_ylabel("lines")
ax.set_title("a  Line shift of 2014-12-04 that best matches 2014-12-28, per block of 50 lines (F1)", loc="left")
ax.legend(frameon=False, loc="lower left", ncol=2)
ax = fig.add_subplot(gs[1], sharex=ax)
ax.plot(mid, O["C"][:, 2], color=CO, lw=1.2, label="ours, unshifted")
ax.plot(mid, L["C"][:, 2], color=CL, lw=1.2, label="the lab's, unshifted")
ax.plot(mid, np.nanmax(L["C"], axis=1), color=CL, lw=1.6, ls=(0, (2, 2)), label="the lab's, shifted 1 line")
ax.set_ylim(0, 1); ax.set_ylabel("coherence p90")
ax.set_title("b  Coherence with 2014-12-28 per block", loc="left"); ax.legend(frameon=False, loc="lower left", ncol=3)
c0, c1 = int(L["c0"]), int(L["c1"])
for i, (z, name, col) in enumerate(((L, "the lab's image", CL), (O, "our re-cut", CO))):
    ax = fig.add_subplot(gs[2 + i])
    ph = z["ph"]
    ax.imshow(ph.T, cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest",
              extent=[c0, c0 + ph.shape[0] * 6, ph.shape[1] * 24, 0])
    ax.axvline(23239, color=col, lw=0.8, ls=(0, (4, 3)))
    ax.set_ylabel("range sample\n(2048, mid-swath)")
    ax.set_title(f"{'cd'[i]}  2014-12-04 x 2014-12-28, interferogram phase, {name} (6 x 24 looks)", loc="left")
ax = fig.add_subplot(gs[4])
dph = np.angle(np.exp(1j * (L["ph"] - O["ph"])))
ax.imshow(dph.T, cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect="auto", interpolation="nearest",
          extent=[c0, c0 + dph.shape[0] * 6, dph.shape[1] * 24, 0])
ax.axvline(23239, color=CL, lw=0.8, ls=(0, (4, 3)))
ax.set_ylabel("range sample\n(2048, mid-swath)"); ax.set_xlabel("line (super-master grid)")
ax.set_title("e  c minus d: identical (0) until line ~22,870, then the lab image's one-line slip to the end "
             "(~5,400 lines; TOPS turns it into fringes)", loc="left")
fig.suptitle("GMTSAR A4: bursts placed by their ANX time land a line early after an epoch step between a date's scenes   "
             "(T173 2014-12-04, F1: the lab's image vs our re-cut)", x=0.125, ha="left", color=INK, fontsize=11)
out = "/Users/erapaportbruck/Documents/IGML/notebooks/figures/a4_burst_epoch_t173.png"
fig.savefig(out, dpi=150, bbox_inches="tight", facecolor="white"); print(out)
