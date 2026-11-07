"""Our interferogram of a pair against the lab's, per subswath (2026-11-01): the phase difference ours - lab (wrapped),
its step across each super-master burst switch (A5: the lab's alignment left one, ours removes it), and the coherence
on the switch line against 20-40 lines away (A1: the lab's band of decorrelated lines). Grids in each subswath's radar
coordinates (y = SLC line). usage: python intf_vs_lab.py DIR cutsF1 cutsF2 cutsF3 [OUT.png]  (cuts comma-separated)"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/erapaportbruck/Documents/IGML")
from igml.real.grd import read_grd

d = Path(sys.argv[1])
cuts = [[int(c) for c in a.split(",")] for a in sys.argv[2:5]]
out = sys.argv[5] if len(sys.argv) > 5 else None
rows = []
for n in (1, 2, 3):
    po, pl = read_grd(d / f"ours_F{n}_phasefilt.grd"), read_grd(d / f"lab_F{n}_phasefilt.grd")
    co, cl = read_grd(d / f"ours_F{n}_corr.grd"), read_grd(d / f"lab_F{n}_corr.grd")
    y = po.y
    ok = np.isfinite(po.z) & np.isfinite(pl.z) & (po.z != 0) & (pl.z != 0)
    dz = np.where(ok, np.exp(1j * (po.z - pl.z)), 0)
    st, dips = [], []
    for c in cuts[n - 1]:
        above = (y > c - 60) & (y < c - 6); below = (y > c + 6) & (y < c + 60)
        if above.sum() < 3 or below.sum() < 3:
            continue
        a, b = dz[above].sum(), dz[below].sum()
        st.append(np.degrees(np.angle(b * np.conj(a))))
        on = np.abs(y - c) <= 4
        near = (np.abs(y - c) >= 20) & (np.abs(y - c) <= 40)
        f = lambda z: (np.nanmedian(z[on]), np.nanmedian(z[near]))
        dips.append((*f(co.z), *f(cl.z)))
    st, dips = np.array(st), np.array(dips)
    rows.append((n, st, dips))
    print(f"F{n}: ours - lab phase: overall {np.degrees(np.angle(dz.sum())):+.1f} deg, step across the {len(st)} burst switches "
          f"median {np.median(st):+.1f} (range {st.min():+.1f}..{st.max():+.1f}) deg")
    print(f"     coherence on the switch line / 20-40 lines away: ours {np.median(dips[:, 0]):.3f} / {np.median(dips[:, 1]):.3f}, "
          f"lab {np.median(dips[:, 2]):.3f} / {np.median(dips[:, 3]):.3f} (lab's worst line {dips[:, 2].min():.3f}, ours {dips[:, 0].min():.3f})")
