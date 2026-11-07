"""Swath-wide burst-boundary steps of sample pairs, pass 1 (geometry-only, OUT/<T>/superseded_nesd1/raw) against
pass 2 (the network's shifts, the product's raw/), 2026-11-01: qc's records measure one 4096-sample window at mid-swath,
which overstates residuals; here the median over 2048-sample blocks across the swath at 6 of the super-master's
burst-switch lines per subswath (as steps_wide.py). Pairs: qc's own (each date against its qc neighbour) -- the worst
of pass 2 and a spread over years and satellite pairings.
usage (seasat, the venv): python steps_sample.py T [n_worst] [n_spread]"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pt)
T = sys.argv[1]
nw = int(sys.argv[2]) if len(sys.argv) > 2 else 6
ns = int(sys.argv[3]) if len(sys.argv) > 3 else 12
tr = pt.Track(T)
st = tr.state
sat = {k: Path(v[0]).name[:3] for k, v in st.get("stitched", {}).items() if v}
sat[tr.master] = "S1A"


def path(d, n, ps):
    if d == tr.master:
        return pt.slc_src(tr, d, n)
    return (tr.out / f"superseded_nesd1/raw/F{n}" / f"S1_{d}_ALL_F{n}.SLC") if ps == 1 else \
        (tr.outp / f"F{n}/raw" / f"S1_{d}_ALL_F{n}.SLC")


def steps(A, B, centres, H=80):
    P = pt._prm(Path(str(A)[:-4] + ".PRM")) if Path(str(A)[:-4] + ".PRM").exists() else pt._prm(pt.prm_src(tr, tr.master, 1))
    w = int(P["bytes_per_line"]) // 4
    nl = min(Path(A).stat().st_size, Path(B).stat().st_size) // (4 * w)
    ma = np.memmap(A, dtype="<i2", mode="r", shape=(nl, w, 2))
    mb = np.memmap(B, dtype="<i2", mode="r", shape=(nl, w, 2))
    out = []
    for c in centres:
        if not (H < c < nl - H):
            continue
        a = ma[c - H:c + H].astype(np.float32); a = a[..., 0] + 1j * a[..., 1]
        b = mb[c - H:c + H].astype(np.float32); b = b[..., 0] + 1j * b[..., 1]
        ok = (np.abs(a) > 0) & (np.abs(b) > 0)
        cols = np.where(ok.mean(axis=0) > 0.99)[0]
        if cols.size < 4096:
            continue
        blk = []
        for x0 in range(cols.min(), cols.max() - 2048, 2048):
            aa, bb = a[:, x0:x0 + 2048], b[:, x0:x0 + 2048]
            ifg = aa * np.conj(bb)
            F = np.fft.fft2(ifg); ky, kx = np.unravel_index(np.argmax(np.abs(F)), F.shape)
            yy, xx = np.mgrid[0:ifg.shape[0], 0:ifg.shape[1]]
            ifg = ifg * np.exp(-2j * np.pi * (ky * yy / ifg.shape[0] + kx * xx / ifg.shape[1]))
            z = ifg.reshape(ifg.shape[0], -1, 8).sum(axis=2)
            p = np.exp(1j * np.angle(z)).mean(axis=1)
            lines = c - H + np.arange(ifg.shape[0])
            sel = ((lines >= c - 33) & (lines < c - 3)) | ((lines > c + 3) & (lines <= c + 33))
            if np.abs(p[sel]).mean() < 0.2:
                continue
            l, q = lines[sel], np.unwrap(np.angle(p[sel]))
            M = np.c_[np.ones(l.size), l - c, (l > c).astype(float)]
            blk.append(np.degrees(np.linalg.lstsq(M, q, rcond=None)[0][2]))
        if len(blk) >= 3:
            out.append(float(np.median(blk)))
    return out


q2, q1 = st["qc"], st.get("qc_nesd1", {})
rec = []
for d, r in q2.items():
    for n in tr.subs:
        v = r.get(f"F{n}") or {}
        if v.get("steps"):
            rec.append((float(np.median(np.abs(v["steps"]))), d, v["vs"], n))
rec.sort(reverse=True)
pick = [x for x in rec[:nw]]
rng = np.random.default_rng(0)
by = {}
for x in rec:
    by.setdefault((x[1][:4], "-".join(sorted({sat.get(x[1], "?"), sat.get(x[2], "?")}))), []).append(x)
keys = sorted(by)
for k in rng.permutation(len(keys))[:ns]:
    g = by[keys[k]]
    pick.append(g[rng.integers(len(g))])
cuts = {n: pt._master_cut_lines(tr, n) for n in tr.subs}
print(f"{T}: pass 1 vs pass 2, swath-wide |step| (deg): median of 6 boundaries; qc's pass-2 window value for reference")
for qv, d, x, n in pick:
    cs = cuts[n]
    cs = [cs[int(i)] for i in np.linspace(0, len(cs) - 1, min(6, len(cs)))]
    try:
        s1 = steps(path(d, n, 1), path(x, n, 1), cs)
        s2 = steps(path(d, n, 2), path(x, n, 2), cs)
    except Exception as e:  # a missing file etc.
        print(f"  {d} x {x} F{n}: {str(e)[:80]}")
        continue
    f = lambda s: f"{np.median(np.abs(s)):5.1f} (max {np.max(np.abs(s)):5.1f}, n {len(s)})" if s else "  -  "
    print(f"  {d} x {x} F{n} {sat.get(d, '?')}-{sat.get(x, '?')}: pass 1 {f(s1)} | pass 2 {f(s2)} | qc window {qv:5.1f}",
          flush=True)
