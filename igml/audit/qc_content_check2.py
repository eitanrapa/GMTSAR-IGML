"""T071's pass-1 set-asides, second look (2026-11-06; qc_content_check.py: 24 of the 25 'content' dates are open water).
Read only, on the pass-1 images (OUT/T/superseded_nesd1/raw/F<n>/) and the lab's.
  1. The dates that failed on coherence: per subswath, coherence p90 and floor against the 3 nearest dates either side,
     whether or not those passed (qc's own neighbours were the last date that passed and a far re-cut date).
  2. 2015-11-22, the one date whose blocks really hold other content: per block, amplitude correlation with its
     neighbours, ours and the lab's image; where the content went (best azimuth offset of each stretch of the image
     against a neighbour, on 20 x 40-look amplitude).
usage (seasat venv): python qc_content_check2.py T"""
import importlib.util
import sys

import numpy as np

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T = sys.argv[1]
tr = pt.Track(T)
RAW = tr.out / "superseded_nesd1/raw"
q = tr.state["qc_nesd1"]
ALL = sorted(set(q) | {tr.master})
LAB = set(tr.state["lab_dates"])


def files(d, n):
    if d == tr.master:
        p = tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
        return p, pt.prm_src(tr, d, n)
    p = RAW / f"F{n}/S1_{d}_ALL_F{n}.SLC"
    return p, p.with_suffix(".PRM")


def lab_files(d, n):
    p = tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
    return (p, pt.prm_src(tr, d, n)) if p.exists() and p.stat().st_size > 0 else (None, None)


print(f"{T}: 1. dates that failed on coherence; (date) p90 / floor; * = qc set it aside", flush=True)
for d in ("20150514", "20220922", "20250503", "20250515"):
    i = ALL.index(d)
    near = ALL[max(0, i - 3):i] + ALL[i + 1:i + 4]
    print(f"{d} (in the lab's dates: {d in LAB}):", flush=True)
    for n in tr.subs:
        A = pt._blocks(*files(d, n))
        out = []
        for x in near:
            cov, coh, flo = pt._qc_pair(A, pt._blocks(*files(x, n)))
            out.append(f"{x[2:]}{'*' if not q.get(x, {'ok': 1}).get('ok', 1) else ''} {coh:.2f}/{flo:.2f}")
        print(f"  F{n}: " + "  ".join(out), flush=True)

d, nb = "20151122", ["20151029", "20151216"]
nb = [x for x in nb if x in ALL] or ["20151029"]
print(f"\n2. {d}: per block, amplitude correlation; ours vs {nb}, lab's image vs ours / vs the same neighbours (lab's)", flush=True)
for n in tr.subs:
    A = pt._blocks(*files(d, n))
    N = [pt._blocks(*files(x, n)) for x in nb]
    lf, lp = lab_files(d, n)
    L = pt._blocks(lf, lp) if lf else None
    print(f" F{n} lab image {'present' if lf else 'ABSENT'}: block: ours-nb ... | lab-ours | lab-nb", flush=True)
    for k in range(len(A)):
        s = " ".join(f"{pt._amp_corr(A[k], B[k]):+.2f}" for B in N)
        if L:
            s += f" | {pt._amp_corr(L[k], A[k]):+.2f} | " + " ".join(f"{pt._amp_corr(L[k], B[k]):+.2f}" for B in N)
        print(f"   {k:2d}: {s}", flush=True)


def looks(path, prm, ly=20, lx=40):
    p = pt._prm(prm)
    w, nl = int(p["bytes_per_line"]) // 4, int(p["num_lines"])
    m = np.memmap(path, dtype="<i2", mode="r", shape=(nl, w, 2))
    rows = []
    for r in range(0, nl - ly, ly):
        a = np.asarray(m[r:r + ly, ::2, :], dtype=np.float32)
        z = np.hypot(a[..., 0], a[..., 1])
        z = z[:, : z.shape[1] // (lx // 2) * (lx // 2)].reshape(ly, -1, lx // 2).mean(axis=(0, 2))
        rows.append(z)
    return np.array(rows, dtype=np.float32), nl


print(f"\n   where each stretch of {d} sits against {nb[0]} (best azimuth offset in lines, corr), per 1000 lines:", flush=True)
for n in tr.subs:
    a, nl = looks(*files(d, n))
    b, _ = looks(*files(nb[0], n))
    H = min(a.shape[0], b.shape[0]); W = min(a.shape[1], b.shape[1])
    a, b = np.log1p(a[:H, :W]), np.log1p(b[:H, :W])
    step, srch = 50, 200                       # 50 x 20 = 1000 lines; offsets to +-4000 lines
    line = []
    for r in range(0, H - step, step):
        x = a[r:r + step] - a[r:r + step].mean()
        best, off = -2.0, 0
        for o in range(-srch, srch + 1):
            if r + o < 0 or r + o + step > H:
                continue
            y = b[r + o:r + o + step]
            y = y - y.mean()
            den = np.sqrt((x * x).sum() * (y * y).sum())
            c = float((x * y).sum() / den) if den > 0 else -2.0
            if c > best:
                best, off = c, o
        line.append((r * 20, off * 20, best))
    print(f" F{n} ({nl} lines): " + " ".join(f"{l}:{o:+d}/{c:.2f}" for l, o, c in line), flush=True)
