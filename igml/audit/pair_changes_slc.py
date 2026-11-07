"""Which of the rebuild's changes reached given pairs of a replaced track, at the level of the aligned images (2026-11-06;
user: "which of our changes affected the interferograms ..."). Read only. Lab = its originals in lab_originals_20261029,
ours = OUT/<T>/product_<T>/F<n>/raw (the lab's raw/ now holds copies of ours). For each date of each pair:
  the state's records (frames: the lab's or a cut of ours, how many; IPF: A2 needs < 002.43; raw_mix / re-cut: F6, A4);
  the lab's image against ours of the SAME date, per 50-line block (qc's measure): coherence and any line shift (A0, A4, F6);
  the PRM's differing alignment terms and whether the orbit (LED) is the same;
  the network's solved shift per frame (A5, A6).
For each pair, at the super-master's burst-switch lines, from the two images (2048-sample blocks across the swath):
  the phase step across the switch, lab and ours (A5, A6);
  the band of decorrelated lines at the switch, lab and ours: lines within +-12 whose coherence is < 0.6 of that 20-40 lines
  away, and the lowest coherence against it (A1).
usage (seasat venv): python pair_changes_slc.py T PAIR [PAIR ...]   (PAIR = yyyyddd_yyyyddd, day 0 = 1 January)"""
import hashlib
import importlib.util
import json
import statistics
import sys
from datetime import date, timedelta
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, pairs = sys.argv[1], sys.argv[2:]
tr = pt.Track(T)
st = tr.state
ARCH = pt.LAB_ARCHIVE / T / f"product_{T}"


def ymd(c):
    return (date(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d")


def lab(d, n):
    return ARCH / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"


def ours(d, n):
    return tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"


deg, sol = {}, {}
for n in tr.subs:
    N = json.loads((tr.out / f"esd/F{n}/nesd.json").read_text())
    sep = []
    for line in open(tr.out / f"esd/F{n}/pairs.jsonl"):
        r = json.loads(line)
        if r.get("spectral_spectrationXdta"):
            sep.append(r["spectral_spectrationXdta"])
    deg[n] = 360 * statistics.median(sep) / statistics.median(N["factors"].values())
    sol[n] = N["solution"]


def frames(d, n):
    return {k: v for k, v in sol[n].items() if k.split("_")[1] == d}


seen = set()
for pair in pairs:
    a, b = (ymd(c) for c in pair.split("_"))
    print(f"\n=== {pair} = {a} x {b}", flush=True)
    for d in (a, b):
        if d in seen:
            continue
        seen.add(d)
        fr = st.get("stitched", {}).get(d, [])
        print(f" {d}: {len(fr)} frame(s) {[('ours' if str(s).startswith(str(tr.out)) else 'lab') + ':' + Path(s).name[:3] + ':' + Path(s).parent.name for s in fr]}; "
              f"IPF {st.get('ipf', {}).get(d)}; raw_mix {st.get('raw_mix', {}).get(d)}; re-cut {d in st.get('recut', {})}; "
              f"in the lab's merged pairs {d in st.get('lab_dates', [])}", flush=True)
        for n in tr.subs:
            A, B = pt._blocks(lab(d, n), lab(d, n).with_suffix(".PRM")), pt._blocks(ours(d, n), ours(d, n).with_suffix(".PRM"))
            cov, coh, flo, blk = pt._qc_pair(A, B, shifts=True)
            shifted = {k: s for k, (c0, cs, s) in enumerate(blk) if cs > c0 + pt.QC_SHIFT_GAIN}
            pa, pb = pt._prm(lab(d, n).with_suffix(".PRM")), pt._prm(ours(d, n).with_suffix(".PRM"))
            keys = ("rshift", "sub_int_r", "stretch_r", "a_stretch_r", "ashift", "sub_int_a", "stretch_a", "a_stretch_a", "num_lines", "bytes_per_line", "SC_clock_start")
            dif = {k: (pa.get(k), pb.get(k)) for k in keys if pa.get(k) != pb.get(k)}
            led = [hashlib.md5(Path(str(f(d, n))[:-4] + ".LED").read_bytes()).hexdigest()[:8] for f in (lab, ours)]
            fs = frames(d, n)
            print(f"   F{n} lab vs ours, same date: coherence p90 {coh:.3f} (floor {flo:.3f}), blocks with a line shift {shifted or 'none'}; LED {'same' if led[0] == led[1] else 'DIFFER'}; "
                  f"network shift {[round(v['x'] * deg[n], 1) for v in fs.values()]} deg", flush=True)
            if dif:
                print(f"      PRM differs: {dif}", flush=True)


def profile(A, B, c, H=80):
    """Per block across the swath: the deramped interferogram around line c -> (phase step in deg, coherence per line)."""
    P = pt._prm(Path(str(A)[:-4] + ".PRM"))
    w = int(P["bytes_per_line"]) // 4
    nl = min(Path(A).stat().st_size, Path(B).stat().st_size) // (4 * w)
    ma = np.memmap(A, dtype="<i2", mode="r", shape=(nl, w, 2))
    mb = np.memmap(B, dtype="<i2", mode="r", shape=(nl, w, 2))
    a = ma[c - H:c + H].astype(np.float32); a = a[..., 0] + 1j * a[..., 1]
    b = mb[c - H:c + H].astype(np.float32); b = b[..., 0] + 1j * b[..., 1]
    ok = (np.abs(a) > 0) & (np.abs(b) > 0)
    cols = np.where(ok.mean(axis=0) > 0.99)[0]
    steps, prof = [], []
    lines = c - H + np.arange(2 * H)
    sel = ((lines >= c - 33) & (lines < c - 3)) | ((lines > c + 3) & (lines <= c + 33))
    for x0 in range(cols.min(), cols.max() - 2048, 2048):
        ifg = a[:, x0:x0 + 2048] * np.conj(b[:, x0:x0 + 2048])
        F = np.fft.fft2(ifg); ky, kx = np.unravel_index(np.argmax(np.abs(F)), F.shape)
        yy, xx = np.mgrid[0:ifg.shape[0], 0:ifg.shape[1]]
        ifg = ifg * np.exp(-2j * np.pi * (ky * yy / ifg.shape[0] + kx * xx / ifg.shape[1]))
        z = ifg.reshape(ifg.shape[0], -1, 8).sum(axis=2)
        p = np.exp(1j * np.angle(z)).mean(axis=1)
        if np.abs(p[sel]).mean() < 0.2:
            continue
        prof.append(np.abs(p))
        l, q = lines[sel], np.unwrap(np.angle(p[sel]))
        M = np.c_[np.ones(l.size), l - c, (l > c).astype(float)]
        steps.append(np.degrees(np.linalg.lstsq(M, q, rcond=None)[0][2]))
    if len(steps) < 3:
        return None
    return float(np.median(steps)), np.median(prof, axis=0), lines


for pair in pairs:
    a, b = (ymd(c) for c in pair.split("_"))
    print(f"\n=== {pair}: at the super-master's burst switches", flush=True)
    for n in tr.subs:
        cs = pt._master_cut_lines(tr, n)
        cs = [cs[int(i)] for i in np.linspace(0, len(cs) - 1, min(10, len(cs)))]
        out = {}
        for name, f in (("lab", lab), ("ours", ours)):
            st_, nb, lo = [], [], []
            for c in cs:
                if c < 200:
                    continue
                r = profile(f(a, n), f(b, n), c)
                if r is None:
                    continue
                s, prof, lines = r
                st_.append(s)
                side = np.median(prof[(np.abs(lines - c) >= 20) & (np.abs(lines - c) <= 40)])
                near = prof[np.abs(lines - c) <= 12]
                nb.append(int(np.sum(near < 0.6 * side)))
                lo.append(float(near.min() / side))
            out[name] = (st_, nb, lo)
        fa, fb = frames(a, n), frames(b, n)
        xa, xb = [v["x"] for v in fa.values()], [v["x"] for v in fb.values()]
        pred = [round((y - x) * deg[n], 1) for x, y in zip(xa, xb)]
        for name in ("lab", "ours"):
            st_, nb, lo = out[name]
            if st_:
                print(f"  F{n} {name:4s}: step median {np.median(st_):+6.1f} deg (|max| {np.max(np.abs(st_)):5.1f}, {len(st_)} switches); band lines "
                      f"{np.mean(nb):.1f} per switch (max {np.max(nb)}), lowest coherence {np.median(lo):.2f} of the level beside", flush=True)
        print(f"  F{n} the network's correction for this pair (b - a): {pred} deg", flush=True)
