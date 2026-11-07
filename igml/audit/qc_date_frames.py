"""One date's frames and raw scenes, against a good neighbour's (2026-11-06; T071 2015-11-22: the last third of every
subswath holds other content, the lab's image too). Read only: the state's records, the pass-1 burst records
(OUT/T/esd/F<n>/<d>/frames.json), the frame folders, the raw scenes, the lab's pairs with the date, and the content's azimuth offset
against the neighbour with the range pattern removed (each range column's mean subtracted).
usage (seasat venv): python qc_date_frames.py T DATE NEIGHBOUR"""
import glob
import importlib.util
import json
import os
import sys

import numpy as np

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, d, x = sys.argv[1:4]
tr = pt.Track(T)
s = tr.state
for k in ("plan", "stitched", "rescued", "datatake", "ipf", "raw_mix", "burst_check", "assembly_check", "raw_bad", "recut"):
    for dd in (d, x):
        v = s.get(k, {}).get(dd) if isinstance(s.get(k), dict) else None
        if v is not None:
            print(f"{k} {dd}: {str(v)[:420]}")
print("stitch_fixed:", d in s.get("stitch_fixed", []), x in s.get("stitch_fixed", []))
for n in tr.subs:
    for dd in (d, x):
        f = tr.out / f"esd/F{n}/{dd}/frames.json"
        if f.exists():
            fr = json.loads(f.read_text())
            print(f"F{n} {dd} frames:", {k: {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items() if a in ("master", "tmp_da", "esd1_shift", "lines", "width")} for k, v in fr.items()})
print("\nframe folders holding the date:")
for dd in (d, x):
    fs = sorted(glob.glob(str(tr.labd / f"F*_F*/*_{dd}T*.SAFE"))) + sorted(glob.glob(str(tr.out / f"data_{T}/F*_F*/*_{dd}T*.SAFE")))
    print(f" {dd}: {len(fs)} SAFE links")
    for f in fs[:8]:
        print("   ", f.replace(str(tr.labd), "LAB").replace(str(tr.out), "OUT"), "->", os.path.realpath(f).split("/")[-1][:70] if os.path.islink(f) else "")
    fr = sorted(glob.glob(str(tr.labd / f"F*_F*/F1/raw/S1_{dd}_*_F1.PRM"))) + sorted(glob.glob(str(tr.out / f"data_{T}/F*_F*/F1/raw/S1_{dd}_*_F1.PRM")))
    for f in fr[:6]:
        p = pt._prm(f)
        print("   frame", f.replace(str(tr.labd), "LAB").replace(str(tr.out), "OUT").split("/raw/")[0], os.path.basename(f), "lines", p.get("num_lines"), "clock_start", p.get("clock_start"))
y, doy = d[:4], (pt.day(d) - pt.day(d[:4] + "0101")).days
tag = f"{y}{doy:03d}"
lp = sorted(p.name for p in (tr.labp / "merge").iterdir() if tag in p.name)
print(f"\nthe lab's merged pairs with {d} ({tag}): {len(lp)} {lp[:12]}")


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


RAW = tr.out / "superseded_nesd1/raw"
print(f"\ncontent of {d} against {x} (range pattern removed), best azimuth offset per 1000 lines, offset: corr; corr at 0 too:")
for n in tr.subs:
    fa = RAW / f"F{n}/S1_{d}_ALL_F{n}.SLC"
    fb = RAW / f"F{n}/S1_{x}_ALL_F{n}.SLC"
    a, nl = looks(fa, fa.with_suffix(".PRM"))
    b, _ = looks(fb, fb.with_suffix(".PRM"))
    H, W = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
    a, b = np.log1p(a[:H, :W]), np.log1p(b[:H, :W])
    a -= a.mean(axis=0, keepdims=True); b -= b.mean(axis=0, keepdims=True)
    step, srch = 50, 500                       # 1000 lines; offsets to +-10000 lines
    out = []
    for r in range(0, H - step + 1, step):
        u = a[r:r + step]; u = u - u.mean()
        best, off, c0 = -2.0, 0, np.nan
        for o in range(-srch, srch + 1):
            if r + o < 0 or r + o + step > H:
                continue
            v = b[r + o:r + o + step]; v = v - v.mean()
            den = np.sqrt((u * u).sum() * (v * v).sum())
            c = float((u * v).sum() / den) if den > 0 else -2.0
            if o == 0:
                c0 = c
            if c > best:
                best, off = c, o
        out.append(f"{r * 20}:{off * 20:+d}/{best:.2f}({c0:.2f})")
    print(f" F{n} ({nl} lines): " + " ".join(out), flush=True)
