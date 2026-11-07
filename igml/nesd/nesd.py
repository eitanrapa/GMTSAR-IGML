#!/usr/bin/env python3
"""Network enhanced spectral diversity (NESD) for Sentinel-1 TOPS stacks processed with GMTSAR.

GMTSAR's preproc_batch_tops_esd.csh estimates each date's azimuth shift against the super-master alone. Years from the
super-master the burst overlaps have decorrelated and the estimate is poor. Here each date is measured against the
next K dates, where the overlaps are still coherent, and one shift per date frame is solved from all those pairs, with the
super-master fixed at zero. The network idea follows Fattahi, Agram & Simons (2017), IEEE TGRS 55(2), 777-786; the choices
here (K, one constant per frame, Huber weighting, GMTSAR's own spectral_diversity as the measurement) are this code's.

One subswath at a time. Each date is aligned twice:

  1. Pass 1, geometry only, keeping the burst overlaps (in each alignment folder):
       IGML_ESD_COLLECT=1 preproc_batch_tops_esd.csh data.in dem.grd 2 1
       nesd.py collect ALIGN_DIR ESD_DIR --master YYYYMMDD
  2. The network: measure the pairs, solve, write ESD_DIR/shifts.txt ("<frame stem> <shift in lines>"):
       nesd.py run ESD_DIR --master YYYYMMDD [--dates good_dates.txt]
  3. Pass 2, every date aligned again with its shift:
       IGML_ESD_SHIFTS=ESD_DIR/shifts.txt preproc_batch_tops_esd.csh data.in dem.grd 2 1

Make interferograms only from pass 2's images, and only from one network's images: one shift per date cancels around any
triangle, so closure is untouched, but a pair made from another solution's images does not close with these.

Needs numpy and, on PATH, spectral_diversity and gmtsar_sharedir.csh from this fork. preproc_batch_tops_esd.csh's two
variables and the A6 scale correction (the .BB header's fifth number) are in the fork's commits 94beffb and 90b6b68.

ESD_DIR layout: <date>/<stem>.{SLCH,SLCL,BB,sd} + <date>/frames.json (stem -> its super-master frame, the burst offset
tmp_da, the estimate against the super-master); pairs.jsonl (every pair measured, a cache); shifts.txt; nesd.json (the
solution, the residuals, a summary).

This is process_track.py's nesd() and helpers without the pipeline around them. Checked against the production runs: the solve
reproduces T071's and T173's shifts.txt on all three subswaths, and the measurement reproduces their pairs.jsonl records."""
import argparse
import json
import re
import shutil
import statistics
import subprocess
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np

K = 4                    # each date is paired with its next K dates
MIN_ROWS = 1000          # ddphase rows (coherent overlap samples) a pair needs to be used
MIN_PHASE = 0.02         # rad: floor of a pair's phase sigma (the formal one ignores the ionosphere and systematics)
FILTER = "gauss25x7"     # as preproc_batch_tops_esd.csh
OVERLAP_EXT = ("SLCH", "SLCL", "BB")
_BURST = re.compile(r"<burst>\s*<azimuthTime>([^<]+)</azimuthTime>")


def parse_sd(text):
    """spectral_diversity's report (and the collecting run's esd<n>_shift line) as numbers."""
    out = {}
    for line in text.split("\n"):
        m = re.match(r"\s*(residual_phase|spectral_spectrationXdta|residual_shift|esd\d_shift)\s*=\s*([-+0-9.eE]+)", line)
        if m:
            out[m.group(1)] = float(m.group(2))
        m = re.match(r"\s*Image analyzed (\d+)x(\d+)", line)
        if m:
            out["width"], out["lines"] = int(m.group(1)), int(m.group(2))
    return out


def collect(work, esd, master, dates=None):
    """Move what IGML_ESD_COLLECT left in an alignment folder into ESD_DIR: each date's frames, frames.json written last so a
    folder with frames.json is complete; the super-master's overlaps once (copied: every alignment folder has them)."""
    work, esd = Path(work), Path(esd)
    rows = {}
    for line in (work / "nesd_frames.txt").read_text().split("\n"):
        p = line.split()
        if len(p) == 3:
            rows[p[0]] = {"master": p[1], "tmp_da": int(p[2])}
    dm = esd / master
    if not (dm / "frames.json").exists():
        ms = [s.strip() for s in (work / "masterlist").read_text().split("\n") if s.strip()]
        dm.mkdir(parents=True, exist_ok=True)
        for s in ms:
            for ext in OVERLAP_EXT:
                shutil.copy(work / f"{s}.{ext}", dm / f"{s}.{ext}")
        (dm / "frames.json").write_text(json.dumps({s: {"master": s, "tmp_da": 0} for s in ms}, indent=1))
    found = sorted({s.split("_")[1] for s in rows})
    for d in (found if dates is None else dates):
        fr = {s: dict(v) for s, v in rows.items() if s.startswith(f"S1_{d}_")}
        if not fr:
            raise RuntimeError(f"no burst overlaps recorded for {d} in {work / 'nesd_frames.txt'}")
        dd = esd / d
        shutil.rmtree(dd, ignore_errors=True)
        dd.mkdir(parents=True)
        for s, v in fr.items():
            for ext in OVERLAP_EXT + ("sd",):
                shutil.move(str(work / f"{s}.{ext}"), str(dd / f"{s}.{ext}"))
            v.update(parse_sd((dd / f"{s}.sd").read_text()))
        (dd / "frames.json").write_text(json.dumps(fr, indent=1))
    return found


def measure_pair(a, da, b, db, filt, scratch):
    """spectral_diversity between frames a and b (stems; their folders da, db) aligned to one super-master frame: b's azimuth
    shift relative to a [lines], the median of ddphase / (2 pi spec_sep) as preproc_batch_tops_esd.csh's esd_mode 1, and sigma
    from the phases' spread (MAD / sqrt(n / 25): ddphase rows are 10 lines apart, a 25-sample filter along range)."""
    wd = Path(scratch) / f"{a}__{b}"
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    for stem, src in ((a, Path(da)), (b, Path(db))):
        for ext in OVERLAP_EXT:
            (wd / f"{stem}.{ext}").symlink_to(src.resolve() / f"{stem}.{ext}")
    r = subprocess.run(["spectral_diversity", a, b, "0", filt], cwd=wd, capture_output=True, text=True)
    rec = {"a": a, "b": b, **parse_sd(r.stdout)}
    dd, sep = wd / "ddphase", rec.get("spectral_spectrationXdta")
    rec["n"] = 0
    if sep and dd.exists() and dd.stat().st_size > 0:
        x = np.loadtxt(dd, ndmin=2)
        if len(x):
            ph = x[:, 2]
            med = float(np.median(ph))
            mad = 1.4826 * float(np.median(np.abs(ph - med)))
            rec.update(n=int(len(ph)), corr=round(float(np.mean(x[:, 3])), 4), phase_med=med, phase_mad=mad,
                       shift=med / (2 * np.pi * sep),
                       sigma=max(mad / np.sqrt(max(len(ph) / 25.0, 1.0)), MIN_PHASE) / (2 * np.pi * sep))
    shutil.rmtree(wd, ignore_errors=True)
    return rec


def burst_cycle(xml):
    """(linesPerBurst, burst cycle in lines) of a TOPS annotation. Only for overlaps collected before the A6 correction (a .BB
    header of four numbers): spec_sep then spans a whole burst, but a target's two looks are one burst CYCLE apart, so
    every shift comes out cycle / linesPerBurst too small, and the factor to multiply by is linesPerBurst / cycle."""
    t = Path(xml).read_text(errors="replace")
    lpb = int(re.search(r"<linesPerBurst>(\d+)</linesPerBurst>", t).group(1))
    dta = float(re.search(r"<azimuthTimeInterval>([-+0-9.eE]+)</azimuthTimeInterval>", t).group(1))
    bt = [datetime.fromisoformat(m.group(1)) for m in _BURST.finditer(t)]
    return lpb, statistics.median((b - a).total_seconds() for a, b in zip(bt, bt[1:])) / dta


def frame_factor(esd, date, stem, old):
    """1 for a frame whose .BB carries the A6 correction (five header numbers); else the factor given for it."""
    bb = Path(esd) / date / f"{stem}.BB"
    if bb.exists() and len(bb.read_text().split("\n", 1)[0].split()) >= 5:
        return 1.0
    if stem not in old:
        raise RuntimeError(f"{bb} predates the A6 correction: give its factor (--factors; burst_cycle gives linesPerBurst / cycle)")
    return float(old[stem])


def solve(master, obs):
    """Per date frame x [lines], the correction its a.grd takes, with x = 0 at the super-master's frame `master`: weighted
    least squares on the pairs (a, b, y, s): y = x_b - x_a, Huber reweighting (k = 1.345 robust sigma) to convergence. Only
    frames connected to the super-master through pairs are solved. Returns {stem: (x, sx, npairs)}, the residuals
    [(a, b, r, z)] and the robust scale of z."""
    adj = defaultdict(set)
    for a, b, _, _ in obs:
        adj[a].add(b)
        adj[b].add(a)
    seen, todo = {master}, [master]
    while todo:
        for v in adj[todo.pop()]:
            if v not in seen:
                seen.add(v)
                todo.append(v)
    obs = [o for o in obs if o[0] in seen and o[1] in seen]
    nodes = sorted(seen - {master})
    if not nodes or not obs:
        return {}, [], float("nan")
    k = {s: i for i, s in enumerate(nodes)}
    A = np.zeros((len(obs), len(nodes)))
    y, s = np.array([o[2] for o in obs]), np.array([o[3] for o in obs])
    for i, (a, b, _, _) in enumerate(obs):
        if b != master:
            A[i, k[b]] += 1
        if a != master:
            A[i, k[a]] -= 1
    w0, wt, scale = 1 / s ** 2, np.ones(len(obs)), 1.0
    for _ in range(50):
        W = w0 * wt
        N = A.T @ (A * W[:, None])
        x = np.linalg.solve(N, A.T @ (W * y))
        z = (y - A @ x) * np.sqrt(w0)
        scale = max(1.4826 * float(np.median(np.abs(z))), 1e-6)
        new = np.minimum(1.0, 1.345 * scale / np.maximum(np.abs(z), 1e-12))
        if np.max(np.abs(new - wt)) < 1e-4:
            break
        wt = new
    sx = np.sqrt(np.diag(np.linalg.inv(N))) * scale
    npr = np.abs(A).sum(axis=0)
    sol = {v: (float(x[k[v]]), float(sx[k[v]]), int(npr[k[v]])) for v in nodes}
    res = [(a, b, float(r), float(zz)) for (a, b, _, _), r, zz in zip(obs, y - A @ x, z)]
    return sol, res, scale


def load_network(esd, master, dates):
    """Every date's frames.json, and the frame pairs of the network: each date with its next K dates, frames on one super-master
    frame and not a whole burst out of step with it (|tmp_da| < 1000 lines)."""
    esd = Path(esd)
    if dates is None:
        dates = sorted(p.name for p in esd.iterdir() if (p / "frames.json").exists() and p.name != master)
    L = sorted(set(dates) | {master})
    fr = {d: json.loads((esd / d / "frames.json").read_text()) for d in L}
    jobs = []
    for i, a in enumerate(L):
        for b in L[i + 1:i + 1 + K]:
            for sb, vb in fr[b].items():
                sa = [s for s, va in fr[a].items() if va["master"] == vb["master"] and abs(va["tmp_da"]) < 1000]
                if abs(vb["tmp_da"]) < 1000 and len(sa) == 1:
                    jobs.append((a, sa[0], b, sb))
    return list(dates), L, fr, jobs


def read_cache(path):
    have = {}
    if path.exists():
        for line in path.read_text().split("\n"):
            if line.strip():
                r = json.loads(line)
                have[(r["a"], r["b"])] = r
    return have


def measure(esd, master, dates, out, procs=16, limit=0):
    """Measure the network's pairs not yet in OUT/pairs.jsonl."""
    esd, out = Path(esd), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    _, L, _, jobs = load_network(esd, master, dates)
    cache = out / "pairs.jsonl"
    have = read_cache(cache)
    todo = [j for j in jobs if (j[1], j[3]) not in have]
    if limit:
        todo = todo[:limit]
    print(f"{len(jobs)} frame pairs ({len(todo)} to measure) among {len(L)} dates, K = {K}", flush=True)
    share = subprocess.run("gmtsar_sharedir.csh", shell=True, capture_output=True, text=True).stdout.strip()
    filt = f"{share}/filters/{FILTER}"
    lock = threading.Lock()
    scratch = out / "scratch"

    def one(j):
        a, sa, b, sb = j
        r = measure_pair(sa, esd / a, sb, esd / b, filt, scratch)
        with lock:
            with open(cache, "a") as fh:
                fh.write(json.dumps(r) + "\n")
        return r

    with ThreadPoolExecutor(max_workers=max(1, procs)) as ex:
        for i, _ in enumerate(ex.map(one, todo), 1):
            if i % 200 == 0:
                print(f"{i} of {len(todo)} pairs measured", flush=True)
    shutil.rmtree(scratch, ignore_errors=True)


def network(esd, master, dates, out, old_factors=None):
    """Solve per super-master frame and write OUT/shifts.txt and OUT/nesd.json. Frames the network cannot reach keep their
    estimate against the super-master, flagged as fallbacks."""
    esd, out = Path(esd), Path(out)
    dates, L, fr, jobs = load_network(esd, master, dates)
    have = read_cache(out / "pairs.jsonl")
    fac = {s: frame_factor(esd, d, s, old_factors or {}) for d in L for s in fr[d]}
    cur = {(j[1], j[3]) for j in jobs}                 # pairs of this network only (the cache may hold more)
    smaster = {s: v["master"] for d in L for s, v in fr[d].items()}
    sol, resid, fallback, scales = {}, [], {}, {}
    for m in sorted(fr[master]):
        obs = [(r["a"], r["b"], r["shift"] * fac[r["b"]], r["sigma"] * fac[r["b"]]) for k, r in have.items()
               if k in cur and smaster.get(k[0]) == m and r.get("n", 0) >= MIN_ROWS and "shift" in r]
        s, res, scale = solve(m, obs)
        sol.update(s)
        resid += res
        scales[m] = scale
    rows = []
    for d in dates:
        for s, v in fr[d].items():
            if s in sol:
                x, sx, npr = sol[s]
            else:
                x, sx, npr = v.get("esd1_shift", 0.0) * fac[s], None, 0
                fallback[s] = x
            rows.append((s, x, sx, npr, v["esd1_shift"] * fac[s] if "esd1_shift" in v else None))
    (out / "shifts.txt").write_text("".join(f"{s} {x:.9f}\n" for s, x, _, _, _ in rows))
    sep = statistics.median([r["spectral_spectrationXdta"] for r in have.values() if r.get("spectral_spectrationXdta")] or [1.0])
    deg = 360.0 * sep / statistics.median(fac.values())     # lines -> degrees of the burst-boundary step
    rr = [abs(r) for _, _, r, _ in resid]
    dif = [x - e for s, x, sx, npr, e in rows if e is not None and s not in fallback]
    summ = {"pairs": len(cur), "used": len(resid), "frames": len(rows), "fallback": len(fallback), "robust_z_scale": scales,
            "resid_median_deg": round(statistics.median(rr) * deg, 2) if rr else None,
            "resid_p95_deg": round(sorted(rr)[int(0.95 * (len(rr) - 1))] * deg, 2) if rr else None,
            "x_minus_vs_master_median_deg": round(statistics.median(dif) * deg, 2) if dif else None,
            "x_range_deg": [round(min(x for _, x, _, _, _ in rows) * deg, 1), round(max(x for _, x, _, _, _ in rows) * deg, 1)] if rows else None}
    (out / "nesd.json").write_text(json.dumps(
        {"summary": summ, "factors": fac,
         "solution": {s: {"x": x, "sx": sx, "pairs": npr, "vs_master": e} for s, x, sx, npr, e in rows},
         "fallback": fallback, "residuals": [{"a": a, "b": b, "r": r, "z": z} for a, b, r, z in resid]}, indent=1))
    print(json.dumps(summ), flush=True)
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect", help="move what IGML_ESD_COLLECT left in an alignment folder into ESD_DIR")
    c.add_argument("align_dir"); c.add_argument("esd_dir"); c.add_argument("--master", required=True)
    for name, hlp in (("measure", "measure the pairs"), ("solve", "solve and write shifts.txt"), ("run", "measure, then solve")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("esd_dir"); p.add_argument("--master", required=True, help="the super-master date, YYYYMMDD")
        p.add_argument("--dates", help="file of the dates to use, one per line (those that passed your qc); default every date in ESD_DIR")
        p.add_argument("--out", help="where pairs.jsonl, shifts.txt and nesd.json go (default ESD_DIR)")
        p.add_argument("--procs", type=int, default=16)
        p.add_argument("--factors", help="JSON {stem: factor} for overlaps collected before the A6 correction")
        p.add_argument("--limit", type=int, default=0, help="measure at most this many pairs (a test)")
    a = ap.parse_args()
    if a.cmd == "collect":
        print("collected", collect(a.align_dir, a.esd_dir, a.master))
        return
    dates = [d for d in Path(a.dates).read_text().split() if d] if a.dates else None
    out = a.out or a.esd_dir
    old = json.loads(Path(a.factors).read_text()) if a.factors else None
    if a.cmd in ("measure", "run"):
        measure(a.esd_dir, a.master, dates, out, a.procs, a.limit)
    if a.cmd in ("solve", "run"):
        network(a.esd_dir, a.master, dates, out, old)


if __name__ == "__main__":
    main()
