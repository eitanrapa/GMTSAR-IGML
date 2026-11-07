"""Second test of filter.csh's amplitude cache (2026-11-02), on the production path: intf_tops_parallel.csh with the
pipeline's environment (proc_env, IGML_AMP_CACHE in a test folder) makes N pairs sharing a reference date, 2 at a time
-- the first two miss the cache at once (both store it: the rename race), the rest hit -- and every output grid is
compared byte for byte with the same pair made by the running pipeline without the cache (work/<T>/intf/F<n>/intf_all).
usage (seasat, the venv): python ampcache_par_test.py T n REF DATE [DATE ...]   (the patched filter.csh in ./bin)"""
import filecmp, importlib.util, os, shutil, subprocess, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, n, a, reps = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4:]
tr = pt.Track(T)
root = Path(__file__).resolve().parent
w = root / f"F{n}"
shutil.rmtree(w, ignore_errors=True); shutil.rmtree(root / f"cache_F{n}", ignore_errors=True)
(w / "raw").mkdir(parents=True); (w / "topo").mkdir()
for d in (a, *reps):
    for ext in ("SLC", "LED", "PRM"):
        f = tr.outp / f"F{n}/raw" / f"S1_{d}_ALL_F{n}.{ext}"
        (os.link if ext != "PRM" else shutil.copy)(f, w / "raw" / f.name)
shutil.copy(tr.labp / f"F{n}/topo/topo_ra.grd", w / "topo/topo_ra.grd")
shutil.copy(tr.labp / f"F{n}/batch_tops.config", w / "batch_tops.config")
(w / "intf.in").write_text("".join(f"S1_{a}_ALL_F{n}:S1_{x}_ALL_F{n}\n" for x in reps))
env = dict(pt.proc_env(tr, n), IGML_AMP_CACHE=str(root / f"cache_F{n}"))
env["PATH"] = f"{root}/bin:" + env["PATH"]
r = subprocess.run("intf_tops_parallel.csh intf.in batch_tops.config 2 > itp.log 2>&1", shell=True, cwd=w, env=env)
hits = sum(f.read_text().count("from the amplitude cache") for f in w.glob("intf_*.log"))
print(f"intf_tops_parallel.csh exit {r.returncode}; {len(reps)} pairs, cache hits {hits}; cache: "
      f"{sorted(f.name for f in (root / f'cache_F{n}').iterdir())}", flush=True)
bad = 0
for x in reps:
    p = f"{pt.doy0(a)}_{pt.doy0(x)}"
    mine, ref = w / "intf_all" / p, tr.work / f"intf/F{n}/intf_all" / p
    grds = sorted(g.name for g in ref.glob("*.grd"))
    same = [g for g in grds if (mine / g).exists() and filecmp.cmp(mine / g, ref / g, shallow=False)]
    diff = [g for g in grds if g not in same]
    bad += len(diff)
    print(f"{p}: {len(same)} of {len(grds)} grids byte-identical to the pipeline's" + (f"; DIFFER: {diff}" if diff else ""))
print("PASS" if bad == 0 and hits >= len(reps) - 2 else "FAIL")
