"""Test of filter.csh's amplitude cache (2026-11-02): two pairs sharing a reference date, made one after the other with
the pipeline's setup (_intf_run) and IGML_AMP_CACHE (the second takes amp1 from the cache), against the same pairs the
running pipeline made without it (work/<T>/intf/F<n>/intf_all): every output grid compared byte for byte.
usage (seasat, the venv): python ampcache_test.py T n REF DATE_B DATE_C"""
import filecmp, importlib.util, os, shutil, subprocess, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, n, a, b, c = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
tr = pt.Track(T)
root = Path("/zdata/datasets/IGML/work/ampcache_test")
w = root / f"F{n}"
shutil.rmtree(w, ignore_errors=True)
(w / "raw").mkdir(parents=True); (w / "topo").mkdir()
for d in (a, b, c):
    for ext in ("SLC", "LED", "PRM"):
        f = tr.outp / f"F{n}/raw" / f"S1_{d}_ALL_F{n}.{ext}"
        (os.link if ext != "PRM" else shutil.copy)(f, w / "raw" / f.name)
shutil.copy(tr.labp / f"F{n}/topo/topo_ra.grd", w / "topo/topo_ra.grd")
shutil.copy(tr.labp / f"F{n}/batch_tops.config", w / "batch_tops.config")
env = dict(pt.proc_env(tr, n), IGML_AMP_CACHE=str(root / f"cache_F{n}"))
env["PATH"] = f"{root}/bin:" + env["PATH"]
for x in (b, c):
    (w / "intf.in").write_text(f"S1_{a}_ALL_F{n}:S1_{x}_ALL_F{n}\n")
    r = subprocess.run("intf_tops.csh intf.in batch_tops.config > intf_%s.log 2>&1" % x, shell=True, cwd=w, env=env)
    print(f"pair {a} x {x}: exit {r.returncode}; cache hits: "
          f"{(w / f'intf_{x}.log').read_text().count('from the amplitude cache')}", flush=True)
for x in (b, c):
    p = f"{pt.doy0(a)}_{pt.doy0(x)}"
    mine, ref = w / "intf_all" / p, tr.work / f"intf/F{n}/intf_all" / p
    grds = sorted(g.name for g in ref.glob("*.grd"))
    same = [g for g in grds if (mine / g).exists() and filecmp.cmp(mine / g, ref / g, shallow=False)]
    diff = [g for g in grds if g not in same]
    print(f"{p}: {len(same)} of {len(grds)} grids byte-identical to the pipeline's" + (f"; DIFFER: {diff}" if diff else ""))
