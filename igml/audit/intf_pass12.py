"""One pair's interferogram made twice with the pipeline's own setup (_intf_run: raw/, topo/topo_ra.grd, batch_tops.config,
intf.in, intf_tops.csh with proc_env), from pass 1's images (geometry-only, OUT/<T>/superseded_nesd1) and from pass 2's
(network spectral diversity, the product's raw/), into /zdata/datasets/IGML/work/visible_step/<pass>/F<n> (2026-11-02;
user: "I can't really see the difference"). Read only on the pipeline's files (hard links).
usage (seasat, the venv): python intf_pass12.py T DATE_A DATE_B n"""
import importlib.util, os, shutil, subprocess, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T, a, b, n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
tr = pt.Track(T)
jobs = []
for ps in (1, 2):
    w = Path(f"/zdata/datasets/IGML/work/visible_step/pass{ps}/F{n}")
    shutil.rmtree(w, ignore_errors=True)
    (w / "raw").mkdir(parents=True); (w / "topo").mkdir()
    for d in (a, b):
        src = (tr.out / f"superseded_nesd1/raw/F{n}") if ps == 1 else (tr.outp / f"F{n}/raw")
        for ext in ("SLC", "LED", "PRM"):
            f = src / f"S1_{d}_ALL_F{n}.{ext}"
            (os.link if ext != "PRM" else shutil.copy)(f, w / "raw" / f.name)
    shutil.copy(tr.labp / f"F{n}/topo/topo_ra.grd", w / "topo/topo_ra.grd")
    shutil.copy(tr.labp / f"F{n}/batch_tops.config", w / "batch_tops.config")
    (w / "intf.in").write_text(f"S1_{a}_ALL_F{n}:S1_{b}_ALL_F{n}\n")
    jobs.append(subprocess.Popen("intf_tops.csh intf.in batch_tops.config > intf.log 2>&1", shell=True, cwd=w,
                                 env=pt.proc_env(tr, n)))
for j in jobs:
    j.wait()
for ps in (1, 2):
    w = Path(f"/zdata/datasets/IGML/work/visible_step/pass{ps}/F{n}")
    print(ps, sorted(str(p.relative_to(w)) for p in w.glob("intf*/*/phasefilt.grd")), flush=True)
