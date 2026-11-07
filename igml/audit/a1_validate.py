"""A1 / A5 validation (2026-10-30): process_track.align()'s recipe for ONE subswath -- the super-master and a few dates,
their frames and orbits as the pipeline links them -- into /zdata/datasets/IGML/work/a1_test/<T>/<RUN>/F<n>, with a
chosen preproc script, make_s1a_tops build and IGML_A1. Nothing of the pipeline's state or outputs is touched.

usage (seasat, the venv): python a1_validate.py T n RUN SCRIPT BIN A1 DATE [DATE ...]
  SCRIPT  preproc_batch_tops.csh | preproc_batch_tops_esd.csh (mode 2 with esd_mode 1, the median)
  BIN     the folder put first on PATH: bin_eap (production), bin_a1, or an absolute path (the fork's install/bin)
  IGML_A3 and other variables in the environment are passed on
  A1      1 sets IGML_A1
"""
import glob
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pt)

T, n, run, script, bindir, a1 = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6] == "1"
dates = sys.argv[7:]
pt.log = lambda *a: print(*a, flush=True)
tr = pt.Track(T)
tr.save = lambda: None
w = Path(f"/zdata/datasets/IGML/work/a1_test/{T}/{run}/F{n}")
shutil.rmtree(w, ignore_errors=True)
w.mkdir(parents=True)
msafes = pt._one_cut_per_frame(pt.real(glob.glob(str(tr.labd / f"F*_F*/*_{tr.master}T*.SAFE"))))
morb = pt.orbit_for(tr, tr.master, Path(msafes[0]).name[:3])
lines = []
for d, safes, orb in [(tr.master, msafes, morb)] + [(d, tr.state["stitched"][d], pt._orbit_of(tr, d)) for d in dates]:
    rec = []
    for t, x in pt._frames_line(safes, n):
        os.symlink(t, w / t.name)
        if d != tr.master and (pt._ipf(Path(x).parent.parent) or pt.EAP_IPF) < pt.EAP_IPF and pt.EAP_FIX.is_dir():
            pt._eap_xml(x, w / x.name)
        else:
            os.symlink(x, w / x.name)
        rec.append(t.stem)
    if not (w / orb.name).exists():
        os.symlink(orb, w / orb.name)
    lines.append(":".join(rec + [orb.name]))
(w / "data.in").write_text("\n".join(lines) + "\n")
shutil.copy(tr.labp / "topo/dem.grd", w / "dem.grd")
bindir = bindir if bindir.startswith("/") else f"/zdata/datasets/IGML/gmtsar_fix/{bindir}"   # 10-31: or the fork's install
env = dict(os.environ, PATH=f"{bindir}:" +
           (f"{pt.STITCH_FIX}:" if pt.STITCH_FIX.is_dir() else "") + f"{pt.TOOLS}:/usr/local/GMTSAR/bin:{os.environ['PATH']}")
env.pop("IGML_A1", None)
if a1:
    env["IGML_A1"] = "1"
mode2 = f"{script} data.in dem.grd 2" + (" 1" if "esd" in script else "")
print(f"{T} F{n} {run}: {script} with {bindir}, IGML_A1={'1' if a1 else 'unset'}, dates {dates}", flush=True)
r = subprocess.run(f"preproc_batch_tops.csh data.in dem.grd 1 > pbt1.log 2>&1 && {mode2} > pbt2.log 2>&1",
                   shell=True, cwd=w, env=env)
print(f"exit {r.returncode}", flush=True)
for d in [tr.master] + dates:
    f = w / f"S1_{d}_ALL_F{n}.SLC"
    print(d, f.stat().st_size if f.exists() else "missing", flush=True)
