"""Per-overlap ESD phase of date pairs from a track's kept burst overlaps (OUT/<T>/esd/F<n>/<date>/), frame by frame:
does the azimuth misregistration drift within a frame (a constant per frame would leave part of it)? 11-01.
usage (seasat, the venv): python peroverlap.py ESD_DIR n DATE1 DATE2 [DATE1 DATE2 ...]"""
import json, os, re, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

esd, n = Path(sys.argv[1]), int(sys.argv[2])
filt = "/zdata/datasets/IGML/gmtsar_igml/install/share/gmtsar/filters/gauss25x7"
env = dict(os.environ, PATH="/zdata/datasets/IGML/gmtsar_igml/install/bin:" + os.environ["PATH"])
for d1, d2 in zip(sys.argv[3::2], sys.argv[4::2]):
    f1 = json.loads((esd / f"F{n}" / d1 / "frames.json").read_text())
    f2 = json.loads((esd / f"F{n}" / d2 / "frames.json").read_text())
    for a, va in sorted(f1.items()):
        b = [s for s, vb in f2.items() if vb["master"] == va["master"]][0]
        with tempfile.TemporaryDirectory(dir="/zdata/datasets/IGML/work") as td:
            for s, dd in ((a, d1), (b, d2)):
                for e in ("SLCH", "SLCL", "BB"):
                    os.symlink(esd / f"F{n}" / dd / f"{s}.{e}", Path(td) / f"{s}.{e}")
            subprocess.run(["spectral_diversity", a, b, "0", filt], cwd=td, capture_output=True, env=env)
            x = np.loadtxt(Path(td) / "ddphase", ndmin=2)
        y, ph = x[:, 1], x[:, 2]
        uy = np.sort(np.unique(y)); cut = [uy[0]] + [uy[i + 1] for i in np.where(np.diff(uy) > 50)[0]] + [uy[-1] + 1]
        meds = [np.degrees(np.median(ph[(y >= cut[i]) & (y < cut[i + 1])])) for i in range(len(cut) - 1)]
        k = np.arange(len(meds)); slope = np.polyfit(k, meds, 1)[0] if len(meds) > 2 else 0
        print(f"F{n} {d1} x {d2} frame {a[12:18]}: {len(meds)} overlaps, ESD deg " + " ".join(f"{m:+.1f}" for m in meds)
              + f" | spread {max(meds) - min(meds):.1f}, trend {slope * (len(meds) - 1):+.1f} over the frame", flush=True)
