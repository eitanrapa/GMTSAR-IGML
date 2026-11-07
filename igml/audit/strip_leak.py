"""Does phasefilt's zero strip (GMTSAR audit I1) change SNAPHU's answer OUTSIDE the strip? (2026-10-28)
Per pair: copies of the merged phasefilt / corr / mask and the track's landmask_ra / mask_def into
work/chk/strip_leak/<T>/<pair>/{ctrl,test}; test = phasefilt and corr NaN where phasefilt == 0 exactly (the strip: then
interpolated for SNAPHU by nearest_grid and masked in unwrap.grd); both unwrapped with snaphu_interp.csh 0.0001 0 (SMOOTH,
our recipe). Then test - ctrl outside the strip, after removing its median: the share of pixels a cycle or more apart
and how far from the strip they lie. Read only on the sources.

usage (seasat): python strip_leak.py T PAIR_DIR [PAIR_DIR ...]     (runs the pairs one after another)
"""
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

T, pairs = sys.argv[1], [Path(p) for p in sys.argv[2:]]
VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
masks = Path(f"/zdata/datasets/seasatHD{VOL[T]}/InSAR_Processing/{T}/product_{T}/merge")
ENV = {"PATH": "/usr/local/GMTSAR/bin:/usr/local/bin:/usr/bin:/bin", "HOME": str(Path.home())}


def load(path):
    f = subprocess.run(["gmt", "grdinfo", "-C", str(path)], capture_output=True, text=True, check=True, env=ENV).stdout.split()
    nx, ny = int(f[9]), int(f[10])
    raw = subprocess.run(["gmt", "grd2xyz", str(path), "-ZTLf"], capture_output=True, check=True, env=ENV).stdout
    return np.frombuffer(raw, dtype=np.float32).reshape(ny, nx)


for src in pairs:
    base = Path(f"/zdata/datasets/IGML/work/chk/strip_leak/{T}/{src.name}")
    for arm in ("ctrl", "test"):
        w = base / arm
        shutil.rmtree(w, ignore_errors=True)
        w.mkdir(parents=True)
        for f in ("phasefilt.grd", "corr.grd", "mask.grd"):
            shutil.copy(src / f, w / f)
        for f in ("landmask_ra.grd", "mask_def.grd"):
            if (masks / f).exists():
                shutil.copy(masks / f, w / f)
        if arm == "test":                     # + (NaN where phasefilt == 0, else 0); corr first, phasefilt last
            for f in ("corr.grd", "phasefilt.grd"):
                subprocess.run(f"gmt grdmath {f} phasefilt.grd 99 DENAN 0 EQ 1 NAN ADD = tmp.grd && mv tmp.grd {f}",
                               shell=True, cwd=w, check=True, env=ENV, capture_output=True)
        r = subprocess.run("snaphu_interp.csh 0.0001 0 > log_unwrap 2>&1", shell=True, cwd=w, env=ENV)
        if r.returncode or not (w / "unwrap.grd").exists():
            sys.exit(f"{src.name} {arm}: snaphu_interp failed (see {w}/log_unwrap)")
    p0 = load(src / "phasefilt.grd")
    strip = np.isfinite(p0) & (p0 == 0)
    rows = np.flatnonzero(strip.sum(1) > 0.5 * np.isfinite(p0).sum(1).clip(1))     # of the row's FINITE pixels
    uc, ut = load(base / "ctrl/unwrap.grd"), load(base / "test/unwrap.grd")
    ok = np.isfinite(uc) & np.isfinite(ut) & ~strip
    d = (ut - uc)[ok]
    med = np.median(d)
    off = np.abs(d - med) > np.pi
    yy = np.nonzero(ok)[0][off]
    dist = np.min(np.abs(yy[:, None] - rows[None, :]), axis=1) if rows.size and yy.size else np.array([])
    us = load(src / "unwrap.grd")
    same = np.isfinite(us) & np.isfinite(uc)
    print(f"{T} {src.name}: strip rows {rows.min() if rows.size else '-'}..{rows.max() if rows.size else '-'} "
          f"({rows.size}); test - ctrl median {med:+.3f} rad; pixels >= 1 cycle apart outside the strip "
          f"{off.sum()} of {ok.sum()} ({100 * off.mean():.4f} %), rows from the strip: "
          f"{'median %.0f, max %.0f' % (np.median(dist), dist.max()) if dist.size else '-'}; "
          f"ctrl vs the stored unwrap.grd: max |diff| {np.nanmax(np.abs(uc - us)[same]):.3g} rad", flush=True)
