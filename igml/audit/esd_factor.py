"""(1) GMTSAR's ESD conversion: spec_sep = lpb x mean kt (make_s1a_tops mode 2), while the two looks of a target are one
burst CYCLE apart (EL - SL + 1 valid lines without the overlap): the shift is under-estimated by cycle / lpb.
(2) Per-overlap ESD of one pair from its ddphase, against the burst-boundary steps. 10-31."""
import re, subprocess, sys, tempfile, os
from pathlib import Path
import numpy as np

run = Path(sys.argv[1]).resolve(); pairs = [p.split(",") for p in sys.argv[2:]]
filt = "/zdata/datasets/IGML/gmtsar_igml/install/share/gmtsar/filters/gauss25x7"
for a, b in pairs:
    bb = [list(map(float, l.split())) for l in (run / f"{b}.BB").read_text().split("\n") if l.strip()]
    nb, width, spec_sep, dta = int(bb[0][0]), int(bb[0][1]), bb[0][2], bb[0][3]
    B = np.array(bb[1:], int)                        # SL SC SH EL EC EH per burst
    tiff = None
    cyc = B[1:-1, 3] - B[1:-1, 0] + 1                # middle bursts: EL - SL + 1
    kov = B[:-1, 5] - B[:-1, 3]                      # EH - EL: overlap with the next burst
    x = sorted(run.glob(f"*{b[3:11]}*{b[12:18]}*.xml")) or sorted(run.glob(f"*{b[3:11]}*.xml"))
    lpb = None
    for f in x:
        m = re.search(r"<linesPerBurst>(\d+)</linesPerBurst>", f.read_text(errors="replace"))
        if m:
            lpb = int(m.group(1)); break
    print(f"{b}: {nb} bursts, lpb {lpb}, cycle {np.median(cyc)} lines (range {cyc.min()}-{cyc.max()}), overlap {np.median(kov)}, "
          f"spec_sep {spec_sep:.1f} Hz -> factor lpb / cycle = {lpb / np.median(cyc):.4f}")
    with tempfile.TemporaryDirectory(dir="/zdata/datasets/IGML/work/a1_test") as td:
        for s in (a, b):
            for e in ("SLCH", "SLCL", "BB"):
                os.symlink(run / f"{s}.{e}", Path(td) / f"{s}.{e}")
        out = subprocess.run(["spectral_diversity", a, b, "0", filt], cwd=td, capture_output=True, text=True,
                             env=dict(os.environ, PATH="/zdata/datasets/IGML/gmtsar_igml/install/bin:" + os.environ["PATH"])).stdout
        sxd = float(re.search(r"spectral_spectrationXdta = ([-0-9.e]+)", out).group(1))
        d = np.loadtxt(Path(td) / "ddphase", ndmin=2)
    y, ph = d[:, 1], d[:, 2]
    gaps = np.where(np.diff(np.sort(np.unique(y))) > 50)[0]
    uy = np.sort(np.unique(y)); edges = [uy[0]] + [uy[g + 1] for g in gaps]
    print(f"  overall median {np.median(ph):+.4f} rad -> GMTSAR {np.median(ph) / (2 * np.pi * sxd):+.5f} lines, "
          f"cycle-corrected {np.median(ph) / (2 * np.pi * sxd) * lpb / np.median(cyc):+.5f}")
    for i, e0 in enumerate(edges):
        e1 = edges[i + 1] if i + 1 < len(edges) else uy[-1] + 1
        sel = (y >= e0) & (y < e1)
        print(f"   overlap at frame line ~{int(np.median(y[sel]))}: {sel.sum()} rows, ESD phase {np.degrees(np.median(ph[sel])):+6.1f} deg")
