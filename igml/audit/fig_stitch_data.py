"""Data for the stitch_tops figure (2026-10-29; user: "can you show me the stitch_tops problem?"), from the 10-28 stitch
test of T173 F1: work/chk/stitchtest_T173_F1 (GMTSAR's stitch_tops) and stitchtest_T173_F1_fix (ours), every date
aligned to the super-master 2016-10-06. 2026-01-29 is two frames whose overlap correlation found ashift -1.
  self          2026-01-29 GMTSAR x conj(2026-01-29 fixed), looks of 8 lines x 16 samples: the same scenes stitched two
                ways (phase, coherence)
  coh_g, coh_f  coherence of the 12-day pair 2026-01-17 / 2026-01-29, GMTSAR / fixed, in 4 x 4 windows (B_perp 214 m
                puts a flat-earth fringe every ~30 samples: wider range windows decorrelate), averaged onto the looks
  prof_g/_f     their median over range, per 4 lines
  shift_g/_f    per block (NB blocks, 400 lines from each block's middle, samples R0..R1), the 12-day coherence median
                with 2026-01-29 moved by SHIFTS lines: line k + s of 2026-01-29 against line k of 2026-01-17
Read only. Runs on seasat next to the SLCs (numpy only; the Mac's NFS read is too slow for 7 GB).

usage: python fig_stitch_data.py OUT.npz      (IGML_ROOT=/Volumes/seasat2 on the Mac; default /zdata/datasets)
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("IGML_ROOT", "/zdata/datasets")) / "IGML/work/chk"
G, F = ROOT / "stitchtest_T173_F1", ROOT / "stitchtest_T173_F1_fix"
W, NL = 85312 // 4, 28292
LY, LX, CH, NB, NS = 8, 16, 512, 24, 400
SHIFTS = (-2, -1, 0, 1, 2)
R0, R1 = 4000, 16000
warnings.filterwarnings("ignore", category=RuntimeWarning)          # all-NaN medians off the bursts


def slc(p):
    return np.memmap(p, dtype="<i2", mode="r", shape=(NL, W, 2))


def cplx(m, a, b, c0=0, c1=W):
    x = np.asarray(m[a:b, c0:c1], dtype=np.float32)
    return x[..., 0] + 1j * x[..., 1]


def looks(x, ly, lx):
    n, m = x.shape[0] // ly * ly, x.shape[1] // lx * lx
    return x[:n, :m].reshape(n // ly, ly, m // lx, lx).sum(axis=(1, 3))


def coh(a, b, ly=4, lx=4):
    num = looks(a * b.conj(), ly, lx)
    den = np.sqrt(looks(np.abs(a) ** 2, ly, lx) * looks(np.abs(b) ** 2, ly, lx))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, np.abs(num) / den, np.nan)


if __name__ == "__main__":
    s17, sg, sf = slc(G / "S1_20260117_ALL_F1.SLC"), slc(G / "S1_20260129_ALL_F1.SLC"), slc(F / "S1_20260129_ALL_F1.SLC")
    assert np.array_equal(slc(F / "S1_20260117_ALL_F1.SLC")[::997], s17[::997])     # the partner is the same in both
    n = NL // LY * LY
    self_, coh_g, coh_f, prof_g, prof_f = [], [], [], [], []
    for a in range(0, n, CH):
        b = min(a + CH, n)
        r, g, f = cplx(s17, a, b), cplx(sg, a, b), cplx(sf, a, b)
        num = looks(g * f.conj(), LY, LX)
        den = np.sqrt(looks(np.abs(g) ** 2, LY, LX) * looks(np.abs(f) ** 2, LY, LX))
        with np.errstate(invalid="ignore", divide="ignore"):
            self_.append(np.where(den > 0, num / den, np.nan).astype(np.complex64))
        for out, prof, x in ((coh_g, prof_g, g), (coh_f, prof_f, f)):
            c = coh(x, r)
            prof.append(np.nanmedian(c, axis=1))
            m = c[: c.shape[0] // 2 * 2, : c.shape[1] // 4 * 4]
            with np.errstate(invalid="ignore"):
                out.append(np.nanmean(m.reshape(m.shape[0] // 2, 2, m.shape[1] // 4, 4), axis=(1, 3)).astype(np.float32))
        print(f"lines {a}-{b}", flush=True)
    shift_g, shift_f = np.full((NB, len(SHIFTS)), np.nan), np.full((NB, len(SHIFTS)), np.nan)
    mid = (np.arange(NB) + 0.5) * NL / NB
    for i, m in enumerate(mid):
        a = int(m) - NS // 2
        r = cplx(s17, a, a + NS, R0, R1)
        for out, mm in ((shift_g, sg), (shift_f, sf)):
            x = cplx(mm, a - 2, a + NS + 2, R0, R1)
            for j, s in enumerate(SHIFTS):
                out[i, j] = np.nanmedian(coh(x[2 + s: 2 + s + NS], r))
    np.savez_compressed(sys.argv[1], self=np.concatenate(self_), coh_g=np.concatenate(coh_g),
                        coh_f=np.concatenate(coh_f), prof_g=np.concatenate(prof_g), prof_f=np.concatenate(prof_f),
                        shift_g=shift_g, shift_f=shift_f, shift_mid=mid, shifts=np.array(SHIFTS), ly=LY, lx=LX)
    print("wrote", sys.argv[1])
