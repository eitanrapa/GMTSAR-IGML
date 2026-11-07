"""Raw merged grids of the T173 super-master pair: which pixels would get an interpolated phase that survives the
re-mask (phase NaN but mask2 = 1), and how big the mask_def-only area is."""
import sys
import numpy as np

D = sys.argv[1]
sh = (4939, 2129)
ld = lambda n: np.fromfile(f"{D}/{n}.bin", dtype=np.float32).reshape(sh)
ph, co, mk = ld("M_phasefilt"), ld("M_corr"), ld("M_mask")
print("NaN: phase", np.isnan(ph).sum(), "corr", np.isnan(co).sum(), "mask", np.isnan(mk).sum())
print("mask values:", np.unique(mk[np.isfinite(mk)])[:5])
mask2 = np.isfinite(co) & (co >= 1e-4) & np.isfinite(mk)
print("phase NaN & mask2 valid (would leak):", (np.isnan(ph) & mask2).sum())
print("phase valid & mask2 NaN (masked out):", (np.isfinite(ph) & ~mask2).sum())
print("corr > 1:", (co > 1).sum(), " corr < 1e-4:", (co < 1e-4).sum())
try:
    md = ld("mask_def")
    lm = ld("landmask_patch")
    print("mask_def NaN", np.isnan(md).sum(), " landmask_patch NaN", np.isnan(lm).sum())
    only = np.isnan(md) & mask2 & np.isfinite(lm)
    print("mask_def-only pixels (valid data, on land, masked only by mask_def):", only.sum())
    r, c = np.nonzero(only)
    if r.size:
        print("  rows", r.min(), r.max(), " cols", c.min(), c.max())
        # where: histogram over columns in 100-col bins
        h, e = np.histogram(c, bins=np.arange(0, 2200, 100))
        print("  col hist:", dict(zip(e[:-1].astype(int), h)))
        h, e = np.histogram(r, bins=np.arange(0, 5000, 500))
        print("  row hist:", dict(zip(e[:-1].astype(int), h)))
    st = ld("ours_corr")
    print("stored corr NaN", np.isnan(st).sum(), "; NaN set == raw NaN | mask_def NaN ?",
          np.array_equal(np.isnan(st), np.isnan(md) | np.isnan(st)))
except FileNotFoundError as e:
    print("skip", e)
