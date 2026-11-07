"""Our patched copies of GMTSAR's interferogram / merge code (2026-10-29; user: "Ok let's do it" -- fix the processing
itself and replace the lab's stacks, originals archived). The lab's /usr/local/GMTSAR is not touched: src/*.orig.* are
its sources (master 3dd43d3), src/* the patched ones written here, built into bin_proc/ by build_proc.sh.

  conv.c       I2: each output cell's filter may be centred off the block start (IGML_CONV_DI / _DJ, input index units)
               and the output grid's origin shifted (IGML_CONV_Y0 / _X0, coordinate units); the input grid's own
               origin is carried to the output (GMTSAR wrote every grid from 0). No variable set = the old geometry.
  conv2d.c     I3: the filter window stopped at row ni / column nj INCLUSIVE (one past the data: the last column read
               the next row's first sample).
  filter.csh   uses those: the second (6 x 16) stage centres each cell on its block (+3 lines, +8 stage-1 columns =
               +16 samples: the half-cell of I2), the first stage starts a subswath at (IGML_SUB_Y0, IGML_SUB_X0) so its
               32 x 6 cells fall exactly on the merged grid's (M6), and the grids carry that origin.
  phasefilt.c  I1: the Goldstein patches now include one aligned with each far edge, and a cell no patch weights is
               NaN instead of an exact 0 (only the outermost row / column; they were 11-16 rows of zeros).
  merge_swath.c  M6: a subswath grid's origin (header x_min / y_min) is added to its PRM's near_range / clock_start
               before the whole-cell placement, which is then exact (checked, printed); M1: F2 no longer lands one
               column too far in the 3-subswath branch.
"""
from pathlib import Path

F = Path("/zdata/datasets/IGML/gmtsar_fix/src")


def patch(name, reps):
    src = F / name.replace(".", ".orig.", 1) if name.count(".") == 1 else None
    s = src.read_text()
    for old, new in reps:
        n = s.count(old)
        assert n == 1, (name, n, old[:80])
        s = s.replace(old, new)
    (F / name).write_text(s)
    print("patched", name)


patch("conv.c", [
    ("	struct GMT_GRID *In = NULL;  /* Grid structure containing ->header and ->data */\n",
     "	struct GMT_GRID *In = NULL;  /* Grid structure containing ->header and ->data */\n"
     "	/* IGML 2026-10-29: filter centre offsets (input index units) and output origin shift (coordinate units) */\n"
     "	int di = 0, dj = 0;\n"
     "	double x0 = 0.0, y0 = 0.0, xlo_in = 0.0, ylo_in = 0.0;\n"
     "	char *ev = NULL;\n"),
    ("		ymax = In->header->wesn[GMT_YHI];\n",
     "		ymax = In->header->wesn[GMT_YHI];\n"
     "		xlo_in = In->header->wesn[GMT_XLO]; /* IGML: carry the input's origin */\n"
     "		ylo_in = In->header->wesn[GMT_YLO];\n"),
    ("	if (verbose)\n		fprintf(stderr, \" read file %s %d %d \\n\", input_name, xdim, ydim);\n",
     "	if (verbose)\n		fprintf(stderr, \" read file %s %d %d \\n\", input_name, xdim, ydim);\n"
     "	if ((ev = getenv(\"IGML_CONV_DI\")))\n		di = atoi(ev);\n"
     "	if ((ev = getenv(\"IGML_CONV_DJ\")))\n		dj = atoi(ev);\n"
     "	if ((ev = getenv(\"IGML_CONV_X0\")))\n		x0 = atof(ev);\n"
     "	if ((ev = getenv(\"IGML_CONV_Y0\")))\n		y0 = atof(ev);\n"
     "	if (di < 0 || dj < 0)\n		die(\"IGML_CONV_DI / IGML_CONV_DJ must be >= 0\", \"\");\n"),
    ("	inc[GMT_X] = round(xmax / (double)jout);\n"
     "	inc[GMT_Y] = round(ymax / (double)iout);\n"
     "	jout = floor(xmax / inc[GMT_X]);\n"
     "	iout = floor(ymax / inc[GMT_Y]);\n",
     "	inc[GMT_X] = round((xmax - xlo_in) / (double)jout);\n"
     "	inc[GMT_Y] = round((ymax - ylo_in) / (double)iout);\n"
     "	jout = floor((xmax - xlo_in) / inc[GMT_X]);\n"
     "	iout = floor((ymax - ylo_in) / inc[GMT_Y]);\n"),
    ("	wesn[GMT_XLO] = 0.0;\n"
     "	wesn[GMT_XHI] = inc[GMT_X] * jout;\n"
     "	wesn[GMT_YLO] = 0.0;\n"
     "	wesn[GMT_YHI] = inc[GMT_Y] * iout;\n",
     "	wesn[GMT_XLO] = xlo_in + x0;\n"
     "	wesn[GMT_XHI] = wesn[GMT_XLO] + inc[GMT_X] * jout;\n"
     "	wesn[GMT_YLO] = ylo_in + y0;\n"
     "	wesn[GMT_YHI] = wesn[GMT_YLO] + inc[GMT_Y] * iout;\n"),
    ("	for (ic = 0, row = 0; ic < iout * idec; ic = ic + idec) {\n",
     "	for (ic = di, row = 0; ic < di + iout * idec; ic = ic + idec) {\n"),
    ("		for (jc = 0; jc < floor(xmax / inc[GMT_X]) * jdec; jc = jc + jdec) {\n",
     "		for (jc = dj; jc < dj + floor((xmax - xlo_in) / inc[GMT_X]) * jdec; jc = jc + jdec) {\n"),
])

patch("conv2d.c", [
    ("	i1 = min(*ni, *ic + nif2);\n", "	i1 = min(*ni - 1, *ic + nif2); /* IGML 2026-10-29: was *ni, one past the data */\n"),
    ("	j1 = min(*nj, *jc + njf2);\n", "	j1 = min(*nj - 1, *jc + njf2); /* IGML 2026-10-29: was *nj (read the next row) */\n"),
])

patch("phasefilt.c", [
    ("	float *outphase = NULL, *wgt = NULL, *amp = NULL, *corr = NULL, *diff = NULL, *ftmp = NULL;\n",
     "	float *outphase = NULL, *wgt = NULL, *amp = NULL, *corr = NULL, *diff = NULL, *ftmp = NULL;\n"
     "	float *wsum = NULL;                 /* IGML 2026-10-29: total patch weight per cell */\n"
     "	int *si = NULL, *sj = NULL, nsi = 0, nsj = 0, a, b, s;\n"),
    ("	make_wgt(wgt, nxp, nyp);\n\n"
     "	for (ii = 0; ii < (ydim - nyp); ii += nyp / 2) {\n"
     "		for (jj = 0; jj < (xdim - nxp); jj += nxp / 2) {\n",
     "	make_wgt(wgt, nxp, nyp);\n\n"
     "	/* IGML 2026-10-29: patch starts every nyp / 2 as before, plus one aligned with the far edge (the old loop\n"
     "	   stopped short of it and left the last 11-16 rows / columns at exactly 0 phase) */\n"
     "	if ((wsum = calloc(T->header->nm, sizeof(float))) == NULL)\n		die(\"error allocating memory\", \"\");\n"
     "	if ((si = malloc((ydim / (nyp / 2) + 3) * sizeof(int))) == NULL || (sj = malloc((xdim / (nxp / 2) + 3) * sizeof(int))) == NULL)\n"
     "		die(\"error allocating memory\", \"\");\n"
     "	for (s = 0; s < (int)ydim - (int)nyp; s += nyp / 2)\n		si[nsi++] = s;\n"
     "	if ((int)ydim >= (int)nyp)\n		si[nsi++] = (int)ydim - (int)nyp;\n"
     "	for (s = 0; s < (int)xdim - (int)nxp; s += nxp / 2)\n		sj[nsj++] = s;\n"
     "	if ((int)xdim >= (int)nxp)\n		sj[nsj++] = (int)xdim - (int)nxp;\n\n"
     "	for (a = 0; a < nsi; a++) {\n"
     "		ii = si[a];\n"
     "		for (b = 0; b < nsj; b++) {\n"
     "			jj = sj[b];\n"),
    ("					fdata[k1].r = fdata[k1].r + wgt[k2] * patch1[k2].r;\n"
     "					fdata[k1].i = fdata[k1].i + wgt[k2] * patch1[k2].i;\n",
     "					fdata[k1].r = fdata[k1].r + wgt[k2] * patch1[k2].r;\n"
     "					fdata[k1].i = fdata[k1].i + wgt[k2] * patch1[k2].i;\n"
     "					wsum[k1] = wsum[k1] + wgt[k2];\n"),
    ("	for (i = 0; i < T->header->nm; i++)\n"
     "		outphase[i] = atan2f(fdata[i].i, fdata[i].r);\n",
     "	for (i = 0; i < T->header->nm; i++)\n"
     "		outphase[i] = (wsum[i] > 0.0f) ? atan2f(fdata[i].i, fdata[i].r) : NAN; /* IGML: no patch, no phase */\n"),
])

patch("merge_swath.c", [
    ("	/* compute coefficients neede for merging*/\n",
     "	/* IGML 2026-10-29: a subswath grid that starts at (y_min, x_min) [lines, samples] covers its PRM's geometry\n"
     "	   from there: put that into near_range / clock_start, so the whole-cell placement below is exact */\n"
     "	{\n"
     "		double drr = c_speed / prm1.fs / 2.0;\n"
     "		struct GMT_GRID *GG[3] = {G1, G2, G3};\n"
     "		struct PRM *PP[3] = {&prm1, &prm2, &prm3};\n"
     "		int q;\n"
     "		for (q = 0; q < nfile; q++) {\n"
     "			PP[q]->near_range += GG[q]->header->wesn[GMT_XLO] * drr;\n"
     "			PP[q]->clock_start += GG[q]->header->wesn[GMT_YLO] / prm1.prf / 86400.0;\n"
     "			PP[q]->SC_clock_start += GG[q]->header->wesn[GMT_YLO] / prm1.prf / 86400.0;\n"
     "			printf(\"IGML merge_swath: grid %d origin x %g y %g\\n\", q + 1, GG[q]->header->wesn[GMT_XLO], GG[q]->header->wesn[GMT_YLO]);\n"
     "		}\n"
     "	}\n\n"
     "	/* compute coefficients neede for merging*/\n"),
    ("	ovl12 = G1->header->n_columns - (int)round((prm2.near_range - prm1.near_range) / (c_speed / prm1.fs / 2) / incx);\n",
     "	ovl12 = G1->header->n_columns - (int)round((prm2.near_range - prm1.near_range) / (c_speed / prm1.fs / 2) / incx);\n"
     "	printf(\"IGML merge_swath: F2 offset %.4f columns, %.4f rows\\n\", (prm2.near_range - prm1.near_range) / (c_speed / prm1.fs / 2) / incx,\n"
     "	       (prm2.clock_start - prm1.clock_start) * 86400.0 * prm1.prf / incy);\n"
     "	if (nfile == 3)\n"
     "		printf(\"IGML merge_swath: F3 offset %.4f columns, %.4f rows\\n\", (prm3.near_range - prm1.near_range) / (c_speed / prm1.fs / 2) / incx,\n"
     "		       (prm3.clock_start - prm1.clock_start) * 86400.0 * prm1.prf / incy);\n"),
    ("				k = (ii - head2) * G2->header->n_columns + jj - G1->header->n_columns + ovl12 - 1;\n",
     "				k = (ii - head2) * G2->header->n_columns + jj - G1->header->n_columns + ovl12; /* IGML M1: no - 1 */\n"),
])

# filter.csh: the two decimation stages
s = (F / "filter.orig.csh").read_text()
reps = [
    ("  echo \"$filter2 $idec $jdec ($az_lks $dec_rng)\" \n",
     "  echo \"$filter2 $idec $jdec ($az_lks $dec_rng)\" \n"
     "#\n# IGML 2026-10-29: centre each decimated cell on its block (the second stage, in stage-1 index units) and start\n"
     "# this subswath at (IGML_SUB_Y0 lines, IGML_SUB_X0 samples) so its cells fall on the merged grid's (first stage)\n#\n"
     "  set igx0 = 0\n  set igy0 = 0\n"
     "  if ($?IGML_SUB_X0) set igx0 = $IGML_SUB_X0\n  if ($?IGML_SUB_Y0) set igy0 = $IGML_SUB_Y0\n"
     "  @ igci = $idec / 2\n  @ igcj = $jdec / 2\n"
     "  if ($az_lks != 1) then\n    echo \"IGML filter.csh: az_lks $az_lks != 1 -- cell centring needs a one-line first stage\"\n    exit 1\n  endif\n"
     "  set st1 = \"env IGML_CONV_DI=$igy0 IGML_CONV_DJ=$igx0 IGML_CONV_Y0=$igy0 IGML_CONV_X0=$igx0\"\n"
     "  set st2 = \"env IGML_CONV_DI=$igci IGML_CONV_DJ=$igcj\"\n"
     "  echo \"IGML filter.csh: origin $igy0 lines, $igx0 samples; cell centres +$igci / +$igcj (stage 2)\"\n"),
]
for name in ("$1 amp1_tmp.grd=bf", "$2 amp2_tmp.grd=bf", "real.grd=bf real_tmp.grd=bf", "imag.grd=bf imag_tmp.grd=bf"):
    reps.append((f"  conv $az_lks $dec_rng $filter1 {name}\n", f"  $st1 conv $az_lks $dec_rng $filter1 {name}\n"))
for name in ("amp1_tmp.grd=bf amp1.grd", "amp2_tmp.grd=bf amp2.grd", "real_tmp.grd=bf realfilt.grd",
             "imag_tmp.grd=bf imagfilt.grd"):
    reps.append((f"  conv $idec $jdec $filter2 {name}\n", f"  $st2 conv $idec $jdec $filter2 {name}\n"))
for name in ("xt.grd=bf xreal.grd", "yt.grd=bf yreal.grd", "xt.grd=bf ximag.grd", "yt.grd=bf yimag.grd"):
    reps.append((f"    conv $idec $jdec $filter2 {name}\n", f"    $st2 conv $idec $jdec $filter2 {name}\n"))
for old, new in reps:
    assert s.count(old) == 1, ("filter.csh", s.count(old), old[:70])
    s = s.replace(old, new)
(F / "filter.csh").write_text(s)
print("patched filter.csh")
