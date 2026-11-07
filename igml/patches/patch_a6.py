"""GMTSAR audit A6 (2026-11-01): spectral diversity's phase-to-shift conversion over the burst CYCLE, not a whole burst.

Applied in the IGML fork (/zdata/datasets/IGML/gmtsar_igml/src, branch igml) after A5 (NESD). make_s1a_tops (mode 2)
writes spec_sep = linesPerBurst x mean kt x dta into the .BB file -- the Doppler span of a WHOLE burst -- and
spectral_diversity / preproc_batch_tops_esd.csh divide the double-difference phase by 2 pi spec_sep dta. A target's two
looks (end of burst k, start of burst k + 1) are one burst cycle apart (consecutive bursts' azimuthTime), so every ESD
shift came out cycle / linesPerBurst too small (T166 F1: 1342 / 1497, -10.3 %; the network test left -9.0 deg of a
-77.6 deg step, -1.1 deg once scaled).

The fix: spec_sep is taken over the burst cycle (consecutive bursts' azimuthTime over azimuthTimeInterval, averaged),
and the .BB header carries that cycle as a FIFTH number ("nc width spec_sep dta cycle"; spectral_diversity reads the
first four), so a reader knows the correction is in (process_track's nesd() then applies none of its own). A cycle
outside (0.5, 1) x linesPerBurst (a datatake across a year boundary in GMTSAR's yyyyddd times, a single burst) leaves
GMTSAR's value and the four-number header, with a message. Only mode 2 (the overlap images) changes; modes 0 / 1 / 3
are untouched.

usage: python patch_a6.py <fork source tree>   (in place; refuses to run twice)
"""
import sys
from pathlib import Path

root = Path(sys.argv[1])
mk = root / "preproc/S1A_preproc/src_tops/make_s1a_tops.c"
s = mk.read_text()
if "IGML A6" in s:
    sys.exit(f"{mk}: A6 already applied")


def sub(text, old, new):
    n = text.count(old)
    assert n == 1, f"expected 1 x {old[:80]!r}, found {n}"
    return text.replace(old, new)


s = sub(s, "void igml_a1_snap(struct burst_bounds *, int, int *, int *);\n",
        "void igml_a1_snap(struct burst_bounds *, int, int *, int *);\n"
        "double igml_burst_cycle(struct tree *, int, double); /* IGML A6 */\n")
s = sub(s, '/* IGML A1 (2026-10-30): the output stem, for <stem>.cuts */\nchar IGML_STEM[256] = "";\n',
        '/* IGML A1 (2026-10-30): the output stem, for <stem>.cuts */\nchar IGML_STEM[256] = "";\n'
        "/* IGML A6 (2026-11-01): the burst cycle (lines) mode 2's spec_sep was taken over (0: GMTSAR's, a whole burst) */\n"
        "double IGML_CYCLE = 0.0;\n")
s = sub(s, "\t\tspec_sep = spec_sep * dta / count / width;\n\t}\n",
        "\t\tspec_sep = spec_sep * dta / count / width;\n"
        "\t\t/* IGML A6 (2026-11-01): the two looks of a target are one burst CYCLE apart, not lpb lines: spec_sep over the\n"
        "\t\t   cycle (every ESD shift was cycle / lpb too small, T166 F1 -10.3 %) */\n"
        "\t\t{\n"
        "\t\t\tdouble cyc = igml_burst_cycle(xml_tree, count, dta);\n"
        "\t\t\tif (cyc > 0.5 * (double)lpb && cyc < (double)lpb) {\n"
        "\t\t\t\tspec_sep = spec_sep * cyc / (double)lpb;\n"
        "\t\t\t\tIGML_CYCLE = cyc;\n"
        "\t\t\t}\n"
        "\t\t\telse\n"
        "\t\t\t\tfprintf(stderr, \"IGML A6: burst cycle %.1f lines not usable (lpb %d): spec_sep over the whole burst\\n\", cyc, lpb);\n"
        "\t\t}\n"
        "\t}\n")
s = sub(s, '\t\tfprintf(BB, "%d %d %.6f %.12f\\n", nc,\n'
           '\t\t        IGML_W > 0 ? IGML_W : (int)str2double(tmp_str) - (int)str2double(tmp_str) % 4, spec_sep, dta); /* IGML A3 */\n',
        '\t\tif (IGML_CYCLE > 0.0) /* IGML A6: spec_sep is over this burst cycle (lines), the header\'s fifth number */\n'
        '\t\t\tfprintf(BB, "%d %d %.6f %.12f %.3f\\n", nc,\n'
        '\t\t\t        IGML_W > 0 ? IGML_W : (int)str2double(tmp_str) - (int)str2double(tmp_str) % 4, spec_sep, dta, IGML_CYCLE);\n'
        '\t\telse\n'
        '\t\t\tfprintf(BB, "%d %d %.6f %.12f\\n", nc,\n'
        '\t\t\t        IGML_W > 0 ? IGML_W : (int)str2double(tmp_str) - (int)str2double(tmp_str) % 4, spec_sep, dta); /* IGML A3 */\n')
s = s.rstrip("\n") + '''

/* IGML A6 (2026-11-01): the burst cycle in lines -- consecutive bursts' azimuthTime (GMTSAR's yyyyddd.fraction days)
   over azimuthTimeInterval, averaged; 0 with fewer than two bursts. Spectral diversity's Doppler separation is the
   Doppler rate times this interval, the time between a target's two looks. */
double igml_burst_cycle(struct tree *xml_tree, int count, double dta) {
	char tmp_c[1024];
	double t1, t2, cyc = 0.0;
	int kk;

	if (count < 2 || dta <= 0.0)
		return 0.0;
	for (kk = 1; kk < count; kk++) {
		search_tree(xml_tree, "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 2, 4, kk);
		t1 = str2double(tmp_c);
		search_tree(xml_tree, "/product/swathTiming/burstList/burst/azimuthTime/", tmp_c, 2, 4, kk + 1);
		t2 = str2double(tmp_c);
		cyc += (t2 - t1) * 86400.0 / dta;
	}
	return cyc / (double)(count - 1);
}
'''
mk.write_text(s)
print(f"A6 applied in {root}")
