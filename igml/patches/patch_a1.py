"""GMTSAR audit A1 (2026-10-30): every date cut between overlapping bursts at the super-master's lines.

Applied to src/make_s1a_tops.c AFTER patch_eap.py (the bin_eap build). make_s1a_tops writes each burst's lines
[SC, EC] one after another; a date's cut between bursts k-1 and k therefore falls at its OWN cumulative count, a few
lines from the super-master's (each date's burst timing differs by 2-9 lines), and between two dates' cut lines a pair
compares burst k of one with burst k+1 of the other (~4.7 kHz apart in Doppler): coherence ~0.1 over 2-9 lines at
every burst boundary, in every pair (T166 F1 2026-05-29 x 06-11; the lab's T071 2019-11-01 x 11-25 the same).

With IGML_A1 set in the environment:
  - without shift tables (the super-master's frames in preproc_batch_tops(_esd).csh mode 2), the time of each interior
    cut is written to <stem>.cuts (PRM clock_start units: day of year and fraction);
  - with them (every other date), the cut times in ./*.cuts are placed on this frame's output grid through ./tmp.PRM
    (the master reference preproc writes before the call: the super-master's stitched PRM shifted by the frame's
    lines, or, in the ESD script, the matching super-master frame's PRM), and each interior cut moves to the nearest
    one if burst k holds the lines before its SC (or burst k-1 after its EC) that the move takes; otherwise it stays.
Each burst's resampling (al_start, the shift tables) is the date's own placement whatever the cuts, so the image
changes only between the old and the new cut lines. Without IGML_A1 the program behaves exactly as before.

usage: python patch_a1.py src/make_s1a_tops.c   (in place; refuses to run twice)
"""
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()
if "IGML A1" in s:
    sys.exit(f"{p}: A1 already applied")


def sub(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, f"expected {count} x {old[:70]!r}, found {n}"
    s = s.replace(old, new)


sub("#include <string.h>\n",
    "#include <string.h>\n"
    "#include <glob.h>\n"
    "\n"
    "/* IGML A1 (2026-10-30): the output stem, for <stem>.cuts */\n"
    "char IGML_STEM[256] = \"\";\n"
    "void igml_a1_snap(struct burst_bounds *, int, int *, int *);\n")
# the prototype needs struct burst_bounds: move it below the typedef
s = s.replace("void igml_a1_snap(struct burst_bounds *, int, int *, int *);\n", "", 1)
sub("int pop_burst(struct PRM *, struct tree *, struct burst_bounds *, char *, char *);\n",
    "int pop_burst(struct PRM *, struct tree *, struct burst_bounds *, char *, char *);\n"
    "void igml_a1_snap(struct burst_bounds *, int, int *, int *);\n")

sub("\tnull_MEM_STR();\n\timode = atoi(argv[4]);\n",
    "\tnull_MEM_STR();\n\tstrncpy(IGML_STEM, argv[3], 255); /* IGML A1 */\n\timode = atoi(argv[4]);\n")

sub('\tfprintf(stderr, "Working on burst ");\n',
    "\t/* IGML A1 (2026-10-30): c_own[k] = the first output line of burst k as GMTSAR places it (the date's own cuts),\n"
    "\t   o_start[k] = where it is written from (the super-master's cut when IGML_A1 moves it; see patch_a1.py) */\n"
    "\tint a1 = getenv(\"IGML_A1\") != NULL && (imode == 1 || imode == 3);\n"
    "\tint *c_own = (int *)malloc((count + 2) * sizeof(int));\n"
    "\tint *o_start = (int *)malloc((count + 2) * sizeof(int));\n"
    "\tc_own[1] = 0;\n"
    "\tfor (kk = 1; kk <= count; kk++)\n"
    "\t\tc_own[kk + 1] = c_own[kk] + (bb[kk].EC - bb[kk].SC + 1);\n"
    "\tfor (kk = 1; kk <= count + 1; kk++)\n"
    "\t\to_start[kk] = c_own[kk];\n"
    "\tif (a1 && dr_table[0] != '\\0' && da_table[0] != '\\0')\n"
    "\t\tigml_a1_snap(bb, count, c_own, o_start);\n"
    "\n"
    '\tfprintf(stderr, "Working on burst ");\n')

sub("\t\t\tal_start = cl - bb[kk].SC;\n",
    "\t\t\tal_start = c_own[kk] - bb[kk].SC; /* IGML A1: = cl - SC as before; the date's own placement */\n")

sub("ii >= bb[kk].SC && ii <= bb[kk].EC",
    "ii >= bb[kk].SC + o_start[kk] - c_own[kk] && ii <= bb[kk].EC + o_start[kk + 1] - c_own[kk + 1]", count=2)

sub('\tfprintf(stderr, "\\n");\n\n\tif (imode == 2) {\n',
    '\tfprintf(stderr, "\\n");\n'
    "\n"
    "\t/* IGML A1: the super-master's frames (no shift tables) record their cut times */\n"
    "\tif (a1 && dr_table[0] == '\\0' && imode == 1 && IGML_STEM[0] != '\\0') {\n"
    "\t\tchar fn[300];\n"
    "\t\tFILE *fc;\n"
    "\t\tsnprintf(fn, sizeof(fn), \"%s.cuts\", IGML_STEM);\n"
    "\t\tif ((fc = fopen(fn, \"w\")) != NULL) {\n"
    "\t\t\tfor (kk = 2; kk <= count; kk++)\n"
    "\t\t\t\tfprintf(fc, \"%.12f\\n\", prm->clock_start + (double)c_own[kk] / prm->prf / 86400.0);\n"
    "\t\t\tfclose(fc);\n"
    "\t\t}\n"
    "\t}\n"
    "\tfree(c_own);\n"
    "\tfree(o_start);\n"
    "\n"
    "\tif (imode == 2) {\n")

s += r'''
/* IGML A1 (2026-10-30): each interior cut c_own[k] (k = 2 .. count) moved to the nearest super-master cut, placed on
   this frame's output lines through ./tmp.PRM (clock_start, PRF), when burst k holds the lines before its SC (a cut
   moved earlier) or burst k - 1 the lines after its EC (moved later) that the move takes, with 2 lines to spare and at
   most 30 lines; otherwise the cut stays where GMTSAR puts it. */
void igml_a1_snap(struct burst_bounds *bb, int count, int *c_own, int *o_start) {
	FILE *f;
	char line[512], key[128];
	double v, t0 = -1.0, prf = 0.0, t, m[8192], bd, d;
	int nm = 0, k, j, dl, room;
	glob_t g;

	if ((f = fopen("tmp.PRM", "r")) == NULL) {
		fprintf(stderr, "IGML A1: no tmp.PRM -- every cut left as GMTSAR puts it\n");
		return;
	}
	while (fgets(line, sizeof(line), f) != NULL)
		if (sscanf(line, " %127s = %lf", key, &v) == 2) {
			if (strcmp(key, "clock_start") == 0)
				t0 = v;
			else if (strcmp(key, "PRF") == 0)
				prf = v;
		}
	fclose(f);
	if (t0 < 0.0 || prf <= 0.0) {
		fprintf(stderr, "IGML A1: tmp.PRM has no clock_start / PRF -- every cut left as GMTSAR puts it\n");
		return;
	}
	if (glob("*.cuts", 0, NULL, &g) != 0) {
		fprintf(stderr, "IGML A1: no super-master cuts (*.cuts) -- every cut left as GMTSAR puts it\n");
		return;
	}
	for (j = 0; j < (int)g.gl_pathc; j++) {
		if ((f = fopen(g.gl_pathv[j], "r")) == NULL)
			continue;
		while (nm < 8192 && fscanf(f, "%lf", &t) == 1)
			m[nm++] = (t - t0) * 86400.0 * prf;
		fclose(f);
	}
	globfree(&g);
	for (k = 2; k <= count; k++) {
		bd = 1.0e9;
		for (j = 0; j < nm; j++) {
			d = m[j] - (double)c_own[k];
			if (fabs(d) < fabs(bd))
				bd = d;
		}
		dl = (int)floor(bd + 0.5);
		room = dl < 0 ? bb[k].SC - bb[k].SL : bb[k - 1].EH - bb[k - 1].EC;
		if (fabs(bd) <= 30.0 && abs(dl) <= room - 2) {
			o_start[k] = c_own[k] + dl;
			fprintf(stderr, "IGML A1: cut %d at line %d -> %d (%+d; super-master cut %+.3f)\n", k, c_own[k], o_start[k], dl,
			        bd);
		}
		else
			fprintf(stderr, "IGML A1: cut %d kept at line %d (nearest super-master cut %+.1f lines, room %d)\n", k,
			        c_own[k], bd, room);
	}
}
'''
p.write_text(s)
print(f"{p}: A1 applied")
