"""GMTSAR audit A3 (2026-10-31): a resampled date keeps its far-range columns where it is narrower than the super-master.

Applied in the IGML fork (/zdata/datasets/IGML/gmtsar_igml/src, branch igml) after A1. preproc_batch_tops(_esd).csh
builds the shift tables over the date's OWN width (rmax = its num_rng_bins) and make_s1a_tops writes the resampled date
at that width; stitch_tops then pads to the super-master's. Where the date's range shift is negative its data reach
past its own width on the super-master's grid, and those columns are dropped: T144 2025-07-07 (S1C) F3, shift -109:
data to column ~24,036, image ends at 23,927 (64-157 columns on T144's S1C and T064's S1A dates).

With IGML_A3 set in the environment:
  - both scripts take rmax = the larger of the date's and the super-master's (tmp.PRM) width;
  - make_s1a_tops writes a date resampled with tables wider than its tiff at the tables' width: the tiff's lines are
    zero-padded, the PRM's num_rng_bins / bytes_per_line, compute_eap's ramp and mode 2's boundary file take that
    width. Where the widths are equal nothing changes (byte-identical); a date wider than the super-master keeps its
    own width as before.

usage: python patch_a3.py <fork source tree>   (in place; refuses to run twice)
"""
import sys
from pathlib import Path

root = Path(sys.argv[1])
mk = root / "preproc/S1A_preproc/src_tops/make_s1a_tops.c"
s = mk.read_text()
if "IGML A3" in s:
    sys.exit(f"{mk}: A3 already applied")


def sub(text, old, new, count=1):
    n = text.count(old)
    assert n == count, f"expected {count} x {old[:70]!r}, found {n}"
    return text.replace(old, new)


s = sub(s, '/* IGML A1 (2026-10-30): the output stem, for <stem>.cuts */\nchar IGML_STEM[256] = "";\n',
        '/* IGML A1 (2026-10-30): the output stem, for <stem>.cuts */\nchar IGML_STEM[256] = "";\n'
        '/* IGML A3 (2026-10-31): the width a resampled date is written at, when wider than its tiff (else 0) */\n'
        'int IGML_W = 0;\n')
s = sub(s, "\twidth = widthi - widthi % 4;\n\twidth2 = 2 * width;\n",
        "\twidth = widthi - widthi % 4;\n\twidth2 = 2 * width;\n"
        "\t/* IGML A3 (2026-10-31): a resampled date is written at the shift tables' width where that is the larger\n"
        "\t   (preproc's rmax: the super-master's width with IGML_A3), its tiff's lines zero-padded */\n"
        "\tint width_in = width;\n"
        "\tif (getenv(\"IGML_A3\") != NULL && dr_table[0] != '\\0' && da_table[0] != '\\0') {\n"
        "\t\tint wt = (int)floor(R->header->wesn[GMT_XHI] + 0.5);\n"
        "\t\twt -= wt % 4;\n"
        "\t\tif (wt > width) {\n"
        "\t\t\tfprintf(stderr, \"IGML A3: written %d samples wide (its own %d)\\n\", wt, width);\n"
        "\t\t\twidth = wt;\n"
        "\t\t\twidth2 = 2 * width;\n"
        "\t\t\tIGML_W = width;\n"
        "\t\t\tprm->num_rng_bins = width;\n"
        "\t\t\tprm->bytes_per_line = 4 * width;\n"
        "\t\t\tprm->good_bytes = prm->bytes_per_line;\n"
        "\t\t}\n"
        "\t}\n")
s = sub(s, "\t\t\t\tbrst[k] = (short)buf[jj];\n",
        "\t\t\t\tbrst[k] = jj < 2 * width_in ? (short)buf[jj] : 0; /* IGML A3: zero beyond the tiff */\n")
s = sub(s, "\tn_samples = n_samples - n_samples % 4;\n",
        "\tn_samples = n_samples - n_samples % 4;\n"
        "\tif (IGML_W > 0)\n"
        "\t\tn_samples = IGML_W; /* IGML A3: the width the resampled date is written at */\n")
s = sub(s, '\t\tfprintf(BB, "%d %d %.6f %.12f\\n", nc, (int)str2double(tmp_str) - (int)str2double(tmp_str) % 4, spec_sep, dta);\n',
        '\t\tfprintf(BB, "%d %d %.6f %.12f\\n", nc,\n'
        '\t\t        IGML_W > 0 ? IGML_W : (int)str2double(tmp_str) - (int)str2double(tmp_str) % 4, spec_sep, dta); /* IGML A3 */\n')
mk.write_text(s)

old = "          set rmax = `grep num_rng_bins $stem.PRM | awk '{print $3}'`\n"
new = (old +
       "          # IGML A3 (2026-10-31): the shift tables span the super-master's width where it is the larger, so a\n"
       "          # date narrower than it keeps its far-range columns (make_s1a_tops writes at the tables' width)\n"
       "          if ($?IGML_A3) then\n"
       "            set rmaxm = `grep num_rng_bins tmp.PRM | awk '{print $3}'`\n"
       "            if ($rmaxm > $rmax) set rmax = $rmaxm\n"
       "          endif\n")
for name in ("preproc_batch_tops.csh", "preproc_batch_tops_esd.csh"):
    f = root / "gmtsar/csh" / name
    f.write_text(sub(f.read_text(), old, new))
print(f"A3 applied in {root}")
