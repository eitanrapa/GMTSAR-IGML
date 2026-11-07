"""GMTSAR audit A5, network spectral diversity (2026-10-31): preproc_batch_tops_esd.csh collects what a network of
date pairs needs, and later applies the network's shifts.

Applied in the IGML fork (/zdata/datasets/IGML/gmtsar_igml/src, branch igml) after A3. GMTSAR's spectral diversity
estimates each date's residual azimuth shift against the SUPER-MASTER (preproc_batch_tops_esd.csh: make_s1a_tops mode 2
writes the burst-overlap images, spectral_diversity the double-difference phase, a.grd += the median). Against a
super-master 11 years away it half-works (T166 F1 2026-05-29 x 06-11: steps -75..-96 deg -> -53..-64 in frame 1). The
same program pairs any two dates aligned to one super-master frame (the overlaps are matched by their line on that
frame's grid), so the shifts can be estimated between dates close in time and solved per date (process_track nesd).

With IGML_ESD_COLLECT set (pass 1):
  - every frame's overlap images and boundary file (<stem>.SLCH / .SLCL / .BB, the super-master's too) are kept;
  - each date frame's estimate against the super-master is recorded (<stem>.sd: spectral_diversity's output + the
    shift of esd_mode <n>, esd<n>_shift) and NOT applied: the images are the geometry-only alignment;
  - nesd_frames.txt lists each date frame with its super-master frame and burst offset (tmp_da);
  - the super-master's overlap images are written before its image, so mode 1's PRM / LED are the ones kept.
With IGML_ESD_SHIFTS=<file> (pass 2): each date frame's shift is read from the file (lines "stem shift", in lines)
  and added to a.grd as GMTSAR's constant ESD does; no overlap images, no spectral_diversity; a frame missing from the
  file stops the run.
Without either the script is unchanged.

usage: python patch_nesd.py <fork source tree>   (in place; refuses to run twice)
"""
import sys
from pathlib import Path

root = Path(sys.argv[1])
f = root / "gmtsar/csh/preproc_batch_tops_esd.csh"
s = f.read_text()
if "IGML NESD" in s:
    sys.exit(f"{f}: NESD already applied")


def sub(text, old, new):
    n = text.count(old)
    assert n == 1, f"expected 1 x {old[:80]!r}, found {n}"
    return text.replace(old, new)


s = sub(s, "          make_s1a_tops $file.xml $file.tiff $stem 1\n          make_s1a_tops $file.xml $file.tiff $stem 2\n",
        "          # IGML NESD (2026-10-31): a collecting run writes the overlap images before the image (mode 1's PRM /\n"
        "          # LED are then the ones kept, as in preproc_batch_tops.csh); a run with network shifts needs none\n"
        "          if ($?IGML_ESD_COLLECT) then\n"
        "            make_s1a_tops $file.xml $file.tiff $stem 2\n"
        "            make_s1a_tops $file.xml $file.tiff $stem 1\n"
        "          else if ($?IGML_ESD_SHIFTS) then\n"
        "            make_s1a_tops $file.xml $file.tiff $stem 1\n"
        "          else\n"
        "            make_s1a_tops $file.xml $file.tiff $stem 1\n"
        "            make_s1a_tops $file.xml $file.tiff $stem 2\n"
        "          endif\n")

s = sub(s, "          # generate the image with point-by-point shifts\n"
           "          make_s1a_tops $file.xml $file.tiff $stem 2 r.grd a.grd\n",
        "          # IGML NESD (2026-10-31): the network's shift for this frame instead of the estimate against the\n"
        "          # super-master\n"
        "          if ($?IGML_ESD_SHIFTS) then\n"
        "            set res_shift = `awk '$1 == \"'$stem'\" {print $2; exit}' $IGML_ESD_SHIFTS`\n"
        "            if (\"$res_shift\" == \"\") then\n"
        "              echo \"IGML NESD: no network shift for $stem in $IGML_ESD_SHIFTS\"\n"
        "              exit 1\n"
        "            endif\n"
        "            echo \"IGML NESD: updating azimuth shift of $stem with the network's constant ($res_shift)\"\n"
        "            gmt grdmath a.grd $res_shift ADD = tmp.grd\n"
        "            mv tmp.grd a.grd\n"
        "          else\n"
        "          if ($?IGML_ESD_COLLECT) echo \"$stem $tmp_master $tmp_da\" >> nesd_frames.txt\n"
        "          # generate the image with point-by-point shifts\n"
        "          make_s1a_tops $file.xml $file.tiff $stem 2 r.grd a.grd\n")

s = sub(s, "              set res_shift = `grep residual_shift tmp | awk '{print $3}'`\n            endif\n"
           "            echo \"Updating azimuth shift with a constant...(medain $res_shift)\"\n"
           "            gmt grdmath a.grd $res_shift ADD = tmp.grd\n"
           "            mv tmp.grd a.grd\n          endif\n",
        "              set res_shift = `grep residual_shift tmp | awk '{print $3}'`\n            endif\n"
        "            if ($?IGML_ESD_COLLECT) then\n"
        "              # IGML NESD: recorded for the network (its fallback), not applied\n"
        "              cp tmp $stem.sd\n"
        "              echo \"esd${esd}_shift = $res_shift\" >> $stem.sd\n"
        "              echo \"IGML NESD: $stem against $tmp_master $res_shift (recorded, not applied)\"\n"
        "            else\n"
        "            echo \"Updating azimuth shift with a constant...(medain $res_shift)\"\n"
        "            gmt grdmath a.grd $res_shift ADD = tmp.grd\n"
        "            mv tmp.grd a.grd\n"
        "            endif\n          endif\n"
        "          endif\n")

s = sub(s, "  if($mode == 2) rm tmp* topo.llt flt.grd atmp* rtmp* *SLCL *SLCH *BB\n",
        "  # IGML NESD (2026-10-31): a collecting run keeps the overlap images and boundary files (a run with network\n"
        "  # shifts makes none)\n"
        "  if ($mode == 2 && ($?IGML_ESD_COLLECT || $?IGML_ESD_SHIFTS)) then\n"
        "    rm tmp* topo.llt flt.grd atmp* rtmp*\n"
        "  else if ($mode == 2) then\n"
        "    rm tmp* topo.llt flt.grd atmp* rtmp* *SLCL *SLCH *BB\n"
        "  endif\n")
f.write_text(s)
print(f"NESD applied in {root}")
