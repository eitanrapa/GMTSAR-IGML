"""IGML fork (2026-11-02): filter.csh reuses each date's decimated amplitude across its pairs (IGML_AMP_CACHE).

intf_tops.csh -> filter.csh decimates BOTH dates' amplitudes for every pair (conv twice per date, two stages: amp1.grd
from the reference PRM's SLC, amp2.grd from the repeat's): 4.5 of a pair's ~8.5 min per subswath on T173, while every
date is in ~20 pairs.

Only amp1 is cached: the stage-2 filter (gauss_<wavelength>) is made per pair from the REFERENCE PRM (its velocity and
height), so the repeat's amp2 depends on its partner (filters differ by up to 0.12 % between reference dates), while
amp1 -- the reference's own amplitude with its own filter -- is the same in every pair it references (~23 % of a pair).
With IGML_AMP_CACHE=<dir> set, amp1 is looked up in <dir> under a key made of the SLC's name, inode, mtime and size and
every setting it depends on (az_lks, dec_rng, idec, jdec, filter1's name, filter2's md5, the IGML cell origin igx0 /
igy0); a hit is copied to amp1.grd, a miss is computed as before and stored (written beside the entry, then renamed:
parallel pairs never read a partial file). amp2 is always computed, and so is amp1 when the SLC or the filter cannot
be identified (stat / md5sum fail: a key without the SLC's identity could serve another version's amplitude). Without
IGML_AMP_CACHE nothing changes.

usage: python patch_ampcache.py <fork source tree>   (in place; refuses to run twice)
"""
import sys
from pathlib import Path

root = Path(sys.argv[1])
f = root / "gmtsar/csh/filter.csh"
s = f.read_text()
if "IGML_AMP_CACHE" in s:
    sys.exit(f"{f}: amplitude cache already applied")

old = ("  $st1 conv $az_lks $dec_rng $filter1 $1 amp1_tmp.grd=bf\n"
       "  $st2 conv $idec $jdec $filter2 amp1_tmp.grd=bf amp1.grd\n"
       "  rm amp1_tmp.grd\n"
       "  $st1 conv $az_lks $dec_rng $filter1 $2 amp2_tmp.grd=bf\n"
       "  $st2 conv $idec $jdec $filter2 amp2_tmp.grd=bf amp2.grd\n"
       "  rm amp2_tmp.grd\n")
assert s.count(old) == 1, "the amplitude block is not as expected"
new = '''  # IGML 2026-11-02: with IGML_AMP_CACHE, the reference's decimated amplitude (amp1) is computed once per date and
  # setup and copied for its other pairs (key: the SLC's name, inode, mtime, size and every setting the result depends
  # on, filter2 by its md5 -- filter2 comes from the reference PRM, so amp2 depends on its partner and is not cached)
  foreach k (1 2)
    set prm = $argv[$k]
    set amp = "amp${k}.grd"
    set key = ""
    if ($?IGML_AMP_CACHE && $k == 1) then
      set slc = `grep "^SLC_file" $prm | awk '{print $3}'`
      set sid = `stat -L -c "%i_%Y_%s" $slc`
      set fmd = `md5sum $filter2 | cut -c1-16`
      if ("$sid" != "" && "$fmd" != "") then
        set key = "$IGML_AMP_CACHE/$slc:t:r_${sid}_${az_lks}_${dec_rng}_${idec}_${jdec}_$filter1:t_${fmd}_${igx0}_${igy0}.grd"
      endif
    endif
    if ("$key" != "" && -e "$key") then
      cp "$key" $amp
      echo "IGML filter.csh: $amp from the amplitude cache ($key:t)"
    else
      $st1 conv $az_lks $dec_rng $filter1 $prm amp${k}_tmp.grd=bf
      $st2 conv $idec $jdec $filter2 amp${k}_tmp.grd=bf $amp
      rm amp${k}_tmp.grd
      if ("$key" != "") then
        mkdir -p $IGML_AMP_CACHE
        cp $amp "$key.$$" && mv "$key.$$" "$key"
      endif
    endif
  end
'''
s = s.replace(old, new)
f.write_text(s)
print(f"amplitude cache applied in {root}")
