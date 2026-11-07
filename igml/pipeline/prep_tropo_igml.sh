#!/bin/bash
# The lab's product_<T>/tropo/prep_tropo.sh for OUR zenith grids (2026-10-27): per <date>.<HHMM>.igml.grd in this folder
# (OUR product_<T>/tropo), proj_ll2ra.csh on the lab's merge/trans.dat (read in place), gmt surface onto the lab's
# tropo/unwrap.grd grid, / the lab's tropo/up.grd x 4 pi / 0.055165 (the lab's wavelength; Sentinel-1's is 0.0554658)
# -> <date>_<doy>_los_rad_igml.grd here. One temporary folder per date (proj_ll2ra writes fixed-name files), 6 at once.
# usage (seasat, in OUR product_<T>/tropo): prep_tropo_igml.sh <the lab's product_<T> folder>
LAB=$1
export LAB
export PATH=/usr/local/GMTSAR/bin:$PATH     # a non-interactive shell on seasat lacks it
one() {
  g=$1
  datestr=${g%%.*}
  doy=$(date -d "${datestr:0:4}-${datestr:4:2}-${datestr:6:2}" +%j)
  doy=$((10#$doy - 1))
  out=${datestr}_${doy}_los_rad_igml.grd
  [ -f "$out" ] && return
  t=$(mktemp -d tmp_${datestr}_XXXX)
  (cd "$t" && proj_ll2ra.csh "$LAB/merge/trans.dat" "../$g" junk1 > log 2>&1 &&
   gmt grd2xyz junk1 > junk2 && gmt surface junk2 -R"$LAB/tropo/unwrap.grd" -Gjunk3 >> log 2>&1 &&
   gmt grdmath junk3 "$LAB/tropo/up.grd" DIV 4 MUL pi MUL 0.055165 DIV = "../$out.part" >> log 2>&1 &&
   mv "../$out.part" "../$out") || echo "FAILED $g (see $t/log)"
  [ -f "$out" ] && rm -rf "$t"
}
export -f one
ls *.igml.grd 2>/dev/null | parallel -j 6 one {}
echo "$(ls *_los_rad_igml.grd 2>/dev/null | wc -l) LOS grids, $(ls *.igml.grd 2>/dev/null | wc -l) zenith grids"
