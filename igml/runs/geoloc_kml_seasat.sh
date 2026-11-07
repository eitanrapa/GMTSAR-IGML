#!/bin/bash
# 2026-11-06 (user: "Google Earth files where I can see the geolocation problem, in the merged and conv stages"): what the Google Earth
# files are drawn from, T173's super-master pair, stock GMTSAR (the lab's originals, archived) against the fixed conv / merge_swath (ours):
#   merged_<build>/: corr and display_amp through each build's own merge_swath (the stored corr.grd is masked), checked against the stored one;
#   every stage / build / subswath: GMTSAR's own proj_ra2ll.csh with the product's trans.dat (100-m pixels, its default for a 400-m filter),
#   which also leaves GMTSAR's lookup grids raln.grd / ralt.grd beside them.
# Nothing outside /zdata/datasets/IGML/work/geoloc_kml is written.
set -u
P=2016279_2016303
A=/zdata/datasets/IGML/lab_originals_20261029/T173/product_T173
L=/zdata/datasets/seasatHD3/InSAR_Processing/T173/product_T173
W=/zdata/datasets/IGML/work/geoloc_kml
export GMTSAR=/usr/local/GMTSAR PATH=/usr/local/GMTSAR/bin:/usr/local/bin:$PATH
mkdir -p $W
geo() {   # dir trans grd name
  cd "$1" || return 1
  proj_ra2ll.csh "$2" "$3" "$4"_ll.grd 400 > log_$4 2>&1
  echo "$(date +%T) $1 $4: $(gmt grdinfo -C "$4"_ll.grd 2>&1 | cut -f2-5,8-11)"
}
for b in stock fixed; do
  if [ $b = stock ]; then R=$A; BIN=/usr/local/GMTSAR/bin/merge_swath; else R=$L; BIN=/zdata/datasets/IGML/gmtsar_igml/install/bin/merge_swath; fi
  D=$W/merged_$b; mkdir -p $D; cd $D
  for g in corr display_amp; do
    for n in 1 2 3; do echo "$R/F$n/intf_all/$P/tmp.PRM:$R/F$n/intf_all/$P/$g.grd"; done > list_$g
    $BIN list_$g ${g}_ra.grd > log_merge_$g 2>&1
  done
  gmt grdmath corr_ra.grd $R/merge/$P/corr.grd SUB = diff_corr.grd
  echo "merged_$b: fresh corr minus stored: $(gmt grdinfo -C diff_corr.grd | cut -f6-7,10-11)"
  grep -hE "IGML|ovl12|Stitching" log_merge_corr | cut -c1-90
done
N=0
run() { "$@" & N=$((N+1)); if [ $N -ge 5 ]; then wait -n; fi; }
for b in stock fixed; do
  if [ $b = stock ]; then R=$A; else R=$L; fi
  ( geo $W/merged_$b $L/merge/trans.dat $W/merged_$b/corr_ra.grd corr; geo $W/merged_$b $L/merge/trans.dat $W/merged_$b/display_amp_ra.grd display_amp ) &
  for n in 1 2 3; do
    D=$W/conv_${b}_F$n; mkdir -p $D
    ( geo $D $L/F$n/topo/trans.dat $R/F$n/intf_all/$P/corr.grd corr; geo $D $L/F$n/topo/trans.dat $R/F$n/intf_all/$P/display_amp.grd display_amp ) &
  done
  wait
done
echo "$(date +%T) done"; du -sh $W
