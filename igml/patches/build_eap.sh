#!/bin/bash
# 2026-10-29: our make_s1a_tops (patch_eap.py) into bin_eap/, with src_tops/Makefile's flags; the unpatched source first
# (must equal /usr/local/GMTSAR/bin/make_s1a_tops byte for byte)
set -e
cd /zdata/datasets/IGML/gmtsar_fix
G=/usr/local/GMTSAR; S=$G/preproc/S1A_preproc
CF="-O2 -fcommon -fPIC -O2 -Wall -fPIC -fno-strict-aliasing -std=gnu99"
INC="-I$S/include -I$G/gmtsar -I/usr/local/include -I/usr/local/include/gmt"
LD="-s -Wl,-rpath,/usr/local/lib64"
LIBS="-L$S/lib -lxmlC -L$G/gmtsar -lgmtsar -L/usr/local/lib64 -lgmt -lopenblas -L/usr/local/lib -ltiff -lm"
mkdir -p build bin_eap
for v in orig patched; do
  src=src/make_s1a_tops.c; out=bin_eap/make_s1a_tops; [ $v = orig ] && { src=src/make_s1a_tops.orig.c; out=build/make_s1a_tops.orig; }
  gcc $CF $INC -c $src -o build/make_s1a_tops.$v.o
  gcc $CF $LD -o $out build/make_s1a_tops.$v.o $LIBS
done
cmp -s build/make_s1a_tops.orig $G/bin/make_s1a_tops && echo "make_s1a_tops: our build of the lab source = the lab binary" || echo "make_s1a_tops: our build of the lab source DIFFERS"
ls -la bin_eap
