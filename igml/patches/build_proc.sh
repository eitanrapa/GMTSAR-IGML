#!/bin/bash
# 2026-10-29: build our patched GMTSAR interferogram / merge programs (patch_proc.py) into bin_proc/, with the lab's
# own flags (/usr/local/GMTSAR/config.mk, gmtsar/Makefile, src_stitch/Makefile). First the UNPATCHED sources, which must
# reproduce /usr/local/GMTSAR/bin byte for byte (the flags are right), then the patched ones. filter.csh is copied.
set -e
cd /zdata/datasets/IGML/gmtsar_fix
G=/usr/local/GMTSAR
CF="-O2 -fcommon -fPIC -O2 -Wall -fPIC -fno-strict-aliasing -std=gnu99"
INC="-I/usr/local/include/gmt -I$G/gmtsar -I/usr/local/include"
LD="-s -Wl,-rpath,/usr/local/lib64"
LIBS="-L$G/gmtsar -lgmtsar -L/usr/local/lib64 -lgmt -lopenblas -lm"
mkdir -p build bin_proc
b() {  # out, sources...; gmtsar/Makefile link line: LDFLAGS obj GMTSAR GMT_LIB LALIBS LIBS -ltiff -lm
  out=$1; shift; objs=""
  for s in "$@"; do o=build/$(basename ${s%.c}).o; gcc $CF $INC -c $s -o $o; objs="$objs $o"; done
  gcc $LD $objs $LIBS -L/usr/local/lib -ltiff -lm -o $out
}
bs() {  # src_stitch/Makefile: CFLAGS LDFLAGS -o out obj CLIBS
  out=$1; s=$2; o=build/$(basename ${s%.c}).o; gcc $CF -I$G/gmtsar -I/usr/local/include/gmt -c $s -o $o
  gcc $CF $LD -o $out $o -L$G/gmtsar -lgmtsar -L/usr/local/lib64 -lgmt -lopenblas -lm
}
b build/conv.orig src/conv.orig.c
b build/phasefilt.orig src/phasefilt.orig.c
bs build/merge_swath.orig src/merge_swath.orig.c
for p in conv phasefilt merge_swath; do
  if cmp -s build/$p.orig $G/bin/$p; then echo "$p: our build of the lab source = the lab binary"; else echo "$p: our build of the lab source DIFFERS from $G/bin/$p"; fi
done
b bin_proc/conv src/conv.c src/conv2d.c
b bin_proc/phasefilt src/phasefilt.c
bs bin_proc/merge_swath src/merge_swath.c
cp src/filter.csh bin_proc/filter.csh; chmod +x bin_proc/filter.csh
ls -la bin_proc
