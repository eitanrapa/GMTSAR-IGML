#!/bin/bash
# Installed as /zdata/datasets/IGML/gmtsar_fix/bin_proc/phasediff (2026-10-29). A slot gate in front of GMTSAR's own
# phasediff, which is unchanged (this only waits, so every output is the same): at most IGML_PHASEDIFF_SLOTS (16) run at
# once on seasat. phasediff holds ~6 GB (F2: both SLCs), and at a batch's first pair all 18 interferogram processes of a
# track start with it -- up to seven tracks reprocess at once (250 GB, ZFS's cache ~155 GB). A slot is a flock on
# locks/phasediff.<i>, held by phasediff itself (exec keeps the descriptor) until it exits.
K=${IGML_PHASEDIFF_SLOTS:-16}
D=/zdata/datasets/IGML/gmtsar_fix/locks
mkdir -p "$D"
while :; do
  for i in $(seq 1 "$K"); do
    exec 9>>"$D/phasediff.$i"
    flock -n 9 && exec /usr/local/GMTSAR/bin/phasediff "$@"
    exec 9>&-
  done
  sleep $((RANDOM % 5 + 2))
done
