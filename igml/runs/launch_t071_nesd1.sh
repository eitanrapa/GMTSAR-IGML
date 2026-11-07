#!/bin/bash
# 11-05: T071, network spectral diversity PASS 1, as T173 10-31 (user 11-05: "T071, then T064, T166", after T173): every
# date aligned again with the IGML fork (A1, antenna pattern, complete frames only), geometry-only, its burst overlaps
# kept (OUT/T071/esd); qc. Twelve chunks at once. Interferograms wait for the network (nesd) and pass 2.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T071_20261105.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
if pgrep -u $U -f "process_track.py T071" > /dev/null; then echo "T071 already running"; exit 1; fi
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): NESD pass 1 (plan,stitch,rescue,align,qc), ALIGN_PAR 12" >> $LOG
ALIGN_PAR=12 setsid -f nohup nice -n 10 $PY process_track.py T071 --stages plan,stitch,rescue,align,qc >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date "+%F %T") launched T071 NESD pass 1" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "process_track.py T071"
