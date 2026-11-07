#!/bin/bash
# 10-31: T173, network spectral diversity PASS 1 (user 10-31: "go ahead on decision 2. And we can begin"; "We do want
# T173 and all its outputs delivered first"): every date aligned again with the IGML fork (A1, antenna pattern, complete
# frames only), geometry-only, its burst overlaps kept (OUT/T173/esd); qc. Four chunks of 12 dates at once.
# Interferograms wait for the network (nesd) and pass 2.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
if pgrep -u $U -f "process_track.py T173" > /dev/null; then echo "T173 already running"; exit 1; fi
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): NESD pass 1 (plan,stitch,rescue,align,qc), ALIGN_PAR 4" >> $LOG
ALIGN_PAR=4 setsid -f nohup nice -n 10 $PY process_track.py T173 --stages plan,stitch,rescue,align,qc >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date "+%F %T") launched T173 NESD pass 1" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "process_track.py T173"
