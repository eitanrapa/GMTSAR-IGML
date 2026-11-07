#!/bin/bash
# 11-01: T173's pairs after pass 2 (network spectral diversity): every interferogram made again with the fork --
# intf / merge / unwrap / store / combine, in windows of <= 60 dates; INTF_NCORES 12 interferograms per subswath at once
# (phasediff, the read-heavy step, behind its 16-slot gate). start_proc sets the proc1 pairs aside first. The lab's
# files are replaced only by the explicit replace stage, after review.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
if pgrep -u $U -f "^$PY process_track.py T173" > /dev/null; then echo "T173 already running"; exit 1; fi
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): pairs (intf,merge,unwrap,store,combine --batched 60), INTF_NCORES 12" >> $LOG
INTF_NCORES=12 setsid -f nohup nice -n 10 $PY process_track.py T173 --stages intf,merge,unwrap,store,combine --batched 60 >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date "+%F %T") launched T173 pairs" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "process_track.py T173"
