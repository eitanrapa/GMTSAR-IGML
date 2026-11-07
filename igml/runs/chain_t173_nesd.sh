#!/bin/bash
# 10-31: once T173's NESD pass 1 has ended cleanly, the network (process_track nesd: every pair through
# spectral_diversity, the solution, then the state turns to pass 2). Pass 2 is launched by hand after the solution is
# reviewed.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
while pgrep -u $U -f "^$PY process_track.py T173 --stages plan,stitch,rescue,align,qc$" > /dev/null; do sleep 120; done
if ! tail -3 $LOG | grep -q "stages \[.plan., .stitch., .rescue., .align., .qc.\] done"; then
    echo "$(date "+%F %T") T173 pass 1 did not end cleanly: nesd not started" >> /zdata/datasets/IGML/logs/queue.log
    exit 1
fi
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
echo "=== chain $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): nesd (pairs + solution), NESD_PROCS 40" >> $LOG
NESD_PROCS=40 nice -n 10 $PY process_track.py T173 --stages nesd >> $LOG 2>&1 < /dev/null
echo "$(date "+%F %T") T173 nesd ended (rc $?)" >> /zdata/datasets/IGML/logs/queue.log
