#!/bin/bash
# 11-06: T071 after the qc fix (AUDIT 11-06), state rolled back to the end of pass 1 (rollback_t071_pass1.py):
#   1. pass 1 closed: 2015-11-22 cut again from raw (its 4 scenes mix linesPerBurst; the lab's frame holds other content from
#      about line 19,000), aligned; qc judges again every date the old qc set aside;
#   2. the network again with those dates (pairs.jsonl keeps the 1,586 pairs measured per subswath);
#   3. chain_t071_rest.sh: gates, pass 2, gate, pairs.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
LOG=/zdata/datasets/IGML/logs/T071_20261105.log
Q=/zdata/datasets/IGML/logs/queue.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
echo "=== chain $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): pass 1 closed (--recut 20151122, qc again), ALIGN_PAR 12" >> $LOG
ALIGN_PAR=12 nice -n 10 $PY process_track.py T071 --stages plan,stitch,rescue,align,qc --recut 20151122 >> $LOG 2>&1 < /dev/null
RC=$?
echo "$(date "+%F %T") T071 pass 1 closed (rc $RC): $(grep "T071 qc: " $LOG | tail -1 | cut -c21-260)" >> $Q
tail -3 $LOG | grep -q "stages \[.plan., .stitch., .rescue., .align., .qc.\] done" || { echo "$(date "+%F %T") T071 pass 1 did not end cleanly: nesd not started" >> $Q; exit 1; }
echo "=== chain $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): nesd (pairs + solution), NESD_PROCS 40" >> $LOG
NESD_PROCS=40 nice -n 10 $PY process_track.py T071 --stages nesd >> $LOG 2>&1 < /dev/null
echo "$(date "+%F %T") T071 nesd ended (rc $?)" >> $Q
/zdata/datasets/IGML/work/chain_t071_rest.sh
