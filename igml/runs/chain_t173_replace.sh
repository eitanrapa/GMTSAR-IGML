#!/bin/bash
# 11-05 (user: "Yes, replace"): once T173 three super-master pairs are stored, the lab T173 stack is replaced by ours
# (process_track replace: the lab pairs / aligned images / masks MOVED to lab_originals_20261029/T173, ours installed,
# every move in its manifest.jsonl) -- only if a fresh dry run finds every lab pair made (none archived without one).
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
Q=/zdata/datasets/IGML/logs/queue.log
while pgrep -u $U -f "^$PY process_track.py T173 --stages intf,merge,unwrap,store,combine --batched 60$" > /dev/null; do sleep 60; done
if ! tail -3 $LOG | grep -q "stages \[.intf., .merge., .unwrap., .store., .combine.\] done"; then
    echo "$(date "+%F %T") T173 pairs did not end cleanly: replace not run" >> $Q; exit 1
fi
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
DRY=$($PY process_track.py T173 --stages replace --dry-run 2>&1 | grep "replace (dry run)")
echo "$DRY" >> $LOG
case "$DRY" in
    *", 0 archived without one)"*) ;;
    *) echo "$(date "+%F %T") T173 replace NOT run, dry run not clean: $DRY" >> $Q; exit 1;;
esac
cp -p /zdata/datasets/IGML/InSAR_Processing/T173/state.json /zdata/datasets/IGML/InSAR_Processing/T173/state.json.pre_replace_$(date +%Y%m%d_%H%M%S)
echo "=== replace $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8)" >> $LOG
nice -n 10 $PY process_track.py T173 --stages replace >> $LOG 2>&1 < /dev/null
echo "$(date "+%F %T") T173 replace ended (rc $?)" >> $Q
