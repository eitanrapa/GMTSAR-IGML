#!/bin/bash
# 2026-10-29 (user: "Ok let's do it" -- fix all the processing and the processed data): every track made again with
# the patched GMTSAR (process_track.py PROC proc1, every pair, the lab's too), all stages, --batched 60, nice 10, at
# most MAX process_track.py runs at once (the old T137 / T064 align runs count). T137 and T064 relaunch the moment
# stop_before_intf.sh has stopped their old runs (align done). T173's batch 1 is the validated test run (10-29 03:52).
# The lab's files are NOT touched here: replace is an explicit stage, per track, after a dry run.
# Log: /zdata/datasets/IGML/logs/queue.log; per track <T>_20261029.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
L=/zdata/datasets/IGML/logs
MAX=${MAX:-6}
PY=/zdata/datasets/IGML/venv/bin/python
n_running() { pgrep -fc "^$PY process_track.py "; }
launch() {
  f=$L/$1_20261029.log
  echo "=== queue launch $(date '+%F %T') (PROC proc1: all stages, --batched 60)" >> $f
  setsid nohup nice -n 10 $PY process_track.py $1 --batched 60 >> $f 2>&1 < /dev/null &
  echo "$(date '+%F %T') launched $1 (proc1) -> $f" >> $L/queue.log
}
relaunch="T137:2934692 T064:3986456"
queue="T173 T071 T166 T144 T115 T042 T035"
echo "$(date '+%F %T') proc1 queue started (MAX $MAX): relaunch $relaunch; then $queue" >> $L/queue.log
while [ -n "$relaunch$queue" ]; do
  left=""
  for item in $relaunch; do
    if kill -0 ${item##*:} 2>/dev/null; then left="$left $item"; else launch ${item%%:*}; fi
  done
  relaunch=$(echo $left)
  set -- $queue
  if [ $# -gt 0 ] && [ "$(n_running)" -lt $MAX ]; then
    launch $1; shift; queue="$*"; sleep 120; continue
  fi
  sleep 300
done
echo "$(date '+%F %T') proc1 queue empty" >> $L/queue.log
