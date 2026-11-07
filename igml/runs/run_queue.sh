#!/bin/bash
# 2026-10-28: the gap-fill queue on seasat (every stage, --batched 60, nice 10), at most MAX process_track.py runs at
# once. First T071 and T064: when their 10-27 17:16 runs (plan..qc) end, every stage again on the synced code (qc v2,
# intf..combine). Then the six tracks not started yet, SoCal -> Central -> North, as slots free up. T173 relaunches on
# its own (t173_followup.sh). Log: /zdata/datasets/IGML/logs/queue.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
L=/zdata/datasets/IGML/logs
MAX=${MAX:-6}
PY=/zdata/datasets/IGML/venv/bin/python
n_running() { pgrep -fc "^$PY process_track.py "; }
launch() {
  f=$L/$1_20261027.log; [ -f $f ] || f=$L/$1_20261028.log
  echo "=== queue launch $(date '+%F %T') (all stages, --batched 60)" >> $f
  setsid nohup nice -n 10 $PY process_track.py $1 --batched 60 >> $f 2>&1 < /dev/null &
  echo "$(date '+%F %T') launched $1 -> $f" >> $L/queue.log
}
relaunch="T071:2381529 T064:2381527"
queue="T166 T144 T137 T115 T042 T035"
echo "$(date '+%F %T') queue started (MAX $MAX): relaunch $relaunch; then $queue" >> $L/queue.log
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
echo "$(date '+%F %T') queue empty" >> $L/queue.log
