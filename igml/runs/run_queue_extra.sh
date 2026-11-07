#!/bin/bash
# 2026-10-29 (user: "Could you also queue the other tracks that were added"): the extra tracks the lab never processed
# (T100, T013, T108; process_track.py NEW: our own base stage, the California subswaths only) after run_queue_proc1.sh:
# waits until that queue (PID $MAIN_QUEUE, default 2036315; "none": no wait) has launched its last track and exited,
# then launches each
# of these when fewer than MAX process_track.py runs are going, like the main queue (all stages, --batched 60, nice 10)
# -- and never while another process_track.py run of the same track is going (a test run, a base being built: a second
# run would remove its half-made OUT/<T>/base.part).
# Log: /zdata/datasets/IGML/logs/queue.log; per track <T>_<date of launch>.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
L=/zdata/datasets/IGML/logs
MAX=${MAX:-6}
MAIN_QUEUE=${MAIN_QUEUE:-2036315}
PY=/zdata/datasets/IGML/venv/bin/python
n_running() { pgrep -fc "^$PY process_track.py "; }
busy() { pgrep -f "process_track.py $1( |$)" > /dev/null; }
launch() {
  f=$L/$1_$(date +%Y%m%d).log
  echo "=== queue launch $(date '+%F %T') (PROC proc1, NEW track: base stage first; all stages, --batched 60)" >> $f
  setsid nohup nice -n 10 $PY process_track.py $1 --batched 60 >> $f 2>&1 < /dev/null &
  echo "$(date '+%F %T') launched $1 (extra track) -> $f" >> $L/queue.log
}
queue="T100 T013 T108"
echo "$(date '+%F %T') extra-track queue started (MAX $MAX, after queue $MAIN_QUEUE): $queue" >> $L/queue.log
[ "$MAIN_QUEUE" = none ] || while kill -0 $MAIN_QUEUE 2>/dev/null; do sleep 300; done     # MAIN_QUEUE=none: start now
while [ -n "$queue" ]; do
  set -- $queue
  if busy $1; then
    sleep 300; continue
  fi
  if [ "$(n_running)" -lt $MAX ]; then
    launch $1; shift; queue="$*"; sleep 120; continue
  fi
  sleep 300
done
echo "$(date '+%F %T') extra-track queue empty" >> $L/queue.log
