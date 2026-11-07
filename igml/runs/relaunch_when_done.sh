#!/bin/bash
# 2026-10-27: T064 / T071 run the pre-17:00 code (empty lab images, old super-master check). When each ends, relaunch
# it with the synced code unless it finished cleanly. Logs: /zdata/datasets/IGML/logs/<T>_20261027.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
P=/zdata/datasets/IGML/venv/bin/python; L=/zdata/datasets/IGML/logs
watch() {
  local T=$1 PID=$2
  while kill -0 $PID 2>/dev/null; do sleep 30; done
  if tail -3 $L/${T}_20261027.log | grep -q "stages .* done"; then echo "$(date '+%F %T') $T finished cleanly: not relaunched" >> $L/relaunch.log; return; fi
  echo "=== relaunch $(date '+%F %T') (empty lab images handled, T064 lab dates re-aligned, rescue)" >> $L/${T}_20261027.log
  setsid nohup nice -n 10 $P process_track.py $T --stages plan,stitch,rescue,align,qc >> $L/${T}_20261027.log 2>&1 < /dev/null &
  echo "$(date '+%F %T') relaunched $T" >> $L/relaunch.log
}
watch T064 2245033 &
watch T071 2245030 &
wait
