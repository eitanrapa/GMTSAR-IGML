#!/bin/bash
# 2026-10-28 (GMTSAR audit): the runs launched before the audit's fixes (T071 03:03, T137 / T144 / T166 01:07-01:11) are
# relaunched, every stage, when each ends, so that they pick up datatake_rescue (T071's four S1C dates, T166 2016...),
# the raw-size and burst checks and the merge / store guards. At most MAX process_track runs at once, like run_queue.sh.
# T173 relaunches through t173_followup.sh, T064 was relaunched by hand 10:50. Log: /zdata/datasets/IGML/logs/queue.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
L=/zdata/datasets/IGML/logs
MAX=${MAX:-6}
PY=/zdata/datasets/IGML/venv/bin/python
n_running() { pgrep -fc "^$PY process_track.py "; }
pending=""
for T in ${@:-T071 T137 T144 T166}; do
  p=$(pgrep -f "^$PY process_track.py $T " | head -1)
  [ -n "$p" ] && pending="$pending $T:$p"
done
echo "$(date '+%F %T') followup_20261028 started: relaunch$pending when each ends (MAX $MAX)" >> $L/queue.log
while [ -n "$pending" ]; do
  left=""
  for item in $pending; do
    T=${item%%:*}; p=${item##*:}
    if kill -0 $p 2>/dev/null || [ "$(n_running)" -ge $MAX ]; then left="$left $item"; continue; fi
    f=$L/${T}_20261027.log; [ -f $f ] || f=$L/${T}_20261028.log
    echo "=== followup relaunch $(date '+%F %T') (all stages, --batched 60; the audit's fixes)" >> $f
    setsid nohup nice -n 10 $PY process_track.py $T --batched 60 >> $f 2>&1 < /dev/null &
    echo "$(date '+%F %T') followup relaunched $T -> $f" >> $L/queue.log
    sleep 120
  done
  pending=$(echo $left)
  [ -n "$pending" ] && sleep 300
done
echo "$(date '+%F %T') followup_20261028 done" >> $L/queue.log
