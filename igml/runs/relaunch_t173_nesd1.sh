#!/bin/bash
# 10-31: stop T173's NESD pass-1 run (its whole session: the run and its GMTSAR jobs) and launch it again with ALIGN_PAR=$1
# (default 10). Chunks already finished stay (state "aligned"); the running ones start again. Refuses unless the PID is
# ours, runs exactly the pass-1 command and leads its own session. state.json is copied first and checked after.
set -u
N=${1:-10}
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
ST=/zdata/datasets/IGML/InSAR_Processing/T173/state.json
P=$(pgrep -u $U -f "^$PY process_track.py T173 --stages plan,stitch,rescue,align,qc$")
[ -n "$P" ] || { echo "no pass-1 run found"; exit 1; }
[ "$(echo "$P" | wc -l)" = 1 ] || { echo "several runs: $P"; exit 1; }
sid=$(ps -o sid= -p $P | tr -d " ")
[ "$sid" = "$P" ] || { echo "$P: not a session leader (sid $sid)"; exit 1; }
cp -p $ST $ST.bak_$(date +%Y%m%d%H%M%S)
n=$(pgrep -u $U -s $sid | wc -l)
pkill -TERM -u $U -s $sid
for i in $(seq 1 30); do pgrep -u $U -s $sid > /dev/null || break; sleep 1; done
pkill -KILL -u $U -s $sid 2>/dev/null; sleep 1
if pgrep -u $U -s $sid > /dev/null; then echo "session $sid still has processes"; exit 1; fi
$PY -c "import json; s = json.load(open(\"$ST\")); print(\"state.json parses:\", len(s.get(\"aligned\", [])), \"aligned\")" || { echo "state.json DAMAGED: restore the .bak"; exit 1; }
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
echo "=== relaunch $(date "+%F %T") (stopped pid $P and its $((n - 1)) jobs; code $(md5sum process_track.py | cut -c1-8)): NESD pass 1, ALIGN_PAR $N" >> $LOG
ALIGN_PAR=$N setsid -f nohup nice -n 10 $PY process_track.py T173 --stages plan,stitch,rescue,align,qc >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date "+%F %T") relaunched T173 NESD pass 1 (ALIGN_PAR $N, was $P)" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "process_track.py T173"
