#!/bin/bash
# 10-30: stop one process_track.py run (its whole session: the run and its GMTSAR jobs) and launch it again with the
# code on disk, in one go, so the proc1 queue (MAX 6 runs) never sees a free slot. Refuses unless PID is ours, is
# exactly "python process_track.py <T> --batched 60" and leads its own session.
# usage: relaunch_track.sh T PID [extra process_track.py args]
set -u
T=$1; P=$2; shift 2
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
cmd=$(tr '\0' ' ' < /proc/$P/cmdline 2>/dev/null)
[ "$(ps -o uid= -p $P | tr -d ' ')" = "$U" ] || { echo "$P: not ours"; exit 1; }
[ "$cmd" = "$PY process_track.py $T --batched 60 " ] || { echo "$P: command line is '$cmd'"; exit 1; }
sid=$(ps -o sid= -p $P | tr -d ' ')
[ "$sid" = "$P" ] || { echo "$P: not a session leader (sid $sid)"; exit 1; }
n=$(pgrep -u $U -s $sid | wc -l)
pkill -TERM -u $U -s $sid
for i in $(seq 1 30); do pgrep -u $U -s $sid > /dev/null || break; sleep 1; done
pkill -KILL -u $U -s $sid 2>/dev/null; sleep 1
if pgrep -u $U -s $sid > /dev/null; then echo "session $sid still has processes"; exit 1; fi
LOG=/zdata/datasets/IGML/logs/${T}_20261029.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
echo "=== relaunch $(date '+%F %T') (stopped pid $P and its $((n - 1)) jobs; code $(md5sum process_track.py | cut -c1-8): 10-30 partial frames / re-cut) $*" >> $LOG
setsid -f nohup nice -n 10 $PY process_track.py $T --batched 60 "$@" >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date '+%F %T') relaunched $T (was $P) $*" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "^$PY process_track.py $T "
