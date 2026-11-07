#!/bin/bash
# While process_track.py runs, kill any GMTSAR ext_orb_s1a that has used > 10 min of CPU (normally seconds): an orbit
# GMTSAR cannot use makes it loop forever (10-24, 10-25 on T173). Logs to the track's log.
log=${1:-scratchpad/s1_gapfill/process_T173.log}
while pgrep -f "process_track.py" >/dev/null; do
  for p in $(pgrep -x ext_orb_s1a); do
    t=$(ps -o time= -p $p | awk -F'[:.]' '{print ($1*60+$2)}')
    [ "${t:-0}" -gt 600 ] && kill $p && echo "$(date '+%F %T') guard: killed ext_orb_s1a $p ($t s CPU) in $(lsof -p $p -a -d cwd -Fn 2>/dev/null | sed -n 's/^n//p')" >> "$log"
  done
  sleep 60
done
