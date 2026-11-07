#!/bin/bash
# 2026-10-29 (user: "Ok let's do it" -- reprocess everything with the patched GMTSAR): T137 and T064 finish their
# align stage (unaffected by the patches); each is stopped the moment its intf stage starts ("pairs to make").
L=/zdata/datasets/IGML/logs
declare -A base
for T in T137 T064; do f=$(ls -t $L/${T}_2026102*.log | head -1); base[$T]=$(grep -c "pairs to make" $f); done
left="T137 T064"
while [ -n "$left" ]; do
  next=""
  for T in $left; do
    f=$(ls -t $L/${T}_2026102*.log | head -1); p=$(pgrep -f "^/zdata/datasets/IGML/venv/bin/python process_track.py $T ")
    if [ -z "$p" ]; then echo "$(date "+%F %T") $T: no run left" >> $L/queue.log; continue; fi
    if [ "$(grep -c "pairs to make" $f)" -gt "${base[$T]}" ]; then
      g=$(ps -o pgid= -p $p | tr -d " "); kill -TERM -- -$g
      echo "=== $(date "+%F %T") STOPPED at the start of intf (align done; reprocessing with the patched GMTSAR)" >> $f
      echo "$(date "+%F %T") stopped $T before intf" >> $L/queue.log
    else next="$next $T"; fi
  done
  left=$(echo $next); [ -n "$left" ] && sleep 15
done
