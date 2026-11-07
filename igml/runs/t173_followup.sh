#!/bin/bash
# 2026-10-27: when T173's 15:05 run (pre-17:00 code) ends, run every stage again on the synced code: it adds the
# year-long summer pairs (and anything left). Logs: /zdata/datasets/IGML/logs/T173_20261027.log, relaunch.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
L=/zdata/datasets/IGML/logs
while kill -0 2240192 2>/dev/null; do sleep 60; done
echo "=== relaunch $(date '+%F %T') (year-long June-August pairs, empty-image handling)" >> $L/T173_20261027.log
setsid nohup nice -n 10 /zdata/datasets/IGML/venv/bin/python process_track.py T173 --batched 60 >> $L/T173_20261027.log 2>&1 < /dev/null &
echo "$(date '+%F %T') relaunched T173 (all stages)" >> $L/relaunch.log
