#!/bin/bash
# 11-01: T173, network spectral diversity PASS 2: every date of the network aligned again with its shift
# (IGML_ESD_SHIFTS = OUT/T173/esd/F<n>/shifts.txt, the fork's preproc_batch_tops_esd.csh, A1 on), then qc (its records
# carry the burst steps; pass 1's are kept as qc_nesd1). Twelve chunks at once. The pairs (intf / merge / unwrap /
# store) after the steps are reviewed.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T173_20261031.log
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
if pgrep -u $U -f "process_track.py T173" > /dev/null; then echo "T173 already running"; exit 1; fi
$PY -c "import json; s = json.load(open(\"/zdata/datasets/IGML/InSAR_Processing/T173/state.json\")); n = s.get(\"nesd\", {}); assert n.get(\"pass\") == 2, \"the network is not solved\"; print(\"network solved\", n.get(\"solved\"), len(n[\"dates\"]), \"dates\")" || exit 1
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): NESD pass 2 (align,qc), ALIGN_PAR 12" >> $LOG
ALIGN_PAR=12 setsid -f nohup nice -n 10 $PY process_track.py T173 --stages align,qc >> $LOG 2>&1 < /dev/null
sleep 3
echo "$(date "+%F %T") launched T173 NESD pass 2" >> /zdata/datasets/IGML/logs/queue.log
pgrep -u $U -fa "process_track.py T173"
