#!/bin/bash
# 11-05 (user: "Ok we can delete those as we go"): T071 pass-1 images (OUT/T071/superseded_nesd1, set aside when pass 2
# starts) are deleted once pass 2 has passed its qc gate (chain_t071_rest.sh writes "T071 pass-2 gate passed" to the
# queue log). Nothing else is touched. 11-06: relaunched behind chain_t071_resume.sh, which runs chain_t071_rest.sh.
set -u
Q=/zdata/datasets/IGML/logs/queue.log
until grep -q "T071 pass-2 gate passed" $Q; do
    if ! pgrep -u $(id -u) -f "[c]hain_t071_(resume|rest)\.sh" > /dev/null; then
        grep -q "T071 pass-2 gate passed" $Q || { echo "$(date "+%F %T") T071 cleanup: chain ended without the pass-2 gate -- pass-1 images kept" >> $Q; exit 1; }
    fi
    sleep 300
done
D=/zdata/datasets/IGML/InSAR_Processing/T071/superseded_nesd1
[ -d /zdata/datasets/IGML/InSAR_Processing/T071/superseded_nesd1 ] || { echo "$(date "+%F %T") T071 cleanup: $D not found -- nothing deleted" >> $Q; exit 1; }
N=$(find "$D" -type f | wc -l); S=$(du -sh "$D" | cut -f1)
rm -rf /zdata/datasets/IGML/InSAR_Processing/T071/superseded_nesd1 && echo "$(date "+%F %T") T071 cleanup: pass-1 images deleted ($N files, $S, $D)" >> $Q
