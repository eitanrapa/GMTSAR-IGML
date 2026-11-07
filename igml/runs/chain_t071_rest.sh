#!/bin/bash
# 11-05 (user: "T071, then T064, T166"): T071 after its network, without waiting for a review:
#   gate 1 -- the network as reviewed on T173 (each subswath: >= 98 % of the pairs used, <= 1 % fallback frames,
#             residual median <= 0.2 deg, p95 <= 1.5 deg; T173: 100 %, 0, 0.05-0.08, 0.36-0.52) and
#             <= 3 % of pass 1's dates set aside by qc (11-06; T173 0 of 430)                        -> pass 2 (align,qc)
#   gate 2 -- pass 2 ended cleanly and >= 97 % of the aligned dates pass qc (T173: 430 / 430)       -> the pairs
# A failed gate stops the chain (logs/queue.log says which). 11-06: gate 1's print had v[pairs] unquoted: it stopped the
# chain with a NameError at 00:43; both gates were then run on real state files before this version went live.
set -u
PY=/zdata/datasets/IGML/venv/bin/python
U=$(id -u)
LOG=/zdata/datasets/IGML/logs/T071_20261105.log
Q=/zdata/datasets/IGML/logs/queue.log
S=/zdata/datasets/IGML/InSAR_Processing/T071/state.json
while pgrep -u $U -f "^$PY process_track.py T071 " > /dev/null; do sleep 120; done
cd /zdata/datasets/IGML/code/scratchpad/s1_gapfill
G1=$($PY - "$S" <<"PYEOF" 2>&1
import json, sys
s = json.load(open(sys.argv[1]))
n = s.get("nesd", {})
assert n.get("pass") == 2, "the network is not solved"
for sw, v in sorted(n["subswaths"].items()):
    print(f"F{sw} pairs {v['pairs']} used {v['used']} frames {v['frames']} fallback {v['fallback']} resid median {v['resid_median_deg']} p95 {v['resid_p95_deg']}")
    assert v["used"] >= 0.98 * v["pairs"] and v["fallback"] <= 0.01 * v["frames"], f"F{sw} pairs / fallbacks"
    assert v["resid_median_deg"] <= 0.2 and v["resid_p95_deg"] <= 1.5, f"F{sw} residuals"
a, b = s.get("aligned_nesd1", []), s.get("bad_align_nesd1", {})
print(f"pass 1 set aside {len(b)} of {len(a)}: {sorted(b)}")
assert len(a) and len(b) <= 0.03 * len(a), "too many dates set aside in pass 1"
PYEOF
) || { echo "$(date "+%F %T") T071 network gate FAILED, pass 2 not started: $G1" >> $Q; exit 1; }
echo "$(date "+%F %T") T071 network gate passed: $G1" | tr "\n" ";" >> $Q; echo >> $Q
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): NESD pass 2 (align,qc), ALIGN_PAR 12" >> $LOG
ALIGN_PAR=12 nice -n 10 $PY process_track.py T071 --stages align,qc >> $LOG 2>&1 < /dev/null
tail -3 $LOG | grep -q "stages \[.align., .qc.\] done" || { echo "$(date "+%F %T") T071 pass 2 did not end cleanly: pairs not started" >> $Q; exit 1; }
G2=$($PY - "$S" <<"PYEOF" 2>&1
import json, sys
s = json.load(open(sys.argv[1]))
q, al = s.get("qc", {}), s.get("aligned", [])
ok = sum(1 for d in al if q.get(d, {}).get("ok"))
print(f"pass 2 qc ok {ok} of {len(al)}")
assert len(al) and ok >= 0.97 * len(al), "too many dates fail qc"
PYEOF
) || { echo "$(date "+%F %T") T071 pass-2 gate FAILED, pairs not started: $G2" >> $Q; exit 1; }
echo "$(date "+%F %T") T071 pass-2 gate passed: $G2" >> $Q
echo "=== launch $(date "+%F %T") code $(md5sum process_track.py | cut -c1-8): pairs (intf,merge,unwrap,store,combine --batched 60), INTF_NCORES 12" >> $LOG
INTF_NCORES=12 nice -n 10 $PY process_track.py T071 --stages intf,merge,unwrap,store,combine --batched 60 >> $LOG 2>&1 < /dev/null
echo "$(date "+%F %T") T071 pairs ended (rc $?)" >> $Q
