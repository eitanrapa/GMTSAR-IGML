import json, glob, datetime as dt
from pathlib import Path
def t(s): return dt.datetime.strptime(s, "%Y%m%dT%H%M%S")
tot = 0; odd = []
for f in sorted(glob.glob("/zdata/datasets/IGML/InSAR_Processing/T*/state.json")):
    T = f.split("/")[-2]
    p = json.load(open(f)).get("plan", {})
    for d, r in sorted(p.items()):
        s = [Path(x).name for x in (r.get("safes") or []) if x]
        s = sorted(s, key=lambda n: n[17:32])
        if len(s) < 2: continue
        tot += 1
        ov = [(t(a[33:48]) - t(b[17:32])).total_seconds() for a, b in zip(s, s[1:])]
        du = [(t(n[33:48]) - t(n[17:32])).total_seconds() for n in s]
        sats = {n[:3] for n in s}
        if max(ov) > 3.5 or max(du) > 29 or min(ov) < 0:
            odd.append((T, d, sorted(sats), ov, du))
print("multi-scene dates", tot)
for o in odd: print(*o)
