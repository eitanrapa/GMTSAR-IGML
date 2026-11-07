import json, glob, re
from pathlib import Path
for f in sorted(glob.glob("/zdata/datasets/IGML/InSAR_Processing/T*/state.json")):
    T = f.split("/")[-2]; s = json.load(open(f)); p = s.get("plan", {})
    logs = "".join(open(g, errors="replace").read() for g in glob.glob(f"/zdata/datasets/IGML/InSAR_Processing/{T}/logs/organize_*.log"))
    for d, r in sorted(p.items()):
        n = [Path(x).name for x in (r.get("safes") or []) if x]
        if len(n) < 2: continue
        dts = sorted({x[56:62] for x in n}); modes = sorted({x[4:6] for x in n})
        if len(dts) > 1 or modes != ["IW"]:
            m = re.search(r"SKIP %s, as it ([a-z ]+)" % d, logs)
            byd = {}
            for x in sorted(n, key=lambda x: x[17:32]): byd.setdefault((x[4:6], x[56:62]), []).append(x[17:32][9:] + "-" + x[33:48][9:])
            print(T, d, "lab" if r.get("lab_aligned") else "", "stitched" if d in s.get("stitched", {}) else "", (m.group(1) if m else "-"), byd)
