import json, os, sys, time, statistics
from pathlib import Path
from collections import defaultdict
t0 = time.time(); med = defaultdict(list); rec = []
for T in sys.argv[1:]:
    s = json.load(open(f"/zdata/datasets/IGML/InSAR_Processing/{T}/state.json")); p = s.get("plan", {})
    for d, r in sorted(p.items()):
        for x in (r.get("safes") or []):
            if not x or Path(x).name[4:6] != "IW": continue
            m = Path(x, "measurement")
            try: names = [q for q in os.listdir(m) if "-vv-" in q and q.endswith(".tiff") and not q.startswith("._")]
            except OSError as e: print("ERR", T, d, x, e); continue
            sz = {q.split("-")[1]: os.stat(m / q).st_size for q in names}
            if len(sz) != 3: print("NOT3", T, d, Path(x).name[17:48], sorted(sz)); continue
            rec.append((T, d, Path(x).name, sz, r.get("lab_aligned"), d in s.get("stitched", {})))
            for k, v in sz.items(): med[(T, k)].append(v)
    print("#", T, f"{time.time()-t0:.0f}s", flush=True)
M = {k: statistics.median(v) for k, v in med.items()}
for T, d, nm, sz, lab, st in rec:
    r = {k: v / M[(T, k)] for k, v in sz.items()}
    lo = min(r, key=r.get); others = max(v for k, v in r.items() if k != lo)
    if r[lo] < 0.8 * others:
        print(T, d, nm[17:48], lo, " ".join(f"{k}={v:.2f}" for k, v in sorted(r.items())), "lab" if lab else "", "stitched" if st else "")
