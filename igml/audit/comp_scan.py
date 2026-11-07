import json, glob, re
from pathlib import Path
rc = re.compile(rb"productComposition>([^<]*)<"); ri = re.compile(rb'Sentinel-1 IPF" version="([^"]*)"')
tot = 0; odd = []; ipfmix = []
for f in sorted(glob.glob("/zdata/datasets/IGML/InSAR_Processing/T*/state.json")):
    T = f.split("/")[-2]; p = json.load(open(f)).get("plan", {})
    for d, r in sorted(p.items()):
        s = sorted([x for x in (r.get("safes") or []) if x], key=lambda x: Path(x).name[17:32])
        if len(s) < 2: continue
        comps, ipfs = [], set()
        for x in s:
            b = b""
            with open(Path(x, "manifest.safe"), "rb") as fh:
                while len(b) < 400000:
                    c = fh.read(16384)
                    if not c: break
                    b += c
                    if rc.search(b) and ri.search(b): break
            m = rc.search(b); i = ri.search(b)
            comps.append(m.group(1).decode() if m else "?"); ipfs.add(i.group(1).decode() if i else "?")
            tot += 1
        if any(c != "Slice" for c in comps): odd.append((T, d, comps))
        if len(ipfs) > 1: ipfmix.append((T, d, sorted(ipfs)))
print("scenes read", tot)
print("dates with non-Slice products:", len(odd)); [print("  ", *o) for o in odd]
print("dates mixing IPF versions:", len(ipfmix)); [print("  ", *o) for o in ipfmix]
