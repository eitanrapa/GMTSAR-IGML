import json, glob, os, struct
from pathlib import Path
bad = []; n = 0; missing = []
import sys
for f in [f"/zdata/datasets/IGML/InSAR_Processing/{t}/state.json" for t in sys.argv[1:]]:
    T = f.split("/")[-2]; s = json.load(open(f)); p = s.get("plan", {})
    for d, r in sorted(p.items()):
        for x in (r.get("safes") or []):
            if not x or Path(x).name[4:6] != "IW": continue
            for sw in ("iw1", "iw2", "iw3"):
                t = [q for q in Path(x, "measurement").glob(f"*{sw}*vv*.tiff") if not q.name.startswith("._")]
                if not t: missing.append((T, d, Path(x).name[17:48], sw)); continue
                n += 1
                with open(t[0], "rb") as fh: h = fh.read(16)
                if h[:4] == b"II*\x00":
                    off = struct.unpack("<I", h[4:8])[0]
                elif h[:4] == b"II+\x00":
                    off = struct.unpack("<Q", h[8:16])[0]
                else:
                    bad.append((T, d, Path(x).name[17:48], sw, "bad magic")); continue
                sz = os.path.getsize(t[0])
                if off + 2 > sz:
                    bad.append((T, d, Path(x).name[17:48], sw, f"{sz/1e9:.2f}/{off/1e9:.2f} GB", "lab" if r.get("lab_aligned") else "", "stitched" if d in s.get("stitched", {}) else ""))
print("tiffs checked", n, "truncated", len(bad), "missing", len(missing))
for b in bad: print("  ", *b)
for m in missing[:20]: print("  missing", *m)
