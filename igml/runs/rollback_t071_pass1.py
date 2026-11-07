"""T071 back to the end of pass 1 (2026-11-06). The old qc set aside 29 of 427 dates (28 false: open water, far neighbours;
process_track qc 11-06), and nesd() then solved the network and turned the state to pass 2 without them. Undoes
start_pass2 exactly: the pass-1 images back from OUT/T071/superseded_nesd1/raw/F<n>/ to product_T071/F<n>/raw/, the
record back from aligned_nesd1 / qc_nesd1 / bad_align_nesd1, "nesd" dropped. Every check before anything moves; state.json is
copied first. esd/ stays (the overlaps; pairs.jsonl spares re-measuring pairs). usage (seasat venv, nothing of T071 running)."""
import json
import os
import shutil
import subprocess
from pathlib import Path

T = "T071"
O = Path(f"/zdata/datasets/IGML/InSAR_Processing/{T}")
S = O / "state.json"
assert subprocess.run(["pgrep", "-u", str(os.getuid()), "-f", f"process_track.py {T} "], capture_output=True).returncode == 1, "T071 is running"
s = json.load(open(S))
assert s.get("nesd", {}).get("pass") == 2, "not in pass 2"
assert not s["aligned"] and not s["qc"] and not s["bad_align"], "pass 2 has started"
a, q, b = s["aligned_nesd1"], s["qc_nesd1"], s["bad_align_nesd1"]
assert len(a) == 427 and len(q) == 427 and len(b) == 29, (len(a), len(q), len(b))
plan = []
for n in (1, 2, 3):
    src, dst = O / f"superseded_nesd1/raw/F{n}", O / f"product_{T}/F{n}/raw"
    assert dst.is_dir() and not any(dst.iterdir()), f"{dst} not empty"
    fs = sorted(src.iterdir())
    assert len(fs) == 427 * 3, (n, len(fs))
    names = {f.name for f in fs}
    for d in a:
        for ext in ("SLC", "PRM", "LED"):
            assert f"S1_{d}_ALL_F{n}.{ext}" in names, (d, n, ext)
    plan += [(f, dst / f.name) for f in fs]
for n in (1, 2, 3):
    for d in a:
        assert (O / f"esd/F{n}/{d}/frames.json").exists(), (d, n, "no overlaps")
copy = O / "state.json.pre_rollback_20261106"
assert not copy.exists(), "already rolled back once"
shutil.copy2(S, copy)
for f, t in plan:
    os.replace(f, t)
for k in ("aligned", "qc", "bad_align"):
    s[k] = s.pop(f"{k}_nesd1")
s.pop("nesd")
tmp = S.with_suffix(".json.tmp_rollback")
tmp.write_text(json.dumps(s, indent=1))
os.replace(tmp, S)
print(f"{T}: {len(plan)} pass-1 files back in product_{T}/F<n>/raw; state back at the end of pass 1: {len(s['aligned'])} aligned, "
      f"{sum(1 for v in s['qc'].values() if v.get('ok'))} qc ok, {len(s['bad_align'])} set aside; copy {copy.name}")
