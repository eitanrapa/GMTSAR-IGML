"""After a track's replacement (process_track --stages replace), is the lab's tree what it should be? (2026-11-05)
  - state "replaced": the pairs installed, none archived without one, the images swapped;
  - every lab pair folder holds OUR pair (igml_proc.json, PROC) and its original sits in the archive;
  - every swapped image: the lab's raw/ file is ours (same size as OUT's; md5 on a few), the lab's original archived;
  - the manifest: one "installed" line per pair.
Read only. usage (seasat, the venv): python replace_check.py T [N_MD5]"""
import hashlib
import importlib.util
import json
import random
import sys

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
T = sys.argv[1]
nmd5 = int(sys.argv[2]) if len(sys.argv) > 2 else 3
tr = pt.Track(T)
rep = tr.state.get("replaced", {})
arch = pt.LAB_ARCHIVE / T / f"product_{T}"
names = sorted(p.name for p in (tr.labp / "merge").iterdir() if p.name[:2] == "20" and p.is_dir())
bad_pair = [n for n in names if json.loads((tr.labp / "merge" / n / "igml_proc.json").read_text()).get("proc") != pt.PROC]
missing_orig = [n for n in names if not (arch / "merge" / n).is_dir()]
print(f"{T} replaced: {rep.get('time')} proc {rep.get('proc')}: {rep.get('pairs')} pairs, {len(rep.get('archived_only', []))} "
      f"archived without one, {len(rep.get('images', []))} images")
print(f"lab merge/: {len(names)} pair folders, {len(names) - len(bad_pair)} ours ({pt.PROC}); not ours {bad_pair[:3]}; "
      f"original missing from the archive {missing_orig[:3]} ({len(missing_orig)})")
size_bad, swapped = [], []
for d in rep.get("images", []):
    for n in tr.subs:
        f = f"S1_{d}_ALL_F{n}.SLC"
        lab, ours, orig = tr.labp / f"F{n}/raw/{f}", tr.outp / f"F{n}/raw/{f}", arch / f"F{n}/raw/{f}"
        if not ours.exists():
            continue
        swapped.append((lab, ours))
        if not lab.exists() or lab.stat().st_size != ours.stat().st_size or lab.stat().st_ino == ours.stat().st_ino:
            size_bad.append(str(lab))
        if not orig.exists():
            size_bad.append(f"no archived original: {orig}")
print(f"images: {len(swapped)} SLCs swapped; problems {len(size_bad)} {size_bad[:3]}")


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


random.seed(0)
for lab, ours in random.sample(swapped, min(nmd5, len(swapped))):
    print(f"  md5 {lab.name}: lab {md5(lab)[:12]} ours {md5(ours)[:12]}")
man = (pt.LAB_ARCHIVE / T / "manifest.jsonl").read_text().splitlines()
n_inst = sum('"installed"' in l for l in man)          # outside the f-string: seasat's python < 3.12
print(f"manifest: {len(man)} lines, {n_inst} installed")
ok = (not bad_pair and not missing_orig and not size_bad and rep.get("pairs") == len(names)
      and not rep.get("archived_only"))
print("REPLACEMENT OK" if ok else "REPLACEMENT HAS PROBLEMS")
