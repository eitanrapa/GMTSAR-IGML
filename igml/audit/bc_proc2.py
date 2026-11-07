"""Dry run of PROC proc2 (10-31) on a deep copy of a track's live state, never saved, nothing moved (_set_aside
recorded): what start_realign would set aside and withdraw, and which dates the launch would stitch / align.
usage (seasat, the venv): python bc_proc2.py NEW_ROOT T173 [...]"""
import copy, importlib.util, json, os, sys
os.environ["IGML_CODE"] = sys.argv[1]
spec = importlib.util.spec_from_file_location("pt", os.path.join(sys.argv[1], "scratchpad/s1_gapfill/process_track.py"))
pt = importlib.util.module_from_spec(spec); spec.loader.exec_module(pt)
for T in sys.argv[2:]:
    live = json.load(open(f"/zdata/datasets/IGML/InSAR_Processing/{T}/state.json"))
    moved, msgs = [], []
    pt.log = lambda *a: msgs.append(" ".join(map(str, a)))
    pt._set_aside = lambda src, d: moved.append((str(src), str(d)))
    tr = pt.Track(T); tr.state = copy.deepcopy(live); tr.save = lambda: None
    print(f"== {T}: PROC {pt.PROC}, REALIGN_ALL {pt.REALIGN_ALL}, USE_FORK {pt.USE_FORK}, A3_ON {pt.A3_ON}, ESD {pt.ESD_MODE!r}")
    pt.start_realign(tr)
    imgs = [m for m in moved if "/raw/" in m[0]]
    print(f"   start_realign: {len(imgs)} aligned-image files set aside, {len(moved) - len(imgs)} partial frames withdrawn")
    print("  ", [m for m in msgs if "every date aligned again" in m][0][:300])
    ready = pt.plan(tr)
    pl = tr.state.get("plan", {})
    ours = [d for d in ready if not (pl.get(d, {}).get("lab_aligned") and not pl[d].get("lab_aligned_bad") and tr.lab_aligned(d))]
    realign = [d for d in pt.realign_dates(tr) if d not in ours]
    print(f"   plan: {len(ready)} ready (ours {len(ours)}); the lab's dates to align again: {len(realign)}; "
          f"lab-aligned images used: {sum(tr.lab_aligned(d) for d in set(tr.lab_dates()))}")
    print(f"   dates to stitch / align: {len(set(ours) | set(realign))}; master {tr.master}")
