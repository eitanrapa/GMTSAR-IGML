"""10-30: one combined process_track.py change (user 10-29: "over 60% sounds good. Let's try all the ones you say are
possible"): partial frames down to 60 %, gap dates from their longest contiguous run, one copy per slice (re-downloads
included), the F6 / A4 frame check on every stitched frame (the lab's too) and the re-cut of flagged dates, qc's
coverage scaled by the share of the pin line a partial frame keeps. Exact replacements, each anchor asserted once."""
import shutil, sys, time
from pathlib import Path

P = Path(sys.argv[1] if len(sys.argv) > 1 else "/Users/erapaportbruck/Documents/IGML/scratchpad/s1_gapfill/process_track.py")
s = P.read_text()
EDITS = []


def edit(old, new):
    EDITS.append((old, new))


# ---------------------------------------------------------------- imports, constants
edit("from datetime import datetime, timedelta\n", "from datetime import datetime, timedelta, timezone\n")

edit('''# 10-27 (user: "We can try to recover them"): dates whose scenes miss one end of the lab's frame (GMTSAR: "SKIP ..., as
# it does not have enough scenes") are stitched with the outer pin slid inside the date's scenes. Only dates whose scenes
# cover >= RESCUE_MIN_COVER of the swath (ASF footprints): the qc stage then needs data in >= 90 % of a neighbour's
# windows (blocks at 2-98 % of the lines), which ~11 % missing at one end still meets.
RESCUE_MIN_COVER = 0.88
''', '''# 10-27 (user: "We can try to recover them"): dates whose scenes miss one end of the lab's frame (GMTSAR: "SKIP ..., as
# it does not have enough scenes") are stitched with the outer pin slid inside the date's scenes. Only dates whose scenes
# cover >= RESCUE_MIN_COVER of the swath (ASF footprints) and keep as much of the lab's pin line; qc then needs data in
# QC_COVER x that share of a neighbour's windows. 10-29 (user: "over 60% sounds good"): 0.88 -> 0.60, and dates with a
# gap keep their longest contiguous run of scenes (datatake_rescue).
RESCUE_MIN_COVER = 0.60
''')

# ---------------------------------------------------------------- plan: one copy per slice
edit('''        for stem in [] if tr.new else sorted(scenes.get(d, [])):
            hit = real(glob.glob(f"{ARCHIVE}/{stem[:3]}/CA/{d[:4]}/T*/{stem}_*.SAFE"))   # the archive, 10-27 layout
            paths.append(sorted(hit)[0] if hit else None)
''', '''        slices = []                         # 10-30: one slice per start time (within 3 s): T042 2025-05-01's last
        for stem in [] if tr.new else sorted(scenes.get(d, []), key=lambda x: x[17:32]):  # is listed twice, its
            t0 = datetime.strptime(stem[17:32], "%Y%m%dT%H%M%S")                          # stops 2 s apart
            if slices and (t0 - slices[-1][0]).total_seconds() <= 3:
                slices[-1][1].append(stem)
            else:
                slices.append((t0, [stem]))
        for _, stems in slices:
            paths.append(_pick_copy([h for st in stems for h in                                  # the archive,
                                     glob.glob(f"{ARCHIVE}/{st[:3]}/CA/{d[:4]}/T*/{st}_*.SAFE")    # 10-27 layout
                                     + glob.glob(f"{REDOWNLOAD}/{st}_*.SAFE")]))
''')
edit('''        la = tr.lab_aligned(d)            # 10-27: the lab aligns new dates before it merges them (its live stacks)
        plan_[d] = {"safes": paths, "ready": la or (bool(paths) and all(paths) and orb is not None),
                    "orbit": str(orb) if orb else None, "lab_aligned": la}
''', '''        la = tr.lab_aligned(d)            # 10-27: the lab aligns new dates before it merges them (its live stacks)
        plan_[d] = {"safes": paths, "ready": la or (bool(paths) and all(paths) and orb is not None),
                    "orbit": str(orb) if orb else None, "lab_aligned": la}
        if tr.state.get("plan", {}).get(d, {}).get("lab_aligned_bad"):     # 10-30: kept (retry_failed, recut)
            plan_[d]["lab_aligned_bad"] = True
''')

edit('''def _usable(tr, d, safes=None):''', '''REDOWNLOAD = "/palsar/IGML/redownload_*"      # complete copies of scenes the archive holds broken (10-28)


def _processed(safe):
    """When the scene was processed (manifest.safe's first processing start; "" if unreadable)."""
    try:
        m = re.search(r'<safe:processing [^>]*start="([^"]+)"', (Path(safe) / "manifest.safe").read_text(errors="replace"))
    except OSError:
        return ""
    return m.group(1) if m else ""


def _pick_copy(hits):
    """One copy of a slice (10-30), None if there is none: an intact one (_raw_bad), then the latest processing, then the
    archive's before a re-download of the same product. T071 2019-11-13 holds the 2019-11-14 and the 2020-01-14
    processing of every slice, their ANX times 24.9 ms apart: picked by name, the plan mixed the two (GMTSAR audit A4).
    T173 2021-01-01's 4AA9 is complete only in redownload_20261028."""
    hits = real(hits)
    if len(hits) <= 1:
        return hits[0] if hits else None
    return max(hits, key=lambda x: (not _raw_bad(x), _processed(x), x.startswith(f"{ARCHIVE}/S1")))


def _usable(tr, d, safes=None):''')

# ---------------------------------------------------------------- stitch: never the lab's frames of a re-cut date
edit('''    reused = []
    for d in list(todo):
        lab = _one_cut_per_frame(real(glob.glob(str(tr.labd / f"F*_F*/*_{d}T*.SAFE"))))''', '''    reused = []
    for d in list(todo):
        if d in tr.state.get("recut", {}):      # 10-30: cut from raw (recut), not from the lab's frames again
            continue
        lab = _one_cut_per_frame(real(glob.glob(str(tr.labd / f"F*_F*/*_{d}T*.SAFE"))))''')

# ---------------------------------------------------------------- rescue: longest run, the 60 % floor
edit('''def clipped_pins(tr, d, w):
    """The lab's pins (in time order) with the outer ones slid along the pin line until they lie RESCUE_MARGIN lines
    inside the date's first / last scene in all three subswaths; the interior pins (frame boundaries) are kept.
    -> (pins, share of the lab's pin line kept) or (None, why)."""
    pins = [l.split() for l in open(tr.pins) if l.strip()]
    if any(len(p) != 2 for p in pins):
        return None, "per-subswath pins: not handled"
    pins = [(float(x), float(y)) for x, y in pins]
    safes = _usable(tr, d)
    if not safes:''', '''def clipped_pins(tr, d, w, safes=None):
    """The lab's pins (in time order) with the outer ones slid along the pin line until they lie RESCUE_MARGIN lines
    inside the date's first / last scene (of `safes`, default its usable ones) in all three subswaths; the interior pins
    (frame boundaries) are kept. -> (pins, share of the lab's pin line kept) or (None, why)."""
    pins = [l.split() for l in open(tr.pins) if l.strip()]
    if any(len(p) != 2 for p in pins):
        return None, "per-subswath pins: not handled"
    pins = [(float(x), float(y)) for x, y in pins]
    safes = _usable(tr, d) if safes is None else safes
    if not safes:''')

edit('''def skipped_before(tr):''', '''def _runs(safes):
    """A date's scenes (time order) in contiguous runs, one version per slice (10-30): a new run where a scene starts
    > 1 s after the previous one stops or belongs to another datatake; a scene starting within 3 s of the previous one
    is another version of the same slice (T137 2015-05-31: 01:58:12 and 01:58:13), the later processing kept."""
    t = lambda x, i: datetime.strptime(Path(x).name[i:i + 15], "%Y%m%dT%H%M%S")
    runs = []
    for x in sorted(safes, key=lambda x: Path(x).name[17:32]):
        prev = runs[-1][-1] if runs else None
        if prev and abs((t(x, 17) - t(prev, 17)).total_seconds()) <= 3:
            if _processed(x) > _processed(prev):
                runs[-1][-1] = x
        elif prev and Path(x).name[56:62] == Path(prev).name[56:62] and (t(x, 17) - t(prev, 33)).total_seconds() <= 1:
            runs[-1].append(x)
        else:
            runs.append([x])
    return runs


def skipped_before(tr):''')

edit('''def rescue(tr, dates=None):
    """Stitch rescue_candidates with clipped_pins, one organize run per date (state "rescued": pins, share kept)."""''',
     '''def rescue(tr, dates=None):
    """Stitch rescue_candidates with clipped_pins, one organize run per date (state "rescued": pins, share kept), from
    the date's longest contiguous run of scenes (_runs; all of them unless it has a gap or two versions of a slice);
    >= RESCUE_MIN_COVER of the lab's pin line kept (10-30)."""''')
edit('''        pins, kept = clipped_pins(tr, d, w)
        if pins is None:
            res[d] = {"ok": False, "why": kept}
            log(f"{tr.T} rescue {d}: {kept}")
        else:
            (w / "pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\\n" for x, y in pins))
            (w / "safelist").write_text("\\n".join(_usable(tr, d)) + "\\n")''', '''        safes = max(_runs(_usable(tr, d)), key=len, default=[])
        pins, kept = clipped_pins(tr, d, w, safes) if safes else (None, "no usable raw scene")
        if pins is None:
            res[d] = {"ok": False, "why": kept}
            log(f"{tr.T} rescue {d}: {kept}")
        elif kept < RESCUE_MIN_COVER:
            res[d] = {"ok": False, "kept": kept, "why": f"{kept * 100:.0f} % of the lab's pin line inside its scenes"}
            log(f"{tr.T} rescue {d}: {res[d]['why']} (< {RESCUE_MIN_COVER * 100:.0f} %)")
        else:
            (w / "pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\\n" for x, y in pins))
            (w / "safelist").write_text("\\n".join(safes) + "\\n")''')

edit('''def datatake_rescue(tr):
    """Dates GMTSAR skipped for a gap ("stopped observation in the middle") whose scenes belong to more than one
    datatake: stitched again from the datatake with the most scenes, the lab's pins, one organize run per date (state
    "datatake"). organize_files_tops.csh tests gaps per date across datatakes (GMTSAR audit 10-28): T071's S1C
    2026-01-10 / 01-22 / 02-03 / 03-11 each have one scene of an earlier datatake 18 s before the frame's own three."""
    res = tr.state.setdefault("datatake", {})
    done, plan_ = tr.state.setdefault("stitched", {}), tr.state.get("plan", {})
    todo = [d for d, why in sorted(skipped_before(tr).items()) if why == "gap" and d not in res and d not in done
            and plan_.get(d, {}).get("ready") and not plan_[d].get("lab_aligned")]
    for d in todo:
        groups = defaultdict(list)
        for s in _usable(tr, d):
            groups[Path(s).name[56:62]].append(s)
        if len(groups) < 2:
            res[d] = {"ok": False, "why": "one datatake"}
            continue
        take, safes = max(groups.items(), key=lambda kv: len(kv[1]))
        w = tr.work / f"datatake_{d}"
        shutil.rmtree(w, ignore_errors=True)
        w.mkdir(parents=True)
        shutil.copy(tr.pins, w / "pins.ll")
        (w / "safelist").write_text("\\n".join(safes) + "\\n")
        orb = Path(plan_[d]["orbit"])
        os.symlink(orb, w / orb.name)
        run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log", check=False)''', '''def datatake_rescue(tr):
    """Dates GMTSAR skipped for a gap ("stopped observation in the middle"): stitched again from their longest
    contiguous run of scenes (_runs: one datatake, one version per slice) with the lab's pins clipped to it as in
    rescue (>= RESCUE_MIN_COVER of the pin line kept; all of it, i.e. the lab's pins, where the run spans them), one
    organize run per date (state "datatake", v 2). organize_files_tops.csh tests gaps per date across datatakes (GMTSAR
    audit 10-28): T071's S1C 2026-01-10 / 01-22 / 02-03 / 03-11 each have one scene of an earlier datatake 18 s before
    the frame's own three. 10-30 (user: "over 60% sounds good"): dates of ONE datatake with a gap (T166 2016-08-31,
    T115 2017-02-11 ...) keep their longest piece; the earlier version's failures are tried again once."""
    res = tr.state.setdefault("datatake", {})
    done, plan_ = tr.state.setdefault("stitched", {}), tr.state.get("plan", {})
    todo = [d for d, why in sorted(skipped_before(tr).items()) if why == "gap" and d not in done
            and (d not in res or (not res[d].get("ok") and res[d].get("v", 1) < 2))
            and plan_.get(d, {}).get("ready") and not plan_[d].get("lab_aligned")]
    for d in todo:
        usable = _usable(tr, d)
        safes = max(_runs(usable), key=len, default=[])
        w = tr.work / f"datatake_{d}"
        shutil.rmtree(w, ignore_errors=True)
        w.mkdir(parents=True)
        pins, kept = clipped_pins(tr, d, w, safes) if safes else (None, "no usable raw scene")
        if pins is None or kept < RESCUE_MIN_COVER:
            why = kept if pins is None else f"its longest run keeps {kept * 100:.0f} % of the lab's pin line"
            res[d] = {"ok": False, "v": 2, "why": why, "run": [Path(s).name for s in safes]}
            log(f"{tr.T} datatake {d}: {why}")
            tr.save()
            shutil.rmtree(w, ignore_errors=True)
            continue
        take = Path(safes[0]).name[56:62]
        (w / "pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\\n" for x, y in pins))
        (w / "safelist").write_text("\\n".join(safes) + "\\n")
        run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log", check=False)''')
edit('''        res[d] = {"ok": bool(frames), "datatake": take,
                  "left_out": sorted(Path(s).name for k, v in groups.items() if k != take for s in v)}
        log(f"{tr.T} datatake {d}: {len(frames)} frames from datatake {take} ({len(safes)} of "
            f"{sum(map(len, groups.values()))} scenes)")''', '''        res[d] = {"ok": bool(frames), "v": 2, "datatake": take, "kept": kept,
                  "pins": [[round(x, 6), round(y, 6)] for x, y in pins],
                  "left_out": sorted(Path(s).name for s in usable if s not in safes),
                  **({} if frames else {"why": "organize cut no frame from the run"})}
        log(f"{tr.T} datatake {d}: {len(frames)} frames from datatake {take} ({len(safes)} of {len(usable)} scenes, "
            f"{kept * 100:.0f} % of the lab's pin line)")''')

# ---------------------------------------------------------------- the F6 / A4 frame check
edit('''def _burst_check(tr, safe):''', '''_utc = lambda v: datetime.strptime(v.strip()[:26], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc).timestamp()
_BURST = re.compile(r"<burst>\\s*<azimuthTime>([^<]+)</azimuthTime>\\s*<azimuthAnxTime>([^<]+)</azimuthAnxTime>")
SCENE_TIMING = Path(__file__).resolve().parent / "size_audit/scene_timing.json"    # every raw scene of 10-29
_TIMING = {}


def _timing(x):
    """(linesPerBurst, ANX epoch = the first burst's azimuthTime - azimuthAnxTime [s]) of an annotation, read up to its
    first burst."""
    lpb, data = None, ""
    with open(x, errors="replace") as f:
        while chunk := f.read(1 << 18):
            data += chunk
            if lpb is None and (m := re.search(r"<linesPerBurst>(\\d+)<", data)):
                lpb = int(m.group(1))
            if lpb is not None and (m := _BURST.search(data)):
                return lpb, _utc(m.group(1)) - float(m.group(2))
    return lpb, None


def _scene_timing(safe):
    """[linesPerBurst iw1, iw2, iw3, ANX epoch] of a raw scene: from SCENE_TIMING (the 10-29 surveys of all 19,549),
    else from its annotations."""
    if not _TIMING:
        _TIMING.update(json.load(open(SCENE_TIMING))["scenes"] if SCENE_TIMING.exists() else {"": None})
    name = Path(safe).name
    if _TIMING.get(name) is None:
        row = [None] * 4
        for n in (1, 2, 3):
            x = real(Path(safe).glob(f"annotation/*iw{n}*vv*.xml"))
            if x:
                row[n - 1], e = _timing(x[0])
                row[3] = row[3] if n > 1 else e
        _TIMING[name] = row
    return _TIMING[name]


def _raw_mix(tr, d):
    """Do the date's raw scenes (every copy in the archive) differ in linesPerBurst in a processed subswath, or in ANX
    epoch by > 0.01 ms (GMTSAR audit F6 / A4)? Cached per date in state "raw_mix": only such dates' frames are read by
    _assembly_defects (theirs sit on the busy /zdata)."""
    cache = tr.state.setdefault("raw_mix", {})
    if d not in cache:
        rows = [_scene_timing(x) for x in real(glob.glob(f"{ARCHIVE}/S1?/CA/{d[:4]}/{tr.T}/S1?_IW_SLC__1S??_{d}T*.SAFE"))]
        ep = [r[3] for r in rows if r[3] is not None]
        cache[d] = any(len({r[n - 1] for r in rows if r[n - 1] is not None}) > 1 for n in tr.subs) or \\
            (len(ep) > 1 and max(ep) - min(ep) > 1e-5)
    return cache[d]


def _assembly_defects(tr, d, safes):
    """How make_s1a_tops would mis-read the date's frames ("" if it would not), per processed subswath: valid-sample
    lists not linesPerBurst long (GMTSAR audit F6: scenes of one datatake differ in linesPerBurst -> a heap overflow or
    each burst a line further off) or bursts its ANX-time placement puts on another line than their absolute times say
    (A4: scenes whose ANX epochs differ, e.g. two processings of the slices mixed). Only dates whose raw scenes mix
    (_raw_mix); cached per frame in state "assembly_check". Frames our assemble_tops cut (bin_proc, live since 10-29
    21:41) pass both; the lab's do not always (T173 2014-12-04, T144 2016-06-30, T071 2019-11-13)."""
    if not _raw_mix(tr, d):
        return ""
    cache, why = tr.state.setdefault("assembly_check", {}), []
    for s in safes:
        if s not in cache:
            w = []
            for n in tr.subs:
                x = real(Path(s).glob(f"annotation/*iw{n}*vv*.xml"))
                if not x:
                    continue
                t = Path(x[0]).read_text(errors="replace")
                lpb = int(re.search(r"<linesPerBurst>(\\d+)<", t).group(1))
                prf = 1 / float(re.search(r"<azimuthTimeInterval>([^<]+)<", t).group(1))
                b = _BURST.findall(t)
                ab, ax = [_utc(a) for a, _ in b], [float(v) for _, v in b]
                cnt = sorted({int(c) for c in re.findall(r'<firstValidSample count="(\\d+)"', t)})
                if cnt != [lpb]:
                    w.append(f"iw{n} valid-sample lists of {cnt} lines under linesPerBurst {lpb}")
                off = [k for k in range(len(b)) if int((ax[k] - ax[0]) * prf + 0.1) != round((ab[k] - ab[0]) * prf)]
                if off:
                    w.append(f"iw{n} bursts {off} placed off their times")
            cache[s] = f"{Path(s).parent.name}: {', '.join(w)}" if w else ""
        if cache[s]:
            why.append(cache[s])
    return "; ".join(why)


def _burst_check(tr, safe):''')

edit('''        why += [w for s in tr.state["stitched"][d] if s.startswith(str(tr.outd)) for w in [_burst_check(tr, s)] if w]
''', '''        why += [w for s in tr.state["stitched"][d] if s.startswith(str(tr.outd)) for w in [_burst_check(tr, s)] if w]
        why += [w for w in [_assembly_defects(tr, d, tr.state["stitched"][d])] if w]     # 10-30: the lab's frames too
''')

# ---------------------------------------------------------------- qc: coverage need of a partial frame
edit('''    res = tr.state.setdefault("qc", {})
    bad = tr.state.setdefault("bad_align", {})
    todo = [d for d in sorted(set(tr.state.get("aligned", [])) | set(lab_aligned_planned(tr)))
            if d not in res or (res[d].get("ok") and res[d].get("v", 1) < QC_VERSION)]''', '''    res = tr.state.setdefault("qc", {})
    bad = tr.state.setdefault("bad_align", {})
    part = _kept(tr)                         # 10-30: a partial frame needs QC_COVER x its share (its record: "need")
    todo = [d for d in sorted(set(tr.state.get("aligned", [])) | set(lab_aligned_planned(tr)))
            if d not in res or (res[d].get("ok") and res[d].get("v", 1) < QC_VERSION)
            or (not res[d].get("ok") and d in part and "need" not in res[d])]
    todo.sort(key=lambda d: d not in tr.state.get("recut", {}))   # re-cut dates first: neighbours recut() reopened
                                                                  # are then judged next to them (stable: time order)''')
edit('''        rec, ok = {"v": QC_VERSION}, bool(nb)
        for n in tr.subs:''', '''        rec, ok = {"v": QC_VERSION}, bool(nb)
        need = QC_COVER * part.get(d, 1)
        if d in part:
            rec["need"] = round(need, 3)
        for n in tr.subs:''')
edit('''                good = cov >= QC_COVER and coh >= flo + QC_MARGIN
''', '''                good = cov >= need and coh >= flo + QC_MARGIN
''')

edit('''def reject(tr):
    """Our stored pairs with a date qc has set aside since (10-28: the block test) -> OUT/<T>/rejected/<pair>, out of
    the state's "stored" list and the combined product (where its link points at ours)."""
    bad = {doy0(d) for d in tr.state.get("bad_align", {})}
    stored = tr.state.get("stored", [])
    gone = [p for p in stored if p[:7] in bad or p[8:] in bad]
    if not gone:
        return
    rj = tr.out / "rejected"
    rj.mkdir(exist_ok=True)
    for p in gone:
        src = tr.outp / "merge" / p
        link = tr.out / "combined" / tr.T / f"product_{tr.T}" / "merge" / p
        if link.is_symlink() and os.readlink(link) == str(src):
            link.unlink()
        if src.exists() and not (rj / p).exists():
            shutil.move(str(src), str(rj / p))
    tr.state["stored"] = [p for p in stored if p not in gone]
    tr.save()
    log(f"{tr.T} rejected {len(gone)} stored pairs of dates qc set aside -> {rj}: {gone[:6]}{' ...' if len(gone) > 6 else ''}")
''', '''def _kept(tr):
    """{date: share of the lab's pin line its frames keep} for the dates cut with clipped pins (rescue,
    datatake_rescue) short of the whole line."""
    return {d: v["kept"] for k in ("rescued", "datatake") for d, v in tr.state.get(k, {}).items()
            if v.get("ok") and isinstance(v.get("kept"), (int, float)) and v["kept"] < 1}


def _set_aside(src, dest_dir):
    """Move src into dest_dir, as <name>.<n> where that name is taken (nothing is deleted)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest, k = dest_dir / src.name, 1
    while dest.exists() or dest.is_symlink():
        dest, k = dest_dir / f"{src.name}.{k}", k + 1
    shutil.move(str(src), str(dest))
    return dest


def _withdraw(tr, dates, why):
    """Our stored pairs with any of `dates` -> OUT/<T>/rejected/<pair> (<pair>.<n> where a rejected one of that name is
    there already), out of the state's "stored" list and the combined product (where its link points at ours)."""
    bad = {doy0(d) for d in dates}
    stored = tr.state.get("stored", [])
    gone = [p for p in stored if p[:7] in bad or p[8:] in bad]
    if not gone:
        return []
    rj = tr.out / "rejected"
    for p in gone:
        src = tr.outp / "merge" / p
        link = tr.out / "combined" / tr.T / f"product_{tr.T}" / "merge" / p
        if link.is_symlink() and os.readlink(link) == str(src):
            link.unlink()
        if src.exists():
            _set_aside(src, rj)
    tr.state["stored"] = [p for p in stored if p not in gone]
    tr.save()
    log(f"{tr.T} rejected {len(gone)} stored pairs of dates {why} -> {rj}: {gone[:6]}{' ...' if len(gone) > 6 else ''}")
    return gone


def reject(tr):
    """Our stored pairs with a date qc has set aside since (10-28: the block test) -> OUT/<T>/rejected/<pair>, out of
    the state's "stored" list and the combined product (where its link points at ours)."""
    _withdraw(tr, tr.state.get("bad_align", {}), "qc set aside")


def recut(tr, asked=()):
    """Dates cut again from raw, once each (state "recut": {date: why, when, frames}; 10-30): planned dates whose frames
    _assembly_defects flags (the lab's included: T173 2014-12-04, the T144 four, T071 2019-11-13, T115 2015-03-06) or
    whose frames are all the lab's and fail any other of bad_stitches' checks (T115 2022-01-04: the lab's second frame
    holds iw1 only, its raw scenes are whole) -- our own cut wrong otherwise would come out the same -- and the dates
    asked for (--recut: T071 2026-08-15 with its standard slices, T173 2021-01-01 with its re-downloaded scene). Our
    frames -> OUT/<T>/rejected/frames/<frame>/, our aligned images -> rejected/aligned/, its stored pairs -> rejected/;
    the date leaves every stage's state, so stitch (never from the lab's frames again) or rescue / datatake_rescue cut
    it with our assemble_tops, and align / qc / its pairs follow in the same run. A flagged date we do not plan (a lab
    pair's partner) is only skipped by align (bad_stitches). A date that failed qc on coherence alone against a
    neighbour farther away than a re-cut date is judged again (qc takes the re-cut dates first)."""
    rc = tr.state.setdefault("recut", {})
    st, plan_ = tr.state.setdefault("stitched", {}), tr.state.get("plan", {})
    todo, bad = {}, bad_stitches(tr)
    for d in sorted(st):
        if d in rc or d not in bad:
            continue
        ad = _assembly_defects(tr, d, st[d])
        if not ad and any(s.startswith(str(tr.outd)) for s in st[d]):
            continue                                       # our own cut, wrong otherwise: a re-cut would repeat it
        if d in plan_:
            todo[d] = ad or f"the lab's frames: {bad[d]}"
        else:
            log(f"{tr.T} recut: {d}'s frames are cut wrong ({bad[d]}), not a date of ours: align skips it")
    for d in asked:
        if d in rc:
            log(f"{tr.T} recut: {d} was cut again already ({rc[d]['when']}), not again")
        elif d in plan_:
            todo[d] = "asked (--recut)"
        else:
            log(f"{tr.T} recut: {d} is not a date of ours (plan)")
    for d, why in todo.items():
        frames = st.pop(d, [])
        for s in frames:
            if s.startswith(str(tr.outd)) and Path(s).exists():
                _set_aside(Path(s), tr.out / "rejected" / "frames" / Path(s).parent.name)
            for k in ("burst_check", "assembly_check"):
                tr.state.get(k, {}).pop(s, None)
        for n in tr.subs:
            for ext in ("SLC", "PRM", "LED"):
                f = tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.{ext}"
                if f.exists():
                    _set_aside(f, tr.out / "rejected" / "aligned")
        for k in ("qc", "bad_align", "align_failed", "bad_stitch", "rescued", "datatake"):
            tr.state.get(k, {}).pop(d, None)
        for k in ("aligned", "retried", "stitch_fixed"):
            if d in tr.state.get(k, []):
                tr.state[k].remove(d)
        if plan_[d].get("lab_aligned"):
            plan_[d]["lab_aligned_bad"] = True                     # the lab's image came from those frames
        w = tr.work / "intf"                  # the batch scratch: its copy of the image set aside, its interferograms
        for n in tr.subs:
            for f in (w / f"F{n}/raw").glob(f"S1_{d}_ALL_F{n}.*"):
                f.unlink()
            for x in [*(w / f"F{n}/intf_all").glob(f"{doy0(d)}_*"), *(w / f"F{n}/intf_all").glob(f"*_{doy0(d)}")]:
                shutil.rmtree(x, ignore_errors=True)
        for x in [*(w / "merge").glob(f"{doy0(d)}_*"), *(w / "merge").glob(f"*_{doy0(d)}")]:
            shutil.rmtree(x, ignore_errors=True)
        rc[d] = {"why": why, "when": time.strftime("%F %T"), "frames": [Path(s).name for s in frames]}
        log(f"{tr.T} recut {d}: {why}")
    if todo:
        _withdraw(tr, todo, "cut again")
    q, ba = tr.state.get("qc", {}), tr.state.get("bad_align", {})
    for e, rec in sorted(q.items()) if todo else []:     # T144 2024-04-13: coherence 0.27 vs 2024-03-20, 24 days away,
        fs = [v for k, v in rec.items() if k.startswith("F") and v]          # its 12-day neighbour 2024-04-01 re-cut
        if rec.get("ok") or e in todo or not fs or any("shifted" in v for v in fs):
            continue                                  # passed, re-cut itself, or misregistered (its own frames)
        far = min(abs((day(e) - day(v["vs"])).days) for v in fs)
        near = [d for d in sorted(todo) if abs((day(e) - day(d)).days) < far]
        if near:
            q.pop(e)
            ba.pop(e, None)
            log(f"{tr.T} recut: {e} failed qc on coherence alone, against a date {far} days away: judged again "
                f"next to {near}")
    tr.save()
    return sorted(todo)
''')

# ---------------------------------------------------------------- run_batched: a relaunch finishes the batch it broke off
edit('''    todo = [p for p in all_pairs(tr) if f"{doy0(p[0])}_{doy0(p[1])}" not in tr.state.get("stored", [])]
    B = batches(todo, max_dates)
    mdn = ''', '''    todo = [p for p in all_pairs(tr) if f"{doy0(p[0])}_{doy0(p[1])}" not in tr.state.get("stored", [])]
    B = _resume_first(tr, todo, max_dates)
    mdn = ''')
edit('''def run_batched(tr, max_dates, ncores=6):''', '''def _resume_first(tr, todo, max_dates):
    """batches(), but the pairs whose dates' SLCs are all still copied in the work directory go first (10-30): a run
    killed mid-batch and launched again (new dates from recut / rescue shift every later boundary) finishes that batch
    instead of freeing its copies (~60 dates x 3 subswaths) only to copy them again batches later."""
    have = {f.name[3:11] for f in (tr.work / f"intf/F{tr.subs[0]}/raw").glob("S1_*_ALL_F*.SLC")}
    first = [p for p in todo if set(p) <= have]
    if not first or len({x for p in first for x in p}) > max_dates + 1:
        return batches(todo, max_dates)
    keep = {tuple(p) for p in first}
    return [sorted(first, key=lambda p: (p[1], p[0]))] + batches([p for p in todo if tuple(p) not in keep], max_dates)


def run_batched(tr, max_dates, ncores=6):''')

# ---------------------------------------------------------------- main
edit('''    ap.add_argument("--dry-run", action="store_true", help="replace: list what would move, move nothing")
''', '''    ap.add_argument("--dry-run", action="store_true", help="replace: list what would move, move nothing")
    ap.add_argument("--recut", default="", help="dates (comma-separated) to cut again from raw: their frames, alignment "
                    "and pairs set aside (OUT/<T>/rejected), once each")
''')
edit('''    ready = plan(tr, only) if "plan" in st else [d for d, v in tr.state.get("plan", {}).items() if v["ready"]]
    if only:
        ready = [d for d in ready if d in only]
''', '''    ready = plan(tr, only) if "plan" in st else [d for d, v in tr.state.get("plan", {}).items() if v["ready"]]
    if only:
        ready = [d for d in ready if d in only]
    if "stitch" in st:                          # 10-30: before `ours`, as a re-cut date of the lab's is ours to align
        recut(tr, [d for d in a.recut.split(",") if d])
''')

missing = [o[:80] for o, _ in EDITS if s.count(o) != 1]
if missing:
    sys.exit("anchors not found exactly once:\n  " + "\n  ".join(repr(m) for m in missing))
for o, n in EDITS:
    s = s.replace(o, n)
bak = P.with_name(P.name + time.strftime(".bak_%Y%m%d_%H%M%S"))
shutil.copy(P, bak)
P.write_text(s)
print(f"{len(EDITS)} edits applied to {P} (backup {bak.name})")
