"""Process the dates the lab's stack lacks, the lab's way (2026-10-16; user: "download and process the gaps ... plus make
the longer pairs up to 84 days"; the chain validated against the lab: val_t173.sh, stitching test).

RUNS ON SEASAT since 2026-10-27 (user: "only raw data on palsar", then "Seasat"): the lab's LIVE stacks are read in
/zdata/datasets/seasatHD<n>/InSAR_Processing/<T> (never written), ours go to /zdata/datasets/IGML/InSAR_Processing/<T>
laid out like the lab's (data_<T>/<frames>/ stitched, product_<T>/F<n>/raw aligned, product_<T>/merge/<pair>), scratch
in /zdata/datasets/IGML/work/<T>; code + the inventory files mirrored to /zdata/datasets/IGML/code (sync_seasat.sh),
python = /zdata/datasets/IGML/venv (numpy, for qc). Orbits are only READ there (the lab's data_<T>, /palsar/IGML/orbits):
fetch missing ones on the Mac (prefetch_orbits.py; the Earthdata login stays on the Mac). The Mac version:
process_track_mac_20261027.py. Dates the lab has ALIGNED but not merged (its stacks were being extended, Sept-Oct 2026)
use the lab's aligned SLCs: qc'd and paired, not re-aligned; combine prefers the lab's pair where both exist.

Per track, stages (each resumable; state in OUT/state.json):
  plan     dates = ASF acquisitions (>= half the swath) the lab stack lacks, whose every scene is on disk
           (/palsar/S1<A-D>/CA/<year>/T<track>/, the lab's and ours), with an orbit (the lab's data_<T>, or
           /palsar/IGML/orbits)
  stitch   organize_files_tops.csh with the lab's pins (orbit download stubbed) -> OUT/data_<T>/F*_F*/; since 10-28 (the
           GMTSAR audit, notes/gapfill_plan.md) raw scenes whose VV tiffs are not their manifest's size are left out
           (state "raw_bad"), frames we cut must be the super-master's (bad_stitches; the T064 pins trap) and their
           subswaths' burst spans iw1's (_burst_check); rescue also takes dates GMTSAR skipped for a gap to ANOTHER
           datatake (datatake_rescue: the datatake with the most scenes)
  align    preproc_batch_tops.csh modes 1, 2 per subswath against the lab's super-master (re-made each batch; checked
           bit-identical to the lab's) -> OUT/product_<T>/F<n>/raw/S1_<date>_ALL_F<n>.{PRM,LED,SLC}; since 10-28 with
           our patched stitch_tops for every date but the super-master (STITCH_FIX: GMTSAR applies the last frame's
           line offset twice)
  qc       every aligned date against its nearest aligned neighbours (the lab's, or ours that passed): image coverage
           and coherence over 12 blocks per subswath; a failing date is set aside (state "bad_align"), not paired
           (10-27: the pilot's 2014-12-04 had no data in a third of iw1 and nothing had flagged it); since 10-28 also
           no block whose coherence peaks at a 1-2 line shift (the lab's 2014-12-04 image: last slice ~1 line off);
           stored pairs of a date that fails later are moved to OUT/<T>/rejected (reject); a date set aside is aligned
           once more in the same run (retry_failed: ours misregistered before the fix, the lab's aligned-only images)
  intf     every pair <= MAXDT days with at least one new date, + the lab's super-master pair (merge needs it first);
           SLCs copied to local disk, intf_tops_parallel.csh per subswath (the lab's batch_tops.config, topo_ra.grd)
  merge    merge_batch.csh (the lab's dem / trans.dat / masks / config, copied)
  unwrap   the track's variant (unwrap_intf.csh): T064 / T071 snaphu_interp.csh 0.0001 60, T035 / T137
           snaphu_ocean_interp.csh 0.0001 0 (MATLAB grid_nearest), others snaphu_interp.csh 0.0001 0
  store    merge/<pair>/{unwrap,corr,conncomp,phasefilt,mask}.grd + supermaster.PRM -> OUT/product_<T>/merge/<pair>/
  combine  OUT/combined/<T>/product_<T>/: merge/ = links to every complete lab pair + ours, topo/ and the rest linked
           to the lab's -- the root igml.real converts
The lab's volumes are only READ: files GMTSAR could rewrite are copied, read-only inputs symlinked.

usage (seasat): /zdata/datasets/IGML/venv/bin/python /zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py
       <TRACK> [--stages plan,stitch,...] [--dates d1,d2] [--batch N]
       [--batched N_DATES]   (intf..store in time windows: local disk holds at most N_DATES dates' SLCs, ~8 GB each)
"""
import argparse
import csv
import glob
import json
import os
import re
import shutil
import signal
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

SEASAT = Path("/zdata/datasets/seasatHD1").is_dir()
REPO = Path(os.environ.get("IGML_CODE", "/zdata/datasets/IGML/code" if SEASAT else "/Users/erapaportbruck/Documents/IGML"))
LAB_ROOT = Path("/zdata/datasets" if SEASAT else "/Volumes/seasat2")    # the lab's live volumes (read only)
OUT_ROOT = Path("/zdata/datasets/IGML/InSAR_Processing")
WORK_ROOT = Path("/zdata/datasets/IGML/work")
VOL = {"T064": 1, "T071": 1, "T137": 2, "T144": 2, "T166": 3, "T173": 3, "T035": 4, "T042": 4, "T115": 4}
# 10-29 (user: "queue the other tracks that were added"): tracks the lab never processed, only the subswaths with
# California in them (extra_tracks/subswaths.json, 10-27). What the lab's stack gives the nine tracks (pins, the
# super-master, DEM, topography) is ours, made by the "base" stage under OUT/<T>/base/ in the lab's layout.
NEW = {"T013": (1, 2), "T100": (3,), "T108": (3,)}
EXTRA = "scratchpad/s1_inventory_20261015/extra_tracks/subswaths.json"
DEM_SRC = Path("/zdata/datasets/IGML/dem")        # <T>/dem.grd: SRTM 3" + EGM96 (make_dem.csh on the Mac; seasat's GMT
                                                  # cannot read the remote earth_relief tiles)
CONFIG_TEMPLATE = LAB_ROOT / "seasatHD3/InSAR_Processing/T173/product_T173/F1/batch_tops.config"   # the lab's recipe
UNWRAP = {"T064": "snaphu_interp.csh 0.0001 60", "T071": "snaphu_interp.csh 0.0001 60",   # the LAB's recipe: for
          "T035": "snaphu_ocean_interp.csh 0.0001 0", "T137": "snaphu_ocean_interp.csh 0.0001 0"}  # the parity tests
# OUR pairs (10-26; user: "let's do smooth except for coseismic pairs"): SNAPHU's SMOOTH mode, and DEFO (defomax 60, the
# lab's T064 / T071 value) only for a pair spanning an earthquake of coseismic_events.json (igml.catalog's rule: >= 50 mm
# in the track's swath). The lab ran DEFO on every T064 / T071 pair; on identical inputs DEFO added whole-cycle steps
# (T173 median 0.24 -> 0.29 % of pixels, 6-10 % on pairs with a steep feature).
UNWRAP_SCRIPT = {"T035": "snaphu_ocean_interp.csh", "T137": "snaphu_ocean_interp.csh"}
DEFOMAX = 60
COSEISMIC = json.load(open(Path(__file__).resolve().parent / "coseismic_events.json"))["tracks"]
TOOLS = REPO / "scratchpad/s1_gapfill/gmtsar_tools"
# 10-28: GMTSAR's stitch_tops applies the last frame's cross-correlation line offset twice (n_start -= line_offset_total,
# then -= line_offset; upstream master too): whenever the overlap shows a shift, the date's last frame lands that many
# lines off (T173 2026-01-29: blocks 6-11 one line off in all three subswaths, coherence peaking at -1). Seasat holds
# our patched copy (gmtsar_fix/bin) and a wrapper (gmtsar_fix/wrap/stitch_tops) that keeps the lab's binary for the
# super-master, so it stays the lab's; the lab's GMTSAR is not touched. Off where the folder is missing.
STITCH_FIX = Path("/zdata/datasets/IGML/gmtsar_fix/wrap")
# 10-29 (user: "Ok let's do it": fix the processing itself and replace the lab's stacks, originals archived): our patched
# GMTSAR (gmtsar_fix/patch_proc.py, patch_eap.py; the lab's /usr/local/GMTSAR untouched; our builds of the lab's own
# sources are byte-identical to its binaries). bin_proc, for interferograms and merges: conv centres each 32 x 6 cell on
# its block (I2, +16 samples / +3 lines) and starts each subswath on the merged grid's cells (IGML_SUB_X0 / _Y0, from
# the super-master's subswath PRMs), conv2d's bounds (I3), phasefilt's zero edges (I1, now NaN), merge_swath's
# placement exact (M6) and without M1's extra F2 column. bin_eap, for alignment: make_s1a_tops whose compute_eap removes
# the elevation antenna pattern of IPF < 002.43 scenes (A2; its swath lookup, weights, bounds fixed) when the AUX_CAL and
# the manifest are appended to the annotation (align does that). Validated on T071 at the Rosamond reflectors
# (audit/fixtest_check.py: declared - true -0.01 samples / +0.02 lines, subswaths exactly where trans.dat puts them, no
# zero strip). The merged grids gain up to a row and a column, so the lab's merged-frame masks are padded (pad_mask.py).
# PROC is stored with every pair (igml_proc.json): pairs from older code are made again; REPLACE swaps the lab's pairs
# for ours per track, its originals MOVED (same file system) to LAB_ARCHIVE with a manifest -- nothing deleted.
PROC_FIX = Path("/zdata/datasets/IGML/gmtsar_fix/bin_proc")
EAP_FIX = Path("/zdata/datasets/IGML/gmtsar_fix/bin_eap")
AUX_CAL = Path("/zdata/datasets/IGML/gmtsar_fix/aux/S1A_AUX_CAL_V20140915T100000_G20151125T103928.SAFE/data/"
               "s1a-aux-cal.xml")
PAD_MASK = Path("/zdata/datasets/IGML/gmtsar_fix/pad_mask.py")
# 10-30 afternoon (user: "I think we should fix A1. My philosophy is the data should be as good as it gets ... We do want
# T173 and all its outputs delivered first"; "if spectral diversity and A3 work, fold them in"): every date aligned again
# by us -- no image of the lab's is used (Track.lab_aligned) -- with make_s1a_tops from patch_eap.py + patch_a1.py
# (bin_a1: IGML_A1 cuts each date between bursts at the super-master's lines, audit A1) and, once validated, GMTSAR's
# spectral-diversity alignment (preproc_batch_tops_esd.csh, ESD_MODE: audit A5). start_realign sets each track's old
# alignment aside once (OUT/<T>/superseded_<old>/raw, nothing deleted).
PROC = "proc2"
REALIGN_ALL = os.environ.get("REALIGN_ALL", "1") == "1"
A1_FIX = Path("/zdata/datasets/IGML/gmtsar_fix/bin_a1")
# 10-31 (user: "Let's try to fix A3. If it is easier for you, we can fork GMTSAR and maintain our own version"): the IGML
# FORK of GMTSAR, /zdata/datasets/IGML/gmtsar_igml/src (git, branch igml = upstream master 3dd43d3 + one commit per fix:
# A0, I1/I2/I3/M1/M6, A2, F6/A4, A1, A3; its unmodified build reproduces the lab's 51 programs and 100 scripts byte for
# byte, and each fix's build our validated one), installed in gmtsar_igml/install. With it every GMTSAR call (align,
# interferograms, merges, unwrapping) goes to the fork -- /usr/local/GMTSAR leaves the PATH -- with IGML_A1 and, once
# validated, IGML_A3 (A3_ON); the phasediff slot gate (FORK_WRAP) sits in front of the fork's phasediff. Our tools
# (orbit-download stubs, the lab's unwrap script) stay ahead of it.
GMTSAR_FORK = Path("/zdata/datasets/IGML/gmtsar_igml/install/bin")
FORK_WRAP = Path("/zdata/datasets/IGML/gmtsar_igml/wrap")
USE_FORK = GMTSAR_FORK.is_dir() and os.environ.get("USE_FORK", "1") == "1"
A3_ON = os.environ.get("A3_ON", "0") == "1"
# 10-31 (user: "go ahead on decision 2"): NETWORK SPECTRAL DIVERSITY (audit A5: geometry-only alignment leaves a phase
# step at every burst boundary; GMTSAR's ESD against an 11-year-old super-master only half-removes it). Pass 1 aligns
# every date geometry-only and keeps each frame's burst overlaps (the fork's preproc_batch_tops_esd.csh with
# IGML_ESD_COLLECT); nesd() measures the residual azimuth shift between each date and its next NESD_K dates (GMTSAR's
# spectral_diversity on the pair) and solves it per date frame against the super-master; pass 2 aligns every date again
# with those shifts (IGML_ESD_SHIFTS). Interferograms only from pass 2.
NESD = USE_FORK and os.environ.get("NESD", "1") == "1"
NESD_K = int(os.environ.get("NESD_K", "4"))
NESD_PROCS = int(os.environ.get("NESD_PROCS", "16"))
NESD_MIN_ROWS = 1000          # ddphase rows (coherent overlap samples) a pair needs
NESD_MIN_PHASE = 0.02         # rad: floor of a pair's phase sigma (the formal one ignores the ionosphere, systematics)
SPECDIV_FILTER = "gauss25x7"  # as preproc_batch_tops_esd.csh
ALIGN_PAR = int(os.environ.get("ALIGN_PAR", "3"))      # chunks of dates aligned at once (10-31: seasat's 64 cores)
# 11-01: interferograms per subswath at once (intf_tops_parallel.csh's Ncores; phasediff, the read-heavy step, stays
# behind its slot gate, IGML_PHASEDIFF_SLOTS = 16 in all)
INTF_NCORES = int(os.environ.get("INTF_NCORES", "6"))
# 11-02: the fork's filter.csh computes each reference date's decimated amplitude (amp1.grd) once per subswath and copies
# it for the date's other pairs (IGML_AMP_CACHE; ~2 of a pair's ~8.5 min; byte-identical: audit/ampcache_test.py, _par_test.py)
AMP_CACHE = os.environ.get("AMP_CACHE", "1") == "1"
ESD_MODE = os.environ.get("ESD_MODE", "")      # preproc_batch_tops_esd.csh's esd_mode (1 = median); "" = the plain one
EAP_IPF = "002.43"            # scenes from an older IPF still carry the elevation antenna pattern phase
LAB_ARCHIVE = Path("/zdata/datasets/IGML/lab_originals_20261029")
REPROCESS_LAB = os.environ.get("REPROCESS_LAB", "1") == "1"     # the lab's pairs too (and every date they need)
GS = LAB_ROOT / "seasatHD1/Files_for_InSAR_Processing/GeneralScripts"
ORB = Path("/palsar/IGML/orbits")
DOWNLOAD_ORBITS = not SEASAT         # ASF s1qc needs the Earthdata login, which stays on the Mac (prefetch_orbits.py)
ARCHIVE = "/palsar"     # raw scenes, the lab's and ours, since the 10-27 reorganisation: <ARCHIVE>/S1<A-D>/CA/<year>/T<track>/
                        # (was /palsar/IGML/S1A/CA: the lab's archive mirrored by hard links + our downloads)
# 10-27: GMTSAR's organize_files_tops.csh is the macOS script (BSD `date -jf` / `date -v-1d`); on Linux those calls fail,
# the orbit name comes out empty and ext_orb_s1a loops (seasat: 3 runs x 1 h 40 min). The _linux twin differs only in
# the date syntax, numeric comparisons and the orbit-download stub's name (download_sentinel_orbits_linux.csh, stubbed).
ORGANIZE = "organize_files_tops_linux.csh" if SEASAT else "organize_files_tops.csh"
MAXDT = int(os.environ.get("MAXDT", "84"))
# 10-27 (user: "let's add year-long interferograms from June-August for all stacks and all years"): every pair of
# June-August dates in consecutive years 365 +- ANNUAL_TOL days apart, between any two dates with a usable aligned image
# (the lab's too: its stacks have no pairs this long). ~2,000 pairs over the nine stacks at 12 days.
ANNUAL_TOL = int(os.environ.get("ANNUAL_TOL", "12"))
SUMMER = ("06", "07", "08")
STAGES = ["base", "plan", "stitch", "rescue", "align", "qc", "nesd", "intf", "merge", "unwrap", "store", "combine"]
# 10-27 (user: "We can try to recover them"): dates whose scenes miss one end of the lab's frame (GMTSAR: "SKIP ..., as
# it does not have enough scenes") are stitched with the outer pin slid inside the date's scenes. Only dates whose scenes
# cover >= RESCUE_MIN_COVER of the swath (ASF footprints) and keep as much of the lab's pin line; qc then needs data in
# QC_COVER x that share of a neighbour's windows. 10-29 (user: "over 60% sounds good"): 0.88 -> 0.60, and dates with a
# gap keep their longest contiguous run of scenes (datatake_rescue). 10-30 afternoon (user: "We only want complete
# interferograms, not partial ones. We can revert to the old rule"): 1.0 -- a frame must keep the lab's whole pin line;
# the partial frames cut since 10-27 are withdrawn (start_realign).
RESCUE_MIN_COVER = 1.0
RESCUE_MARGIN = 300           # lines inside the first / last scene, in every subswath
day = lambda s: datetime.strptime(s, "%Y%m%d")


def doy0(d):
    t = day(d)
    return f"{t.year}{(t - datetime(t.year, 1, 1)).days:03d}"


def ymd0(s):
    """GMTSAR's yyyyddd (day 0 = 1 January) -> YYYYMMDD; the inverse of doy0."""
    return (datetime(int(s[:4]), 1, 1) + timedelta(days=int(s[4:]))).strftime("%Y%m%d")


def unwrap_cmd(T, a, b):
    """(SNAPHU recipe, coseismic?) of a pair (a, b), YYYYMMDD: DEFO when the pair overlaps one of the track's earthquakes
    -- a merged sequence from its first to its last member that qualifies on its own (Ridgecrest 2019-07-04 .. 07-06;
    T166's 2019-07-05 sits between them, 10-29) -- by date, both ends included; SMOOTH otherwise (user 10-26 / 10-29:
    "smooth unless there is an earthquake", the lab's pairs too)."""
    co = any(a <= e.get("end_date", e["date"]) and e["date"] <= b for e in COSEISMIC.get(T, []))
    return f"{UNWRAP_SCRIPT.get(T, 'snaphu_interp.csh')} 0.0001 {DEFOMAX if co else 0}", co


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def guard_extorb(limit=600, every=60):
    """Kill any ext_orb_s1a of this run (same session) past ``limit`` s of CPU; it normally needs ~1 s. GMTSAR's orbit
    interpolation loops forever on an orbit that does not cover the scene (10-24 / 10-25 on the Mac: guard_extorb.sh;
    10-27 on seasat: an empty orbit name). The date then gets no frame and the run goes on. Linux /proc only."""
    sid, tck = os.getsid(0), os.sysconf("SC_CLK_TCK")
    while True:
        time.sleep(every)
        for p in Path("/proc").glob("[0-9]*"):
            try:
                s = (p / "stat").read_text()
                f = s[s.rindex(")") + 2:].split()          # state ppid pgrp session ... utime (11) stime (12)
                cpu = (int(f[11]) + int(f[12])) / tck
                if s[s.index("(") + 1:s.rindex(")")] == "ext_orb_s1a" and int(f[3]) == sid and cpu > limit:
                    cwd = os.readlink(p / "cwd")
                    os.kill(int(p.name), signal.SIGTERM)
                    log(f"guard: killed ext_orb_s1a {p.name} ({cpu:.0f} s CPU) in {cwd}")
            except (OSError, ValueError, IndexError):
                continue


def _holds_data(f):
    """f exists and at least half of its length is stored on disk (an emptied SLC is a sparse file of zeros)."""
    try:
        st = os.stat(f)
    except OSError:
        return False
    return st.st_size > 0 and st.st_blocks * 512 >= 0.5 * st.st_size


def gmt_env():
    if USE_FORK:                                                       # 10-31: the IGML fork of GMTSAR
        env = dict(os.environ, PATH=f"{FORK_WRAP}:{TOOLS}:{GMTSAR_FORK}:{os.environ['PATH']}", IGML_A1="1")
        if A3_ON:
            env["IGML_A3"] = "1"
        return env
    fix = f"{PROC_FIX}:" if PROC_FIX.is_dir() else ""                 # merge_batch -> our merge_swath (10-29)
    return dict(os.environ, PATH=f"{fix}{TOOLS}:/usr/local/GMTSAR/bin:/opt/homebrew/bin:{os.environ['PATH']}")


def proc_env(tr, n):
    """intf_tops_parallel.csh's environment for subswath n: our filter.csh / conv / phasefilt first, and the
    subswath's origin on the merged grid (10-29)."""
    x0, y0 = tr.origins[n]
    env = dict(gmt_env(), IGML_SUB_X0=str(x0), IGML_SUB_Y0=str(y0))
    if AMP_CACHE and USE_FORK:
        env["IGML_AMP_CACHE"] = str(tr.work / f"intf/F{n}/ampcache")
    return env


def _prune_amp_cache(tr, keep):
    """11-02: drop the cached reference amplitudes (IGML_AMP_CACHE, ~13 MB each) of dates outside keep."""
    for n in tr.subs:
        for f in (tr.work / f"intf/F{n}/ampcache").glob("S1_*"):
            if f.name[3:11] not in keep:
                f.unlink(missing_ok=True)


def subswath_origins(labp, subs=(1, 2, 3)):
    """{n: (x0 samples, y0 lines)} where subswath n's decimated cells must start so that they fall on the merged grid's
    32 x 6 cells: its near-range and start-time offsets from the first subswath (F1 on the nine stacks; whole samples
    and lines on every stack; checked) taken up to the next multiple (10-29). One subswath: (0, 0)."""
    p = {n: _prm(labp / f"F{n}/topo/master.PRM") for n in subs}
    f = subs[0]
    dr = 299792458.0 / (2 * float(p[f]["rng_samp_rate"]))
    prf = float(p[f]["PRF"])
    out = {}
    for n in subs:
        c = (float(p[n]["near_range"]) - float(p[f]["near_range"])) / dr
        e = (float(p[n]["clock_start"]) - float(p[f]["clock_start"])) * 86400 * prf
        if abs(c - round(c)) > 0.01 or abs(e - round(e)) > 0.01:
            raise RuntimeError(f"{labp}: F{n} is {c:.4f} samples / {e:.4f} lines from F1, not whole: no exact merge")
        out[n] = ((-round(c)) % 32, (-round(e)) % 6)
    return out


def run(cmd, cwd, logfile=None, check=True):
    env = gmt_env()
    with open(logfile or os.devnull, "a") as fh:
        r = subprocess.run(cmd, shell=True, cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if check and r.returncode != 0:
        raise RuntimeError(f"{cmd} in {cwd}: exit {r.returncode} (log {logfile})")
    return r.returncode


class Track:
    def __init__(self, T):
        self.T = T
        if not SEASAT:            # this Mac's NFS writes to seasat hang; the outputs live there since 2026-10-27
            sys.exit("process_track.py runs on seasat since 2026-10-27 (the Mac version: process_track_mac_20261027.py)")
        self.new = T in NEW                 # 10-29: no lab stack; the base stage makes ours in its place
        self.subs = NEW.get(T, (1, 2, 3))   # the subswaths processed
        self.lab = OUT_ROOT / T / "base" if self.new else LAB_ROOT / f"seasatHD{VOL[T]}/InSAR_Processing/{T}"
        self.labp, self.labd = self.lab / f"product_{T}", self.lab / f"data_{T}"
        self.out = OUT_ROOT / T
        self.outp, self.outd = self.out / f"product_{T}", self.out / f"data_{T}"
        self.work = WORK_ROOT / T
        self.out.mkdir(parents=True, exist_ok=True)
        self.work.mkdir(parents=True, exist_ok=True)
        self.state_path = self.out / "state.json"
        self.state = json.load(open(self.state_path)) if self.state_path.exists() else {}
        self.master = self.master_pair = self.origins = self.pins = None
        if not self.new or (self.lab / "base.json").exists():
            self.load_base()

    def load_base(self):
        """The super-master, its pair, the subswath origins and the pins: from the lab's stack, or from ours (NEW)."""
        cfg = (self.labp / f"F{self.subs[0]}/batch_tops.config").read_text()
        self.master = re.search(rf"master_image\s*=\s*S1_(\d{{8}})_ALL_F{self.subs[0]}", cfg).group(1)
        mm = (self.labp / "merge/master_merge").read_text().split(",")[0]
        a, b = re.search(r"intf_all/(\d{7})_(\d{7})/", mm).groups()
        conv = lambda c: (datetime(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d")
        self.master_pair = (conv(a), conv(b))
        self.origins = subswath_origins(self.labp, self.subs)
        T = self.T
        # 10-28: T064's data_ also holds old_T64_pins.ll (north pin 35.54 N, not 35.87), which the glob gave first: our
        # T064 frames were cut ~2 bursts short (F0557_F0568, the lab's F0557_F0574)
        pins = [p for p in sorted(self.labd.glob("*pins.ll")) if not p.name.lower().startswith("old")]
        if len(pins) != 1:
            raise RuntimeError(f"{T}: {len(pins)} current pins files in {self.labd}: {[p.name for p in pins]}")
        self.pins = pins[0]

    def save(self):
        # 10-31: written whole beside it, then renamed over it -- a run stopped (or dying) mid-save no longer leaves a
        # truncated state.json
        tmp = Path(f"{self.state_path}.tmp")
        with open(tmp, "w") as fh:
            json.dump(self.state, fh, indent=1)
        os.replace(tmp, self.state_path)

    def lab_dates(self):
        """Dates with a complete merged pair in the lab's stack (unwrap.grd and corr.grd), cached in the state with the
        mtime of the lab's merge/ (10-27: the live stack had grown from the cached 366 dates to 375)."""
        mt = (self.labp / "merge").stat().st_mtime
        if "lab_dates" in self.state and self.state.get("lab_merge_mtime") == mt:
            return set(self.state["lab_dates"])
        self.state["lab_merge_mtime"] = mt
        out = set()
        for d in (self.labp / "merge").iterdir():
            m = re.match(r"^(\d{7})_(\d{7})$", d.name)
            if m and (d / "unwrap.grd").exists() and (d / "corr.grd").exists():
                for c in m.groups():
                    out.add((datetime(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d"))
        self.state["lab_dates"] = sorted(out)
        self.save()
        return out

    def lab_aligned(self, d):
        """The lab has an aligned image of d we can use: all three subswaths, holding data (10-27: the lab had emptied
        its T064 SLCs -- all 295 -- and T137's / T042's 2025 ones: full-length sparse files with nothing stored).
        A NEW track's base holds only the super-master. 10-30 afternoon: none with REALIGN_ALL -- every date is ours."""
        if REALIGN_ALL:
            return False
        return all(_holds_data(self.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC") for n in self.subs)


# ---------------------------------------------------------------------------------------------------- orbits
_LIST = {}
# 10-25: OFF. GMTSAR's organize_files_tops.csh (line 132) takes only an EOF whose NAME holds the day before and the day
# after (a precise orbit's window); with a restituted one it finds nothing, downloads nothing, and ext_orb_s1a loops on
# an empty orbit (T173 2026-10-09 again, 4 h, after the 10-24 _covers fix). Dates without a POEORB wait (~20 days).
USE_RESORB = False


def _covers(name, when, margin=600):
    """Does an EOF's validity window (..._V<start>_<stop>.EOF) cover ``when`` (YYYYMMDDTHHMMSS) +- margin seconds?"""
    m = re.search(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$", name)
    if not m or not when:
        return False
    t0, t1 = (datetime.strptime(x, "%Y%m%dT%H%M%S") for x in m.groups())
    t = datetime.strptime(when, "%Y%m%dT%H%M%S")
    return t0 <= t - timedelta(seconds=margin) and t + timedelta(seconds=margin) <= t1


def orbit_for(tr, d, sat, when=None):
    """A precise orbit (else restituted) covering date d: the lab's data_<T>, /palsar/IGML/orbits, or ASF s1qc.
    A restituted orbit spans ~3.5 h: it must cover the acquisition time ``when`` (10-24: one that only started on the
    right day sent GMTSAR's ext_orb_s1a into an endless loop on T173 2026-10-09 / 10-15)."""
    n1, n2 = (day(d) - timedelta(days=1)).strftime("%Y%m%d"), (day(d) + timedelta(days=1)).strftime("%Y%m%d")
    for root in (tr.labd, ORB):
        hit = sorted(root.glob(f"{sat}_OPER_AUX_POEORB_*_V{n1}T*_{n2}T*.EOF"))
        if hit:
            return hit[-1]
    if not DOWNLOAD_ORBITS:
        return None
    if "poeorb" not in _LIST:
        html = subprocess.run(["curl", "-s", "-m", "120", "https://s1qc.asf.alaska.edu/aux_poeorb/"],
                              capture_output=True, text=True).stdout
        _LIST["poeorb"] = sorted(set(re.findall(r"S1[ABCD]_OPER_AUX_POEORB_[^\"<>]*?\.EOF", html)))
    hit = [n for n in _LIST["poeorb"] if n.startswith(sat) and f"_V{n1}T" in n and f"_{n2}T" in n]
    kind = "aux_poeorb"
    if not hit and USE_RESORB:                             # recent date: restituted orbits cover a few hours each
        if "resorb" not in _LIST:
            html = subprocess.run(["curl", "-s", "-m", "300", "https://s1qc.asf.alaska.edu/aux_resorb/"],
                                  capture_output=True, text=True).stdout
            _LIST["resorb"] = sorted(set(re.findall(r"S1[ABCD]_OPER_AUX_RESORB_[^\"<>]*?\.EOF", html)))
        hit = [n for n in _LIST["resorb"] if n.startswith(sat) and _covers(n, when)]
        kind = "aux_resorb"
    if not hit:
        return None
    ORB.mkdir(parents=True, exist_ok=True)
    f = ORB / hit[-1]
    if not f.exists():
        subprocess.run(["curl", "-s", "-n", "-L", "-c", str(Path.home() / ".urs_cookies"), "-b",
                        str(Path.home() / ".urs_cookies"), "-o", str(f), f"https://s1qc.asf.alaska.edu/{kind}/{hit[-1]}"],
                       check=True)
    return f


# ---------------------------------------------------------------------------------------------------- NEW tracks
def _footprint(safe):
    """A scene's footprint, [(lon, lat)] from manifest.safe's gml:coordinates (lat,lon pairs); None if unreadable."""
    try:
        m = re.search(r"<gml:coordinates>([^<]+)</gml:coordinates>",
                      (Path(safe) / "manifest.safe").read_text(errors="replace"))
    except OSError:
        return None
    if not m:
        return None
    return [(float(p.split(",")[1]), float(p.split(",")[0])) for p in m.group(1).split()]


def _inside(x, y, poly):
    """Is (x, y) inside the polygon [(x, y), ...] (ray casting)?"""
    c = False
    for (x0, y0), (x1, y1) in zip(poly, poly[1:] + poly[:1]):
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            c = not c
    return c


def _pass(safe):
    m = re.search(r"<s1:pass>(\w+)</s1:pass>", (Path(safe) / "manifest.safe").read_text(errors="replace"))
    return m.group(1) if m else None


def target_line(tr, descending):
    """The target's pin line (10-29): its outline (extra_tracks/subswaths.json, the processed subswaths' California
    part) along its long axis, through the outline's middle, from the first vertex the satellite reaches to the last
    (descending: north to south). [(lon, lat), (lon, lat)] in flight order."""
    import numpy as np
    o = np.array(json.load(open(REPO / EXTRA))[tr.T]["outline"][:-1], dtype=float)
    lat0 = o[:, 1].mean()
    xy = np.c_[(o[:, 0] - o[:, 0].mean()) * np.cos(np.radians(lat0)), o[:, 1] - lat0]
    ax = np.linalg.svd(xy, full_matrices=False)[2][0]            # the long axis
    if (ax[1] < 0) != descending:                                # pointing the way the satellite flies
        ax = -ax
    t = xy @ ax
    ends = [t.min() * ax, t.max() * ax]
    return [(round(float(o[:, 0].mean() + e[0] / np.cos(np.radians(lat0))), 6), round(float(lat0 + e[1]), 6))
            for e in ends]


PIN_SHARE, PIN_SAMPLE = 0.9, 6     # pin_line: the share of the dates an end must reach; dates measured where it is not


def pin_line(tr, by, fp, desc):
    """A NEW track's pins (10-29): the target's line (target_line), an end that fewer than PIN_SHARE of the dates
    reaching the line cover (by the scenes' footprints) moved inside, to where PIN_SAMPLE of the dates short there
    (those that reach the first point PIN_SHARE of the dates cover; spread in time) have it RESCUE_MARGIN lines inside
    every subswath of their first / last scene (_clip_line: GMTSAR's geometry; the footprints are ~10 km off -- T108's
    S1A scenes of 2021-12 .. 2025 start at 37.95 N by their footprints, iw1 at 38.06 N, the target at 37.85 N: IW3 off
    Point Reyes, ocean). T100 / T013: unchanged. ``by`` = {date: [SAFE, ...]}, ``fp`` = {SAFE: footprint}. Once the
    base exists, its pins."""
    import math
    if (tr.lab / "base.json").exists():
        return [tuple(p) for p in json.load(open(tr.lab / "base.json"))["pins"]]
    ends = target_line(tr, desc)
    (x0, y0), (x1, y1) = ends
    K = 200
    line = [(x0 + k / K * (x1 - x0), y0 + k / K * (y1 - y0)) for k in range(K + 1)]
    cov = {d: [any(fp[p] and _inside(x, y, fp[p]) for p in ps) for x, y in line] for d, ps in by.items()}
    cov = {d: c for d, c in cov.items() if any(c)}
    share = [sum(c[k] for c in cov.values()) / len(cov) for k in range(K + 1)]
    ok = [k for k in range(K + 1) if share[k] >= PIN_SHARE]
    if not ok:
        raise RuntimeError(f"{tr.T}: no point of the target's line is covered by {PIN_SHARE:.0%} of the dates")
    km = math.hypot((x1 - x0) * 111.32 * math.cos(math.radians((y0 + y1) / 2)), (y1 - y0) * 110.57) / K
    t = lambda p: ((p[0] - x0) * (x1 - x0) + (p[1] - y0) * (y1 - y0)) / ((x1 - x0) ** 2 + (y1 - y0) ** 2)   # 0 .. 1
    out = list(ends)
    for e, (k_end, k_ok) in enumerate(((0, ok[0]), (K, ok[-1]))):
        if k_end == k_ok:
            continue
        short = sorted(d for d, c in cov.items() if c[k_ok] and not c[k_end] and not any(_raw_bad(s) for s in by[d]))
        sample = [short[round(i * (len(short) - 1) / max(PIN_SAMPLE - 1, 1))] for i in range(min(PIN_SAMPLE, len(short)))]
        got = []
        for d in dict.fromkeys(sample):
            f = Path(by[d][0]).name
            orb = orbit_for(tr, d, f[:3], f[17:32])
            if orb is None:
                continue
            w = tr.work / f"pins_{d}"
            shutil.rmtree(w, ignore_errors=True)
            w.mkdir(parents=True)
            new, kept = _clip_line(ends, sorted(by[d], key=lambda s: Path(s).name[17:32]), orb, w)
            shutil.rmtree(w, ignore_errors=True)
            if new:
                got.append(t(new[0] if e == 0 else new[-1]))
        if not got:
            raise RuntimeError(f"{tr.T}: no date short at the line's {'start' if e == 0 else 'end'} could be measured")
        f = max(got) if e == 0 else min(got)                     # the most inside of the sample
        out[e] = (x0 + f * (x1 - x0), y0 + f * (y1 - y0))
        log(f"{tr.T}: the pin line's {'start' if e == 0 else 'end'} moved {(f if e == 0 else 1 - f) * K * km:.1f} km "
            f"inside the target: {share[k_end]:.0%} of the dates reach it ({len(short)} short; {len(got)} measured, "
            f"{RESCUE_MARGIN} lines inside every subswath)")
    return [(round(x, 6), round(y, 6)) for x, y in out]


def new_scenes(tr):
    """{date: {"safes": [...], "cover": share of the pin line inside the scenes' footprints}} for a NEW track (10-29):
    every scene of the track in the archive (one copy of each, by start / stop time), the dates' scenes that reach the
    pin line (pin_line), in time order. Cached in the state with the archive's listing and the line."""
    names = real(glob.glob(f"{ARCHIVE}/S1?/CA/20*/{tr.T}/*.SAFE"))
    key = f"v3:{len(names)}:{max((os.path.getmtime(p) for p in names), default=0):.0f}"
    if tr.state.get("new_scenes", {}).get("key") == key:
        return tr.state["new_scenes"]["dates"]
    one = {}
    for p in names:
        one.setdefault(Path(p).name[17:48], p)                   # the same scene twice: one copy
    by = defaultdict(list)
    for p in sorted(one.values(), key=lambda p: Path(p).name[17:32]):
        by[Path(p).name[17:25]].append(p)
    desc = _pass(next(iter(one.values()))) == "DESCENDING"
    fp = {p: _footprint(p) for ps in by.values() for p in ps}
    pins = pin_line(tr, by, fp, desc)
    (x0, y0), (x1, y1) = pins
    line = [(x0 + k / 60 * (x1 - x0), y0 + k / 60 * (y1 - y0)) for k in range(61)]
    out = {}
    for d, ps in sorted(by.items()):
        hit = [any(fp[p] and _inside(x, y, fp[p]) for p in ps) for x, y in line]
        use = [p for p in ps if fp[p] and any(_inside(x, y, fp[p]) for x, y in line)]
        if use:
            out[d] = {"safes": use, "cover": round(sum(hit) / len(hit), 3)}
    tr.state["new_scenes"] = {"key": key, "pass": "DESCENDING" if desc else "ASCENDING", "line": [list(p) for p in pins],
                              "dates": out}
    tr.save()
    return out


def _template_config(n, master):
    """The lab's batch_tops.config (T173's: every stack's is the same but for master_image) for subswath n."""
    return re.sub(r"master_image\s*=\s*\S+", f"master_image = S1_{master}_ALL_F{n}", CONFIG_TEMPLATE.read_text())


def build_base(tr):
    """The base of a NEW track (10-29), what the lab's stack holds for the nine, under OUT/<T>/base/ in the lab's layout:
    data_<T>/<T>_pins.ll (pin_line) with the super-master's stitched frame(s) and orbit;
    product_<T>/topo/dem.grd (SRTM 3" + EGM96, DEM_SRC/<T>); per processed subswath n, F<n>/batch_tops.config (the lab's
    recipe), raw/ the super-master SLC / PRM / LED (preproc_batch_tops.csh modes 1 and 2, as align re-makes it), topo/
    master.PRM, topo_ra.grd, trans.dat (dem2topo_ra.csh); merge/batch_tops.config, master_merge (the super-master
    pair) and, for one subswath, its trans.dat (two: the first merge makes it, merge copies it back); base.json.
    The super-master: among the dates whose scenes cover the whole pin line, processed by IPF >= EAP_IPF (no antenna
    pattern) with a precise orbit and every scene whole, the one nearest the middle of the record (S1's orbital tube is
    narrow: time is the spread that matters); its partner: the next such date."""
    if (tr.lab / "base.json").exists():
        return
    sc = new_scenes(tr)
    desc = tr.state["new_scenes"]["pass"] == "DESCENDING"
    pins = [tuple(p) for p in tr.state["new_scenes"]["line"]]
    full = []
    for d, v in sorted(sc.items()):
        if v["cover"] < 1 or any(_raw_bad(s) for s in v["safes"]) or min(_ipf(s) or "0" for s in v["safes"]) < EAP_IPF:
            continue
        f = Path(v["safes"][0]).name
        if orbit_for(tr, d, f[:3], f[17:32]):
            full.append(d)
    if len(full) < 2:
        raise RuntimeError(f"{tr.T} base: {len(full)} dates cover the whole target with an orbit -- no super-master")
    t0, t1 = day(min(sc)), day(max(sc))
    mid = t0 + (t1 - t0) / 2
    master = min(full, key=lambda d: abs(day(d) - mid))
    partner = next((d for d in full if d > master), None) or max(d for d in full if d < master)
    pair = tuple(sorted((master, partner)))
    dem = DEM_SRC / tr.T / "dem.grd"
    if not dem.exists():
        raise RuntimeError(f"{tr.T} base: no {dem} (make it on the Mac: make_dem.csh W E S N 2, then scp)")
    log(f"{tr.T} base: super-master {master} (nearest {mid:%Y-%m-%d}, the middle of {len(sc)} dates; {len(full)} "
        f"cover the whole target), its pair {pair}, pins {pins}, subswaths {tr.subs}")
    tmp = tr.lab.with_name("base.part")                      # built aside, renamed at the end: no half base
    shutil.rmtree(tmp, ignore_errors=True)
    bd, bp = tmp / f"data_{tr.T}", tmp / f"product_{tr.T}"
    bd.mkdir(parents=True)
    (bd / f"{tr.T}_pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\n" for x, y in pins))
    f = Path(sc[master]["safes"][0]).name
    orb = orbit_for(tr, master, f[:3], f[17:32])
    os.symlink(orb, bd / orb.name)
    w = tr.work / "base_stitch"
    shutil.rmtree(w, ignore_errors=True)
    w.mkdir(parents=True)
    shutil.copy(bd / f"{tr.T}_pins.ll", w / "pins.ll")
    (w / "safelist").write_text("\n".join(sc[master]["safes"]) + "\n")
    os.symlink(orb, w / orb.name)
    run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log")
    frames = real(w.glob("F*_F*/*.SAFE"))
    if not frames:
        raise RuntimeError(f"{tr.T} base: organize cut no frame of {master} (see {w}/organize.log)")
    for s in frames:
        shutil.copytree(s, bd / s.parent.name / s.name, symlinks=False)
    msafes = real(bd.glob(f"F*_F*/*_{master}T*.SAFE"))
    (bp / "topo").mkdir(parents=True)
    shutil.copy(dem, bp / "topo/dem.grd")
    env = dict(os.environ, PATH=(f"{EAP_FIX}:" if EAP_FIX.is_dir() else "") + (f"{STITCH_FIX}:" if STITCH_FIX.is_dir()
               else "") + f"{TOOLS}:/usr/local/GMTSAR/bin:{os.environ['PATH']}")
    for n in tr.subs:
        a = tr.work / f"base_align/F{n}"
        shutil.rmtree(a, ignore_errors=True)
        a.mkdir(parents=True)
        rec = []
        for t, x in _frames_line(msafes, n):
            os.symlink(t, a / t.name)
            os.symlink(x, a / x.name)
            rec.append(t.stem)
        os.symlink(orb, a / orb.name)
        (a / "data.in").write_text(":".join(rec + [orb.name]) + "\n")
        shutil.copy(dem, a / "dem.grd")
        r = subprocess.run("preproc_batch_tops.csh data.in dem.grd 1 > pbt1.log 2>&1 && preproc_batch_tops.csh data.in "
                           "dem.grd 2 > pbt2.log 2>&1", shell=True, cwd=a, env=env)
        raw, topo = bp / f"F{n}/raw", bp / f"F{n}/topo"
        raw.mkdir(parents=True)
        topo.mkdir(parents=True)
        for ext in ("PRM", "LED", "SLC"):
            g = a / f"S1_{master}_ALL_F{n}.{ext}"
            if r.returncode or not g.exists():
                raise RuntimeError(f"{tr.T} base F{n}: {g.name} was not made (see {a}/pbt2.log)")
            shutil.copy(g, raw / g.name)
        shutil.copy(raw / f"S1_{master}_ALL_F{n}.PRM", topo / "master.PRM")
        os.symlink(f"../raw/S1_{master}_ALL_F{n}.LED", topo / f"S1_{master}_ALL_F{n}.LED")
        os.symlink("../../topo/dem.grd", topo / "dem.grd")
        run("dem2topo_ra.csh master.PRM dem.grd", topo, topo / "dem2topo_ra.log")
        if not (topo / "topo_ra.grd").exists() or not (topo / "trans.dat").exists():
            raise RuntimeError(f"{tr.T} base F{n}: dem2topo_ra.csh made no topo_ra.grd / trans.dat (see its log)")
        (bp / f"F{n}/batch_tops.config").write_text(_template_config(n, master))
        shutil.rmtree(a, ignore_errors=True)
    (bp / "merge").mkdir()
    (bp / "merge/batch_tops.config").write_text(_template_config(tr.subs[0], master))
    dn = f"{doy0(pair[0])}_{doy0(pair[1])}"
    (bp / "merge/master_merge").write_text(",".join(f"../F{n}/intf_all/{dn}/:S1_{pair[0]}_ALL_F{n}.PRM:S1_{pair[1]}_ALL_F{n}"
                                                    f".PRM" for n in tr.subs) + "\n")
    if len(tr.subs) == 1:                                    # one subswath: its own geometry is the frame's
        shutil.copy(bp / f"F{tr.subs[0]}/topo/trans.dat", bp / "merge/trans.dat")
    (tmp / "base.json").write_text(json.dumps({
        "time": time.strftime("%F %T"), "master": master, "pair": list(pair), "subswaths": list(tr.subs),
        "pins": pins, "pass": "DESCENDING" if desc else "ASCENDING", "dates": len(sc), "full_dates": len(full),
        "middle": f"{mid:%Y%m%d}", "dem": str(dem), "config": str(CONFIG_TEMPLATE),
        "frames": sorted(s.parent.name for s in frames)}, indent=1) + "\n")
    os.replace(tmp, tr.lab)
    shutil.rmtree(w, ignore_errors=True)
    log(f"{tr.T} base made: {tr.lab} (frames {sorted(s.parent.name for s in frames)})")


# ---------------------------------------------------------------------------------------------------- stages
def plan(tr, only=None):
    plan_ = {}
    if tr.new:                       # 10-29: the archive's scenes that reach the target, every date (no lab dates)
        sc = new_scenes(tr)
        want = sorted(d for d, v in sc.items() if v["cover"] >= 0.5)
        if only:
            want = [d for d in want if d in only]
    else:
        inv = json.load(open(REPO / "scratchpad/s1_inventory_20261015/s1_acquisitions.json"))
        lab = tr.lab_dates()
        want = sorted(r["date"] for r in inv["rows"]
                      if r["track"] == tr.T and r["coverage"] >= 0.5 and r["date"] not in lab)
        if only:
            want = [d for d in want if d in only]
        # the frames that cover the swath (size_audit/needed_frames.py: the ones on palsar + the ones downloaded), not
        # every frame ASF lists -- edge frames under 1 % of the swath were never downloaded (10-23)
        scenes = {}
        for o in json.load(open(REPO / "scratchpad/s1_gapfill/size_audit/needed_frames.json")):
            if o["track"] == tr.T:
                scenes[o["date"]] = {n[:-5] for n in o["on_palsar"]} | {a["scene"][:-5] for a in o["add"]}
    for d in want:
        paths = []
        slices = []                         # 10-30: one slice per start time (within 3 s): T042 2025-05-01's last
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
        if tr.new:
            paths = list(sc[d]["safes"])
        first = Path(next((p for p in paths if p), "S1A")).name if any(paths) else ""
        sat, when = (first[:3], first[17:32]) if first else (None, None)
        orb = orbit_for(tr, d, sat, when) if paths and all(paths) else None
        la = tr.lab_aligned(d)            # 10-27: the lab aligns new dates before it merges them (its live stacks)
        plan_[d] = {"safes": paths, "ready": la or (bool(paths) and all(paths) and orb is not None),
                    "orbit": str(orb) if orb else None, "lab_aligned": la}
        if tr.state.get("plan", {}).get(d, {}).get("lab_aligned_bad"):     # 10-30: kept (retry_failed, recut)
            plan_[d]["lab_aligned_bad"] = True
        if tr.new:
            plan_[d]["cover"] = sc[d]["cover"]
    tr.state["plan"] = plan_
    tr.save()
    ready = [d for d, v in plan_.items() if v["ready"]]
    nl = sum(v["lab_aligned"] for v in plan_.values())
    log(f"{tr.T}: {len(want)} dates without a merged pair in the lab's stack: {nl} aligned by the lab (its SLCs are "
        f"used), {len(ready) - nl} to align (scenes + orbit), {len(want) - len(ready)} waiting" +
        (" (no lab stack: the base's super-master counts as the lab's)" if tr.new else ""))
    return ready


def _raw_bad(safe):
    """Why a raw scene cannot be used ("" if it can): a VV tiff missing or not the size its manifest.safe declares, or a
    VV annotation missing (GMTSAR audit 10-28: T173 2021-01-01's iw2 tiff is 0.52 of 1.36 GB; make_s1a_tops fails on
    it inside create_frame_tops.csh, which goes on and hands assemble_tops the next scenes' stems as burst limits --
    that subswath silently becomes all bursts of the later scenes)."""
    s = Path(safe)
    try:
        man = (s / "manifest.safe").read_text(errors="replace")
    except OSError:
        return "no manifest.safe"
    size = {m.group(2): int(m.group(1)) for m in re.finditer(
        r'<byteStream[^>]*size="(\d+)"[^>]*>\s*<fileLocation[^>]*href="\./measurement/([^"]+)"', man)
        if "-vv-" in m.group(2)}
    if len(size) != 3:
        return f"{len(size)} vv tiffs in manifest.safe"
    for n, z in size.items():
        t = s / "measurement" / n
        if not t.exists():
            return f"{n} missing"
        if t.stat().st_size != z:
            return f"{n} is {t.stat().st_size / 1e9:.2f} of {z / 1e9:.2f} GB"
    if len(list(s.glob("annotation/*-vv-*.xml"))) != 3:
        return "a vv annotation missing"
    return ""


REDOWNLOAD = "/palsar/IGML/redownload_*"      # complete copies of scenes the archive holds broken (10-28)


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


def _usable(tr, d, safes=None):
    """The date's scenes (the plan's, or `safes`) that pass _raw_bad, in time order; the others go to state "raw_bad"
    {scene: why} (to download again)."""
    rb = tr.state.setdefault("raw_bad", {})
    out = []
    for s in safes if safes is not None else tr.state["plan"][d]["safes"]:
        why = _raw_bad(s)
        if why:
            if rb.get(Path(s).name) != why:
                log(f"{tr.T} {d}: raw scene {Path(s).name} left out: {why}")
            rb[Path(s).name] = why
        else:
            rb.pop(Path(s).name, None)
            out.append(s)
    return sorted(out, key=lambda s: Path(s).name[17:32])


def _one_cut_per_frame(safes):
    """One SAFE per frame directory: the lab's data_<T>/F*_F*/ can hold a date cut twice into one frame (T144 2014-11-08:
    a 2017 and a 2025 cut, names a second apart, the same bursts), which bad_stitches then flags ("an image name in both
    frames") and align skips. The newest cut (the lab's latest) is kept."""
    by = {}
    for s in safes:
        k = Path(s).parent.name
        if k not in by or os.path.getmtime(s) > os.path.getmtime(by[k]):
            by[k] = s
    return sorted(by.values())


def stitch(tr, dates):
    done = tr.state.setdefault("stitched", {})
    fixed = []                                  # 10-29: dates reused with two cuts of a frame before
    for d, v in done.items():
        if v and all(str(x).startswith(str(tr.labd)) for x in v) and len(_one_cut_per_frame(v)) < len(v):
            done[d] = _one_cut_per_frame(v)
            fixed.append(d)
    if fixed:
        log(f"{tr.T} stitch: one cut per frame of the lab's for {fixed} (it holds two)")
        tr.save()
    todo = list(dict.fromkeys(d for d in dates if d not in done))
    # 10-27: where the lab stitched a date (data_<T>/F*_F*/, intact) but emptied its aligned image, align from the lab's
    # frames (read only: align links them) instead of stitching again
    reused = []
    for d in list(todo):
        if d in tr.state.get("recut", {}):      # 10-30: cut from raw (recut), not from the lab's frames again
            continue
        lab = _one_cut_per_frame(real(glob.glob(str(tr.labd / f"F*_F*/*_{d}T*.SAFE"))))
        if lab and all(_holds_data(t) for s in lab for t in real(Path(s).glob("measurement/*vv*.tiff"))):
            done[d] = lab
            todo.remove(d)
            reused.append(d)
    if reused:
        log(f"{tr.T} stitch: the lab's own frames for {len(reused)} dates (its aligned images are empty or failed qc)")
        tr.save()
    before = skipped_before(tr)
    if any(d in before for d in todo):
        log(f"{tr.T} stitch: {sum(d in before for d in todo)} dates GMTSAR skipped before are left to rescue")
    todo = [d for d in todo if d not in before]
    lost = [d for d in todo if d not in tr.state.get("plan", {})]
    if lost:
        log(f"{tr.T} stitch: no usable frames of the lab's for {lost}: not aligned")
    todo = [d for d in todo if d in tr.state.get("plan", {})]
    if not todo:
        return
    w = tr.work / f"stitch_{time.strftime('%Y%m%d%H%M%S')}"
    w.mkdir(parents=True)
    shutil.copy(tr.pins, w / "pins.ll")
    # time order (organize_files_tops groups consecutive lines by date and satellite): by the start time in the
    # name, not the path -- a path sort put S1A 2020-02-12 before S1B 2020-02-06
    safes = sorted((s for d in todo for s in _usable(tr, d)), key=lambda s: Path(s).name[17:32])
    (w / "safelist").write_text("\n".join(safes) + "\n")
    for d in todo:
        o = Path(tr.state["plan"][d]["orbit"])
        if not (w / o.name).exists():
            os.symlink(o, w / o.name)
    log(f"{tr.T} stitch {len(todo)} dates in {w}")
    run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log")
    for f in real(w.glob("F*_F*/*.SAFE")):
        d = re.search(r"_(\d{8})T", f.name).group(1)
        dest = tr.outd / f.parent.name / f.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            shutil.copytree(f, dest, symlinks=False)
            for ad in dest.rglob("._*"):                       # the copy's AppleDouble companions on palsar
                ad.unlink()
        done.setdefault(d, [])
        if str(dest) not in done[d]:
            done[d].append(str(dest))
    missing = [d for d in todo if d not in done]
    nf = tr.state.setdefault("no_frame", {})                      # GMTSAR cut no frame (partial swath, gaps)
    nf.update({d: time.strftime("%F") for d in missing})
    for d in done:
        nf.pop(d, None)
    keep = tr.out / "logs" / f"organize_{w.name}.log"             # the work folder goes; GMTSAR's log stays
    keep.parent.mkdir(exist_ok=True)
    shutil.copy(w / "organize.log", keep)
    if missing:
        log(f"{tr.T} stitch: no frame for {missing} (see {keep}: not enough scenes / gaps)")
    tr.save()
    shutil.rmtree(w, ignore_errors=True)


def _azimuths(prm, pts):
    """Azimuth line of each lon / lat point in a PRM (GMTSAR SAT_llt2rat with its precise orbit; outside the image too)."""
    out = subprocess.run(["SAT_llt2rat", prm.name, "1"], input="".join(f"{x:.6f} {y:.6f} 0\n" for x, y in pts),
                         capture_output=True, text=True, cwd=prm.parent, env=gmt_env()).stdout.split("\n")
    return [float(l.split()[1]) for l in out if l.strip()]


def clipped_pins(tr, d, w, safes=None):
    """The lab's pins (in time order) with the outer ones slid along the pin line until they lie RESCUE_MARGIN lines
    inside the date's first / last scene (of `safes`, default its usable ones) in all three subswaths; the interior pins
    (frame boundaries) are kept. -> (pins, share of the lab's pin line kept) or (None, why)."""
    pins = [l.split() for l in open(tr.pins) if l.strip()]
    if any(len(p) != 2 for p in pins):
        return None, "per-subswath pins: not handled"
    pins = [(float(x), float(y)) for x, y in pins]
    safes = _usable(tr, d) if safes is None else safes
    if not safes:
        return None, "no usable raw scene"
    return _clip_line(pins, safes, Path(tr.state["plan"][d]["orbit"]), w)


def _clip_line(pins, safes, orb, w):
    """clipped_pins for any pins, scenes (time order) and orbit, in the work folder w. All three subswaths whatever
    a track processes: GMTSAR's organize tests the pins in iw1, create_frame_tops.csh cuts all three."""
    os.symlink(orb, w / orb.name)
    prm = {}
    for tag, s in (("first", safes[0]), ("last", safes[-1])):
        for n in (1, 2, 3):
            x = real(Path(s).glob(f"annotation/*iw{n}*vv*.xml"))[0]
            t = real(Path(s).glob(f"measurement/*iw{n}*vv*.tiff"))[0]
            run(f"make_s1a_tops {x} {t} {tag}{n} 0 && ext_orb_s1a {tag}{n}.PRM {orb.name} {tag}{n}", w, w / "pins.log")
            prm[(tag, n)] = w / f"{tag}{n}.PRM"
    K = 200                                                   # candidate points per pin interval
    path = [(x0 + k / K * (x1 - x0), y0 + k / K * (y1 - y0))
            for (x0, y0), (x1, y1) in zip(pins, pins[1:]) for k in range(K)] + [pins[-1]]
    a0 = [_azimuths(prm[("first", n)], path) for n in (1, 2, 3)]
    a1 = [_azimuths(prm[("last", n)], path) for n in (1, 2, 3)]
    nl = [int(_prm(prm[("last", n)])["num_lines"]) for n in (1, 2, 3)]
    if any(len(a) != len(path) for a in a0 + a1):
        return None, "SAT_llt2rat failed"
    start = [i for i in range(len(path)) if all(a[i] >= RESCUE_MARGIN for a in a0)]
    end = [i for i in range(len(path)) if all(a[i] <= m - RESCUE_MARGIN for a, m in zip(a1, nl))]
    if not start or not end or max(end) <= min(start):
        return None, "no part of the pin line inside the scenes"
    i0, i1 = min(start), max(end)
    new = [path[i0]] + [pins[k] for k in range(1, len(pins) - 1) if i0 < k * K < i1] + [path[i1]]
    return new, round((i1 - i0) / (len(path) - 1), 3)


def _runs(safes):
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


def skipped_before(tr):
    """{date: "pins" | "gap"}: dates GMTSAR's organize skipped in an earlier stitch (its logs, kept in OUT/<T>/logs since
    10-27 15:05, and the state's no_frame) -- not re-tried by stitch; rescue takes the "pins" ones."""
    out = {d: "pins" for d in tr.state.get("no_frame", {})}
    for f in sorted((tr.out / "logs").glob("organize_stitch_*.log")):
        for d, why in re.findall(r"^SKIP (\d{8}), as it (does not have enough scenes|stopped observation)",
                                 f.read_text(errors="replace"), re.M):
            out[d] = "pins" if why.startswith("does") else "gap"
    return out


def rescue_candidates(tr):
    """Dates GMTSAR skipped for want of scenes at one end of the frame (its organize logs, kept since 10-27 15:05; the
    state's no_frame), not stitched, covering >= RESCUE_MIN_COVER of the swath. Gaps between frames are not rescued."""
    skipped = set(tr.state.get("no_frame", {}))
    for f in (tr.out / "logs").glob("organize_*.log"):
        skipped |= set(re.findall(r"^SKIP (\d{8}), as it does not have enough scenes", f.read_text(errors="replace"), re.M))
    gaps = set()
    for f in (tr.out / "logs").glob("organize_*.log"):
        gaps |= set(re.findall(r"^SKIP (\d{8}), as it stopped observation", f.read_text(errors="replace"), re.M))
    if tr.new:                                 # 10-29: the share of the target's pin line the scenes cover (plan)
        cov = {d: v.get("cover", 0) for d, v in tr.state.get("plan", {}).items()}
    else:
        cov = {r["date"]: r["coverage"] for r in json.load(open(
            REPO / "scratchpad/s1_inventory_20261015/s1_acquisitions.json"))["rows"] if r["track"] == tr.T}
    plan_ = tr.state.get("plan", {})
    return [d for d in sorted(skipped - gaps) if d in plan_ and plan_[d]["ready"] and not plan_[d].get("lab_aligned")
            and d not in tr.state.get("stitched", {}) and cov.get(d, 0) >= RESCUE_MIN_COVER]


def rescue(tr, dates=None):
    """Stitch rescue_candidates with clipped_pins, one organize run per date (state "rescued": pins, share kept), from
    the date's longest contiguous run of scenes (_runs; all of them unless it has a gap or two versions of a slice);
    >= RESCUE_MIN_COVER of the lab's pin line kept (10-30)."""
    res = tr.state.setdefault("rescued", {})
    todo = [d for d in (dates or rescue_candidates(tr)) if d not in res]
    log(f"{tr.T} rescue {len(todo)} dates: {todo}")
    done = tr.state.setdefault("stitched", {})
    for d in todo:
        w = tr.work / f"rescue_{d}"
        shutil.rmtree(w, ignore_errors=True)
        w.mkdir(parents=True)
        safes = max(_runs(_usable(tr, d)), key=len, default=[])
        pins, kept = clipped_pins(tr, d, w, safes) if safes else (None, "no usable raw scene")
        if pins is None:
            res[d] = {"ok": False, "why": kept}
            log(f"{tr.T} rescue {d}: {kept}")
        elif kept < RESCUE_MIN_COVER:
            res[d] = {"ok": False, "kept": kept, "why": f"{kept * 100:.0f} % of the lab's pin line inside its scenes"}
            log(f"{tr.T} rescue {d}: {res[d]['why']} (< {RESCUE_MIN_COVER * 100:.0f} %)")
        else:
            (w / "pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\n" for x, y in pins))
            (w / "safelist").write_text("\n".join(safes) + "\n")
            run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log", check=False)
            frames = real(w.glob("F*_F*/*.SAFE"))
            for f in frames:
                dest = tr.outd / f.parent.name / f.name
                dest.parent.mkdir(parents=True, exist_ok=True)
                if not dest.exists():
                    shutil.copytree(f, dest, symlinks=False)
                done.setdefault(d, [])
                if str(dest) not in done[d]:
                    done[d].append(str(dest))
            res[d] = {"ok": bool(frames), "pins": [[round(x, 6), round(y, 6)] for x, y in pins], "kept": kept,
                      **({} if frames else {"why": "organize cut no frame with the clipped pins"})}
            log(f"{tr.T} rescue {d}: {len(frames)} frames, {kept * 100:.0f} % of the lab's pin line kept")
            (tr.out / "logs").mkdir(exist_ok=True)
            shutil.copy(w / "organize.log", tr.out / "logs" / f"organize_rescue_{d}.log")
        tr.save()
        shutil.rmtree(w, ignore_errors=True)


def datatake_rescue(tr):
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
        (w / "pins.ll").write_text("".join(f"{x:.6f} {y:.6f}\n" for x, y in pins))
        (w / "safelist").write_text("\n".join(safes) + "\n")
        run(f"{ORGANIZE} safelist pins.ll 2", w, w / "organize.log", check=False)
        frames = real(w.glob("F*_F*/*.SAFE"))
        for f in frames:
            dest = tr.outd / f.parent.name / f.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                shutil.copytree(f, dest, symlinks=False)
            done.setdefault(d, [])
            if str(dest) not in done[d]:
                done[d].append(str(dest))
        res[d] = {"ok": bool(frames), "v": 2, "datatake": take, "kept": kept,
                  "pins": [[round(x, 6), round(y, 6)] for x, y in pins],
                  "left_out": sorted(Path(s).name for s in usable if s not in safes),
                  **({} if frames else {"why": "organize cut no frame from the run"})}
        log(f"{tr.T} datatake {d}: {len(frames)} frames from datatake {take} ({len(safes)} of {len(usable)} scenes, "
            f"{kept * 100:.0f} % of the lab's pin line)")
        (tr.out / "logs").mkdir(exist_ok=True)
        shutil.copy(w / "organize.log", tr.out / "logs" / f"organize_datatake_{d}.log")
        tr.save()
        shutil.rmtree(w, ignore_errors=True)


def real(paths):
    """Drop macOS AppleDouble companions (._name): every file this Mac writes on palsar gets one, and a glob for
    *iw1*vv*.tiff matches them too (10-25: the first align batch fed make_s1a_tops a 4-KB ._ file for 9 h)."""
    return sorted(p for p in paths if not Path(p).name.startswith("._"))


def _frames_line(safes, n):
    """data.in record of one date: every stitched frame's iw<n> stem, in time order."""
    stems = []
    for s in sorted(safes, key=lambda x: Path(x).name):
        t = real(Path(s).glob(f"measurement/*iw{n}*vv*.tiff"))[0]
        stems.append((t, real(Path(s).glob(f"annotation/*iw{n}*vv*.xml"))[0]))
    return stems


_utc = lambda v: datetime.strptime(v.strip()[:26], "%Y-%m-%dT%H:%M:%S.%f").replace(tzinfo=timezone.utc).timestamp()
_BURST = re.compile(r"<burst>\s*<azimuthTime>([^<]+)</azimuthTime>\s*<azimuthAnxTime>([^<]+)</azimuthAnxTime>")
SCENE_TIMING = Path(__file__).resolve().parent / "size_audit/scene_timing.json"    # every raw scene of 10-29
_TIMING = {}


def _timing(x):
    """(linesPerBurst, ANX epoch = the first burst's azimuthTime - azimuthAnxTime [s]) of an annotation, read up to its
    first burst."""
    lpb, data = None, ""
    with open(x, errors="replace") as f:
        while chunk := f.read(1 << 18):
            data += chunk
            if lpb is None and (m := re.search(r"<linesPerBurst>(\d+)<", data)):
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
        cache[d] = any(len({r[n - 1] for r in rows if r[n - 1] is not None}) > 1 for n in tr.subs) or \
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
                lpb = int(re.search(r"<linesPerBurst>(\d+)<", t).group(1))
                prf = 1 / float(re.search(r"<azimuthTimeInterval>([^<]+)<", t).group(1))
                b = _BURST.findall(t)
                ab, ax = [_utc(a) for a, _ in b], [float(v) for _, v in b]
                cnt = sorted({int(c) for c in re.findall(r'<firstValidSample count="(\d+)"', t)})
                if cnt != [lpb]:
                    w.append(f"iw{n} valid-sample lists of {cnt} lines under linesPerBurst {lpb}")
                off = [k for k in range(len(b)) if int((ax[k] - ax[0]) * prf + 0.1) != round((ab[k] - ab[0]) * prf)]
                if off:
                    w.append(f"iw{n} bursts {off} placed off their times")
            cache[s] = f"{Path(s).parent.name}: {', '.join(w)}" if w else ""
        if cache[s]:
            why.append(cache[s])
    return "; ".join(why)


def _burst_check(tr, safe):
    """Why a frame we cut is assembled wrong ("" if not), from its three annotations' burst times (cached in state
    "burst_check"): a subswath whose bursts start or end > 5 s from iw1's (create_frame_tops.csh hands assemble_tops
    shifted arguments when make_s1a_tops fails on a subswath; GMTSAR audit 10-28), or burst times running backwards
    (overlapping non-slice products: T071 2026-08-15). Exact repeats (0 s) are harmless."""
    cache = tr.state.setdefault("burst_check", {})
    if safe in cache:
        return cache[safe]
    bt = {}
    for n in (1, 2, 3):                 # every subswath of the frame, processed or not: the shift shows against iw1
        x = real(Path(safe).glob(f"annotation/*iw{n}*vv*.xml"))
        t = re.findall(r"<burst>\s*<azimuthTime>([^<]+)</azimuthTime>", x[0].read_text(errors="replace")) if x else []
        bt[n] = [datetime.fromisoformat(v.strip()[:26]).timestamp() for v in t]
    why = []
    if not all(bt.values()):
        why.append("no burst list")
    else:
        if any(b - a < -0.1 for v in bt.values() for a, b in zip(v, v[1:])):
            why.append("burst times run backwards")
        far = [n for n in (2, 3) if abs(bt[n][0] - bt[1][0]) > 5 or abs(bt[n][-1] - bt[1][-1]) > 5]
        if far:
            why.append(f"iw{far} bursts span {bt[far[0]][0] - bt[1][0]:+.0f} / {bt[far[0]][-1] - bt[1][-1]:+.0f} s "
                       f"from iw1's")
    cache[safe] = f"{Path(safe).parent.name}: {', '.join(why)}" if why else ""
    return cache[safe]


def bad_stitches(tr):
    """Stitched dates whose frames GMTSAR cut wrong: an image name in both frames, an image outside 0.5-1.5x the median
    of its frame and subswath over the track's stitched dates, or not three images per frame (10-26: T173 2021-01-01's
    iw2 was the whole 4.1-GB image in BOTH frames -- a lab bad date too; align crashed on the repeated name)."""
    size, rec = defaultdict(list), {}
    for d, safes in tr.state.get("stitched", {}).items():
        rec[d] = []
        for s in safes:
            for t in real(Path(s).glob("measurement/*.tiff")):
                k = (Path(s).parent.name, t.name.split("-")[1])
                rec[d].append((k, t.name, t.stat().st_size))
                size[k].append(rec[d][-1][2])
    med = {k: statistics.median(v) for k, v in size.items()}
    # 10-28: frames we cut with the lab's pins (not rescued) must be the super-master's frames, give or take the
    # ANX-second rounding of the names (the lab's own F0518_F0529 / F0518_F0530) -- the T064 pins trap
    ref = sorted({Path(s).parent.name for s in real(glob.glob(str(tr.labd / f"F*_F*/*_{tr.master}T*.SAFE")))})
    num = lambda name: [int(x) for x in re.findall(r"F(\d+)", name)]
    rescued = {d for d, v in tr.state.get("rescued", {}).items() if v.get("ok")}
    out = {}
    for d, r in rec.items():
        why = [f"{k[0]} {k[1]} {b / 1e9:.2f} GB" for k, _, b in r if not 0.5 * med[k] <= b <= 1.5 * med[k]]
        if len({f for _, f, _ in r}) != len(r):
            why.append("an image name in both frames")
        if len(r) != 3 * len(tr.state["stitched"][d]):
            why.append(f"{len(r)} images")
        ours = [Path(s).parent.name for s in tr.state["stitched"][d] if s.startswith(str(tr.outd))]
        off = [f for f in ours if not any(all(abs(a - b) <= 1 for a, b in zip(num(f), num(m))) for m in ref)]
        if ref and d not in rescued and off:
            why.append(f"frames {off}, the super-master's {ref}")
        why += [w for s in tr.state["stitched"][d] if s.startswith(str(tr.outd)) for w in [_burst_check(tr, s)] if w]
        why += [w for w in [_assembly_defects(tr, d, tr.state["stitched"][d])] if w]     # 10-30: the lab's frames too
        if why:
            out[d] = "; ".join(why)
    return out


MASTER_KEYS = ("num_lines", "num_rng_bins", "bytes_per_line", "SC_clock_start", "SC_clock_stop", "PRF", "near_range",
               "rng_samp_rate", "radar_wavelength", "first_line", "first_sample")


def _master_differs(ours, lab, blocks=3, block=1 << 20, image=True):
    """Why our re-made super-master is not the lab's, or "" (same). Same geometry in the PRM (the lab's carries keys
    later steps append, SC_vel and the like), both SLCs hold num_lines x bytes_per_line, and blocks of the image agree
    to <= 0.1 % of their bytes (seasat's build: ~300 bytes per 2.4 GB differ). 10-27: T064's lab SLC is 3 lines shorter
    than ours beyond num_lines -- a file-size check had stopped the align stage -- and EMPTY (a sparse file of zeros):
    image=False then compares the geometry only."""
    p, q = _prm(ours.with_suffix(".PRM")), _prm(lab.with_suffix(".PRM"))

    def num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return v
    bad = [k for k in MASTER_KEYS if num(p.get(k)) != num(q.get(k))]
    if bad:
        return "PRM " + ", ".join(f"{k} {p.get(k)} vs {q.get(k)}" for k in bad)
    size = int(p["num_lines"]) * int(p["bytes_per_line"])
    if not image:
        return "" if ours.stat().st_size >= size else f"our SLC is shorter than num_lines x bytes_per_line ({size})"
    if min(ours.stat().st_size, lab.stat().st_size) < size:
        return f"an SLC shorter than num_lines x bytes_per_line ({ours.stat().st_size}, {lab.stat().st_size} < {size})"
    with open(ours, "rb") as a, open(lab, "rb") as b:
        for k in range(blocks):
            off = (size - block) * (2 * k + 1) // (2 * blocks)
            a.seek(off)
            b.seek(off)
            x, y = a.read(block), b.read(block)
            n = sum(1 for i, j in zip(x, y) if i != j)
            if n > block // 1000:
                return f"{n} of {block} bytes differ at offset {off}"
    return ""


def _ipf(safe):
    """The IPF version that processed a (frame) SAFE, from its manifest.safe; None if unreadable."""
    try:
        m = re.search(r'software name="Sentinel-1 IPF" version="([0-9.]+)"',
                      (Path(safe) / "manifest.safe").read_text(errors="replace"))
    except OSError:
        return None
    return m.group(1) if m else None


def eap_dates(tr):
    """Dates whose scenes an IPF older than EAP_IPF processed (S1A to ~2015-03; the version is per SCENE: T115
    2015-02-22 is 002.43 between 002.36 dates): their phase still carries the elevation antenna pattern (audit A2,
    -76 deg at the F1/F2 seam), so we align them again from the lab's frames with the pattern removed (10-29). From
    the lab's frames and ours; cached in the state."""
    cache = tr.state.setdefault("ipf", {})
    ds = {ymd0(c) for name in lab_pair_dirs(tr) for c in (name[:7], name[8:])} | set(tr.state.get("stitched", {}))
    for d in sorted(ds - set(cache)):
        safes = tr.state.get("stitched", {}).get(d) or real(glob.glob(str(tr.labd / f"F*_F*/*_{d}T*.SAFE")))
        v = [_ipf(s) for s in safes]
        cache[d] = min(x for x in v if x) if any(v) else None
    tr.save()
    return sorted(d for d, v in cache.items() if v and v < EAP_IPF)


def _eap_xml(x, out):
    """An annotation with the AUX_CAL and its SAFE's manifest appended (first lines dropped), as GMTSAR's make_s1a_tops
    expects for the EAP correction (its usage note)."""
    with open(out, "w") as fh:
        fh.write(Path(x).read_text(errors="replace"))
        fh.write("".join(AUX_CAL.read_text(errors="replace").splitlines(keepends=True)[1:]))
        fh.write("".join((Path(x).parent.parent / "manifest.safe").read_text(errors="replace")
                         .splitlines(keepends=True)[1:]))


def nesd_pass(tr):
    """0 without the network spectral diversity, else the pass the track is in: 1 (every date aligned geometry-only,
    its burst overlaps collected) until nesd() solves the network, then 2 (every date aligned again with its shift)."""
    if not NESD:
        return 0
    return 2 if tr.state.get("nesd", {}).get("pass") == 2 else 1


def _esd_dir(tr, n, d=None):
    p = tr.out / "esd" / f"F{n}"
    return p / d if d else p


def _sd_parse(text):
    """spectral_diversity's report (+ the collecting run's esd<n>_shift line) as numbers."""
    out = {}
    for line in text.split("\n"):
        m = re.match(r"\s*(residual_phase|spectral_spectrationXdta|residual_shift|esd\d_shift)\s*=\s*([-+0-9.eE]+)", line)
        if m:
            out[m.group(1)] = float(m.group(2))
        m = re.match(r"\s*Image analyzed (\d+)x(\d+)", line)
        if m:
            out["width"], out["lines"] = int(m.group(1)), int(m.group(2))
    return out


def _esd_collect(tr, n, w, chunk, failed, lock):
    """Pass 1 (10-31): each date's frames from the chunk folder w -- overlap images (.SLCH / .SLCL), boundary file
    (.BB), estimate against the super-master (.sd) -- into OUT/<T>/esd/F<n>/<date>/, with frames.json (stem -> its
    super-master frame, burst offset tmp_da, that estimate) written last; the super-master's frames once."""
    rows = {}
    if (w / "nesd_frames.txt").exists():
        for line in (w / "nesd_frames.txt").read_text().split("\n"):
            p = line.split()
            if len(p) == 3:
                rows[p[0]] = {"master": p[1], "tmp_da": int(p[2])}
    with lock:
        dm = _esd_dir(tr, n, tr.master)
        if not (dm / "frames.json").exists():
            ms = [s.strip() for s in (w / "masterlist").read_text().split("\n") if s.strip()]
            dm.mkdir(parents=True, exist_ok=True)
            for s in ms:
                for ext in ("SLCH", "SLCL", "BB"):
                    shutil.copy(w / f"{s}.{ext}", dm / f"{s}.{ext}")
            (dm / "frames.json").write_text(json.dumps({s: {"master": s, "tmp_da": 0} for s in ms}, indent=1))
    for d in chunk:
        if d in failed:
            continue
        fr = {s: dict(v) for s, v in rows.items() if s.startswith(f"S1_{d}_")}
        if not fr:
            raise RuntimeError(f"{tr.T} F{n} align: no burst overlaps recorded for {d} (nesd_frames.txt)")
        dd = _esd_dir(tr, n, d)
        shutil.rmtree(dd, ignore_errors=True)
        dd.mkdir(parents=True)
        for s, v in fr.items():
            for ext in ("SLCH", "SLCL", "BB", "sd"):
                shutil.move(str(w / f"{s}.{ext}"), str(dd / f"{s}.{ext}"))
            v.update(_sd_parse((dd / f"{s}.sd").read_text()))
        (dd / "frames.json").write_text(json.dumps(fr, indent=1))


def align(tr, dates, batch=12):
    done = tr.state.setdefault("aligned", [])
    # once each: a date can be both ours and a lab pair's partner to align again (T144 10-29: EAP dates twice)
    todo = list(dict.fromkeys(d for d in dates if d in tr.state.get("stitched", {}) and d not in done))
    if tr.master in todo:                     # 10-31: the reference of every chunk, never a date of its own (REALIGN_ALL)
        todo.remove(tr.master)
    bad = bad_stitches(tr)
    for d in sorted(set(bad) & set(todo)):
        log(f"{tr.T} align: SKIP {d}, its stitched frames are cut wrong ({bad[d]})")
    tr.state["bad_stitch"] = bad
    # 10-29: and dates GMTSAR could not align before (state "align_failed": T115 2015-03-06, make_s1a_tops aborts on
    # its iw1 annotation, the lab's copy too -- the lab never aligned it)
    failed_before = tr.state.get("align_failed", {})
    for d in sorted(set(failed_before) & set(todo)):
        log(f"{tr.T} align: SKIP {d}, GMTSAR could not align it ({failed_before[d]})")
    todo = [d for d in todo if d not in bad and d not in failed_before]
    # 11-01: a date whose precise orbit is not on disk is left out by itself (state "no_orbit"), not its whole chunk
    # (T173: six S1A dates of 2025-10..12 whose POEORBs ESA made in 2026-05 -- fetch them on the Mac, prefetch_orbits.py)
    for d in list(todo):
        try:
            _orbit_of(tr, d)
        except RuntimeError as e:
            log(f"{tr.T} align: SKIP {d}, {e} (fetch it on the Mac)")
            tr.state.setdefault("no_orbit", {})[d] = time.strftime("%F")
            todo.remove(d)
        else:
            tr.state.get("no_orbit", {}).pop(d, None)
    npass = nesd_pass(tr)
    if npass == 2:                            # 10-31: pass 2 aligns the network's dates only
        net = set(tr.state["nesd"]["dates"])
        out = [d for d in todo if d not in net]
        if out:
            log(f"{tr.T} align (network pass 2): SKIP {len(out)} dates outside the network {out}")
        todo = [d for d in todo if d in net]
    if not todo:
        return
    msafes = _one_cut_per_frame(real(glob.glob(str(tr.labd / f"F*_F*/*_{tr.master}T*.SAFE"))))
    morb = orbit_for(tr, tr.master, Path(msafes[0]).name[:3])
    chunks = [todo[b:b + batch] for b in range(0, len(todo), batch)]
    # 10-31: ALIGN_PAR chunks at once, each in its own folder (align/c<i>/F<n>); the state is written under a lock, the
    # first error stops new chunks and is raised once the running ones end
    lock, stop, errs = threading.Lock(), threading.Event(), []

    def one(ci, chunk):
        if stop.is_set():
            return
        try:
            _align_chunk(tr, ci, chunk, msafes, morb, npass, lock)
        except BaseException as e:  # noqa: B902 -- re-raised below
            stop.set()
            errs.append(e)
            log(f"{tr.T} align: chunk {chunk[0]}..{chunk[-1]} failed: {str(e)[:300]}")

    with ThreadPoolExecutor(max_workers=max(1, min(ALIGN_PAR, len(chunks)))) as ex:
        list(ex.map(lambda a: one(*a), enumerate(chunks)))
    if errs:
        raise errs[0]


def _align_chunk(tr, ci, chunk, msafes, morb, npass, lock):
    done = tr.state["aligned"]
    log(f"{tr.T} align{f' (network pass {npass})' if npass else ''} {chunk}")
    jobs = []
    for n in tr.subs:
        w = tr.work / f"align/c{ci}/F{n}"
        shutil.rmtree(w, ignore_errors=True)
        w.mkdir(parents=True)
        lines = []
        for d, safes, orb in [(tr.master, msafes, morb)] + [(d, tr.state["stitched"][d], _orbit_of(tr, d))
                                                            for d in chunk]:
            rec = []
            for t, x in _frames_line(safes, n):
                os.symlink(t, w / t.name)
                if d != tr.master and (_ipf(Path(x).parent.parent) or EAP_IPF) < EAP_IPF and (EAP_FIX.is_dir() or USE_FORK):
                    _eap_xml(x, w / x.name)              # 10-29: the antenna pattern removed (make_s1a_tops)
                else:
                    os.symlink(x, w / x.name)
                rec.append(t.stem)
            if not (w / orb.name).exists():
                os.symlink(orb, w / orb.name)
            lines.append(":".join(rec + [orb.name]))
        (w / "data.in").write_text("\n".join(lines) + "\n")
        shutil.copy(tr.labp / "topo/dem.grd", w / "dem.grd")
        # 10-30 afternoon: bin_a1 (patch_eap.py + patch_a1.py) with IGML_A1, the spectral-diversity script with ESD_MODE
        if USE_FORK:                     # 10-31: the IGML fork (A1, A3 by IGML_A1 / IGML_A3)
            env = gmt_env()
        else:
            mk = A1_FIX if REALIGN_ALL and A1_FIX.is_dir() else EAP_FIX
            env = dict(os.environ, PATH=(f"{mk}:" if mk.is_dir() else "") +
                       (f"{STITCH_FIX}:" if STITCH_FIX.is_dir() else "") +
                       f"{TOOLS}:/usr/local/GMTSAR/bin:{os.environ['PATH']}")
            env.pop("IGML_A1", None)
            if REALIGN_ALL:
                env["IGML_A1"] = "1"
        env.pop("IGML_ESD_COLLECT", None)
        env.pop("IGML_ESD_SHIFTS", None)
        if npass:                        # 10-31: network spectral diversity (the fork's preproc_batch_tops_esd.csh)
            if npass == 1:
                env["IGML_ESD_COLLECT"] = "1"
            else:
                env["IGML_ESD_SHIFTS"] = str(_esd_dir(tr, n) / "shifts.txt")
            mode2 = "preproc_batch_tops_esd.csh data.in dem.grd 2 1"
        else:
            mode2 = (f"preproc_batch_tops_esd.csh data.in dem.grd 2 {ESD_MODE}" if ESD_MODE
                     else "preproc_batch_tops.csh data.in dem.grd 2")
        jobs.append(subprocess.Popen(
            f"preproc_batch_tops.csh data.in dem.grd 1 > pbt1.log 2>&1 && {mode2} > pbt2.log 2>&1", shell=True, cwd=w,
            env=env))
    for j in jobs:
        j.wait()
    with lock:
        for n in tr.subs:                   # 10-30: the super-master's cut times (IGML_A1), for qc's burst steps
            cuts = sorted((tr.work / f"align/c{ci}/F{n}").glob("*.cuts"))
            if cuts:
                (tr.out / "a1_cuts" / f"F{n}").mkdir(parents=True, exist_ok=True)
                for c in cuts:
                    shutil.copy(c, tr.out / "a1_cuts" / f"F{n}" / c.name)
    # 10-29: a date whose image GMTSAR did not make in some subswath is set aside, not the whole run (most of a
    # chunk failing is still an error); its logs are kept in OUT/<T>/logs
    failed = {}
    for d in chunk:
        for n in tr.subs:
            w = tr.work / f"align/c{ci}/F{n}"
            miss = [e for e in ("PRM", "LED", "SLC") if not (w / f"S1_{d}_ALL_F{n}.{e}").exists()
                    or (e == "SLC" and (w / f"S1_{d}_ALL_F{n}.SLC").stat().st_size == 0)]
            if miss:
                failed[d] = f"F{n}: S1_{d}_ALL_F{n}.{'/'.join(miss)} not made"
                break
    if len(failed) * 2 > len(chunk):
        raise RuntimeError(f"{tr.T} align: {len(failed)} of {len(chunk)} dates not made: {failed}")
    if failed:
        (tr.out / "logs").mkdir(parents=True, exist_ok=True)
        for n in tr.subs:
            for k in (1, 2):
                if (tr.work / f"align/c{ci}/F{n}/pbt{k}.log").exists():
                    shutil.copy(tr.work / f"align/c{ci}/F{n}/pbt{k}.log",
                                tr.out / f"logs/align_{chunk[0]}_F{n}_pbt{k}.log")
        for d, why in failed.items():
            log(f"{tr.T} align: {d} set aside, {why} (logs {tr.out}/logs/align_{chunk[0]}_F*_pbt*.log)")
    for n in tr.subs:
        w = tr.work / f"align/c{ci}/F{n}"
        m_ours, m_lab = w / f"S1_{tr.master}_ALL_F{n}.SLC", tr.labp / f"F{n}/raw/S1_{tr.master}_ALL_F{n}.SLC"
        why = _master_differs(m_ours, m_lab, image=_holds_data(m_lab))
        if why:
            raise RuntimeError(f"{tr.T} F{n}: the re-made super-master differs from the lab's: {why}")
        dest = tr.outp / f"F{n}/raw"
        dest.mkdir(parents=True, exist_ok=True)
        with lock:
            if not _holds_data(m_lab) and not (dest / m_ours.name).exists():    # the lab emptied its image: keep ours
                for ext in ("PRM", "LED", "SLC"):
                    shutil.copy(m_ours.with_suffix(f".{ext}"), dest / m_ours.with_suffix(f".{ext}").name)
        for d in chunk:
            if d in failed:
                continue
            for ext in ("PRM", "LED", "SLC"):
                f = w / f"S1_{d}_ALL_F{n}.{ext}"
                shutil.move(str(f), str(dest / f.name))     # 11-01: a rename on one filesystem (was a 2.8 GB copy)
        if npass == 1:
            _esd_collect(tr, n, w, chunk, failed, lock)
    with lock:
        if failed:
            tr.state.setdefault("align_failed", {}).update(failed)
        done += [d for d in chunk if d not in failed]
        if STITCH_FIX.is_dir():
            tr.state.setdefault("stitch_fixed", []).extend(d for d in chunk if d not in tr.state["stitch_fixed"])
        # 10-29: aligned without the antenna pattern from the start (a new track's dates; 10-31: every date the fork
        # aligns, so an interrupted run does not align them once more)
        if (tr.new and EAP_FIX.is_dir()) or USE_FORK:
            eap = [d for d in chunk if d not in failed
                   and min(_ipf(s) or EAP_IPF for s in tr.state["stitched"][d]) < EAP_IPF]
            tr.state["eap_done"] = sorted(set(tr.state.get("eap_done", [])) | set(eap))
        tr.save()
    for n in tr.subs:
        shutil.rmtree(tr.work / f"align/c{ci}/F{n}", ignore_errors=True)


def _orbit_of(tr, d):
    """A date's precise orbit: the plan's, or (the lab's dates we re-align) the lab's data_<T> / /palsar/IGML/orbits."""
    if d in tr.state.get("plan", {}) and tr.state["plan"][d].get("orbit"):
        return Path(tr.state["plan"][d]["orbit"])
    f = orbit_for(tr, d, Path(tr.state["stitched"][d][0]).name[:3])
    if f is None:
        raise RuntimeError(f"{tr.T} {d}: no precise orbit")
    return f


QC_BLOCKS, QC_ROWS, QC_SHIFT = 12, 50, 25   # blocks per subswath (2-98 % of the lines), rows per block, noise-floor shift
QC_COVER, QC_MARGIN = 0.90, 0.10            # ours must have data in >= 90 % of the neighbour's windows; coherence p90
                                            # >= the noise floor + 0.10 (pilot 10-27: good dates 0.39-0.96, floor ~0.18)
# 10-28: and no block may be MISREGISTERED. T173 2014-12-04 (the lab's aligned image, never merged by the lab) passed
# the whole-subswath test (p90 0.73) though its last slice is ~1 line off: coherence 0.20 vs 0.41 in the neighbours'
# pair there, a steady -0.04 rad / line azimuth phase ramp (39 cycles) in all its interferograms, and per block its
# coherence peaks with the neighbour shifted by -1 line (0.71 vs 0.48 unshifted). Lost coherence alone is NOT the test:
# 2026-10-03 and 2025-11-30 keep 0.2-0.6 of the neighbours' excess coherence over whole subswaths (weather) and peak
# at 0 in every block. A block fails if, against every neighbour, a shift of QC_SHIFTS lines beats 0 by QC_SHIFT_GAIN.
QC_VERSION, QC_SHIFTS, QC_SHIFT_GAIN = 2, (-2, -1, 1, 2), 0.03
QC_PRM = 2      # a date that failed before is judged once more: 10-30 SLCs read through prm_src; 11-06 qc's neighbours,
                # water blocks
QC_NEIGHBOURS = 4   # 11-06: a date is judged against the nearest date each side; if it fails, against the nearest 4


def _prm(path):
    return {a.strip(): b.strip() for a, s, b in (line.partition("=") for line in open(path)) if s}


def _blocks(slc, prm=None):
    """QC_BLOCKS blocks of QC_ROWS lines spread over an aligned SLC (GMTSAR: int16 I, Q per pixel), as complex64; its
    size from ``prm`` (prm_src) or the SLC's own PRM."""
    import numpy as np
    p = _prm(prm or slc.with_suffix(".PRM"))
    w, nl = int(p["bytes_per_line"]) // 4, int(p["num_lines"])
    out = []
    with open(slc, "rb") as fh:
        for r in np.linspace(0.02 * nl, 0.98 * nl - QC_ROWS, QC_BLOCKS).astype(int):
            fh.seek(int(r) * w * 4)
            a = np.fromfile(fh, dtype="<i2", count=QC_ROWS * w * 2).astype(np.float32)
            a = a[: a.size // (2 * w) * 2 * w].reshape(-1, w, 2)
            out.append(a[..., 0] + 1j * a[..., 1])
    return out


def _windows(a, b, ly=5, lx=20):
    """Per 5 x 20-pixel window: has a data, has b data, coherence |sum a b*| / sqrt(sum|a|^2 sum|b|^2) (NaN if either
    is empty)."""
    import numpy as np
    n, m = min(len(a), len(b)) // ly * ly, min(a.shape[1], b.shape[1]) // lx * lx
    a = a[:n, :m].reshape(n // ly, ly, m // lx, lx)
    b = b[:n, :m].reshape(n // ly, ly, m // lx, lx)
    pa, pb = (np.abs(a) ** 2).sum(axis=(1, 3)), (np.abs(b) ** 2).sum(axis=(1, 3))
    both = (pa > 0) & (pb > 0)
    c = np.full(pa.shape, np.nan)
    c[both] = np.abs((a * b.conj()).sum(axis=(1, 3)))[both] / np.sqrt(pa[both] * pb[both])
    return pa > 0, pb > 0, c


def _shifted(a, b, s):
    """a's line k + s against b's line k."""
    return (a[s:], b[:-s]) if s > 0 else (a[:s], b[-s:])


def _qc_pair(ours, theirs, shifts=False):
    """(coverage, coherence p90, noise-floor p90): coverage = share of the neighbour's windows with data where ours
    has data too; the floor = the same coherence with the neighbour shifted QC_SHIFT lines (speckle decorrelated).
    shifts=True adds, per block, (p90 unshifted, the best p90 over the QC_SHIFTS line shifts, that shift)."""
    import numpy as np
    va_vb = vb_n = 0
    coh, flo, blk = [], [], []
    p90 = lambda x: float(np.percentile(x, 90)) if x.size else 0.0
    fin = lambda x: x[np.isfinite(x)]
    for a, b in zip(ours, theirs):
        va, vb, c = _windows(a, b)
        va_vb += int((va & vb).sum())
        vb_n += int(vb.sum())
        coh.append(fin(c))
        flo.append(fin(_windows(a[:-QC_SHIFT], b[QC_SHIFT:])[2]))
        if shifts:
            blk.append((p90(coh[-1]),) + max((p90(fin(_windows(*_shifted(a, b, s))[2])), s) for s in QC_SHIFTS))
    coh, flo = np.concatenate(coh), np.concatenate(flo)
    out = (va_vb / vb_n if vb_n else 0.0), p90(coh), p90(flo)
    return out + (blk,) if shifts else out


QC_CONTENT_MIN, QC_CONTENT_AGREE, QC_CONTENT_COH, QC_CONTENT_SHARE = 0.3, 0.6, 0.5, 0.01


def _amp_corr(a, b, ly=5, lx=20):
    """Correlation of two blocks' amplitude patterns (5 x 20 looks, where both hold data): ~0.9 for the same ground,
    whatever its coherence; ~0 when one image holds other content (10-30: the qc's coherence p90 over a whole subswath
    hid a quarter of a subswath read at the wrong stride)."""
    import numpy as np
    n, m = min(len(a), len(b)) // ly * ly, min(a.shape[1], b.shape[1]) // lx * lx
    if n == 0 or m == 0:
        return float("nan")
    A = np.abs(a[:n, :m]).reshape(n // ly, ly, m // lx, lx).mean(axis=(1, 3))
    B = np.abs(b[:n, :m]).reshape(n // ly, ly, m // lx, lx).mean(axis=(1, 3))
    ok = (A > 0) & (B > 0)
    return float(np.corrcoef(A[ok], B[ok])[0, 1]) if ok.sum() > 50 else float("nan")


def _master_cut_lines(tr, n):
    """The super-master's burst-switch lines in its stitched image, from its frames' cut times (a1_cuts, align)."""
    cuts = [float(v) for f in sorted((tr.out / "a1_cuts" / f"F{n}").glob("*.cuts")) for v in f.read_text().split()]
    if not cuts:
        return []
    p = _prm(prm_src(tr, tr.master, n))
    t0, prf = float(p["clock_start"]), float(p["PRF"])
    return sorted({int(round((t - t0) * 86400.0 * prf)) for t in cuts})


def _burst_steps(tr, d, x, n, k=4, h=40, wd=4096):
    """The interferometric phase step (deg) of d x x across up to k of the super-master's burst-switch lines (audit A5:
    a residual azimuth misregistration steps the phase at every burst boundary; spectral diversity should leave ~0):
    80 lines x 4096 samples at mid-swath, 1 x 8 looks, the window's dominant fringe removed, the mean phase per line
    either side of the line (+-2 lines left out) with one linear gradient. [] without the cut times."""
    import numpy as np
    lines = _master_cut_lines(tr, n)
    p = _prm(prm_src(tr, d, n))
    nl, w = int(p["num_lines"]), int(p["bytes_per_line"]) // 4
    lines = [c for c in lines if h + 2 < c < nl - h - 2]
    if not lines:
        return []

    def rd(slc, r0, xm):
        with open(slc, "rb") as fh:
            fh.seek(r0 * w * 4)
            a = np.fromfile(fh, dtype="<i2", count=2 * h * w * 2).astype(np.float32).reshape(2 * h, w, 2)
        return a[:, xm:xm + wd, 0] + 1j * a[:, xm:xm + wd, 1]

    out = []
    xm = max(0, w // 2 - wd // 2)
    for c in sorted({lines[int(i)] for i in np.linspace(0, len(lines) - 1, min(k, len(lines)))}):
        ifg = rd(slc_src(tr, d, n), c - h, xm) * np.conj(rd(slc_src(tr, x, n), c - h, xm))
        z = ifg.reshape(2 * h, -1, 8).sum(axis=2)
        F = np.fft.fft2(z)
        ky, kx = np.unravel_index(np.argmax(np.abs(F)), F.shape)
        yy, xx = np.mgrid[0:z.shape[0], 0:z.shape[1]]
        z = z * np.exp(-2j * np.pi * (ky * yy / z.shape[0] + kx * xx / z.shape[1]))
        m = np.exp(1j * np.angle(z)).mean(axis=1)
        L = np.arange(2 * h) + c - h
        sel = np.abs(L - c + 0.5) > 2.5
        q, l = np.unwrap(np.angle(m[sel])), L[sel]
        A = np.c_[np.ones(l.size), l - c, (l >= c).astype(float)]
        out.append(round(float(np.degrees(np.linalg.lstsq(A, q, rcond=None)[0][2])), 1))
    return out


def _content_blocks(mine, nbs):
    """Blocks whose amplitude pattern matches none of the neighbours (< QC_CONTENT_MIN) where the neighbours match each
    other (the median over their pairs >= QC_CONTENT_AGREE) and no neighbour finds it coherent (> QC_CONTENT_COH in
    QC_CONTENT_SHARE of its windows): content that is not this ground. Needs two neighbours.
    11-06: open water matches nothing only on average. Over 50 lines its amplitude is one range trend, so two dates correlate
    anywhere in -0.6..+0.6 (T071 F2 block 9, 5-95 %; land 0.82-0.97), and one neighbour pair reaching 0.6 set 24 dates
    aside for nothing. With two neighbours the agreement is that one pair: qc then asks again with four."""
    import numpy as np
    if len(nbs) < 2:
        return []
    out = []
    for k, a in enumerate(mine):
        agree = [_amp_corr(x[k], y[k]) for i, x in enumerate(nbs) for y in nbs[i + 1:]]
        if not np.nanmedian(agree) >= QC_CONTENT_AGREE or not np.nanmax([_amp_corr(a, x[k]) for x in nbs]) < QC_CONTENT_MIN:
            continue
        coh = [_windows(a, x[k])[2] for x in nbs]
        if not any(np.mean(c[np.isfinite(c)] > QC_CONTENT_COH) >= QC_CONTENT_SHARE for c in coh if np.isfinite(c).any()):
            out.append(k)
    return out


def _qc_judge(tr, d, nb, need):
    """Date d against the neighbours nb, per subswath: (the record, ok). ok = coverage >= need and coherence p90 >= floor +
    QC_MARGIN with at least one of them, no misregistered block, no block of other content. A block is misregistered
    if a shift of QC_SHIFTS lines beats 0 by QC_SHIFT_GAIN against every neighbour that sees its ground (11-06: best p90
    >= floor + QC_MARGIN; open water has no say, nor has a neighbour that decorrelated, which used to cancel the test)."""
    rec, ok = {"v": QC_VERSION, "pv": QC_PRM}, bool(nb)
    for n in tr.subs:
        mine = _blocks(slc_src(tr, d, n), prm_src(tr, d, n))
        best, votes, theirs = None, {}, []
        for x in nb:
            theirs.append(_blocks(slc_src(tr, x, n), prm_src(tr, x, n)))
            cov, coh, flo, blk = _qc_pair(mine, theirs[-1], shifts=True)
            good = cov >= need and coh >= flo + QC_MARGIN
            r = {"vs": x, "cover": round(cov, 3), "coh": round(coh, 3), "floor": round(flo, 3), "ok": good}
            if best is None or (good, coh) > (best["ok"], best["coh"]):
                best = r
            for k, (c0, cs, s) in enumerate(blk):
                if max(c0, cs) >= flo + QC_MARGIN:
                    votes.setdefault(k, []).append(s if cs > c0 + QC_SHIFT_GAIN else None)
        off = {k: v[-1] for k, v in votes.items() if all(s is not None for s in v)}
        if best and off:
            best.update(ok=False, shifted={str(k): s for k, s in sorted(off.items())})
        bad_c = _content_blocks(mine, theirs)          # 10-30: content that is not this ground
        if best and bad_c:
            best.update(ok=False, content=bad_c)
        if best:                                       # 10-30: recorded, not judged yet (A5)
            try:
                best["steps"] = _burst_steps(tr, d, best["vs"], n)
            except Exception as e:  # never stops the qc
                best["steps_err"] = str(e)[:120]
        rec[f"F{n}"] = best
        ok = ok and bool(best and best["ok"])
    rec["ok"] = ok
    return rec, ok


def qc(tr):
    """Check each new date (ours aligned, or the lab's aligned but never merged), in time order, against the nearest
    aligned date before and after it -- the lab's merged ones, new ones that passed, or ones still to judge -- per
    subswath (_qc_judge); if it fails, against the nearest QC_NEIGHBOURS. 11-06: dates still to judge count. With every
    date re-aligned (REALIGN_ALL) the only date "after" that had passed was a re-cut one years away, so a date had one
    near neighbour: T071 2022-09-22 failed beside a rainy 09-10 (0.30; 12 days on 0.48-0.57), the first date against 2019,
    the S1C pair 2025-05-03 / 05-15 (0.91 with each other, ~0.30 with S1A) against S1A only.
    A failing date goes to state "bad_align" and is never paired; dates that passed an older version are checked again
    (reject() then withdraws their stored pairs if they fail)."""
    res = tr.state.setdefault("qc", {})
    bad = tr.state.setdefault("bad_align", {})
    part = _kept(tr)                         # 10-30: a partial frame needs QC_COVER x its share (its record: "need")
    todo = [d for d in sorted(set(tr.state.get("aligned", [])) | set(lab_aligned_planned(tr)))
            if d not in res or (res[d].get("ok") and res[d].get("v", 1) < QC_VERSION)
            or (not res[d].get("ok") and d in part and "need" not in res[d])
            or (not res[d].get("ok") and res[d].get("pv", 0) < QC_PRM)]   # 10-30: failures judged before prm_src
    todo.sort(key=lambda d: d not in tr.state.get("recut", {}))   # re-cut dates first: neighbours recut() reopened
                                                                  # are then judged next to them (stable: time order)
    if not todo:
        return
    lab = sorted(d for d in tr.lab_dates() if tr.lab_aligned(d))
    pool = sorted(set(tr.state.get("aligned", [])) | {tr.master})    # T064: no lab image left -- our own aligned dates
    held = lambda x: slc_src(tr, x, tr.subs[0]).exists() and _holds_data(slc_src(tr, x, tr.subs[0]))
    pending = set(todo)
    for d in todo:
        pending.discard(d)
        ref = sorted(set(lab) | {x for x, v in res.items() if v.get("ok") and x != d} | pending)
        before, after = [x for x in ref if x < d], [x for x in ref if x > d]
        if not before:
            before = [x for x in pool if x < d and held(x)]
        if not after:
            after = [x for x in pool if x > d and held(x)]
        near = sorted(before[-3:] + after[:3], key=lambda x: abs((day(x) - day(d)).days))
        nb = before[-1:] + after[:1]
        nb += [x for x in near if x not in nb][:2 - len(nb)]
        need = QC_COVER * part.get(d, 1)
        rec, ok = _qc_judge(tr, d, nb, need)
        spare = [x for x in near if x not in nb][:QC_NEIGHBOURS - len(nb)]
        if not ok and spare:
            rec, ok = _qc_judge(tr, d, nb + spare, need)
        if d in part:
            rec["need"] = round(need, 3)
        res[d] = rec
        if ok:
            bad.pop(d, None)
        else:
            bad[d] = "; ".join(f"{k} vs {v['vs']}: cover {v['cover']}, coh {v['coh']} (floor {v['floor']})" +
                               (f", misregistered blocks {v['shifted']} (block: lines)" if "shifted" in v else "") +
                               (f", content off in blocks {v['content']}" if "content" in v else "")
                               for k, v in rec.items() if k.startswith("F") and v and not v["ok"]) or "no neighbour"
        tr.save()
        log(f"{tr.T} qc {d}: {'ok' if ok else 'FAILED'} " +
            " ".join(f"{k} {v['cover']:.2f}/{v['coh']:.2f}/{v['floor']:.2f}" +
                     (f" shifted {v['shifted']}" if "shifted" in v else "") +
                     (f" content off {v['content']}" if "content" in v else "")
                     for k, v in rec.items() if k.startswith("F") and v))
    log(f"{tr.T} qc: {sum(v['ok'] for v in res.values())} ok, {len(bad)} set aside {sorted(bad)}")


def _kept(tr):
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


def retry_failed(tr):
    """Dates qc set aside that another alignment can cure, once each (10-28): our alignments that were misregistered
    before the stitch fix (their files -> OUT/<T>/rejected/aligned/, aligned again), and the lab's aligned-only images
    that failed (plan "lab_aligned_bad": stitched and aligned by us instead). main() then stitches / aligns / checks
    them again in the same run."""
    if not STITCH_FIX.is_dir():
        return []
    bad, q, plan = tr.state.get("bad_align", {}), tr.state.get("qc", {}), tr.state.get("plan", {})
    aligned, fixed = tr.state.setdefault("aligned", []), set(tr.state.get("stitch_fixed", []))
    tried = tr.state.setdefault("retried", [])
    out = []
    for d in sorted(bad):
        if d in tried:
            continue
        shifted = any("shifted" in (v or {}) for k, v in q.get(d, {}).items() if k.startswith("F"))
        if d in aligned and shifted and d not in fixed:
            rj = tr.out / "rejected" / "aligned"
            rj.mkdir(parents=True, exist_ok=True)
            for n in tr.subs:
                for ext in ("SLC", "PRM", "LED"):
                    f = tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.{ext}"
                    if f.exists():
                        shutil.move(str(f), str(rj / f.name))
            aligned.remove(d)
        elif d not in aligned and plan.get(d, {}).get("lab_aligned") and not plan[d].get("lab_aligned_bad"):
            plan[d]["lab_aligned_bad"] = True
        else:
            continue
        q.pop(d, None)
        bad.pop(d, None)
        tried.append(d)
        out.append(d)
    if out:
        tr.save()
        log(f"{tr.T}: {len(out)} dates qc set aside go through align again (stitch fix / our own alignment): {out}")
    return out


def lab_aligned_planned(tr):
    """Dates the lab has aligned (product_<T>/F1/raw) but not merged: theirs are used, not re-aligned (if they still hold
    data: the lab emptied T064's, T137's and T042's 2025 ones)."""
    return [d for d, v in tr.state.get("plan", {}).items()
            if v.get("lab_aligned") and not v.get("lab_aligned_bad") and tr.lab_aligned(d)]


def realign_dates(tr):
    """The lab's dates (complete pairs in its stack) whose aligned image it emptied, within MAXDT days of a date we add
    (or the super-master pair's): we align them again from the lab's own stitched frames, as partners for the new pairs
    only (10-27, user: "we can re-align ourselves"; T064: all 295 emptied)."""
    new = [d for d, v in tr.state.get("plan", {}).items() if v.get("ready")]
    near = lambda d: REPROCESS_LAB or any(abs((day(d) - day(x)).days) <= MAXDT for x in new) or \
        d in tr.master_pair or (ANNUAL_TOL >= 0 and d[4:6] in SUMMER)
    labds = {ymd0(c) for name in lab_pair_dirs(tr) for c in (name[:7], name[8:])} if REPROCESS_LAB else set()
    out = []
    for d in sorted(set(tr.lab_dates()) | labds):            # 10-29: every date of the lab's pairs (T064: all 295)
        if d not in tr.state.get("plan", {}) and not tr.lab_aligned(d) and near(d) \
                and real(glob.glob(str(tr.labd / f"F*_F*/*_{d}T*.SAFE"))):
            out.append(d)
    return out


def good_aligned(tr):
    """The new dates (our plan's) that passed qc, ours and the lab's aligned-only ones (a date qc has not looked at is not
    used). The lab's dates we re-aligned are partners only: pairs_for."""
    q, plan = tr.state.get("qc", {}), tr.state.get("plan", {})
    return [d for d in sorted(set(tr.state.get("aligned", [])) | set(lab_aligned_planned(tr)))
            if d in plan and q.get(d, {}).get("ok")]


def pairs_for(tr, new):
    """Every pair <= MAXDT days with at least one new date; partners: the lab's dates with a usable image, and those we
    re-aligned that passed qc."""
    q, plan = tr.state.get("qc", {}), tr.state.get("plan", {})
    have = set(new) | {d for d in tr.lab_dates() if tr.lab_aligned(d)} | \
        {d for d in tr.state.get("aligned", []) if d not in plan and q.get(d, {}).get("ok")}
    ds = sorted(have)
    out = [(a, b) for i, a in enumerate(ds) for b in ds[i + 1:]
           if (day(b) - day(a)).days <= MAXDT and (a in new or b in new)]
    return out


def annual_pairs(tr):
    """The year-long summer pairs (ANNUAL_TOL): June-August dates of consecutive years, 365 +- ANNUAL_TOL days apart,
    among the dates with a usable aligned image -- the lab's (complete pairs and an image holding data), and ours or the
    lab's aligned-only ones that passed qc. ANNUAL_TOL=-1 turns them off."""
    if ANNUAL_TOL < 0:
        return []
    q = tr.state.get("qc", {})
    have = {d for d in tr.lab_dates() if tr.lab_aligned(d)} | \
        {d for d in set(tr.state.get("aligned", [])) | set(lab_aligned_planned(tr)) if q.get(d, {}).get("ok")}
    sm = sorted(d for d in have if d[4:6] in SUMMER)
    return [(a, b) for i, a in enumerate(sm) for b in sm[i + 1:] if abs((day(b) - day(a)).days - 365) <= ANNUAL_TOL]


def same_orbit_trap(tr, a, b):
    """SAT_baseline.c refines the baseline only when the two images' clock_start, in seconds within the year, differ by
    more than 10 s; otherwise B_offset stays -99999.9999 and the topographic phase is garbage (GMTSAR audit 10-28).
    Two dates on the same day of year within 10 s of each other -- possible for a 366-day pair across 29 February."""
    try:
        t = [float(_prm(prm_src(tr, d, tr.subs[0]))["SC_clock_start"]) for d in (a, b)]
    except (OSError, KeyError, ValueError):
        return False
    return int(t[0]) % 1000 == int(t[1]) % 1000 and abs(t[0] % 1 - t[1] % 1) * 86400 <= 10


def lab_pair_dirs(tr):
    """Names (yyyyddd_yyyyddd) of the lab's merged pair folders, complete or not (T042's 2025-26 ones are all-NaN)."""
    return sorted(d.name for d in (tr.labp / "merge").iterdir() if re.match(r"^\d{7}_\d{7}$", d.name))


def lab_pairs(tr):
    """The lab's merged pairs, to be made again with PROC and to replace the lab's (10-29), where both dates have an
    image we can use: the lab's (holding data) or ours (aligned, passed qc; the EAP dates are ours once re-aligned)."""
    q = tr.state.get("qc", {})
    ours = {d for d in tr.state.get("aligned", []) if q.get(d, {}).get("ok")}
    out, left = [], 0
    for name in lab_pair_dirs(tr):
        a, b = ymd0(name[:7]), ymd0(name[8:])
        # 11-05: the super-master is never re-aligned (it is the reference; its image is the lab's), so with REALIGN_ALL
        # lab_aligned() said no to it and the lab's pairs ENDING on it were never made (T173: three, 2016-07..09 x 10-06)
        if all(d == tr.master or d in ours or tr.lab_aligned(d) for d in (a, b)):
            out.append((a, b))
        else:
            left += 1
    if left:
        log(f"{tr.T}: {left} of the lab's pairs have a date without a usable image yet")
    return out


def all_pairs(tr):
    """Every pair this track gets from us: <= MAXDT days with a new date, and the year-long summer pairs (but none
    SAT_baseline would take for one orbit); since 10-29 also every pair of the lab's (REPROCESS_LAB)."""
    prs = set(pairs_for(tr, good_aligned(tr))) | set(annual_pairs(tr))
    if REPROCESS_LAB:
        prs |= set(lab_pairs(tr))
    trap = {p for p in prs if doy0(p[0])[4:] == doy0(p[1])[4:] and same_orbit_trap(tr, *p)}
    if trap:
        log(f"{tr.T}: {len(trap)} pairs left out, same day of year within 10 s (SAT_baseline): {sorted(trap)}")
    return sorted(prs - trap, key=lambda p: (p[1], p[0]))


def slc_src(tr, d, n):
    ours = tr.outp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"
    return ours if ours.exists() else tr.labp / f"F{n}/raw/S1_{d}_ALL_F{n}.SLC"


def _prm_whole(prm, slc, fit=True):
    """A PRM that describes its SLC: bytes_per_line x num_lines = the SLC's size and (fit) the second rshift ... block
    that preproc_batch_tops.csh mode 2 appends when it aligns a date."""
    try:
        txt = prm.read_text()
        p = _prm(prm)
        size = int(p["bytes_per_line"]) * int(p["num_lines"]) == slc.stat().st_size
    except (OSError, KeyError, ValueError):
        return False
    return size and (not fit or len(re.findall(r"^rshift\s*=", txt, re.M)) > 1)


_PRM_SRC = {}


def prm_src(tr, d, n):
    """The PRM of slc_src(tr, d, n): its own, unless that no longer describes the lab's aligned SLC -- then the aligned
    copy the lab's own pairs kept (F<n>/intf_all/<pair>/, used only if every copy that fits is the same). 10-30: on
    T071 the lab's PRMs of 66 dates (F3 of 2017-01-21 and of every date 2019-08-15 .. 2021-12-14, F2 too from
    2021-08-04) are preproc's raw ones, written over the aligned ones after the lab made its pairs: no alignment fit,
    the scene's own width and line count (24156 x 28300 for a 24420 x 28296 SLC). GMTSAR's phasediff refuses those
    pairs, and the qc read the SLC at the wrong stride (2019-11-13 set aside: both its neighbours were such dates).
    The date's own PRM when no copy fits."""
    key = (tr.T, d, n)
    if key not in _PRM_SRC:
        slc = slc_src(tr, d, n)
        own = out = slc.with_suffix(".PRM")
        if not _prm_whole(own, slc, fit=d != tr.master) and str(slc).startswith(f"{tr.labp}/"):
            good = [c for c in sorted(tr.labp.glob(f"F{n}/intf_all/*{doy0(d)}*/{own.name}")) if _prm_whole(c, slc)]
            texts = {c.read_text() for c in good}
            if len(texts) == 1:
                out = good[0]
                log(f"{tr.T} F{n} {d}: the lab's raw/{own.name} does not describe its SLC -- using the aligned copy "
                    f"its pairs kept ({out.parent.name}, {len(good)} identical copies)")
            else:
                log(f"{tr.T} F{n} {d}: the lab's raw/{own.name} does not describe its SLC and "
                    f"{'no pair kept a copy that does' if not texts else f'its pairs kept {len(texts)} versions'}")
        _PRM_SRC[key] = out
    return _PRM_SRC[key]


def _copy_prm(s, t):
    """A PRM into the work tree when its content differs (a fixed PRM can have the stale copy's size)."""
    if not t.exists() or t.read_bytes() != s.read_bytes():
        shutil.copy(s, t)


def _copy_whole(s, t, link=False):
    """t = s in the work tree: a hard link where `link` (our own images in OUT, one filesystem with the work tree; 11-01),
    else a copy through t.part + rename, so that a run killed mid-copy leaves no short file a later run would take as
    done (10-28). An existing t is reused only if it IS s (same inode) or a copy of this very s (same size and mtime:
    copy2). 11-01: the size alone had kept proc1's copies of every date after pass 2 aligned them again, so T173's first
    pass-2 interferograms were made from the old images (stopped, set aside)."""
    ss = s.stat()
    if t.exists():
        ts = t.stat()
        if (ts.st_dev, ts.st_ino) == (ss.st_dev, ss.st_ino) or \
                (ts.st_size == ss.st_size and abs(ts.st_mtime - ss.st_mtime) < 1.0):
            return
        t.unlink()
    if link:
        try:
            os.link(s, t)
            return
        except OSError:
            pass
    tmp = t.with_name(t.name + ".part")
    shutil.copy2(s, tmp)
    os.replace(tmp, t)


def intf(tr, ncores=6):
    prs = [p for p in all_pairs(tr) if f"{doy0(p[0])}_{doy0(p[1])}" not in tr.state.get("stored", [])]
    if not prs:
        log(f"{tr.T} intf: nothing new")
        return []
    prs = [tr.master_pair] + [p for p in prs if tuple(p) != tuple(tr.master_pair)]   # a pair twice: parallel clobbers
    w = tr.work / "intf"
    log(f"{tr.T} intf: {len(prs) - 1} pairs (+ the super-master pair) in {w}")
    jobs = []
    for n in tr.subs:
        f = w / f"F{n}"
        (f / "raw").mkdir(parents=True, exist_ok=True)
        (f / "topo").mkdir(parents=True, exist_ok=True)
        for d in sorted({x for p in prs for x in p}):
            src = slc_src(tr, d, n)
            for ext in ("SLC", "LED"):
                s, t = src.with_suffix("." + ext), f / "raw" / src.with_suffix("." + ext).name
                _copy_whole(s, t, link=str(src).startswith(f"{tr.out}/"))  # local (11-01: ours hard-linked)
            _copy_prm(prm_src(tr, d, n), f / "raw" / src.with_suffix(".PRM").name)
        for name, src in (("topo/topo_ra.grd", tr.labp / f"F{n}/topo/topo_ra.grd"),
                          ("batch_tops.config", tr.labp / f"F{n}/batch_tops.config")):
            if not (f / name).exists():
                shutil.copy(src, f / name)
        (f / "intf.in").write_text("".join(f"S1_{a}_ALL_F{n}:S1_{b}_ALL_F{n}\n" for a, b in prs))
        jobs.append(subprocess.Popen(f"intf_tops_parallel.csh intf.in batch_tops.config {ncores} > itp.log 2>&1",
                                     shell=True, cwd=f,
                                     env=proc_env(tr, n)))
    for j in jobs:
        j.wait()
    tr.state["intf_pairs"] = [list(p) for p in prs]
    tr.save()
    return prs


def start_realign(tr):
    """Every date aligned again with PROC's GMTSAR (10-30 afternoon, REALIGN_ALL): once per track and version, our
    aligned images move to OUT/<T>/superseded_<old>/raw/F<n>/ (nothing deleted), the alignment record is cleared (aligned,
    qc, bad_align, align_failed, retried, stitch_fixed, eap_done), and the partial frames cut since 10-27 (rescue,
    datatake_rescue keeping < the whole pin line) are withdrawn to OUT/<T>/rejected/partial_frames (user: "We only want
    complete interferograms, not partial ones"); their records stay, marked, so they are not cut again. The frames
    themselves (the lab's, ours, the re-cut ones) are kept. start_proc then sets the old pairs aside."""
    if tr.state.get("realign") == PROC:
        return
    old = tr.state.get("proc", "proc0")
    moved = 0
    for n in tr.subs:
        src = tr.outp / f"F{n}/raw"
        if src.is_dir():
            for f in sorted(src.glob("S1_*_ALL_F*.*")):
                if f.suffix in (".SLC", ".PRM", ".LED"):
                    _set_aside(f, tr.out / f"superseded_{old}/raw/F{n}")
                    moved += 1
    st = tr.state.get("stitched", {})
    partial = sorted({d for k in ("rescued", "datatake") for d, v in tr.state.get(k, {}).items()
                      if v.get("ok") and isinstance(v.get("kept"), (int, float)) and v["kept"] < 1})
    for d in partial:
        for s in st.pop(d, []):
            if str(s).startswith(f"{tr.outd}/") and Path(s).exists():
                _set_aside(Path(s), tr.out / "rejected/partial_frames")
        for k in ("rescued", "datatake"):
            if d in tr.state.get(k, {}):
                tr.state[k][d].update(ok=False, v=3, why="partial frame: complete frames only since 10-30")
    for k in ("aligned", "retried", "stitch_fixed", "eap_done"):
        tr.state[k] = []
    for k in ("qc", "bad_align", "align_failed"):
        tr.state[k] = {}
    tr.state["realign"] = PROC
    tr.save()
    log(f"{tr.T}: every date aligned again ({PROC}): {moved} files of our {old} alignment set aside "
        f"(OUT/{tr.T}/superseded_{old}/raw), {len(partial)} partial dates withdrawn {partial}")


def _sd_pair(tr, n, a, da, b, db, filt):
    """GMTSAR's spectral_diversity between frames a and b (stems; their folders da, db) aligned to one super-master
    frame, in a folder of its own: the double-difference phase of their burst overlaps as b's azimuth shift relative
    to a [lines] -- the median of ddphase / (2 pi spec_sep), as preproc_batch_tops_esd.csh's esd_mode 1 (GMTSAR's
    conversion: nesd() multiplies by linesPerBurst / burst cycle, audit A6), and the
    coherent sum (residual_shift); sigma from the phases' spread (MAD / sqrt(n / 25): ddphase rows are 10 lines apart,
    a 25-sample filter along range). ddphase is deleted."""
    import numpy as np
    wd = tr.work / "nesd" / f"F{n}" / f"{a}__{b}"
    shutil.rmtree(wd, ignore_errors=True)
    wd.mkdir(parents=True)
    for stem, src in ((a, da), (b, db)):
        for ext in ("SLCH", "SLCL", "BB"):
            os.symlink(src / f"{stem}.{ext}", wd / f"{stem}.{ext}")
    r = subprocess.run(["spectral_diversity", a, b, "0", filt], cwd=wd, env=gmt_env(), capture_output=True, text=True)
    rec = {"a": a, "b": b, **_sd_parse(r.stdout)}
    dd, sep = wd / "ddphase", rec.get("spectral_spectrationXdta")
    rec["n"] = 0
    if sep and dd.exists() and dd.stat().st_size > 0:
        x = np.loadtxt(dd, ndmin=2)
        if len(x):
            ph = x[:, 2]
            med = float(np.median(ph))
            mad = 1.4826 * float(np.median(np.abs(ph - med)))
            rec.update(n=int(len(ph)), corr=round(float(np.mean(x[:, 3])), 4), phase_med=med, phase_mad=mad,
                       shift=med / (2 * np.pi * sep),
                       sigma=max(mad / np.sqrt(max(len(ph) / 25.0, 1.0)), NESD_MIN_PHASE) / (2 * np.pi * sep))
    shutil.rmtree(wd, ignore_errors=True)
    return rec


def _burst_cycle(xml):
    """(linesPerBurst, burst cycle in lines) of a TOPS annotation: the cycle from consecutive bursts' azimuthTime over
    azimuthTimeInterval. 10-31, audit A6: GMTSAR's spectral diversity converts the double-difference phase with the
    Doppler span of a whole burst (make_s1a_tops mode 2: spec_sep = linesPerBurst x mean kt x dta), but the two looks of
    a target are one burst CYCLE apart, so every ESD shift comes out cycle / linesPerBurst too small (T166 F1: 1343 /
    1497, -10.3 %)."""
    t = Path(xml).read_text(errors="replace")
    lpb = int(re.search(r"<linesPerBurst>(\d+)</linesPerBurst>", t).group(1))
    dta = float(re.search(r"<azimuthTimeInterval>([-+0-9.eE]+)</azimuthTimeInterval>", t).group(1))
    bt = [datetime.fromisoformat(m.group(1)) for m in _BURST.finditer(t)]
    return lpb, statistics.median((b - a).total_seconds() for a, b in zip(bt, bt[1:])) / dta


def _esd_factors(tr, d, n):
    """{frame stem: linesPerBurst / burst cycle} of date d's frames in subswath n (the super-master: the lab's frames);
    1 for a frame whose .BB the fork wrote with the A6 fix (five numbers in its header)."""
    safes = (_one_cut_per_frame(real(glob.glob(str(tr.labd / f"F*_F*/*_{tr.master}T*.SAFE")))) if d == tr.master
             else tr.state["stitched"][d])
    out = {}
    for t, x in _frames_line(safes, n):
        stem = f"S1_{t.name[15:23]}_{t.name[24:30]}_F{n}"
        bb = _esd_dir(tr, n, d) / f"{stem}.BB"
        # 11-01: the fork's make_s1a_tops takes spec_sep over the burst cycle itself (A6) and writes that cycle as the
        # .BB header's fifth number: no correction of ours then (T173's pass-1 overlaps predate it: four numbers)
        if bb.exists() and len(bb.read_text().split("\n", 1)[0].split()) >= 5:
            out[stem] = 1.0
            continue
        lpb, cyc = _burst_cycle(x)
        out[stem] = lpb / cyc
    return out


def _nesd_solve(master, obs):
    """Per date frame x [lines], the correction its a.grd takes, with x = 0 at the super-master's frame `master`:
    weighted least squares on the pairs (a, b, y, s): y = x_b - x_a, Huber reweighting (k = 1.345 robust sigma) to
    convergence. Only frames connected to the super-master through pairs are solved. Returns {stem: (x, sx, npairs)},
    the residuals [(a, b, r, z)] and the robust scale of z."""
    import numpy as np
    adj = defaultdict(set)
    for a, b, _, _ in obs:
        adj[a].add(b)
        adj[b].add(a)
    seen, todo = {master}, [master]
    while todo:
        for v in adj[todo.pop()]:
            if v not in seen:
                seen.add(v)
                todo.append(v)
    obs = [o for o in obs if o[0] in seen and o[1] in seen]
    nodes = sorted(seen - {master})
    if not nodes or not obs:
        return {}, [], float("nan")
    k = {s: i for i, s in enumerate(nodes)}
    A = np.zeros((len(obs), len(nodes)))
    y, s = np.array([o[2] for o in obs]), np.array([o[3] for o in obs])
    for i, (a, b, _, _) in enumerate(obs):
        if b != master:
            A[i, k[b]] += 1
        if a != master:
            A[i, k[a]] -= 1
    w0, wt, scale = 1 / s ** 2, np.ones(len(obs)), 1.0
    for _ in range(50):
        W = w0 * wt
        N = A.T @ (A * W[:, None])
        x = np.linalg.solve(N, A.T @ (W * y))
        z = (y - A @ x) * np.sqrt(w0)
        scale = max(1.4826 * float(np.median(np.abs(z))), 1e-6)
        new = np.minimum(1.0, 1.345 * scale / np.maximum(np.abs(z), 1e-12))
        if np.max(np.abs(new - wt)) < 1e-4:
            break
        wt = new
    sx = np.sqrt(np.diag(np.linalg.inv(N))) * scale
    npr = np.abs(A).sum(axis=0)
    sol = {v: (float(x[k[v]]), float(sx[k[v]]), int(npr[k[v]])) for v in nodes}
    res = [(a, b, float(r), float(zz)) for (a, b, _, _), r, zz in zip(obs, y - A @ x, z)]
    return sol, res, scale


def nesd(tr, expected):
    """Network spectral diversity (10-31; user: "go ahead on decision 2"): once pass 1 has aligned and checked every
    date, the residual azimuth shift between each date and its next NESD_K dates (and the super-master's neighbours)
    from GMTSAR's spectral_diversity on the burst overlaps pass 1 kept (OUT/<T>/esd/F<n>/<date>/), solved per date
    frame with the super-master at 0 (_nesd_solve), written per subswath to OUT/<T>/esd/F<n>/shifts.txt (frames the
    network cannot reach: their estimate against the super-master, flagged) with nesd.json (every pair, residual,
    solution); then start_pass2. Pairs are cached in pairs.jsonl. Returns True when the network was solved now."""
    if nesd_pass(tr) != 1:
        return False
    q, st = tr.state.get("qc", {}), tr.state.get("stitched", {})
    aligned = set(tr.state.get("aligned", []))
    aside = set(tr.state.get("bad_stitch", {})) | set(tr.state.get("align_failed", {})) | set(tr.state.get("no_orbit", {}))
    # every date pass 1 should hold, whatever stages this run was given (main's lists are empty without stitch / align)
    exp = (set(expected) | set(realign_dates(tr)) | {d for d, v in tr.state.get("plan", {}).items() if v.get("ready")})
    exp = {d for d in exp if d in st and d != tr.master}
    pending = sorted(d for d in exp if d not in aligned and d not in aside)
    unchecked = sorted(d for d in aligned if d not in q)
    if pending or unchecked:
        log(f"{tr.T} nesd: waits for pass 1 -- {len(pending)} dates not aligned {pending[:8]}, {len(unchecked)} not "
            f"checked {unchecked[:8]}")
        return False
    dates = sorted(d for d in aligned if q.get(d, {}).get("ok"))
    L = sorted(set(dates) | {tr.master})
    info = {"K": NESD_K, "subswaths": {}}
    share = subprocess.run("gmtsar_sharedir.csh", shell=True, capture_output=True, text=True, env=gmt_env()).stdout.strip()
    filt = f"{share}/filters/{SPECDIV_FILTER}"
    for n in tr.subs:
        fr = {}
        for d in L:
            f = _esd_dir(tr, n, d) / "frames.json"
            if not f.exists():
                raise RuntimeError(f"{tr.T} F{n} nesd: no burst overlaps of {d} ({f}) -- align it again in pass 1")
            fr[d] = json.loads(f.read_text())
        jobs = []
        for i, a in enumerate(L):
            for b in L[i + 1:i + 1 + NESD_K]:
                for sb, vb in fr[b].items():
                    sa = [s for s, va in fr[a].items() if va["master"] == vb["master"] and abs(va["tmp_da"]) < 1000]
                    if abs(vb["tmp_da"]) < 1000 and len(sa) == 1:
                        jobs.append((a, sa[0], b, sb))
        cache = _esd_dir(tr, n) / "pairs.jsonl"
        have = {}
        if cache.exists():
            for line in cache.read_text().split("\n"):
                if line.strip():
                    r = json.loads(line)
                    have[(r["a"], r["b"])] = r
        todo = [j for j in jobs if (j[1], j[3]) not in have]
        log(f"{tr.T} F{n} nesd: {len(jobs)} frame pairs ({len(todo)} to measure) among {len(L)} dates, K = {NESD_K}")
        lock = threading.Lock()

        def one(j):
            a, sa, b, sb = j
            r = _sd_pair(tr, n, sa, _esd_dir(tr, n, a), sb, _esd_dir(tr, n, b), filt)
            with lock:
                with open(cache, "a") as fh:
                    fh.write(json.dumps(r) + "\n")
                have[(sa, sb)] = r
            return r

        with ThreadPoolExecutor(max_workers=max(1, NESD_PROCS)) as ex:
            for i, _ in enumerate(ex.map(one, todo), 1):
                if i % 200 == 0:
                    log(f"{tr.T} F{n} nesd: {i} of {len(todo)} pairs measured")
        fac = {}                                       # A6: each frame's linesPerBurst / burst cycle
        for d in L:
            fac.update(_esd_factors(tr, d, n))
        miss = sorted({s for d in L for s in fr[d]} - set(fac))
        if miss:
            raise RuntimeError(f"{tr.T} F{n} nesd: no annotation for frames {miss[:5]}")
        cur = {(j[1], j[3]) for j in jobs}             # pairs of this network only (the cache may hold more)
        smaster = {s: v["master"] for d in L for s, v in fr[d].items()}
        sol, resid, fallback, scales = {}, [], {}, {}
        for m in sorted(fr[tr.master]):
            obs = [(r["a"], r["b"], r["shift"] * fac[r["b"]], r["sigma"] * fac[r["b"]]) for k, r in have.items()
                   if k in cur and smaster.get(k[0]) == m and r.get("n", 0) >= NESD_MIN_ROWS and "shift" in r]
            s, res, scale = _nesd_solve(m, obs)
            sol.update(s)
            resid += res
            scales[m] = scale
        rows = []
        for d in dates:
            for s, v in fr[d].items():
                if s in sol:
                    x, sx, npr = sol[s]
                else:
                    x, sx, npr = v.get("esd1_shift", 0.0) * fac[s], None, 0
                    fallback[s] = x
                rows.append((s, x, sx, npr, v["esd1_shift"] * fac[s] if "esd1_shift" in v else None))
        sh = _esd_dir(tr, n) / "shifts.txt"
        sh.write_text("".join(f"{s} {x:.9f}\n" for s, x, _, _, _ in rows))
        sep = statistics.median([r["spectral_spectrationXdta"] for r in have.values() if r.get("spectral_spectrationXdta")]
                                or [1.0])
        deg = 360.0 * sep / statistics.median(fac.values())    # lines -> degrees of the burst-boundary step (A6)
        rr = [abs(r) for _, _, r, _ in resid]
        dif = [x - e for s, x, sx, npr, e in rows if e is not None and s not in fallback]
        summ = {"pairs": len(cur), "used": len(resid), "frames": len(rows), "fallback": len(fallback),
                "robust_z_scale": scales, "resid_median_deg": round(statistics.median(rr) * deg, 2) if rr else None,
                "resid_p95_deg": round(sorted(rr)[int(0.95 * (len(rr) - 1))] * deg, 2) if rr else None,
                "x_minus_vs_master_median_deg": round(statistics.median(dif) * deg, 2) if dif else None,
                "x_range_deg": [round(min(x for _, x, _, _, _ in rows) * deg, 1),
                                round(max(x for _, x, _, _, _ in rows) * deg, 1)] if rows else None}
        (_esd_dir(tr, n) / "nesd.json").write_text(json.dumps(
            {"summary": summ, "factors": fac,
             "solution": {s: {"x": x, "sx": sx, "pairs": npr, "vs_master": e} for s, x, sx, npr, e in rows},
             "fallback": fallback, "residuals": [{"a": a, "b": b, "r": r, "z": z} for a, b, r, z in resid]}, indent=1))
        info["subswaths"][str(n)] = {"shifts": str(sh), **summ}
        log(f"{tr.T} F{n} nesd: {summ}")
    start_pass2(tr, dates, info)
    return True


def start_pass2(tr, dates, info):
    """Pass 2 of the network spectral diversity (10-31): the pass-1 images (geometry-only) move to
    OUT/<T>/superseded_nesd1/raw/F<n>/ (nothing deleted; the super-master's stay), pass 1's record is kept under
    aligned_nesd1 / qc_nesd1 / bad_align_nesd1 and cleared, and align() then aligns the network's `dates` again with
    its shifts (IGML_ESD_SHIFTS). Interferograms wait for pass 2 (_nesd_ready)."""
    moved = 0
    for n in tr.subs:
        src = tr.outp / f"F{n}/raw"
        if src.is_dir():
            for f in sorted(src.glob("S1_*_ALL_F*.*")):
                if f.suffix in (".SLC", ".PRM", ".LED") and not f.name.startswith(f"S1_{tr.master}_"):
                    _set_aside(f, tr.out / f"superseded_nesd1/raw/F{n}")
                    moved += 1
    for k in ("aligned", "qc", "bad_align"):
        tr.state[f"{k}_nesd1"] = tr.state.get(k)
    tr.state["aligned"], tr.state["qc"], tr.state["bad_align"] = [], {}, {}
    tr.state["nesd"] = dict(info, dates=dates, solved=time.strftime("%F %T"))
    tr.state["nesd"]["pass"] = 2
    tr.save()
    log(f"{tr.T}: network spectral diversity solved for {len(dates)} dates; pass 2 aligns them again with its shifts "
        f"({moved} pass-1 files set aside in OUT/{tr.T}/superseded_nesd1/raw)")


def _nesd_ready(tr):
    """Interferograms may be made: no network spectral diversity, or pass 2 has aligned (or set aside) every date of the
    network and qc has checked them."""
    if not NESD:
        return True
    if nesd_pass(tr) != 2:
        return False
    aligned, q = set(tr.state.get("aligned", [])), tr.state.get("qc", {})
    aside = set(tr.state.get("bad_stitch", {})) | set(tr.state.get("align_failed", {})) | set(tr.state.get("no_orbit", {}))
    return all(d in aligned or d in aside for d in tr.state["nesd"]["dates"]) and all(d in q for d in aligned)


def start_proc(tr):
    """Pairs stored by code older than PROC are made again (10-29): ours set aside (product_<T>/superseded_<old>/), the
    work copies of every interferogram and merge cleared (run_batched reuses a pair whose phasefilt.grd exists; the
    aligned SLC copies stay), the padded masks re-made. Once per track and version."""
    old = tr.state.get("proc", "proc0")
    if old == PROC:
        return
    if (tr.outp / "merge").exists():
        (tr.outp / f"superseded_{old}").mkdir(parents=True, exist_ok=True)
        os.replace(tr.outp / "merge", tr.outp / f"superseded_{old}/merge")
    tr.state[f"stored_{old}"] = tr.state.pop("stored", [])
    w = tr.work / "intf"
    for n in tr.subs:
        shutil.rmtree(w / f"F{n}/intf_all", ignore_errors=True)
    if (w / "merge").exists():
        for d in (w / "merge").iterdir():
            if re.match(r"^\d{7}_\d{7}$", d.name):
                shutil.rmtree(d, ignore_errors=True)
        for f in ("landmask_ra.grd", "mask_def.grd"):
            (w / "merge" / f).unlink(missing_ok=True)
    tr.state.pop("intf_pairs", None)
    tr.state["proc"] = PROC
    tr.save()
    log(f"{tr.T}: processing {old} -> {PROC}: {len(tr.state[f'stored_{old}'])} stored pairs set aside "
        f"(product_{tr.T}/superseded_{old}), work interferograms cleared; every pair is made again")


def batches(prs, max_dates):
    """The pairs in time order, cut so that a batch needs at most ``max_dates`` dates' SLCs on local disk."""
    out, cur, ds = [], [], set()
    for p in sorted(prs, key=lambda p: (p[1], p[0])):
        if cur and len(ds | set(p)) > max_dates:
            out.append(cur)
            cur, ds = [], set()
        cur.append(p)
        ds |= set(p)
    return out + ([cur] if cur else [])


def _resume_first(tr, todo, max_dates):
    """batches(), but the pairs whose dates' SLCs are all still copied in the work directory go first (10-30): a run
    killed mid-batch and launched again (new dates from recut / rescue shift every later boundary) finishes that batch
    instead of freeing its copies (~60 dates x 3 subswaths) only to copy them again batches later."""
    have = {f.name[3:11] for f in (tr.work / f"intf/F{tr.subs[0]}/raw").glob("S1_*_ALL_F*.SLC")}
    first = [p for p in todo if set(p) <= have]
    if not first or len({x for p in first for x in p}) > max_dates + 1:
        return batches(todo, max_dates)
    keep = {tuple(p) for p in first}
    return [sorted(first, key=lambda p: (p[1], p[0]))] + batches([p for p in todo if tuple(p) not in keep], max_dates)


def run_batched(tr, max_dates, ncores=6):
    """intf -> merge -> unwrap -> store, a time window at a time, clearing the local copies after each (10-23: T173's
    1,249 pairs need 230 dates' SLCs, 1.9 TB if copied at once). The super-master pair is made once and kept: every
    merge_list starts with it."""
    todo = [p for p in all_pairs(tr) if f"{doy0(p[0])}_{doy0(p[1])}" not in tr.state.get("stored", [])]
    B = _resume_first(tr, todo, max_dates)
    mdn = f"{doy0(tr.master_pair[0])}_{doy0(tr.master_pair[1])}"
    log(f"{tr.T}: {len(todo)} pairs to make, in {len(B)} batches of <= {max_dates} dates")
    limit = int(os.environ.get("BATCH_LIMIT", "0"))          # 10-29: a test run stops after this many batches
    for k, chunk in enumerate(B, 1):
        if limit and k > limit:
            log(f"{tr.T}: BATCH_LIMIT {limit} reached -- stopping before batch {k}")
            break
        log(f"{tr.T} batch {k}/{len(B)}: {len(chunk)} pairs, {chunk[0][1]} .. {chunk[-1][1]}")
        prs = [tr.master_pair] + [p for p in chunk if tuple(p) != tuple(tr.master_pair)]
        w = tr.work / "intf"
        need = {x for p in prs for x in p}
        for n in tr.subs:
            raw = w / f"F{n}/raw"
            if raw.exists():                                  # SLCs the batch does not need: free the disk
                for f in raw.glob("S1_*_ALL_F*.SLC"):
                    if f.name[3:11] not in need:
                        f.unlink()
        mine = [p for p in prs if not all((w / f"F{n}/intf_all/{doy0(p[0])}_{doy0(p[1])}/phasefilt.grd").exists()
                                          for n in tr.subs)]
        if mine:
            _intf_run(tr, mine, ncores)
        tr.state["intf_pairs"] = [list(p) for p in prs]
        tr.save()
        merge(tr)
        unwrap(tr)
        store(tr)
        stored = set(tr.state.get("stored", []))
        for a, b in chunk:                                    # this batch's scratch products, now in product_<T>/merge
            dn = f"{doy0(a)}_{doy0(b)}"
            if dn in stored and dn != mdn:
                shutil.rmtree(w / f"merge/{dn}", ignore_errors=True)
                for n in tr.subs:
                    shutil.rmtree(w / f"F{n}/intf_all/{dn}", ignore_errors=True)
        _prune_amp_cache(tr, {tr.master_pair[0]} | {a for ch in B[k:] for a, _ in ch})   # references still to come


def _intf_run(tr, prs, ncores):
    w = tr.work / "intf"
    jobs = []
    for n in tr.subs:
        f = w / f"F{n}"
        # 10-30: only the pairs this subswath lacks (a run launched again mid-batch made F1 for more pairs than F3)
        mine = [p for p in prs if not (f / f"intf_all/{doy0(p[0])}_{doy0(p[1])}/phasefilt.grd").exists()]
        if not mine:
            continue
        (f / "raw").mkdir(parents=True, exist_ok=True)
        (f / "topo").mkdir(parents=True, exist_ok=True)
        for d in sorted({x for p in mine for x in p}):
            src = slc_src(tr, d, n)
            for ext in ("SLC", "LED"):
                s, t = src.with_suffix("." + ext), f / "raw" / src.with_suffix("." + ext).name
                _copy_whole(s, t, link=str(src).startswith(f"{tr.out}/"))     # 11-01: our images hard-linked
            _copy_prm(prm_src(tr, d, n), f / "raw" / src.with_suffix(".PRM").name)
        for name, src in (("topo/topo_ra.grd", tr.labp / f"F{n}/topo/topo_ra.grd"),
                          ("batch_tops.config", tr.labp / f"F{n}/batch_tops.config")):
            if not (f / name).exists():
                shutil.copy(src, f / name)
        (f / "intf.in").write_text("".join(f"S1_{a}_ALL_F{n}:S1_{b}_ALL_F{n}\n" for a, b in mine))
        jobs.append(subprocess.Popen(f"intf_tops_parallel.csh intf.in batch_tops.config {ncores} > itp.log 2>&1",
                                     shell=True, cwd=f,
                                     env=proc_env(tr, n)))
    for j in jobs:
        j.wait()


def merge(tr):
    prs = [tuple(p) for p in tr.state.get("intf_pairs", [])]
    w = tr.work / "intf/merge"
    w.mkdir(parents=True, exist_ok=True)
    for name in ("dem.grd", "trans.dat", "batch_tops.config"):     # the masks: _pad_masks, after merging (10-29)
        src = tr.labp / ("topo/dem.grd" if name == "dem.grd" else f"merge/{name}")
        if src.exists() and not (w / name).exists():
            shutil.copy(src, w / name)
    lines = []
    for a, b in prs:
        dn = f"{doy0(a)}_{doy0(b)}"
        if all((tr.work / f"intf/F{n}/intf_all/{dn}/phasefilt.grd").exists() for n in tr.subs):
            lines.append(",".join(f"../F{n}/intf_all/{dn}/:S1_{a}_ALL_F{n}.PRM:S1_{b}_ALL_F{n}.PRM" for n in tr.subs))
        else:
            log(f"{tr.T} merge: {dn} missing a subswath's phasefilt.grd -- skipped")
    # merge_batch.csh takes the FIRST line's reference PRM as every pair's supermaster.PRM: without the super-master
    # pair there, the whole batch would be merged on another date's geometry, silently (audit 10-28)
    sm = f"{doy0(tr.master_pair[0])}_{doy0(tr.master_pair[1])}"
    if not lines or f"/intf_all/{sm}/:" not in lines[0].split(",")[0]:
        raise RuntimeError(f"{tr.T} merge: the super-master pair {sm} is not the first line of merge_list")
    (w / "merge_list").write_text("\n".join(lines) + "\n")
    log(f"{tr.T} merge {len(lines)} pairs")
    if len(tr.subs) == 1:                     # 10-29: one subswath (NEW tracks): merge_swath takes two or three
        _merge_single(tr, lines, w)
    else:
        run("merge_batch.csh merge_list batch_tops.config", w, w / "merge.log")
    if tr.new:
        _new_masks(tr, w)
    else:
        _pad_masks(tr, w)


def _merge_single(tr, lines, w):
    """merge_batch.csh for a single subswath (10-29): its grids are the frame's. Per pair merge/<pair>/ gets
    phasefilt.grd, corr.grd and mask.grd as they are, and supermaster.PRM as merge_batch.csh writes it (the first
    line's reference PRM with this pair's rshift, and its first_sample if larger); trans.dat is the subswath's own
    (dem2topo_ra.csh's, in the base: the same SAT_llt2rat projection merge_unwrap_geocode_tops.csh makes)."""
    first = lines[0].split(":")
    sm_prm = w / first[0] / first[1]
    for line in lines:
        pth, ref, _ = line.split(":")
        src = w / pth
        dn = Path(pth).name
        d = w / dn
        d.mkdir(exist_ok=True)
        for g in ("phasefilt.grd", "corr.grd", "mask.grd"):
            if not (d / g).exists():
                shutil.copy(src / g, d / g)
        p = _prm(src / ref)
        txt = sm_prm.read_text()
        sm = _prm(sm_prm)
        txt = re.sub(r"(?m)^(rshift\s*=\s*).*$", lambda m: m.group(1) + p.get("rshift", "0"), txt)
        if float(p.get("first_sample", 0)) > float(sm.get("first_sample", 0)):
            txt = re.sub(r"(?m)^(first_sample\s*=\s*).*$", lambda m: m.group(1) + p["first_sample"], txt)
        (d / "supermaster.PRM").write_text(txt)
        for f in ("trans.dat", "dem.grd", "batch_tops.config"):
            if (w / f).exists() and not (d / f).exists():
                os.symlink(f"../{f}", d / f)


def _new_masks(tr, w):
    """A NEW track's merged-frame land mask (10-29): the lab's tracks bring landmask_ra.grd (padded to our grid,
    _pad_masks); ours is made once on our own frame by GMTSAR's landmask.csh (GSHHG full coastline, the super-master
    pair's region) and kept in the base. No mask_def.grd (the lab's own). Two subswaths: the merged trans.dat the
    first merge made goes back to the base too. The merged frame (the super-master pair's supermaster.PRM) is kept
    there with it and every later merge must give the same, as _pad_masks requires the lab's of the nine (trans.dat is
    built on it); one subswath: it must be the subswath's own (topo/master.PRM, which dem2topo_ra.csh's trans.dat is
    built on)."""
    base = tr.labp / "merge"
    mdn = f"{doy0(tr.master_pair[0])}_{doy0(tr.master_pair[1])}"
    sm = w / mdn / "supermaster.PRM"
    ref = base / "supermaster.PRM" if (base / "supermaster.PRM").exists() else (
        tr.labp / f"F{tr.subs[0]}/topo/master.PRM" if len(tr.subs) == 1 else None)
    if ref is not None:
        a, b = _prm(sm), _prm(ref)
        if abs(float(a["clock_start"]) - float(b["clock_start"])) * 86400 * float(b["PRF"]) > 0.01 or \
                abs(float(a["near_range"]) - float(b["near_range"])) > 0.01:
            raise RuntimeError(f"{tr.T} merge: the merged frame differs from {ref} (clock_start {a['clock_start']} vs "
                               f"{b['clock_start']}, near_range {a['near_range']} vs {b['near_range']}): trans.dat "
                               f"would not fit")
    if not (base / "supermaster.PRM").exists():
        shutil.copy(sm, base / "supermaster.PRM")
    if (w / "trans.dat").exists() and not (base / "trans.dat").exists():
        shutil.copy(w / "trans.dat", base / "trans.dat")
    if (base / "landmask_ra.grd").exists():
        if not (w / "landmask_ra.grd").exists():
            shutil.copy(base / "landmask_ra.grd", w / "landmask_ra.grd")
        return
    ref = w / mdn / "phasefilt.grd"
    if not ref.exists():
        raise RuntimeError(f"{tr.T} merge: the super-master pair has no merged phasefilt.grd -- no land mask")
    region = subprocess.run(["gmt", "grdinfo", "-I-", str(ref)], capture_output=True, text=True,
                            env=gmt_env()).stdout.strip()[2:]
    run(f"landmask.csh {region}", w, w / "landmask.log")
    if not (w / "landmask_ra.grd").exists():
        raise RuntimeError(f"{tr.T} merge: landmask.csh made no landmask_ra.grd (see {w}/landmask.log)")
    shutil.copy(w / "landmask_ra.grd", base / "landmask_ra.grd")
    log(f"{tr.T} merge: land mask made on our frame ({region}) and kept in the base")


def _grid_region(g):
    out = subprocess.run(["gmt", "grdinfo", "-C", str(g)], capture_output=True, text=True, env=gmt_env()).stdout.split()
    return out[1:5] + out[7:9] if len(out) > 8 else None


def _pad_masks(tr, w):
    """The lab's merged-frame masks (landmask_ra.grd, mask_def.grd) on OUR merged grid (10-29): with each subswath
    placed exactly the merged grid gains up to a row and a column, and snaphu_interp.csh multiplies every pair's grids
    by these (a size mismatch only prints an error and unwraps unmasked). pad_mask.py copies the last row / column into
    the new cells and leaves every other value as it was (checked there)."""
    mdn = f"{doy0(tr.master_pair[0])}_{doy0(tr.master_pair[1])}"
    ref = w / mdn / "phasefilt.grd"
    if not ref.exists():
        raise RuntimeError(f"{tr.T} merge: the super-master pair has no merged phasefilt.grd -- no masks")
    # the merged frame must stay the lab's (trans.dat is built on it): exact placement takes each subswath's head up to
    # a whole row, which must not move the earliest subswath's (merge_swath's minh; true on all nine stacks, 10-29)
    ours, lab = _prm(w / mdn / "supermaster.PRM"), _prm(tr.labp / f"merge/{mdn}/supermaster.PRM")
    if abs(float(ours["clock_start"]) - float(lab["clock_start"])) * 86400 * float(lab["PRF"]) > 0.01 or \
            abs(float(ours["near_range"]) - float(lab["near_range"])) > 0.01:
        raise RuntimeError(f"{tr.T} merge: our merged frame differs from the lab's (clock_start {ours['clock_start']} vs "
                           f"{lab['clock_start']}, near_range {ours['near_range']} vs {lab['near_range']}): trans.dat "
                           f"would not fit")
    for name in ("landmask_ra", "mask_def"):
        src, out = tr.labp / f"merge/{name}.grd", w / f"{name}.grd"
        if src.exists() and (not out.exists() or _grid_region(out) != _grid_region(ref)):
            out.unlink(missing_ok=True)
            run(f"{sys.executable} {PAD_MASK} {src} {ref} {out}", w, w / "pad_mask.log")


def unwrap(tr, njobs=8):
    w = tr.work / "intf/merge"
    todo = [d for d in sorted(w.iterdir()) if re.match(r"^\d{7}_\d{7}$", d.name) and (d / "phasefilt.grd").exists()
            and not (d / "unwrap.grd").exists()]
    cmds = {d: unwrap_cmd(tr.T, ymd0(d.name[:7]), ymd0(d.name[8:])) for d in todo}
    log(f"{tr.T} unwrap {len(todo)} pairs: SNAPHU SMOOTH, DEFO (defomax {DEFOMAX}) for the "
        f"{sum(co for _, co in cmds.values())} spanning an earthquake ({UNWRAP_SCRIPT.get(tr.T, 'snaphu_interp.csh')})")
    if tr.T in UNWRAP_SCRIPT:
        for f in ("grid_nearest.m",):
            for d in todo:
                if not (d / f).exists():
                    shutil.copy(GS / f, d / f)
    (w / "unwrap.cmd").write_text("".join(
        f"cd {d} && ln -sf ../landmask_ra.grd . && ln -sf ../mask_def.grd . && {cmds[d][0]} > log_unwrap 2>&1\n"
        for d in todo))
    run(f"parallel --jobs {njobs} < unwrap.cmd", w, w / "unwrap_parallel.log", check=False)


def store(tr):
    w = tr.work / "intf/merge"
    stored = tr.state.setdefault("stored", [])
    n = 0
    matlab_bad = []
    for d in sorted(w.iterdir()):
        if re.match(r"^\d{7}_\d{7}$", d.name) and (d / "unwrap.grd").exists() and d.name not in stored:
            # the lab's snaphu_ocean_interp.csh goes on without its MATLAB interpolation if MATLAB fails (audit 10-28)
            if tr.T in UNWRAP_SCRIPT and (d / "log_unwrap").exists() and re.search(
                    r"License checkout failed|Error using|Undefined function|matlab: [Cc]ommand not found|Out of memory",
                    (d / "log_unwrap").read_text(errors="replace")):
                matlab_bad.append(d.name)
                continue
            dest = tr.outp / "merge" / d.name
            dest.mkdir(parents=True, exist_ok=True)
            # merge_log*: merge_swath's seams (n1 / n2 vary by pair), needed to undo its F2 column shift (audit 10-28)
            for f in ("unwrap.grd", "corr.grd", "conncomp.grd", "phasefilt.grd", "mask.grd", "supermaster.PRM",
                      "merge_log", "merge_log_corr", "merge_log_mask"):
                if (d / f).exists():
                    shutil.copy(d / f, dest / f)
            # 10-29: the processing that made it (igml.real applies no conversion fixes to these) and, for a lab pair,
            # its subswath products too, for replace()
            a, b = ymd0(d.name[:7]), ymd0(d.name[8:])
            (dest / "igml_proc.json").write_text(json.dumps({
                "proc": PROC, "subswath_origins": tr.origins, "gmtsar_fix": str(PROC_FIX.parent),
                "eap_removed": sorted(x for x in (a, b) if x in tr.state.get("eap_done", []))}) + "\n")
            if REPROCESS_LAB and (tr.labp / "merge" / d.name).is_dir():
                mdn = f"{doy0(tr.master_pair[0])}_{doy0(tr.master_pair[1])}"
                for sw in tr.subs:                 # not n: that counts the stored pairs
                    src, dst = tr.work / f"intf/F{sw}/intf_all/{d.name}", tr.outp / f"F{sw}/intf_all/{d.name}"
                    if src.is_dir():
                        shutil.rmtree(dst, ignore_errors=True)
                        dst.parent.mkdir(parents=True, exist_ok=True)
                        if d.name == mdn:            # every batch's merge needs it: copied (links kept relative)
                            shutil.copytree(src, dst, symlinks=True)
                        else:                        # run_batched deletes it next: moved (same file system)
                            os.replace(src, dst)
            stored.append(d.name)
            n += 1
    tr.save()
    if matlab_bad:
        log(f"{tr.T} store: {len(matlab_bad)} pairs NOT stored, MATLAB failed in snaphu_ocean_interp.csh: {matlab_bad}")
    log(f"{tr.T} stored {n} pairs (total {len(stored)})")


def combine(tr):
    """A product directory igml.real can open: the lab's pairs and ours side by side, the rest the lab's."""
    root = tr.out / "combined" / tr.T / f"product_{tr.T}"
    m = root / "merge"
    m.mkdir(parents=True, exist_ok=True)
    for name in ("topo", "F1", "F2", "F3", "tropo"):
        if (tr.labp / name).exists() and not (root / name).exists():
            os.symlink(tr.labp / name, root / name)
    for name in ("trans.dat", "supermaster.PRM"):
        if (tr.labp / "merge" / name).exists() and not (m / name).exists():
            os.symlink(tr.labp / "merge" / name, m / name)
    n_lab = n_ours = 0
    for d in (tr.labp / "merge").iterdir():
        if re.match(r"^\d{7}_\d{7}$", d.name) and (d / "unwrap.grd").exists() and not (m / d.name).exists():
            os.symlink(d, m / d.name)
            n_lab += 1
    for d in (tr.outp / "merge").iterdir() if (tr.outp / "merge").exists() else []:
        if re.match(r"^\d{7}_\d{7}$", d.name) and not (m / d.name).exists():
            os.symlink(d, m / d.name)
            n_ours += 1
    log(f"{tr.T} combined: +{n_lab} lab pairs, +{n_ours} ours -> {root.parent}")


def replace(tr, dry=False):
    """Our pairs into the lab's stack (10-29; user: "we are good to replace the lab's files, as long as they are stored
    somewhere else and not deleted"). Only once every lab pair we can make is stored with PROC, so the lab's tree never
    mixes the old and the new merged-grid sizes. MOVED to LAB_ARCHIVE/<T>/<the same path> (one file system: renames,
    nothing deleted, nothing copied): each lab pair folder (merge/<pair>, F<n>/intf_all/<pair>), the merged-frame masks,
    and the aligned images we re-made (the EAP dates; the images the lab had emptied). Ours then move (pairs) or are
    copied (masks, images: our pipeline keeps reading them) into their place, and the lab's pair folders get back the
    links they had to the merge level. Every move is a line of LAB_ARCHIVE/<T>/manifest.jsonl. A lab pair we could not
    make is archived without a replacement and listed. Explicit stage only (--stages replace; --dry-run lists)."""
    if tr.new:
        log(f"{tr.T} replace: no lab stack (the base is ours) -- nothing to replace")
        return
    names = lab_pair_dirs(tr)
    stored = set(tr.state.get("stored", []))
    makeable = {f"{doy0(a)}_{doy0(b)}" for a, b in lab_pairs(tr)}
    missing = sorted(makeable - stored)
    if missing:
        log(f"{tr.T} replace: {len(missing)} of the lab's pairs not made yet (e.g. {missing[:3]}) -- nothing replaced")
        return
    ours = {n for n in names if (tr.outp / "merge" / n / "igml_proc.json").exists()}
    bad = [n for n in ours if json.loads((tr.outp / "merge" / n / "igml_proc.json").read_text()).get("proc") != PROC]
    if bad:
        raise RuntimeError(f"{tr.T} replace: {len(bad)} of our pairs are not {PROC}: {bad[:3]}")
    arch = LAB_ARCHIVE / tr.T / f"product_{tr.T}"
    q = tr.state.get("qc", {})
    # 11-05: every lab image of a date we re-aligned (qc ok) is swapped, the dates the lab aligned but never paired too
    # (T173: 7), so the lab's raw/ holds the images its pairs were made from; the super-master stays the lab's
    lab_raw = {m.group(1) for f in os.listdir(tr.labp / f"F{tr.subs[0]}/raw")
               for m in [re.match(rf"S1_(\d{{8}})_ALL_F{tr.subs[0]}\.SLC$", f)] if m}
    images = sorted(({d for d in tr.state.get("eap_done", [])} |
                     {d for d in tr.state.get("aligned", []) if q.get(d, {}).get("ok") and not tr.lab_aligned(d)
                      and (d in lab_raw or d in {ymd0(c) for nm in names for c in (nm[:7], nm[8:])})}) - {tr.master})
    log(f"{tr.T} replace{' (dry run)' if dry else ''}: {len(names)} lab pairs ({len(ours)} ours to put in, "
        f"{len(names) - len(ours)} archived without one), {len(images)} aligned images, the merged-frame masks "
        f"-> archive {arch}")
    if dry:
        return
    arch.mkdir(parents=True, exist_ok=True)
    man = open(LAB_ARCHIVE / tr.T / "manifest.jsonl", "a")

    def move(src, dst, why):
        if not (src.exists() or src.is_symlink()):
            return False
        if dst.exists() or dst.is_symlink():
            raise RuntimeError(f"replace: {dst} exists already -- not overwriting the archive")
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.replace(src, dst)
        man.write(json.dumps({"time": time.strftime("%F %T"), "from": str(src), "to": str(dst), "why": why}) + "\n")
        man.flush()
        return True

    for name in ("landmask_ra.grd", "mask_def.grd"):          # the merged-frame masks, padded to our grid
        if move(tr.labp / "merge" / name, arch / "merge" / name, "lab mask, superseded by the padded copy"):
            shutil.copy(tr.work / "intf/merge" / name, tr.labp / "merge" / name)
    left = []
    for name in names:
        move(tr.labp / "merge" / name, arch / "merge" / name, "lab pair (GMTSAR master 3dd43d3)")
        for n in (1, 2, 3):
            move(tr.labp / f"F{n}/intf_all/{name}", arch / f"F{n}/intf_all/{name}", "lab subswath pair")
        if name not in ours:
            left.append(name)
            continue
        os.replace(tr.outp / "merge" / name, tr.labp / "merge" / name)
        for f in ("dem.grd", "trans.dat", "landmask_ra.grd", "mask_def.grd", "batch_tops.config"):
            if (tr.labp / "merge" / f).exists() and not (tr.labp / "merge" / name / f).exists():
                os.symlink(f"../{f}", tr.labp / "merge" / name / f)
        for n in (1, 2, 3):
            if (tr.outp / f"F{n}/intf_all/{name}").is_dir():
                os.replace(tr.outp / f"F{n}/intf_all/{name}", tr.labp / f"F{n}/intf_all/{name}")
        man.write(json.dumps({"time": time.strftime("%F %T"), "installed": name, "proc": PROC}) + "\n")
    for d in images:
        for n in (1, 2, 3):
            for ext in ("SLC", "PRM", "LED"):
                f = f"S1_{d}_ALL_F{n}.{ext}"
                if (tr.outp / f"F{n}/raw/{f}").exists() and \
                        move(tr.labp / f"F{n}/raw/{f}", arch / f"F{n}/raw/{f}", "lab aligned image (EAP / emptied)"):
                    shutil.copy(tr.outp / f"F{n}/raw/{f}", tr.labp / f"F{n}/raw/{f}")
    man.close()
    tr.state["replaced"] = {"time": time.strftime("%F %T"), "proc": PROC, "pairs": len(names) - len(left),
                            "archived_only": left, "images": images}
    tr.save()
    log(f"{tr.T} replaced: {len(names) - len(left)} pairs, {len(images)} images; {len(left)} lab pairs archived "
        f"without a replacement {left[:5]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("track")
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--dates", default="")
    ap.add_argument("--batch", type=int, default=12)
    ap.add_argument("--batched", type=int, default=0,
                    help="run intf / merge / unwrap / store in time windows of at most this many dates")
    ap.add_argument("--dry-run", action="store_true", help="replace: list what would move, move nothing")
    ap.add_argument("--recut", default="", help="dates (comma-separated) to cut again from raw: their frames, alignment "
                    "and pairs set aside (OUT/<T>/rejected), once each")
    a = ap.parse_args()
    tr = Track(a.track)
    if Path("/proc/self/stat").exists():
        threading.Thread(target=guard_extorb, daemon=True).start()
    only = set(a.dates.split(",")) if a.dates else None
    st = a.stages.split(",")
    if tr.new and tr.master is None:            # 10-29: a track the lab never processed: our base first
        if "base" not in st:
            raise RuntimeError(f"{tr.T}: no base yet (OUT/{tr.T}/base): run the base stage first")
        build_base(tr)
        tr.load_base()
    log(f"{tr.T}: super-master {tr.master}, super-master pair {tr.master_pair}, out {tr.out}")
    if REALIGN_ALL and ("stitch" in st or "align" in st):     # 10-30 afternoon: once per track and PROC
        start_realign(tr)
    ready = plan(tr, only) if "plan" in st else [d for d, v in tr.state.get("plan", {}).items() if v["ready"]]
    if only:
        ready = [d for d in ready if d in only]
    if "stitch" in st:                          # 10-30: before `ours`, as a re-cut date of the lab's is ours to align
        recut(tr, [d for d in a.recut.split(",") if d])
    pl = tr.state.get("plan", {})
    ours = [d for d in ready if not (pl.get(d, {}).get("lab_aligned") and not pl[d].get("lab_aligned_bad")
                                     and tr.lab_aligned(d))]
    realign = realign_dates(tr) if ("stitch" in st or "align" in st) else []
    if realign:
        log(f"{tr.T}: {len(realign)} of the lab's dates to align again (its images are empty), as pair partners")
    # 10-29: the dates an IPF < EAP_IPF processed are aligned again (once) with the antenna pattern removed
    eap_new = [d for d in eap_dates(tr) if d not in tr.state.get("eap_done", [])] if "align" in st else []
    if eap_new:
        log(f"{tr.T}: {len(eap_new)} dates processed by IPF < {EAP_IPF}, aligned again without the elevation antenna "
            f"pattern: {eap_new}")
        tr.state["aligned"] = [d for d in tr.state.get("aligned", []) if d not in eap_new]
        for d in eap_new:
            tr.state.get("qc", {}).pop(d, None)
        tr.save()
        realign = sorted(set(realign) | set(eap_new))
    realign = [d for d in realign if d not in ours]     # a date in both lists once (T144 10-29: align crashed)
    if "stitch" in st:
        stitch(tr, ours + realign)
    if "rescue" in st:
        rescue(tr)
        datatake_rescue(tr)
    if "align" in st:
        align(tr, ours + realign, a.batch)
        if eap_new:
            tr.state["eap_done"] = sorted(set(tr.state.get("eap_done", [])) |
                                          (set(eap_new) & set(tr.state.get("aligned", []))))
            tr.save()
    if "qc" in st:
        qc(tr)
        redo = retry_failed(tr) if "align" in st else []
        if redo:                                  # once per date (state "retried"): align again, check again
            stitch(tr, redo)
            align(tr, redo, a.batch)
            qc(tr)
        reject(tr)
    if "nesd" in st and nesd(tr, ours + realign):     # 10-31: the network solved now -> its pass 2 in this run
        if "align" in st:
            align(tr, ours + realign, a.batch)
        if "qc" in st:
            qc(tr)
            reject(tr)
    if NESD and not _nesd_ready(tr):                  # 10-31: interferograms only from pass 2
        held = [s for s in ("intf", "merge", "unwrap", "store", "combine") if s in st]
        if held:
            log(f"{tr.T}: {held} wait for the network spectral diversity (pass {nesd_pass(tr)})")
            st = [s for s in st if s not in held]
    if "intf" in st:
        start_proc(tr)                        # 10-29: pairs from older code are made again
    if a.batched and "intf" in st:
        run_batched(tr, a.batched, ncores=INTF_NCORES)
    else:
        if "intf" in st:
            intf(tr, ncores=INTF_NCORES)
        if "merge" in st:
            merge(tr)
        if "unwrap" in st:
            unwrap(tr)
        if "store" in st:
            store(tr)
    if "combine" in st:
        combine(tr)
    if "replace" in st:                       # explicit only: not in STAGES
        replace(tr, dry=a.dry_run)
    log(f"{tr.T}: stages {st} done")


if __name__ == "__main__":
    main()
