"""Does network spectral diversity break phase closure? (2026-11-06; user: "Do you think this NESD will break the phase closure that SBAS relies on?") One shift per date, applied
once to that date's image, is a per-date phase: it should cancel around any triangle, so closure should not change and the lab's burst steps should not show in closure either. Read
only, over the Mac's mount, T173's lab triplets, the lab's originals and ours, merged filtered phase where every leg's coherence > 0.3:
  wrapped closure: arg exp i(ab + bc - ac): circular rms [deg], share beyond 90 deg, and per subswath the along-track profile's spread (a burst sawtooth would show here);
  unwrapped closure: the share of pixels a whole cycle off, as unw_check.py;
  MIXED legs, on purpose: ab from the lab, bc and ac ours: what closure looks like if pre- and post-NESD pairs share a network (the along-track profile should show the sawtooth).
Triplets: those with the largest network correction on a leg (mixed satellites), and a spread of S1A-only ones.
usage (insar_sim env, PYTHONPATH=IGML): python nesd_closure.py OUT.npz [N_MIXED N_S1A]"""
import json
import re
import sys
from collections import defaultdict
from datetime import date, timedelta
from pathlib import Path

import numpy as np

from igml.real.grd import read_grd

out = sys.argv[1]
NM, NS = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (8, 8)
ARCH = Path("/seasat/IGML/lab_originals_20261029/T173/product_T173/merge")
OURS = Path("/seasat/seasatHD3/InSAR_Processing/T173/product_T173/merge")
OUT = Path("/seasat/IGML/InSAR_Processing/T173")
SEAMS = (636, 1381)
s = json.load(open(OUT / "state.json"))
sat = {d: Path(fr[0]).name[:3] for d, fr in s["stitched"].items() if fr}


def ymd(c):
    return (date(int(c[:4]), 1, 1) + timedelta(days=int(c[4:]))).strftime("%Y%m%d")


pairs = sorted(p.name for p in ARCH.iterdir() if re.match(r"^\d{7}_\d{7}$", p.name) and p.is_dir())
S = set(pairs)
by = defaultdict(list)
for p in pairs:
    a, b = p.split("_")
    by[a].append(b)
trip = [(a, b, c) for a in by for b in by[a] for c in by.get(b, []) if f"{a}_{c}" in S]
N = json.loads((OUT / "esd/F1/nesd.json").read_text())
sep = [json.loads(l).get("spectral_spectrationXdta") for l in open(OUT / "esd/F1/pairs.jsonl")]
sep = sorted(x for x in sep if x)
deg = 360 * sep[len(sep) // 2] / sorted(N["factors"].values())[len(N["factors"]) // 2]
x = defaultdict(list)
for k, v in N["solution"].items():
    x[k.split("_")[1]].append(v["x"] * deg)


def corr_of(a, b):
    A, B = x.get(ymd(a)), x.get(ymd(b))
    return abs(np.mean(B) - np.mean(A)) if A and B else 0.0


def sats(t):
    return {sat.get(ymd(c), "?") for c in t}


rank = sorted(((max(corr_of(t[0], t[1]), corr_of(t[1], t[2]), corr_of(t[0], t[2])), t) for t in trip), key=lambda r: -r[0])
days = lambda t: (date(int(ymd(t[2])[:4]), int(ymd(t[2])[4:6]), int(ymd(t[2])[6:])) - date(int(ymd(t[0])[:4]), int(ymd(t[0])[4:6]), int(ymd(t[0])[6:]))).days
mixed = [t for c, t in rank if len(sats(t)) > 1 and days(t) <= 60][:NM]
s1a = [t for t in trip if sats(t) == {"S1A"} and days(t) <= 48]
pick = [s1a[int(i)] for i in np.linspace(0, len(s1a) - 1, NS)]
print(f"{len(pairs)} lab pairs, {len(trip)} triplets; {NM} mixed-satellite (largest correction, <= 60 days), {NS} S1A-only", flush=True)


def load(root, p, name):
    z = read_grd(root / p / name).z.astype(np.float32)
    return z


def closure(ab, bc, ac):
    r, c = min(ab.shape[0], bc.shape[0], ac.shape[0]), min(ab.shape[1], bc.shape[1], ac.shape[1])
    return np.exp(1j * (ab[:r, :c] + bc[:r, :c] - ac[:r, :c]))


def stats(e, cm):
    r, c = min(e.shape[0], cm.shape[0]), min(e.shape[1], cm.shape[1])
    e, cm = e[:r, :c], cm[:r, :c]
    ok = np.isfinite(e.real) & cm & (e.real != 0)
    m = e[ok].mean()
    ang = np.degrees(np.angle(e[ok] * np.conj(m) / np.abs(m)))
    prof = []
    for c0, c1 in ((20, SEAMS[0] - 20), (SEAMS[0] + 20, SEAMS[1] - 20), (SEAMS[1] + 20, c - 20)):
        g = ok[:, c0:c1]
        z = np.where(g, e[:, c0:c1], 0).sum(axis=1)
        n = g.sum(axis=1)
        p = np.full(r, np.nan)
        p[n > 30] = np.degrees(np.angle(z[n > 30] * np.conj(z[n > 30].sum())))
        prof.append(p)
    return float(np.sqrt(np.mean(ang ** 2))), float(np.mean(np.abs(ang) > 90) * 100), [float(np.nanstd(p)) for p in prof], prof


rows = []
keep = {}
for kind, L in (("mixed", mixed), ("S1A", pick)):
    for a, b, c in L:
        ab, bc, ac = f"{a}_{b}", f"{b}_{c}", f"{a}_{c}"
        Pl = [load(ARCH, p, "phasefilt.grd") for p in (ab, bc, ac)]
        Po = [load(OURS, p, "phasefilt.grd") for p in (ab, bc, ac)]
        Cl = [load(ARCH, p, "corr.grd") for p in (ab, bc, ac)]
        Co = [load(OURS, p, "corr.grd") for p in (ab, bc, ac)]
        r0, c0 = min(z.shape[0] for z in Pl + Po), min(z.shape[1] for z in Pl + Po)
        ml = np.all([z[:r0, :c0] > 0.3 for z in Cl], axis=0)
        mo = np.all([z[:r0, :c0] > 0.3 for z in Co], axis=0)
        el, eo = closure(*Pl), closure(*Po)
        mix = closure(Pl[0], Po[1], Po[2])
        rl, bl, pl, prl = stats(el, ml)
        ro, bo, po, pro = stats(eo, mo)
        rm, bm, pm, prm = stats(mix, mo)
        Ul = [load(ARCH, p, "unwrap.grd").astype(np.float64) for p in (ab, bc, ac)]
        Uo = [load(OURS, p, "unwrap.grd").astype(np.float64) for p in (ab, bc, ac)]

        def cyc(U):
            r, cc = min(z.shape[0] for z in U), min(z.shape[1] for z in U)
            k = U[0][:r, :cc] + U[1][:r, :cc] - U[2][:r, :cc]
            k[(U[0][:r, :cc] == 0) | (U[1][:r, :cc] == 0) | (U[2][:r, :cc] == 0)] = np.nan
            v = k[np.isfinite(k)]
            return float(np.mean(np.abs(v - np.median(v)) > np.pi) * 100)
        ul, uo = cyc(Ul), cyc(Uo)
        name = f"{ymd(a)}({sat.get(ymd(a))}) {ymd(b)}({sat.get(ymd(b))}) {ymd(c)}({sat.get(ymd(c))})"
        cc = [round(corr_of(a, b)), round(corr_of(b, c)), round(corr_of(a, c))]
        print(f"{kind:5s} {name} correction on legs {cc} deg | wrapped rms lab {rl:5.1f} ours {ro:5.1f} mixed legs {rm:5.1f} | >90 deg lab {bl:5.2f} ours {bo:5.2f} mixed {bm:5.2f} % | "
              f"profile sd F1/F2/F3 lab {pl[0]:.1f}/{pl[1]:.1f}/{pl[2]:.1f} ours {po[0]:.1f}/{po[1]:.1f}/{po[2]:.1f} mixed {pm[0]:.1f}/{pm[1]:.1f}/{pm[2]:.1f} | a cycle off lab {ul:.3f} ours {uo:.3f} %", flush=True)
        rows.append({"kind": kind, "triplet": name, "corr": cc, "rms": [rl, ro, rm], "gt90": [bl, bo, bm], "prof": [pl, po, pm], "cyc": [ul, uo]})
        if kind == "mixed" and "lab" not in keep:
            keep.update(lab=np.array(prl), ours=np.array(pro), mix=np.array(prm), name=name)
for kind in ("mixed", "S1A"):
    R = [r for r in rows if r["kind"] == kind]
    print(f"{kind}: median wrapped rms lab {np.median([r['rms'][0] for r in R]):.1f}, ours {np.median([r['rms'][1] for r in R]):.1f}, mixed legs {np.median([r['rms'][2] for r in R]):.1f} deg; "
          f"F1 profile sd lab {np.median([r['prof'][0][0] for r in R]):.1f}, ours {np.median([r['prof'][1][0] for r in R]):.1f}, mixed {np.median([r['prof'][2][0] for r in R]):.1f}; "
          f"a cycle off lab {np.median([r['cyc'][0] for r in R]):.3f}, ours {np.median([r['cyc'][1] for r in R]):.3f} %; ours lower on {sum(r['cyc'][1] < r['cyc'][0] for r in R)} of {len(R)}", flush=True)
np.savez_compressed(out, rows=json.dumps(rows), **keep)
print("->", out, flush=True)
