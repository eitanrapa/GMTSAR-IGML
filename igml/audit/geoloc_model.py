"""The geolocation error of every lab stack's merged products, per subswath, from the parts MEASURED on 10-29
(2026-10-29; user: "let's dig into the geolocation problem and find out what is really happening before doing
anything"):

  1. Sentinel-1 + GMTSAR's SLC geometry: right (Rosamond reflectors at +0.4 samples / -0.2 lines, cr_peak.py)
  2. conv's decimation: every 32 x 6 cell holds the filter centred on its FIRST sample / line, declared at the cell
     centre -> declared - true = +16 samples, +3 lines (reflectors in amp1 / amp2 / display_amp; geoloc_conv2.py
     on T173 F1-F3)
  3. merge_swath places each subswath at a whole column / row: F_n's column k at merged column P_n + k, its row s at
     merged row H_n + s (P = 0, N1 - ovl12 + 1, N1 + N2 - ovl12 - ovl23; H from round(start-time difference / 6
     lines), shifted so the earliest is 0), where the continuous geometry wants C_n = (near_n - near_1) / dr samples
     and D_n = (start_n - start_merged) * PRF lines (value-for-value match on T071: 624 / 77, 1369 / 155)
  4. trans.dat = SAT_llt2rat on the merged supermaster.PRM -- tested here for its rshift (merge_swath copies the
     pair reference's rshift into that PRM; SAT_llt2rat subtracts it)

so merged declared - true = (16 + 32 P_n - C_n, 3 + 6 H_n - D_n), and what the geocoder sees adds trans.dat's own
shift. Ground metres through SAT_llt2rat's Jacobian at the middle of each subswath. Our F2 column fix moves F2 by
-32 samples. Read only; writes only the work directory.

usage (seasat venv, GMTSAR on PATH): python geoloc_model.py WORKDIR T:PAIR [T:PAIR ...]
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

C = 299792458.0


def prm(p):
    return {k.strip(): v.strip() for k, s, v in (l.partition("=") for l in open(p)) if s}


def llt2rat(prm_path, llh, cwd):
    out = subprocess.run(["SAT_llt2rat", str(prm_path), "1"], input="\n".join(" ".join(f"{v:.9f}" for v in p) for p in llh),
                         capture_output=True, text=True, cwd=cwd, check=True).stdout
    return np.array([[float(v) for v in l.split()[:2]] for l in out.strip().splitlines()])


if __name__ == "__main__":
    work = Path(sys.argv[1])
    for arg in sys.argv[2:]:
        T, pair = arg.split(":")                                  # T:PAIR (listing merge/ is slow on a busy disk)
        prod = next(Path("/zdata/datasets").glob(f"seasatHD*/InSAR_Processing/{T}/product_{T}"))
        sub = {n: prm(prod / f"F{n}/topo/master.PRM") for n in (1, 2, 3)}
        pdir = prod / "merge" / pair
        log = (pdir / "merge_log").read_text()
        ovl12, ovl23 = (int(v) for v in log.split("ovl12,23:")[1].split("\n")[0].split(","))
        M = prm(pdir / "supermaster.PRM")
        fs, prf = float(sub[1]["rng_samp_rate"]), float(sub[1]["PRF"])
        dr = C / (2 * fs)
        N = {n: int(float(sub[n]["num_rng_bins"]) // 32) for n in (1, 2, 3)}
        nr = {n: float(sub[n]["near_range"]) for n in (1, 2, 3)}
        cs = {n: float(sub[n]["clock_start"]) for n in (1, 2, 3)}
        head = {n: int(round((cs[n] - cs[1]) * 86400 * prf / 6)) for n in (1, 2, 3)}
        minh = min(head.values())
        H = {n: head[n] - minh for n in (1, 2, 3)}
        cs_c = cs[1] + minh * 6 / prf / 86400                     # what merge_swath's heads imply
        cs_m = float(M["clock_start"])                            # what the merged PRM (and so trans.dat) uses
        P = {1: 0, 2: N[1] - ovl12 + 1, 3: N[1] + N[2] - ovl12 - ovl23}
        Cn = {n: (nr[n] - nr[1]) / dr for n in (1, 2, 3)}
        Dn = {n: (cs[n] - cs_m) * 86400 * prf for n in (1, 2, 3)}
        # trans.dat's own frame: SAT_llt2rat on the merged PRM with rshift as written vs 0, on 200 nodes
        wd = work / T
        wd.mkdir(parents=True, exist_ok=True)
        led = M["led_file"]
        for src in (prod / "merge" / led, prod / f"F1/topo/{led}", prod / f"F1/raw/{led}"):
            if src.exists():
                (wd / led).write_bytes(src.read_bytes())
                break
        t = np.memmap(prod / "merge/trans.dat", dtype="<f8", mode="r").reshape(-1, 5)
        idx = np.sort(np.random.default_rng(0).choice(len(t), 200, replace=False))
        s = np.array(t[idx])
        np.savetxt(wd / "n.txt", s[:, [3, 4]], fmt="%.9f")
        h = subprocess.run(["gmt", "grdtrack", str(wd / "n.txt"), f"-G{prod / 'merge/dem.grd'}", "-nl"], capture_output=True,
                           text=True, check=True).stdout
        llh = np.array([[float(v) for v in l.split()] for l in h.strip().splitlines() if len(l.split()) == 3])
        keep = np.isin(np.round(s[:, 3], 7), np.round(llh[:, 0], 7))
        rs = float(M.get("rshift", 0) or 0)
        frame = {}
        for tag, val in (("as written", rs), ("0", 0.0)):
            txt = "".join(f"rshift = {val:g}\n" if l.split("=")[0].strip() == "rshift" else l for l in open(pdir / "supermaster.PRM"))
            (wd / f"m_{tag[0]}.PRM").write_text(txt)
            ra = llt2rat(wd / f"m_{tag[0]}.PRM", llh, wd)
            frame[tag] = np.median(s[keep][:, 0] - ra[:, 0]), np.median(s[keep][:, 1] - ra[:, 1])
        print(f"\n{T}: pair {pdir.name}, merged PRM rshift {rs:g}, clock_start {float(M['clock_start']):.12f} (from the "
              f"subswaths' heads {cs_c:.12f}, {(cs_m - cs_c) * 86400 * prf:+.2f} lines); trans.dat - SAT_llt2rat: rshift as written {frame['as written'][0]:+.2f} samples, "
              f"rshift 0 {frame['0'][0]:+.2f} samples", flush=True)
        tr = frame["0"][0]                                         # trans.dat range minus the rshift-free geometry
        # a node at continuous range R has trans r = R + tr and gets the cell declared there, whose content sits at
        # r - (declared - true): features move by (declared - true) - tr
        # Jacobian at the middle of each subswath: a trans.dat node there, offset 100 m east / north
        m0 = (wd / "m_0.PRM")
        for n in (1, 2, 3):
            rmid = Cn[n] + 32 * N[n] / 2
            k = np.argmin(np.abs(s[:, 0] - rmid) + 1e9 * ~keep)
            lon, lat = s[k, 3], s[k, 4]
            hh = llh[np.argmin(np.abs(llh[:, 0] - lon) + np.abs(llh[:, 1] - lat)), 2]
            de, dn = 100 / (111320 * np.cos(np.radians(lat))), 100 / 110950
            ra = llt2rat(m0, [(lon, lat, hh), (lon + de, lat, hh), (lon, lat + dn, hh)], wd)
            J = np.array([(ra[1] - ra[0]) / 100, (ra[2] - ra[0]) / 100]).T      # d(r, a) / d(E, N)
            for tag, fix in (("lab", 0), ("ours", -32 if n == 2 else 0)):
                d = np.array([16 + 32 * P[n] - Cn[n] + fix - tr, 3 + 6 * H[n] - Dn[n] - frame["0"][1]])
                e, nn = np.linalg.solve(J, d)
                print(f"  F{n} {tag:4s}: declared - true = range {d[0]:+6.1f} samples, azimuth {d[1]:+5.2f} lines -> "
                      f"product - truth = east {e:+6.1f} m, north {nn:+6.1f} m ({np.hypot(e, nn):5.0f} m)   "
                      f"[P {P[n]} vs {Cn[n] / 32:.3f} cells, H {H[n]} vs {Dn[n] / 6:.3f} rows]", flush=True)
