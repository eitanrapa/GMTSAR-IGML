"""Old (pre-d235ea2) vs new (3dd43d3) SAT_baseline alpha, re-implemented, on T173 F1 aligned PRM/LED files.

Mirrors SAT_baseline.c: get_seconds, coarse search on the repeat orbit (dt = (t12-t11)/(ns-1), ns2 = 2 ns for S1),
poly_interp refinement (100 samples, +-5 ms, quadratic fit of dist^2), get_sign, and the two bh definitions:
  old: bh = sign * sqrt(b^2 - bv^2)                (b = refined min distance, bv = d . r_hat(ref))
  new: bh = sign * |d . cross|, cross = r_hat x along, along = repeat velocity projected normal to r_hat
Then phasediff's flat-earth drho = sqrt(rho^2 + B^2 - 2 rho B sin(theta - alpha) - Bx^2) - rho with B = baseline_start
is evaluated across the swath for both alphas; the phase difference (mean removed) is reported.
Orbit interpolation: cubic Hermite on the 10-s state vectors (GMTSAR uses a 6-point Hermite; mm-level agreement).
"""
import re
import sys
from datetime import date

import numpy as np

src = sys.argv[1]
txt = open(src).read().split("=== ")[1:]
LED, PRM = {}, {}
for blk in txt:
    name, body = blk.split("\n", 1)
    name = name.strip()
    d = name[3:11]
    if name.endswith(".LED"):
        rows = [l.split() for l in body.strip().splitlines()]
        hdr = rows[0]
        a = np.array([[float(v) for v in r] for r in rows[1:]])
        LED[d] = dict(t=86400.0 * a[:, 1] + a[:, 2], p=a[:, 3:6], v=a[:, 6:9])
    elif name.endswith(".PRM"):
        p = {}
        for l in body.splitlines():
            m = re.match(r"\s*(\S+)\s*=\s*(\S*)", l)
            if m:
                p[m.group(1)] = m.group(2)           # last occurrence wins, as get_sio_struct
        PRM[d] = p


def herm(o, t):
    tt, P, V = o["t"], o["p"], o["v"]
    t = np.asarray(t, float)
    sc = t.ndim == 0
    t = np.atleast_1d(t)
    i = np.clip(np.searchsorted(tt, t) - 1, 0, len(tt) - 2)
    h = (tt[i + 1] - tt[i])[:, None]
    s = ((t - tt[i]) / (tt[i + 1] - tt[i]))[:, None]
    h00, h10, h01, h11 = 2 * s**3 - 3 * s**2 + 1, s**3 - 2 * s**2 + s, -2 * s**3 + 3 * s**2, s**3 - s**2
    r = h00 * P[i] + h10 * h * V[i] + h01 * P[i + 1] + h11 * h * V[i + 1]
    return r[0] if sc else r


def secs(p):
    prf = float(p["PRF"])
    m = int(p["nrows"]) - int(p["num_valid_az"])
    s = 86400.0 * float(p["clock_start"]) + (float(p.get("ashift", 0)) + float(p.get("sub_int_a", 0))) / prf + m / (2 * prf)
    return s, s + int(p["num_patches"]) * int(p["num_valid_az"]) / prf


def baseline(d1, d2):
    p1, p2, o1, o2 = PRM[d1], PRM[d2], LED[d1], LED[d2]
    t11, t12 = secs(p1)
    t13 = (t11 + t12) / 2
    prf = float(p1["PRF"])
    ns = int((t12 - t11) / (0.5 / prf))
    dt = (t12 - t11) / (ns - 1)
    ns2 = ns * 2
    X1 = [herm(o1, t) for t in (t11, t12, t13)]
    t21, t22 = secs(p2)
    ks = np.arange(-ns2, ns + ns2)
    ts = t21 + ks * dt
    XS = herm(o2, ts)
    out = []
    for X in X1:
        dist = np.linalg.norm(XS - X, axis=1)
        m = ks[np.argmin(dist)]
        # poly_interp
        tk = (np.arange(100) - 50 + 0.5) * (0.01 / 100)
        bs = np.sum((herm(o2, t21 + m * dt + tk) - X) ** 2, axis=1)
        c2, c1, c0 = np.polyfit(tk, bs, 2)
        tmin = t21 + m * dt - c1 / (2 * c2)
        b = np.sqrt(c0 - c1 * c1 / 4 / c2)
        X2 = herm(o2, tmin)
        dvec = X2 - X
        rad = X / np.linalg.norm(X)
        al = herm(o2, tmin + 0.1) - herm(o2, tmin - 0.1)
        al = al - (al @ rad) * rad
        al /= np.linalg.norm(al)
        cr = np.cross(rad, al)
        cr /= np.linalg.norm(cr)
        bv, bx, bhr = dvec @ rad, dvec @ al, dvec @ cr
        sign = 1
        if p1["orbdir"].startswith("D"):
            sign = -sign
        if p1["lookdir"].startswith("L"):
            sign = -sign
        if np.arctan2(X2[1], X2[0]) < np.arctan2(X[1], X[0]):
            sign = -sign
        arg = b * b - bv * bv
        bh_old = sign * np.sqrt(arg) if arg > 0 else np.nan
        bh_new = sign * abs(bhr)
        out.append(dict(b=b, bv=bv, bx=bx, bh_raw=bhr, bh_old=bh_old, bh_new=bh_new,
                        a_old=np.degrees(np.arctan2(bv, bh_old)), a_new=np.degrees(np.arctan2(bv, bh_new)),
                        Boff=(tmin - t21) * float(p2["SC_vel"])))
    return out


def ramp(d1, r, rho_lo, rho_hi):
    """Phase (rad) of phasediff's flat-earth term, new alpha minus old, across [rho_lo, rho_hi], mean removed (start)."""
    p1 = PRM[d1]
    lam = float(p1["radar_wavelength"])
    RE = float(p1["earth_radius"])
    c = RE + float(p1["SC_height_start"])
    rho = np.linspace(rho_lo, rho_hi, 200)
    cost = (rho**2 + c * c - RE * RE) / (2 * rho * c)
    sint = np.sqrt(1 - cost**2)
    B = r["b"]
    out = []
    for a in (r["a_new"], r["a_old"]):
        al = np.radians(a)
        drho = -rho + np.sqrt(rho**2 + B * B - 2 * rho * B * (sint * np.cos(al) - cost * np.sin(al)))
        out.append(-4 * np.pi / lam * drho)
    d = out[0] - out[1]
    return d - d.mean()


def day(s):
    return date(int(s[:4]), int(s[4:6]), int(s[6:]))


ds = sorted(set(LED) & set(PRM))
pairs = [(a, b) for i, a in enumerate(ds) for b in ds[i + 1:]
         if (day(b) - day(a)).days <= 84 or (a[4:6] in ("06", "07", "08") and b[4:6] in ("06", "07", "08")
                                              and abs((day(b) - day(a)).days - 365) <= 12)]
near = float(PRM[ds[0]]["near_range"])
rho_f1 = (799876.0, 799876.0 + 21328 * 2.3295)
rho_iw = (799876.0, 955000.0)          # approximate IW1 near .. IW3 far slant range
rows = []
for a, b in pairs:
    try:
        r = baseline(a, b)[0]
    except Exception as e:
        print("fail", a, b, e)
        continue
    f1 = ramp(a, r, *rho_f1)
    iw = ramp(a, r, *rho_iw)
    rows.append((a, b, r["b"], r["bv"], r["bx"], r["bh_raw"], r["bh_old"], r["bh_new"], r["a_old"], r["a_new"],
                 np.ptp(f1), np.ptp(iw), r["Boff"]))
A = np.array([x[2:] for x in rows], float)
print(f"{len(rows)} pairs; |bx| median {np.median(abs(A[:,2])):.3f} m, max {abs(A[:,2]).max():.3f} m; "
      f"|bx/bv| median {np.median(abs(A[:,2]/A[:,1])):.2e}")
dBh = A[:, 5] - A[:, 4]
print(f"|Bh_new - Bh_old| median {np.median(abs(dBh))*1e3:.3f} mm, 95% {np.percentile(abs(dBh),95)*1e3:.2f} mm, "
      f"max {abs(dBh).max()*1e3:.1f} mm")
pf1, piw = A[:, 8], A[:, 9]
for thr in (0.05, 0.1, 0.3, 1.0):
    print(f"pairs with new-old ramp > {thr} rad: across F1 {np.sum(pf1 > thr)}, across IW1-IW3 {np.sum(piw > thr)}")
print("worst 10 (ref rep b bv bx bh_new  ptp_F1 ptp_IW rad):")
for k in np.argsort(-piw)[:10]:
    x = rows[k]
    print(f"  {x[0]} {x[1]} b={x[2]:7.2f} bv={x[3]:7.2f} bx={x[4]:6.3f} bh={x[7]:7.3f}  F1 {x[10]:.3f}  IW {x[11]:.3f}")
print(f"B_offset_start: median {np.median(abs(A[:,10])):.3f} m, max {abs(A[:,10]).max():.3f} m")
