"""Estimate the elevation-antenna-pattern (EAP) phase that make_s1a_tops.c compute_eap would remove
from IPF 2.36 data (it does so only when manifest + aux-cal are concatenated to the annotation xml).
Re-implements compute_eap (make_s1a_tops.c:1001-1204 at 3dd43d3) in numpy; prints the phase that is
left uncorrected in our runs, per subswath.  Usage: python -I eap_est.py aux.xml iw1.xml [iw2.xml ...]
"""
import re
import sys
from datetime import datetime

import numpy as np


def tag(x, path_re, n=1):
    m = re.findall(path_re, x, flags=re.S)
    return m[n - 1] if n > 0 else m


def dayfrac(s):
    t = datetime.fromisoformat(s)
    return (t - t.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds() / 86400.0


def aux_eap(aux, swath, pol):
    for blk in re.findall(r"<calibrationParams>(.*?)</calibrationParams>", aux, flags=re.S):
        if f"<swath>{swath}</swath>" in blk and f"<polarisation>{pol}</polarisation>" in blk:
            ep = re.search(r"<elevationAntennaPattern>(.*?)</elevationAntennaPattern>", blk, flags=re.S).group(1)
            inc = float(re.search(r"<elevationAngleIncrement>([^<]+)<", ep).group(1))
            v = np.array(re.search(r"<values[^>]*>([^<]+)<", ep).group(1).split(), float)
            return inc, v[0::2] + 1j * v[1::2]
    raise SystemExit(f"no aux entry for {swath} {pol}")


def main():
    aux = open(sys.argv[1]).read()
    for f in sys.argv[2:]:
        x = open(f).read()
        swath = tag(x, r"<swath>([^<]+)</swath>")
        pol = tag(x, r"<polarisation>([^<]+)</polarisation>")
        ns = int(tag(x, r"<numberOfSamples>([^<]+)<"))
        ns -= ns % 4
        fs = float(tag(x, r"<rangeSamplingRate>([^<]+)<"))
        lpb = int(tag(x, r"<linesPerBurst>([^<]+)<"))
        dta = float(tag(x, r"<azimuthTimeInterval>([^<]+)<"))
        tau0 = float(tag(x, r"<slantRangeTime>([^<]+)</slantRangeTime>"))
        bl = x[x.index("<burstList"):x.index("</burstList>")]
        bt = re.findall(r"<azimuthTime>([^<]+)</azimuthTime>", bl)
        nb = len(bt) // 2 + 1
        t_brst = dayfrac(bt[nb - 1]) + dta * lpb / 2 / 86400.0
        ap = x[x.index("<antennaPatternList"):x.index("</antennaPatternList>")]
        ents = re.findall(r"<antennaPattern>(.*?)</antennaPattern>", ap, flags=re.S)
        at = np.array([dayfrac(re.search(r"<azimuthTime>([^<]+)<", e).group(1)) for e in ents])
        e = ents[int(np.argmin(np.abs(at - t_brst)))]
        tau_sub = np.array(re.search(r"<slantRangeTime[^>]*>([^<]+)<", e).group(1).split(), float)
        th_sub = np.array(re.search(r"<elevationAngle[^>]*>([^<]+)<", e).group(1).split(), float)
        anx = dayfrac(tag(x, r"<ascendingNodeTime>([^<]+)<")) * 86400.0
        t1 = t_brst * 86400.0 - anx
        H = [707714.8, 8351.5, 8947.0, 23.32, 11.74]
        ph = [0, 3.1495, -1.5655, -3.1297, 4.7222]
        W = 2 * np.pi / (12.0 * 86400.0 / 175.0)
        h = H[0] + sum(H[k] * np.sin(W * k * t1 + ph[k]) for k in range(1, 5))
        th_off = 29.450 - 0.0566 * (h / 1000.0 - 711.7)
        inc, G = aux_eap(aux, swath, pol)
        n = len(G)
        th_eap = (-(n - 1) // 2 + np.arange(n)) * inc + th_off
        tau = tau0 + (np.arange(ns) - 1) / fs
        theta = np.interp(tau, tau_sub, th_sub)          # correct weights (GMTSAR swaps them)
        p = np.angle(np.interp(theta, th_eap, G.real) + 1j * np.interp(theta, th_eap, G.imag))
        p = np.unwrap(p)
        d = np.degrees(p)
        print(f"{f.split('/')[-1]} {swath} {pol} burst {nb}: EAP phase near/mid/far = "
              f"{d[0]:.2f} / {d[ns // 2]:.2f} / {d[-1]:.2f} deg, range of variation {d.max() - d.min():.2f} deg "
              f"(look angle {theta[0]:.2f}-{theta[-1]:.2f}, boresight model {th_off:.2f})")


if __name__ == "__main__":
    main()


def residual(f, aux):
    """Same as main() for one file, returning (look angle, phase deg) for a linear-fit check."""
    x = open(f).read()
    swath = tag(x, r"<swath>([^<]+)</swath>"); pol = tag(x, r"<polarisation>([^<]+)</polarisation>")
    ns = int(tag(x, r"<numberOfSamples>([^<]+)<")); ns -= ns % 4
    fs = float(tag(x, r"<rangeSamplingRate>([^<]+)<")); lpb = int(tag(x, r"<linesPerBurst>([^<]+)<"))
    dta = float(tag(x, r"<azimuthTimeInterval>([^<]+)<")); tau0 = float(tag(x, r"<slantRangeTime>([^<]+)</slantRangeTime>"))
    bl = x[x.index("<burstList"):x.index("</burstList>")]
    bt = re.findall(r"<azimuthTime>([^<]+)</azimuthTime>", bl)
    out = []
    for nb in (1, len(bt) // 2 + 1, len(bt)):
        t_brst = dayfrac(bt[nb - 1]) + dta * lpb / 2 / 86400.0
        ap = x[x.index("<antennaPatternList"):x.index("</antennaPatternList>")]
        ents = re.findall(r"<antennaPattern>(.*?)</antennaPattern>", ap, flags=re.S)
        at = np.array([dayfrac(re.search(r"<azimuthTime>([^<]+)<", e).group(1)) for e in ents])
        e = ents[int(np.argmin(np.abs(at - t_brst)))]
        tau_sub = np.array(re.search(r"<slantRangeTime[^>]*>([^<]+)<", e).group(1).split(), float)
        th_sub = np.array(re.search(r"<elevationAngle[^>]*>([^<]+)<", e).group(1).split(), float)
        anx = dayfrac(tag(x, r"<ascendingNodeTime>([^<]+)<")) * 86400.0
        t1 = t_brst * 86400.0 - anx
        H = [707714.8, 8351.5, 8947.0, 23.32, 11.74]; ph = [0, 3.1495, -1.5655, -3.1297, 4.7222]
        W = 2 * np.pi / (12.0 * 86400.0 / 175.0)
        h = H[0] + sum(H[k] * np.sin(W * k * t1 + ph[k]) for k in range(1, 5))
        th_off = 29.450 - 0.0566 * (h / 1000.0 - 711.7)
        inc, G = aux_eap(aux, swath, pol); n = len(G)
        th_eap = (-(n - 1) // 2 + np.arange(n)) * inc + th_off
        tau = tau0 + (np.arange(ns) - 1) / fs
        theta = np.interp(tau, tau_sub, th_sub)
        p = np.degrees(np.unwrap(np.angle(np.interp(theta, th_eap, G.real) + 1j * np.interp(theta, th_eap, G.imag))))
        j = np.arange(ns); c = np.polyfit(j, p, 1); r = p - np.polyval(c, j)
        out.append((nb, p[0], p[-1], p.max() - p.min(), r.max() - r.min(), np.std(r)))
    print(swath, " ".join("burst %d: near %.1f far %.1f span %.1f, after linear fit p2p %.1f rms %.1f deg;" % o for o in out))


if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[-1] == "--resid":
    pass
