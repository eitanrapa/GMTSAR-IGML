import sys, numpy as np
sys.argv = ["x", "t173_f1_raw_led_prm.txt"]
src = open("satbase_oldnew.py").read().split("ds = sorted(set(LED)")[0]
exec(src)
def bparperp(ref, rep, a):
    p0, p1 = PRM[ref], PRM[rep]
    dr = 0.5 * 299792458.0 / float(p0["rng_samp_rate"])
    rc = float(p0["earth_radius"]) + float(p0["SC_height"]); ra = float(p0["earth_radius"])
    nr = float(p1["near_range"]); fr = nr + dr * int(p1["num_rng_bins"])
    a1 = (nr*nr + rc*rc - ra*ra)/(2*nr*rc); a2 = (fr*fr + rc*rc - ra*ra)/(2*fr*rc)
    rl = np.arccos((a1+a2)/2)
    a1 = (-nr*nr + rc*rc + ra*ra)/(2*ra*rc); a2 = (-fr*fr + rc*rc + ra*ra)/(2*ra*rc)
    rl = rl + np.arccos((a1+a2)/2)
    return a["b"]*np.sin(rl - np.radians(a["a_old"])), a["b"]*np.cos(rl - np.radians(a["a_old"])), \
           a["b"]*np.sin(rl - np.radians(a["a_new"])), a["b"]*np.cos(rl - np.radians(a["a_new"]))
for rep, tab in (("20141110", (26.415004285904, 31.173348868250)), ("20141204", (8.085986149380, 8.147352418559)),
                 ("20141228", (-10.769325785831, -1.530035495214))):
    r = baseline("20161006", rep)[0]
    print(rep, "table", tab, "mine old/new (bpar,bperp):", [round(v, 4) for v in bparperp("20161006", rep, r)])
