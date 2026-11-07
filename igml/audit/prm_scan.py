# Every aligned SLC's PRM in the lab's stacks: a subswath without an alignment fit (no second rshift / ashift block,
# stretch_a 0) while the date's other subswaths have one, or a width / near range unlike the super-master's
import glob, os, re, sys
from collections import defaultdict

def prm(p):
    d = defaultdict(list)
    for line in open(p, errors="replace"):
        m = re.match(r"\s*(\S+)\s*=\s*(\S+)", line)
        if m:
            d[m.group(1)].append(m.group(2))
    return d

out = []
for prod in sorted(glob.glob("/zdata/datasets/seasatHD*/InSAR_Processing/T*/product_T*")):
    T = os.path.basename(prod).split("_")[1]
    by_date = defaultdict(dict)
    for p in glob.glob(f"{prod}/F[123]/raw/S1_*_ALL_F[123].PRM"):
        m = re.search(r"S1_(\d{8})_ALL_F(\d)\.PRM$", p)
        if m:
            by_date[m.group(1)][m.group(2)] = p
    if not by_date:
        continue
    fitted = lambda d: len(d.get("rshift", [])) > 1 or any(float(x) != 0 for x in d.get("stretch_a", ["0"]))
    info = {dt: {f: prm(p) for f, p in fs.items()} for dt, fs in by_date.items()}
    # the super-master: the date whose PRMs carry no fit in any subswath, with the most common width
    unf = [dt for dt, fs in info.items() if fs and not any(fitted(v) for v in fs.values())]
    bad = []
    for dt, fs in sorted(info.items()):
        fit = {f: fitted(v) for f, v in fs.items()}
        if any(fit.values()) and not all(fit.values()):
            bad.append(f"{dt} unfitted {''.join('F' + f for f, ok in sorted(fit.items()) if not ok)}")
    print(f"{T} {prod.split('/')[3]}: {len(info)} dates, all-unfitted {unf[:6]}{'...' if len(unf) > 6 else ''} ({len(unf)}), "
          f"partly unfitted {len(bad)}")
    for b in bad:
        print("   ", b)
