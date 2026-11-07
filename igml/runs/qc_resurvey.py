"""A finished track's images under the 11-06 qc (read only; nothing saved). Until 11-06 a re-aligned date had one near neighbour
and a far one, and the far one cancelled the misregistration test: every date is judged here as qc now would, against the
nearest date each side and, if it fails, the nearest 4. Prints what fails and a summary. usage (seasat venv): python qc_resurvey.py T"""
import importlib.util
import sys

spec = importlib.util.spec_from_file_location("pt", "/zdata/datasets/IGML/code/scratchpad/s1_gapfill/process_track.py")
P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
T = sys.argv[1]
tr = P.Track(T)
L = sorted(d for d, v in tr.state["qc"].items() if v.get("ok"))
print(f"{T}: {len(L)} dates that passed qc", flush=True)
fail, second = [], []
for i, d in enumerate(L):
    before, after = L[:i], L[i + 1:]
    near = sorted(before[-3:] + after[:3], key=lambda x: abs((P.day(x) - P.day(d)).days))
    nb = before[-1:] + after[:1]
    nb += [x for x in near if x not in nb][:2 - len(nb)]
    rec, ok = P._qc_judge(tr, d, nb, P.QC_COVER)
    spare = [x for x in near if x not in nb][:P.QC_NEIGHBOURS - len(nb)]
    if not ok and spare:
        second.append(d)
        rec, ok = P._qc_judge(tr, d, nb + spare, P.QC_COVER)
    if not ok:
        fail.append(d)
        print(f" {d} FAILS: " + " ".join(f"{k} vs {v['vs']} {v['coh']:.2f}/{v['floor']:.2f}" + (f" shifted {v['shifted']}" if "shifted" in v else "") +
              (f" content {v['content']}" if "content" in v else "") for k, v in rec.items() if k.startswith("F") and v and not v["ok"]), flush=True)
    if (i + 1) % 100 == 0:
        print(f" {i + 1} judged", flush=True)
print(f"{T}: {len(L)} judged; {len(second)} needed the second judge {second}; {len(fail)} fail {fail}", flush=True)
