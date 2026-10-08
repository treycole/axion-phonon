"""L-point NSCF variants (4 k: (0,0,.5), (0,.5,0), (.5,0,0), (.5,.5,.5) crystal). The four L points are
one m-3m' orbit; (.5,.5,.5) is on Ir1's [111] axis. Prints max |eps_n(L_i) - eps_n(L_111)|, bands 1-156.
Usage: python nscf_L_analyse.py <nscf.out> ..."""
import re, sys
import numpy as np
from pathlib import Path
for path in sys.argv[1:]:
    t = Path(path).read_text()
    done = "JOB DONE" in t
    blocks = re.split(r"\n\s+k =", t.split("End of band structure calculation")[-1])[1:]
    ev = np.array([[float(x) for x in re.findall(r"-?\d+\.\d+", b.split("bands (ev):")[1].split("occupation numbers")[0].split("highest")[0])][:156] for b in blocks])
    d = [1e3 * np.abs(ev[i] - ev[3]).max() for i in range(3)]
    gap = re.findall(r"highest occupied, lowest unoccupied level \(ev\):\s+([-\d.]+)\s+([-\d.]+)", t)
    print(f"{path}{'' if done else '  (NOT DONE)'}\n  L(001),L(010),L(100) vs L(111): " + "  ".join(f"{x:.1f}" for x in d) + f" meV   gap {gap}")
