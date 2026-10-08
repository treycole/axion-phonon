"""eig_symmetry.py for any run directory (k list from .nnkp); prints the worst |d eps| per op class
for the lowest N_BANDS bands. Usage: edit RUN below."""
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict
REPO = Path("/Users/treycole/Repos/axion-phonon")
RUN = Path(sys.argv[1])
N_BANDS = 156
GROUPS = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
text = (RUN / "Y2Ir2O7.nnkp").read_text().splitlines()
def block(name):
    i = next(i for i, l in enumerate(text) if l.strip() == f"begin {name}")
    j = next(j for j in range(i, len(text)) if text[j].strip() == f"end {name}")
    return text[i + 1:j]
recip = np.array([[float(x) for x in l.split()] for l in block("recip_lattice")])
kl = block("kpoints"); kred = np.array([[float(x) for x in l.split()[:3]] for l in kl[1:]])
mp = int(round(len(kred) ** (1 / 3)))
raw = np.loadtxt(RUN / "Y2Ir2O7.eig"); nb, nk = int(raw[:, 0].max()), int(raw[:, 1].max())
eig = raw[:, 2].reshape(nk, nb)
key = lambda k: tuple(np.round(np.mod(k, 1) * mp).astype(int) % mp)
index = {key(k): i for i, k in enumerate(kred)}
def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name
by = defaultdict(list); worst = []
for g, (R, pr) in enumerate(zip(GROUPS["base_ops"], GROUPS["base_priming"])):
    s = -1.0 if str(pr) == "primed" else 1.0
    kt = s * (kred @ recip) @ R.T @ np.linalg.inv(recip)
    perm = np.array([index[key(k)] for k in kt])
    m = np.abs(eig[perm, :N_BANDS] - eig[:, :N_BANDS]).max(); worst.append(m); by[(cls(R), str(pr))].append(m)
print(f"{RUN}\n  {nk} k ({mp}^3), {nb} bands; max over ops {1e3 * max(worst):.3f} meV")
for k, v in sorted(by.items(), key=lambda kv: max(kv[1])):
    print(f"  {k[0]:10s} {k[1]:9s} n={len(v)}  {1e3 * min(v):7.3f} / {1e3 * max(v):7.3f} meV")
