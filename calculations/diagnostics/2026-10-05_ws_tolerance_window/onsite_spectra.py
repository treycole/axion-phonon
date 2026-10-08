"""Spectra of the on-site H(R=0) blocks of symmetry-equivalent atoms (Wannier _hr.dat, base).
Gauge-free: g maps one atom's WF set onto the other's, so the spectra must coincide."""
import sys
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
chk = read_wannier_checkpoint(BASE / "Y2Ir2O7.chk", include_m_matrix=False)
c = np.asarray(chk["wannier_centres"], float); lat = np.asarray(chk["real_lattice"], float)
lines = (BASE / "Y2Ir2O7.win").read_text().splitlines()
i0 = next(i for i, l in enumerate(lines) if l.strip().lower() == "begin atoms_cart")
atoms = []
for l in lines[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    p = l.split(); atoms.append((p[0], np.array(p[1:4], float)))
pos = np.array([a[1] for a in atoms])
img = np.array([(a, b, cc) for a in range(-2, 3) for b in range(-2, 3) for cc in range(-2, 3)], float) @ lat
owner = np.array([np.unravel_index(np.argmin(np.linalg.norm(pos[:, None] + img[None] - c[w], axis=-1)), (len(pos), len(img)))[0] for w in range(len(c))])
# R = 0 block of _hr.dat
with open(BASE / "Y2Ir2O7_hr.dat") as f:
    f.readline(); nw = int(f.readline()); nR = int(f.readline())
    deg = []
    while len(deg) < nR: deg += [int(x) for x in f.readline().split()]
    H0 = None
    for iR in range(nR):
        block = [f.readline() for _ in range(nw * nw)]
        R = tuple(int(x) for x in block[0].split()[:3])
        if R == (0, 0, 0):
            H0 = np.zeros((nw, nw), complex)
            for line in block:
                p = line.split(); H0[int(p[3]) - 1, int(p[4]) - 1] = float(p[5]) + 1j * float(p[6])
            H0 /= deg[iR]
            break
print(f"R=0 found (ndegen {deg[iR]}), hermiticity {np.abs(H0 - H0.conj().T).max():.1e}")
groups = {}
for a, (name, _) in enumerate(atoms):
    wf = np.where(owner == a)[0]
    ev = np.linalg.eigvalsh(H0[np.ix_(wf, wf)])
    groups.setdefault((name.rstrip("1234") if name.startswith("Ir") else name, len(wf)), []).append((a, name, ev))
for (kind, n), members in groups.items():
    # split O into orbits by spectrum similarity to the first member
    ref = members[0][2]
    print(f"\n{kind} ({n} WFs each), {len(members)} atoms")
    for a, name, ev in members:
        print(f"  atom {a + 1:2d} {name:4s}  spread from atom {members[0][0] + 1}: {1e3 * np.abs(ev - ref).max():9.3f} meV   ev[0..2] {np.round(ev[:3], 4)}")
