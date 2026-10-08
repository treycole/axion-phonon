"""Physical tie-tolerance window for YIO base.
Ideal centres = each WF's own atom (nearest atom mod lattice, kept in the WF's cell).
An entry (left, right, R) is a symmetry tie if its nearest/2nd-nearest image distances are
exactly equal (<1e-9 A) for the ideal centres. Physical tolerance must be >= every real-centre
gap of a tie entry, and < every real-centre gap of a non-tie entry."""
import sys
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
B = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
chk = read_wannier_checkpoint(B / "Y2Ir2O7.chk", include_m_matrix=False)
c = np.asarray(chk["wannier_centres"], float); lat = np.asarray(chk["real_lattice"], float)
lines = (B / "Y2Ir2O7.win").read_text().splitlines()
i0 = next(i for i, l in enumerate(lines) if l.strip().lower() == "begin atoms_cart")
atoms = []
for l in lines[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    p = l.split(); atoms.append((p[0], np.array(p[1:4], float)))
pos = np.array([a[1] for a in atoms])
img = np.array([(a, b, cc) for a in range(-2, 3) for b in range(-2, 3) for cc in range(-2, 3)], float) @ lat
ideal = np.empty_like(c); owner = []; off = []
for w in range(len(c)):
    cand = pos[:, None, :] + img[None]                       # (natom, nimg, 3)
    d = np.linalg.norm(cand - c[w], axis=-1)
    a, k = np.unravel_index(np.argmin(d), d.shape)
    ideal[w] = cand[a, k]; owner.append(atoms[a][0]); off.append(d[a, k])
off = np.array(off)
print(f"WF-to-own-atom offset: median {np.median(off):.4f} A, max {off.max():.4f} A")
for name in sorted(set(owner)):
    m = np.array([o == name for o in owner]); print(f"  {name:4s} n={m.sum():3d} offset max {off[m].max():.4f} A")

R = []
with open(B / "Y2Ir2O7_hr.dat") as f:
    f.readline(); nw = int(f.readline()); nR = int(f.readline())
    for _ in range((nR + 14) // 15): f.readline()
    for i in range(nR):
        R.append([int(x) for x in f.readline().split()[:3]])
        for _ in range(nw * nw - 1): f.readline()
Rc = np.array(R, float) @ lat
shift = (np.array([(a, b, cc) for a in (-1, 0, 1) for b in (-1, 0, 1) for cc in (-1, 0, 1)], float) * 8) @ lat

def gaps(centres, l):
    disp = centres[None, :, :] - centres[l]
    d = np.linalg.norm(Rc[:, None, None, :] + shift[None, :, None, :] + disp[:, None], axis=-1)
    order = np.argsort(d, axis=1)
    s = np.take_along_axis(d, order, axis=1)
    return s[:, 1] - s[:, 0], order[:, 0]

real_tie, real_non = [], []
for l in range(nw):
    gi, _ = gaps(ideal, l); gr, _ = gaps(c, l)
    tie = gi < 1e-9
    real_tie.append(gr[tie]); real_non.append(gr[~tie])
real_tie = np.concatenate(real_tie); real_non = np.concatenate(real_non)
print(f"\nsymmetry-tie entries (ideal gap < 1e-9): {real_tie.size}; non-tie: {real_non.size}")
for q in (0.5, 0.9, 0.99, 0.999, 1.0):
    print(f"  tie entries, real gap quantile {q:5.3f}: {np.quantile(real_tie, q):.3e} A")
print(f"  NON-tie entries, smallest real gap: {real_non.min():.3e} A;  "
      f"count < 1e-2: {(real_non < 1e-2).sum()},  < 2e-2: {(real_non < 2e-2).sum()},  < 5e-2: {(real_non < 5e-2).sum()}")
for t in (1e-5, 1e-3, 2e-3, 3e-3, 5e-3, 7e-3, 1e-2):
    print(f"  tol {t:.0e}: ties missed {(real_tie > t).sum():8d}   non-ties wrongly merged {(real_non <= t).sum():8d}")
