"""Set-wise symmetry defect of base's Wannier centres under each of the 48 ops of m-3m'
(origin of the .win frame, as derive_point_groups.py found the basis maps onto itself there).
defect(g) = max_w min_{w', L} |g c_w - c_w' - L|."""
import sys
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
B = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
chk = read_wannier_checkpoint(B / "Y2Ir2O7.chk", include_m_matrix=False)
c = np.asarray(chk["wannier_centres"], float); lat = np.asarray(chk["real_lattice"], float)
ops = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)["base_ops"]
img = np.array([(a, b, cc) for a in range(-2, 3) for b in range(-2, 3) for cc in range(-2, 3)], float) @ lat
lines = (B / "Y2Ir2O7.win").read_text().splitlines()
i0 = next(i for i, l in enumerate(lines) if l.strip().lower() == "begin atoms_cart")
atoms = []
for l in lines[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    p = l.split(); atoms.append(np.array(p[1:4], float))
atoms = np.array(atoms)
import itertools
TAUS = np.array(list(itertools.product([0, .25, .5, .75], repeat=3))) @ lat
def setwise(points, R, tau=np.zeros(3)):
    gp = points @ R.T + tau
    d = np.linalg.norm(gp[:, None, None, :] - points[None, :, None, :] - img[None, None], axis=-1)
    return d.min(axis=(1, 2))
tau_of = []
for R in ops:
    defects = [setwise(atoms, R, t).max() for t in TAUS]
    tau_of.append(TAUS[int(np.argmin(defects))])
atom_def = [setwise(atoms, R, t).max() for R, t in zip(ops, tau_of)]
print("atom-basis defect with best quarter translation, max over ops:", max(atom_def))
per_op = np.array([setwise(c, R, t) for R, t in zip(ops, tau_of)])          # (48, 176)
print(f"centre-set defect over all ops: max {per_op.max():.3e} A, median of per-op max {np.median(per_op.max(1)):.3e} A")
worst = np.argsort(per_op.max(0))[::-1][:8]
print("worst WFs (index: max defect):", [(int(w), f"{per_op[:, w].max():.2e}") for w in worst])
h, e = np.histogram(per_op.max(0), bins=[0, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 1e-1, 10])
print("per-WF max defect histogram:", list(zip([f"{x:.0e}" for x in e[:-1]], h)))
