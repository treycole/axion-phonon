import sys, numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
B = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
chk = read_wannier_checkpoint(B / "Y2Ir2O7.chk", include_m_matrix=False)
c = np.asarray(chk["wannier_centres"], float); lat = np.asarray(chk["real_lattice"], float)
R = []
with open(B / "Y2Ir2O7_hr.dat") as f:
    f.readline(); nw = int(f.readline()); nR = int(f.readline())
    for _ in range((nR + 14) // 15): f.readline()
    for i in range(nR):
        R.append([int(x) for x in f.readline().split()[:3]])
        for _ in range(nw * nw - 1): f.readline()
Rc = np.array(R, float) @ lat
shift = (np.array([(a, b, cc) for a in (-1, 0, 1) for b in (-1, 0, 1) for cc in (-1, 0, 1)], float) * 8) @ lat
inwin = 0; below = 0; lo_max = 0.0; hi_min = np.inf
for l in range(nw):
    d = np.linalg.norm(Rc[:, None, None, :] + shift[None, :, None, :] + (c - c[l])[None, None], axis=-1)
    x = d - d.min(axis=1, keepdims=True)
    x = x[x > 0]
    inwin += ((x > 2.0461e-3) & (x < 2.1247e-2)).sum()
    lo_max = max(lo_max, x[x < 1e-2].max(initial=0)); hi_min = min(hi_min, x[x >= 1e-2].min(initial=np.inf))
print(f"any image (not just 2nd) with excess distance in (2.05e-3, 2.12e-2): {inwin}")
print(f"largest excess < 1e-2: {lo_max:.4e}   smallest excess >= 1e-2: {hi_min:.4e}")
