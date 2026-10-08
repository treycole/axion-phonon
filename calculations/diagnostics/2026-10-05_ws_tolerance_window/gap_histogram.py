"""Distribution of (2nd-nearest - nearest) image distance over all (left, right, R) entries
that the WS fold decides, for YIO base on the .chk centres. Pure geometry."""
import sys
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
B = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
chk = read_wannier_checkpoint(B / "Y2Ir2O7.chk", include_m_matrix=False)
c = np.asarray(chk["wannier_centres"], float)
lat = np.asarray(chk["real_lattice"], float) if "real_lattice" in chk else None
print("chk keys:", sorted(chk))
# R list from _hr.dat header via wsvec? use the H(R) R-set from pythtb is heavy; read R list from _hr.dat quickly
import gzip
R = []
with open(B / "Y2Ir2O7_hr.dat") as f:
    f.readline(); nw = int(f.readline()); nR = int(f.readline())
    deg_lines = (nR + 14) // 15
    for _ in range(deg_lines): f.readline()
    for i in range(nR):
        line = f.readline().split(); R.append(tuple(int(x) for x in line[:3]))
        for _ in range(nw * nw - 1): f.readline()
R = np.array(R, float)
print("nw", nw, "nR", len(R), "lattice\n", lat)
mp = np.array([8, 8, 8], float)
img = np.array([(a, b, cc) for a in (-1, 0, 1) for b in (-1, 0, 1) for cc in (-1, 0, 1)], float)
shift = (img * mp) @ lat
Rc = R @ lat
gaps = []
for l in range(nw):
    disp = c[None, :, :] - c[l]                      # (1, nw, 3) right - left
    d = np.linalg.norm(Rc[:, None, None, :] + shift[None, :, None, :] + disp[:, None], axis=-1)  # (nR, 27, nw)
    s = np.sort(d, axis=1)
    gaps.append((s[:, 1] - s[:, 0]).ravel())
g = np.concatenate(gaps)
np.save(Path(__file__).with_name("gaps_base.npy"), g)
print("entries", g.size)
edges = [0, 1e-12, 1e-9, 1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1, 3, 100]
h, _ = np.histogram(g, bins=edges)
for a, b, n in zip(edges[:-1], edges[1:], h):
    print(f"  [{a:8.0e}, {b:8.0e})  {n:10d}")
