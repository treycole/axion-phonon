"""Per-operation correlation: DFT eigenvalue asymmetry (.eig, bands 1-156) vs curvature residual
(base_full_group_tolerances.npz, tol 5e-3)."""
import numpy as np
from pathlib import Path
HERE = Path(__file__).parent; REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
text = (BASE / "Y2Ir2O7.nnkp").read_text().splitlines()
blk = lambda n: text[text.index(f"begin {n}") + 1:text.index(f"end {n}")]
recip = np.array([[float(x) for x in l.split()] for l in blk("recip_lattice")])
kred = np.array([[float(x) for x in l.split()[:3]] for l in blk("kpoints")[1:]])
raw = np.loadtxt(BASE / "Y2Ir2O7.eig"); eig = raw[:, 2].reshape(int(raw[:, 1].max()), int(raw[:, 0].max()))
key = lambda k: tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8); idx = {key(k): i for i, k in enumerate(kred)}
dft_max, dft_rms = [], []
for R, pr in zip(G["base_ops"], G["base_priming"]):
    s = -1.0 if str(pr) == "primed" else 1.0
    perm = np.array([idx[key(k)] for k in s * (kred @ recip) @ R.T @ np.linalg.inv(recip)])
    d = eig[perm, :156] - eig[:, :156]
    dft_max.append(np.abs(d).max()); dft_rms.append(np.sqrt((d ** 2).mean()))
dft_max, dft_rms = np.array(dft_max)[1:] * 1e3, np.array(dft_rms)[1:] * 1e3
curv = np.load(HERE / "base_full_group_tolerances.npz")["tol_5e-03"][1:] * 100
for j, part in enumerate(["internal", "cross", "external", "total"]):
    print(f"{part:9s} corr with DFT max: {np.corrcoef(dft_max, curv[:, j])[0, 1]:+.3f}   with DFT rms: {np.corrcoef(dft_rms, curv[:, j])[0, 1]:+.3f}")
print(f"\nratio total-curvature% / DFT-rms-meV: median {np.median(curv[:, 3] / dft_rms):.3f}, spread {np.std(curv[:, 3] / dft_rms) / np.median(curv[:, 3] / dft_rms):.2f} (rel. std)")
order = np.argsort(dft_rms)
print("\n op  DFT max  DFT rms (meV)   curvature total / external (%)")
for g in order[::4]:
    print(f" #{g + 1:02d}  {dft_max[g]:6.3f}  {dft_rms[g]:7.4f}        {curv[g, 3]:5.2f} / {curv[g, 2]:5.1f}")
