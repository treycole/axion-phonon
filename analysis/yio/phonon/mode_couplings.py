"""Gamma2- normal modes of Y2Ir2O7 and their axion couplings, for the paper tables (draft/main_draft.tex: yio-charges, yio-modes).

Force constants Phi_ij = -v_i . [F(Q v_j) - F(0)] / Q from the conv_thr 1e-10 SCFs at Q = 0.01 A (amplitude series,
log 2026-10-17), Phi c = w^2 M c with c^T M c = 1 amu, and dtheta/dq_n = sum_i c_in dtheta/dQ_i (Eq. mode-coupling in the
SMODES basis). dtheta/dQ_i: the production values (pipeline 2026-11-02, window f12.422_w22.644_k8, Q -> 0 limits of
a + bQ^2 + cQ^4; c2 ends, mean and Chern-Simons 3-form agree to 2e-6, data/Y2Ir2O7/README.md).
Ported from archive/diagnostics/2026-10-21_yio_mode_couplings/ (same formulas; paths and dtheta/dQ updated).
Axion charge z_kappa = (dtheta/dQ_i) / sum over the moved atoms of |v_kappa| (Y, Ir: 4 x sqrt(3); O 48f: 12 x 1).
"""
import re
from pathlib import Path

import numpy as np

# ---- settings
DATA = Path("/Users/treycole/Repos/axion-phonon/data")
BASE = DATA / "Y2Ir2O7/dft/base/scf/out"                                  # SCF 213998
MODES = [DATA / "Y2Ir2O7/dft/GM2-_p1_Q0.01/scf/out",       # Y, 213999
         DATA / "Y2Ir2O7/dft/GM2-_p2_Q0.01/scf/out",       # Ir, 214003
         DATA / "Y2Ir2O7/dft/GM2-_p3_Q0.01/scf/out"]       # O 48f, 214016
Q = 0.01                                                   # A per Cartesian component (raw SMODES vector, entries +-1)
DTHETA_DQ = np.array([-0.017120, 0.055037, 0.063880])      # rad/A, Q -> 0, pipeline 2026-11-02 (fig5.ipynb prints them)
MASSES = {"Y": 88.90584, "Ir": 192.217, "O": 15.9994}
NAT = 22
# --------------
RYBOHR_EVA = 13.605693 / 0.529177
CM = np.sqrt(1.602176634e-19 / 1e-20 / 1.66053907e-27) / (2 * np.pi * 2.99792458e10)   # sqrt(eV/A^2/amu) -> cm^-1
HBAR, AMU, C_CM = 1.054571817e-34, 1.66053907e-27, 2.99792458e10


def positions(folder):
    t = (folder / "Y2Ir2O7.scf.in").read_text()
    rows = [line.split() for line in t.split("ATOMIC_POSITIONS angstrom")[1].strip().splitlines()[:NAT]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def forces(folder):
    t = (folder / "Y2Ir2O7.scf.out").read_text()
    return np.array(re.findall(r"atom\s+\d+ type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", t), float)[:NAT] * RYBOHR_EVA


species, R0 = positions(BASE)
m = np.array([MASSES[re.sub(r"\d", "", s)] for s in species])
F0 = forces(BASE)
V = [(positions(f)[1] - R0) / Q for f in MODES]
dF = [forces(f) - F0 for f in MODES]
Phi = np.array([[-(V[i] * dF[j]).sum() / Q for j in range(3)] for i in range(3)])
M = np.array([[(m[:, None] * V[i] * V[j]).sum() for j in range(3)] for i in range(3)])
L = np.linalg.cholesky(M)
Li = np.linalg.inv(L)
w2, y = np.linalg.eigh(Li @ ((Phi + Phi.T) / 2) @ Li.T)
C = Li.T @ y                                               # columns c_n, c^T M c = 1 amu
print(f"Phi (eV/A^2):\n{np.round(Phi, 3)}\nasymmetry {np.abs(Phi - Phi.T).max() / np.abs(Phi).max():.2%}")
print("\naxion charges z = (dtheta/dQ_i) / sum_kappa |v_kappa|:")
for i, name in enumerate(["Y (16d)", "Ir (16c)", "O (48f)"]):
    norms = np.linalg.norm(V[i], axis=1)
    print(f"  {name}: {len(norms[norms > 1e-6])} atoms x |v| {norms[norms > 1e-6][0]:.4f}: z = {DTHETA_DQ[i] / norms.sum():+.5f} rad/A")
print("\nomega (cm^-1)  c (Y, Ir, O)                 Y/Ir/O share  dtheta/dq (rad/(A sqrt(amu)))  per A   delta_theta_zp (rad)")
for n in range(3):
    c = C[:, n]
    omega = np.sqrt(w2[n]) * CM
    dr = sum(c[i] * V[i] for i in range(3))                # displacement per unit q
    share = np.array([np.sum((c[i] * V[i]) ** 2) for i in range(3)])
    share /= share.sum()
    g = DTHETA_DQ @ c
    q_zp = np.sqrt(HBAR / (2 * 2 * np.pi * C_CM * omega)) / (1e-10 * np.sqrt(AMU))   # sqrt(hbar / 2 omega) in A sqrt(amu)
    print(f"{omega:9.1f}   {np.round(c, 5)}   {np.round(share, 2)}   {g:+.5f}   {g / np.linalg.norm(dr):+.4f}   {abs(g) * q_zp:.2e}")
