"""Covariance of a QE charge-density.dat (rho, m) under base's 48 m-3m' ops {R|t}.
Symmetric density: rho(RG) = exp(-i RG.t) rho(G);  m(RG) = exp(-i RG.t) s det(R) R m(G), s = -1 if primed.
t recovered from the atom basis (quarter-lattice grid). Usage: python density_symmetry.py <charge-density.dat> ..."""
import sys, itertools
import numpy as np
from pathlib import Path
from scipy.io import FortranFile
REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
lines = (BASE / "Y2Ir2O7.win").read_text().splitlines()
lat = np.array([[float(x) for x in lines[lines.index("begin unit_cell_cart") + 2 + i].split()] for i in range(3)])
i0 = lines.index("begin atoms_cart")
atoms = []
for l in lines[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    atoms.append(np.array(l.split()[1:4], float))
atoms = np.array(atoms)
img = np.array(list(itertools.product(range(-2, 3), repeat=3)), float) @ lat
def defect(R, tau_cart):
    gp = atoms @ R.T + tau_cart
    return np.linalg.norm(gp[:, None, None] - atoms[None, :, None] - img[None, None], axis=-1).min(axis=(1, 2)).max()
TAUS = np.array(list(itertools.product([0, .25, .5, .75], repeat=3)))     # crystal
tau_cryst = []
for R in G["base_ops"]:
    d = [defect(R, t @ lat) for t in TAUS]
    assert min(d) < 1e-8
    tau_cryst.append(TAUS[int(np.argmin(d))])
def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name

for path in sys.argv[1:]:
    f = FortranFile(path, "r")
    gamma_only, ngm, nspin = f.read_ints(np.int32)
    b = f.read_reals(np.float64).reshape(3, 3)
    mill = f.read_ints(np.int32).reshape(ngm, 3)
    comp = np.array([f.read_reals(np.complex128) for _ in range(nspin)])   # rho, mx, my, mz
    index = {tuple(m): i for i, m in enumerate(mill)}
    scale = np.abs(comp).max(axis=1)
    print(f"\n{path}: ngm {ngm}, nspin {nspin}; max |rho| {scale[0]:.3e}, max |m| {scale[1:].max():.3e}")
    rows = []
    for g, (R, pr) in enumerate(zip(G["base_ops"], G["base_priming"])):
        Rc = (b.T @ np.linalg.inv(b.T))                 # identity, keeps shapes explicit
        G_cart = mill @ b
        mill_new = np.rint((G_cart @ R.T) @ np.linalg.inv(b)).astype(int)
        perm = np.array([index.get(tuple(m), -1) for m in mill_new])
        ok = perm >= 0
        phase = np.exp(-2j * np.pi * (mill_new[ok] @ tau_cryst[g]))
        s = -1.0 if str(pr) == "primed" else 1.0
        rho_res = np.abs(comp[0, perm[ok]] - phase * comp[0, ok]).max() / scale[0]
        m_pred = phase[None] * (s * np.linalg.det(R) * R @ comp[1:, ok])
        m_res = np.abs(comp[1:, perm[ok]] - m_pred).max() / scale[1:].max()
        rows.append((cls(R), str(pr), tau_cryst[g], rho_res, m_res, ok.mean()))
    print(f"  {'class':10s} {'prime':9s} {'t (cryst)':18s} {'rho rel':>10s} {'m rel':>10s}")
    for c, p, t, r, m, cov in sorted(rows, key=lambda x: x[4]):
        print(f"  {c:10s} {p:9s} {str(np.round(t, 2)):18s} {r:10.2e} {m:10.2e}" + ("" if cov == 1 else f"  (G coverage {cov:.3f})"))
