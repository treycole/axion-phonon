"""Exact covariance test of the saved Hubbard ns_nc (occup.txt: complex ns(m1,m2,is,na), 5x5x4x22)
under base's 48 ops. ns(g a) must equal  W ns(a) W^dagger,  W = D_d(R) (x) U(R)  (spin index is=2*(s1-1)+s2),
with time reversal (ns -> (i sy) ns* (i sy)^dagger on the spin part) for primed ops.
Conventions (D vs D^T, U vs U^*, spin-block order) are not assumed: every combination is tried and the
one under which the 12 ops of Ir1's site group hold is reported, then applied to all 48."""
import sys, itertools
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
lines = (BASE / "Y2Ir2O7.win").read_text().splitlines()
lat = np.array([[float(x) for x in lines[lines.index("begin unit_cell_cart") + 2 + i].split()] for i in range(3)])
i0 = lines.index("begin atoms_cart"); atoms = []
for l in lines[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    atoms.append(np.array(l.split()[1:4], float))
atoms = np.array(atoms)
img = np.array(list(itertools.product(range(-2, 3), repeat=3)), float) @ lat
def mapping(R, t):
    gp = atoms @ R.T + t
    d = np.linalg.norm(gp[:, None, None] - atoms[None, :, None] - img[None, None], axis=-1).min(axis=2)
    return d.argmin(axis=1), d.min(axis=1).max()
TAUS = np.array(list(itertools.product([0, .25, .5, .75], repeat=3))) @ lat
maps = []
for R in G["base_ops"]:
    best = min((mapping(R, t) for t in TAUS), key=lambda x: x[1]); assert best[1] < 1e-8; maps.append(best[0])

def ylm_d(r):   # QE l=2 real harmonics order: z2, xz, yz, x2-y2, xy (unnormalised; D is fitted)
    x, y, z = r.T; rr = np.sqrt((r ** 2).sum(1))
    x, y, z = x / rr, y / rr, z / rr
    return np.stack([3 * z ** 2 - 1, np.sqrt(3) * 2 * x * z, np.sqrt(3) * 2 * y * z, np.sqrt(3) * (x ** 2 - y ** 2), np.sqrt(3) * 2 * x * y], 1)
pts = np.random.default_rng(0).normal(size=(200, 3))
SIGNS = np.ones(5)              # relative signs of QE's harmonics vs ylm_d (searched below)
def D_d(R):     # phi_m(R^-1 r) = sum_n phi_n(r) D_nm
    Rp = R * np.linalg.det(R)    # d orbitals are even under inversion
    A = ylm_d(pts) * SIGNS; B = ylm_d(pts @ Rp) * SIGNS   # rows: phi(R^-1 r) with R^-1 = R^T -> r @ R
    return np.linalg.lstsq(A, B, rcond=None)[0]
sig = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]
def U_spin(R):
    Rp = R * np.linalg.det(R); w, v = np.linalg.eig(Rp)
    ax = np.real(v[:, np.argmin(abs(w - 1))]); ang = np.arccos(np.clip((np.trace(Rp) - 1) / 2, -1, 1))
    if abs(ang) > 1e-9 and not np.allclose(Rp, np.eye(3)):   # fix the axis orientation so R(ax, ang) == Rp
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        if not np.allclose(np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K, Rp): ax = -ax
    return np.cos(ang / 2) * np.eye(2) - 1j * np.sin(ang / 2) * sum(a * s for a, s in zip(ax, sig))
iSy = 1j * sig[1]

def load(path):
    v = []
    for l in open(path):
        l = l.strip().strip("()")
        if l: re_, im = l.split(","); v.append(complex(float(re_), float(im)))
    ns = np.array(v).reshape(22, 4, 5, 5).transpose(3, 2, 1, 0) if False else np.array(v).reshape(22, 4, 5, 5)  # Fortran (m1,m2,is,na) -> [na,is,m2,m1]
    ns = ns.transpose(0, 1, 3, 2)                      # -> [na, is, m1, m2]
    return ns
def full(nsa, order):     # 10x10, index (s, m)
    M = np.zeros((10, 10), complex)
    for is_ in range(4):
        s1, s2 = order[is_]
        M[s1 * 5:(s1 + 1) * 5, s2 * 5:(s2 + 1) * 5] = nsa[is_]
    return M
ORDERS = {"is=2(s1-1)+s2": [(0, 0), (0, 1), (1, 0), (1, 1)], "is=s1+2(s2-1)": [(0, 0), (1, 0), (0, 1), (1, 1)]}
IR = [4, 5, 6, 7]
good = [g for g in range(48) if all(maps[g][a] == a for a in [4])]        # Ir1's site group
def residuals(ns, order, tD, cU):
    out = []
    for g, (R, pr) in enumerate(zip(G["base_ops"], G["base_priming"])):
        Dd = D_d(R); Dd = Dd.T if tD else Dd
        U = U_spin(R); U = {0: U, 1: U.conj(), 2: U.T, 3: U.conj().T}[cU]
        W = np.kron(U, Dd.T)                    # basis (s, m): spin outer
        worst = 0.0
        for a in IR:
            M = full(ns[a], order)
            if str(pr) == "primed":
                T = np.kron(iSy, np.eye(5)); M = T @ M.conj() @ T.conj().T
            pred = W @ M @ W.conj().T
            worst = max(worst, np.abs(full(ns[maps[g][a]], order) - pred).max())
        out.append(worst)
    return np.array(out)
for path in sys.argv[1:]:
    ns = load(path)
    print(f"\n{path}\n  Ir traces: " + " ".join(f"{np.trace(full(ns[a], ORDERS['is=2(s1-1)+s2'])).real:.5f}" for a in IR) +
          f";  non-Hubbard atoms max |ns| {np.abs(ns[[i for i in range(22) if i not in IR]]).max():.1e}")
    best = None
    for (oname, order), tD, cU in itertools.product(ORDERS.items(), (False, True), range(4)):
        r = residuals(ns, order, tD, cU)
        if best is None or r[good].max() < best[0][good].max(): best = (r, oname, tD, cU)
    r, oname, tD, cU = best
    print(f"  convention: {oname}, D transposed={tD}, U conjugated={cU}")
    print(f"  Ir1 site group ({len(good)} ops): max residual {r[good].max():.2e};  other {48 - len(good)} ops: min {np.delete(r, good).min():.2e} max {np.delete(r, good).max():.2e}  (|ns| ~ {np.abs(ns[4]).max():.2f})")

def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name
if "--detail" in sys.argv[-1:] or True:
    ns = load(sys.argv[1])
    print("\nper-op exact passes (residual < 1e-6) by convention, file", sys.argv[1])
    counts = []
    for sg in itertools.product((1, -1), repeat=4):
        SIGNS[:] = (1,) + sg
        for (oname, order), tD, cU in itertools.product(ORDERS.items(), (False, True), range(4)):
            r = residuals(ns, order, tD, cU); counts.append((int(np.sum(r < 1e-6)), tuple(SIGNS), oname, tD, cU))
    counts.sort(key=lambda x: -x[0])
    for c in counts[:6]:
        print(f"  {c[0]:2d}/48 exact: signs {c[1]} {c[2]} tD={c[3]} U-variant={['U','U*','U^T','U^dag'][c[4]]}")
    n, sg, oname, tD, cU = counts[0]; SIGNS[:] = sg
    r = residuals(ns, ORDERS[oname], tD, cU)
    from collections import defaultdict
    by = defaultdict(list)
    for g, (R, pr) in enumerate(zip(G["base_ops"], G["base_priming"])):
        by[(cls(R), str(pr), "fixes Ir1" if g in good else "moves Ir1")].append(r[g])
    for k, v in sorted(by.items()):
        print(f"  {k[0]:10s} {k[1]:9s} {k[2]:10s} n={len(v)}: " + " ".join(f"{x:.0e}" for x in sorted(v)))
