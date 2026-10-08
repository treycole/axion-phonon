"""Layer 1: are the DFT eigenvalues (.eig, NSCF on the full 8x8x8 mesh, nosym) symmetric under
base's 48 m-3m' ops? eps_n(k) must equal eps_n(g k) (unprimed) or eps_n(-g k) (primed)."""
import sys
import numpy as np
from pathlib import Path
REPO = Path("/Users/treycole/Repos/axion-phonon"); sys.path.insert(0, str(REPO))
from wannier.wannier_io import read_wannier_checkpoint
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
GROUPS = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
chk = read_wannier_checkpoint(BASE / "Y2Ir2O7.chk", include_m_matrix=False)
kred = np.asarray(chk["kpt_red"], float).reshape(-1, 3)
recip = np.asarray(chk["recip_lattice"], float)
raw = np.loadtxt(BASE / "Y2Ir2O7.eig")
nb, nk = int(raw[:, 0].max()), int(raw[:, 1].max())
eig = raw[:, 2].reshape(nk, nb)
print(f"{nk} k-points, {nb} bands; kpt_red shape {kred.shape}")
index = {tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8): i for i, k in enumerate(kred)}
BANDS = {"occupied (1-156... of all 300 in window)": slice(0, 156), "below frozen max 15 eV": None, "all 300": slice(0, nb)}
worst = []
for g, (R, pr) in enumerate(zip(GROUPS["base_ops"], GROUPS["base_priming"])):
    sign = -1.0 if str(pr) == "primed" else 1.0
    kt = sign * (kred @ recip) @ R.T @ np.linalg.inv(recip)
    perm = np.array([index[tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8)] for k in kt])
    d = np.abs(eig[perm] - eig)
    frozen = (eig < 15.0) & (eig[perm] < 15.0)
    worst.append((d[:, :156].max(), d[frozen].max(), d[:, :290].max()))
w = np.array(worst)
print("max |eps(gk) - eps(k)| over ops (eV):  bands 1-156: %.2e   below 15 eV: %.2e   bands 1-290: %.2e" % tuple(w.max(0)))
print("per-op max (bands 1-156), sorted:", np.round(np.sort(w[:, 0])[::-1][:6], 7))

def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name
from collections import defaultdict
by = defaultdict(list)
for g, (R, pr) in enumerate(zip(GROUPS["base_ops"], GROUPS["base_priming"])):
    by[(cls(R), str(pr))].append(w[g, 0])
print("\nper class, max |d eps| bands 1-156 (meV): min / max over the class")
for k, v in sorted(by.items(), key=lambda kv: max(kv[1])):
    print(f"  {k[0]:10s} {k[1]:9s} n={len(v)}  {1e3*min(v):7.3f} / {1e3*max(v):7.3f}")

good = [g for g in range(48) if w[g, 0] < 3e-3]
print(f"\n{len(good)} ops with max |d eps| < 3 meV:")
for g in good:
    R = GROUPS["base_ops"][g]; P = R * np.linalg.det(R)
    wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))]) if round(np.trace(P)) != 3 else np.zeros(3)
    ax = ax / (np.abs(ax).max() or 1)
    print(f"  #{g:02d} {cls(R):10s} {str(GROUPS['base_priming'][g]):9s} axis {np.round(ax, 2)}  {1e3 * w[g, 0]:.3f} meV")
