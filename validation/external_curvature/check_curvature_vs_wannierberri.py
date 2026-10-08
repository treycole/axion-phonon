#!/usr/bin/env python3
"""modules/curvature.py against WannierBerri's non-Abelian curvature.

WannierBerri's ``formula.covariant.Omega.nn(ik, inn, out)`` is the M x M non-Abelian
curvature of the band group ``inn`` with intermediate states ``out`` (three Cartesian
axial components).  ``wb.evaluate_k(..., formula=..., iband=...)`` returns that block at
one k-point, so both codes evaluate the same formula from the *same file*:
SrTiO3 trial04's ``_tb.dat`` for WannierBerri, and its ``_hr.dat`` / ``_r.dat`` (same
Wannier90 run) for us via PythTB.

trial04 is occupied-only (all 38 Wannier bands occupied), so for the full set the
complement is empty and only the external term is exercised.  A proper subset of the
bands (the lowest 30, 20) has a non-empty complement: the internal and cross terms are
then exercised as well.  At a generic k-point the curvature of any band subset is
algebraically well defined, and both codes compute the same formula for it.

Each code fixes its own eigenvector phases and the blocks live in the band basis, so only
gauge-invariant quantities are compared: the eigenvalues of every plane's Hermitian
matrix, and its band trace.  WannierBerri returns Cartesian axial components in Angstrom^2;
they are converted to our reduced two-form components with the same ``two_form`` the
curvature module uses (``pos.two_form``).

    python check_curvature_vs_wannierberri.py [N_K]
"""

import sys
from pathlib import Path

import numpy as np
import wannierberri as wb
from pythtb import W90

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from validation._paths import DATA
from modules.curvature import (
    connection, fields, omega_cross, omega_external, omega_internal, position_terms,
)

RUN = DATA / "SrTiO3/Ti_Q_2A/output/188914bc889/trial_04_Sr_sp_Ti_p_O_sp/proj_gauge/191583bc889"
PREFIX = "SrTiO3"
NUM_WANN = 38
N_K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
BAND_SETS = (38, 30, 20)  # 38 = external only (empty complement); 30, 20 exercise everything

k_points = np.random.default_rng(7).uniform(0.0, 1.0, size=(N_K, 3))

# ---- WannierBerri: from the _tb.dat of the same Wannier90 run -------------------------
system = wb.System_tb(str(RUN / f"{PREFIX}_tb.dat"), berry=True)
Omega = wb.formula.covariant.Omega
formulas = {"total": Omega, "internal": Omega, "external": Omega}
options = {
    "internal": {"external_terms": False},   # WannierBerri: internal only
    "external": {"internal_terms": False},   # WannierBerri: cross + external
}

# ---- ours: the same run through PythTB ----------------------------------------------
w90 = W90(str(RUN), PREFIX)
tau = w90.lattice.orb_vecs.copy()
model = w90.model(orb_vecs=tau)
del w90
pos = position_terms(model, dtype=np.complex128)
f = fields(model, pos, k_points)
two_form = pos.two_form.astype(float)

# ---- sanity: are the two codes looking at the same Hamiltonian? -----------------------
energies_ours = np.linalg.eigvalsh(f.H)
energies_wb = np.array([wb.evaluate_k(system, k=tuple(k), quantities=["energy"]) for k in k_points])
print(f"H(k) check: max |E_ours - E_WannierBerri| over {N_K} k-points and {NUM_WANN} bands "
      f"= {np.abs(energies_ours - energies_wb).max():.2e} eV")

worst = {}


def to_reduced_dual(block):
    """(n, n, 3) Cartesian axial components -> (3, n, n) reduced dual components."""
    return np.einsum("la,mna->lmn", two_form, block)


def compare(name, ours, theirs):
    """ours, theirs: (3, n, n) Hermitian matrices, one per plane.  Returns the worst relative error.

    The overall sign of WannierBerri's curvature is a convention, so the better of +/- is used
    (and reported); the eigenvalues themselves must agree in size.
    """
    worst_error = 0.0
    for plane in range(3):
        a = np.linalg.eigvalsh(ours[plane])
        errors = {s: np.abs(a - np.linalg.eigvalsh(s * theirs[plane])).max() for s in (+1.0, -1.0)}
        sign = min(errors, key=errors.get)
        scale = max(np.abs(a).max(), 1e-300)
        worst_error = max(worst_error, errors[sign] / scale)
        trace_error = abs(np.trace(ours[plane]).real - sign * np.trace(theirs[plane]).real)
        print(f"      {name:9s} plane {plane}: sign {sign:+.0f}  max|eig diff| {errors[sign]:.2e}  "
              f"(scale {scale:.2e}, relative {errors[sign] / scale:.1e})  |trace diff| {trace_error:.1e}")
    return worst_error


for n_occ in BAND_SETS:
    print(f"\nbands 0..{n_occ - 1} as the group (complement = {NUM_WANN - n_occ} bands)")
    c = connection(f, n_occ)
    pieces = {
        "total": omega_internal(c) + omega_cross(c) + omega_external(c),
        "internal": omega_internal(c),
        "external": omega_cross(c) + omega_external(c),
    }
    for ik, k in enumerate(k_points):
        print(f"   k = {np.round(k, 3)}")
        wb_blocks = wb.evaluate_k(
            system, k=tuple(k), formula=formulas, param_formula=options,
            iband=list(range(n_occ)), return_single_as_dict=True,
        )
        for name in ("total", "internal", "external"):
            ours = np.stack([pieces[name][(l + 1) % 3, (l + 2) % 3, ik] for l in range(3)])
            theirs = to_reduced_dual(wb_blocks[name])
            if n_occ == NUM_WANN and name != "total":
                continue  # empty complement: internal is identically zero, external == total
            error = compare(name, ours, theirs)
            worst[(n_occ, name)] = max(worst.get((n_occ, name), 0.0), error)

print("\nworst relative eigenvalue disagreement:")
for (n_occ, name), error in sorted(worst.items(), reverse=True):
    print(f"   bands 0..{n_occ - 1:<2d} {name:9s} {error:.1e}")
