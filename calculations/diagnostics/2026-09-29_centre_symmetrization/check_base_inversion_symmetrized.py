#!/usr/bin/env python3
"""Does symmetrizing base's Wannier centres under inversion fix the A(R) residual?

check_magnetic_symmetry (jaemo route, QE 7.6 pair, 2026-09-29) found base_inversion's
H(R)-only residual fell 2-4x after the new_ns_nc DFT fix (0.045 -> 0.013) but the
TOTAL residual got worse (0.097 -> 0.152), now dominated by the A(R)-derived cross/
external parts (external relative residual 1.6). Hypothesis: transl_inv_full's
per-image phase exp(-i b.R/2) makes tied images differ by (-1)^n, so A(R) is not
covariant when the Wannier centres used to choose images are only *approximately*
inversion-symmetric (Wannierization noise). This script tests that in isolation:
rebuild base's H(R) and A(R) using centres forced to be EXACTLY inversion-symmetric
(the raw .chk centres, symmetrized), and compare the base_inversion covariance
residual to the un-symmetrized one already on record.

Reuses exactly the same base_inversion check as
validation/symmetry/check_magnetic_symmetry.ipynb (same RNG seed, same 7 generic
k-points, same n_occ/batch), so the numbers are directly comparable to
validation/results/magnetic_symmetry_jaemo_qe76/metrics.csv.

Run in the axion conda env: python check_base_inversion_symmetrized.py
"""

import sys
from pathlib import Path

import numpy as np
from pythtb import W90
from scipy.optimize import linear_sum_assignment

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from modules import curvature as cv
from wannier.position_matrix import (
    _position_matrix_for_model,
    apply_orbital_dependent_R_mapping,
    position_matrix_transl_inv_full,
)
from wannier.wannier_io import read_wannier_checkpoint

PREFIX = "Y2Ir2O7"
BASE_DIRECTORY = (
    REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
)
MP_GRID = (8, 8, 8)
N_OCCUPIED = 156
BATCH_SIZE = 4
OPTIONS = dict(zero_energy=0, min_hopping_norm=1e-5, max_distance=125, ignorable_imaginary_part=None)
INVERSION = -np.eye(3)

# ============================================================== 1. symmetrize centres
checkpoint = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
raw_centres = np.asarray(checkpoint["wannier_centres"], dtype=float)  # (176, 3) Cartesian
n = len(raw_centres)
print(f"Loaded {n} raw Wannier centres from {BASE_DIRECTORY / f'{PREFIX}.chk'}")

w90 = W90(str(BASE_DIRECTORY), PREFIX)
lattice = np.asarray(w90.lat, dtype=float)  # (3,3), rows = lattice vectors, Cartesian

shell = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
shift_cart = shell @ lattice  # (27, 3)

target = raw_centres @ INVERSION.T  # (n,3): where inversion sends each centre
# periodic distance from target[s] to raw_centres[t] + shift_cart[l], minimized over l
diff = target[:, None, None, :] - (raw_centres[None, :, None, :] + shift_cart[None, None, :, :])
dist_all = np.linalg.norm(diff, axis=-1)  # (n, n, 27)
best_l = dist_all.argmin(axis=-1)  # (n, n)
cost = dist_all.min(axis=-1)  # (n, n): best periodic distance from target[s] to raw_centres[t]

row, col = linear_sum_assignment(cost)
assert np.array_equal(row, np.arange(n)), "linear_sum_assignment should return sorted rows"
partner = col  # partner[s] = t
L = shift_cart[best_l[np.arange(n), partner]]  # (n,3): -c_s = c_{partner[s]} + L[s] (approx)

match_distance = cost[np.arange(n), partner]
print(f"Inversion partner match: max residual distance {match_distance.max():.4e} A, "
      f"mean {match_distance.mean():.4e} A")

involution_ok = np.array_equal(partner[partner], np.arange(n))
n_self_paired = int((partner == np.arange(n)).sum())
print(f"Matching is an involution: {involution_ok}; self-paired (on-axis) WFs: {n_self_paired}")
if not involution_ok:
    bad = np.where(partner[partner] != np.arange(n))[0]
    raise RuntimeError(f"Inversion matching is not an involution at indices {bad[:10]}...; "
                        "refusing to proceed (degenerate-centre ambiguity or hypothesis wrong).")

r = -raw_centres - L - raw_centres[partner]  # (n,3) defect per WF
symmetrized_centres = raw_centres + r / 2.0

# Verify: INVERSION @ e[s] - L[s] == e[partner[s]] to machine precision
check = symmetrized_centres @ INVERSION.T - L - symmetrized_centres[partner]
print(f"Symmetrized-centre exactness check (should be ~0): max |defect| = "
      f"{np.abs(check).max():.3e} A")
assert np.abs(check).max() < 1e-9, "symmetrization algebra is wrong"

print(f"Centre shift from symmetrizing: max |e-c| = {np.abs(symmetrized_centres - raw_centres).max():.4e} A, "
      f"mean = {np.linalg.norm(symmetrized_centres - raw_centres, axis=1).mean():.4e} A "
      f"(compare to raw inversion defect max {match_distance.max():.4e} A)")


def build_model(centers_cartesian, label):
    """Base model with H(R)/A(R) folded on the given centres (own w90, own .chk)."""
    w = W90(str(BASE_DIRECTORY), PREFIX)
    apply_orbital_dependent_R_mapping(w, MP_GRID, centers_cartesian=centers_cartesian, verbose=True)
    embedding = w.lattice.orb_vecs.copy()  # tau: base's own nominal .win positions, untouched
    model = w.model(orb_vecs=embedding, **OPTIONS)

    ck = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
    ck.pop("v_matrix", None)
    position = position_matrix_transl_inv_full(
        ck, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centers_cartesian, verbose=False
    )
    keys = sorted(position)
    data = {
        "AA": np.stack([position.pop(R).transpose(1, 2, 0) for R in keys]),
        "iRvec": np.asarray(keys, dtype=int),
        "wcc_cart": np.asarray(ck["wannier_centres"]),
    }
    model._pos_r, info = _position_matrix_for_model(model, data)
    print(f"  [{label}] A(R): {len(keys)} R vectors; "
          f"centres_max_error={info.get('centres_max_error', float('nan')):.3e} A")
    return model


# ============================================================== 2. curvature covariance test
def trace_curvature_vector(model, pos, kpoints, batch_size=BATCH_SIZE):
    kpoints = np.asarray(kpoints, dtype=float)
    reciprocal_inverse = np.linalg.inv(np.asarray(model.recip_lat_vecs, dtype=float))
    pieces = {name: [] for name in ("internal", "cross", "external", "total")}
    for start in range(0, len(kpoints), batch_size):
        stop = min(start + batch_size, len(kpoints))
        c = cv.connection(cv.fields(model, pos, kpoints[start:stop]), N_OCCUPIED)
        omega = {
            "internal": cv.omega_internal(c),
            "cross": cv.omega_cross(c),
            "external": cv.omega_external(c),
        }
        omega["total"] = omega["internal"] + omega["cross"] + omega["external"]
        for name in pieces:
            value_cartesian = np.einsum(
                "ia,jb,ab...->ij...", reciprocal_inverse, reciprocal_inverse, omega[name], optimize=True
            )
            trace = np.trace(value_cartesian, axis1=-2, axis2=-1)
            pieces[name].append(np.stack((trace[1, 2], trace[2, 0], trace[0, 1]), axis=-1))
    return {name: np.concatenate(values, axis=0) for name, values in pieces.items()}


rng = np.random.default_rng(20260825)
generic_cartesian = rng.uniform(-0.23, 0.23, size=(7, 3))  # same draw as the notebook


def run_inversion_test(centers_cartesian, label):
    model = build_model(centers_cartesian, label)
    pos = cv.position_terms(model, dtype=np.complex128)
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)
    generic = generic_cartesian @ reciprocal_inverse
    transformed = (generic_cartesian @ INVERSION.T) @ reciprocal_inverse  # antiunitary=False -> no extra sign

    original = trace_curvature_vector(model, pos, generic)
    evaluated = trace_curvature_vector(model, pos, transformed)
    print(f"  [{label}] base_inversion residuals (matching metrics.csv's definition: "
          f"max|residual| / max|[evaluated, predicted]|):")
    for part in ("internal", "cross", "external", "total"):
        predicted = np.linalg.det(INVERSION) * original[part] @ INVERSION.T
        residual = evaluated[part] - predicted
        max_residual = float(np.max(np.abs(residual)))
        scale = float(np.max(np.abs(np.concatenate((evaluated[part], predicted)))))
        rel = max_residual / max(scale, 1e-30)
        print(f"    {part:9s} rel={rel:.6f}  max_residual={max_residual:.6f} A^2  scale={scale:.6f} A^2")
    return model


print("\n=== raw (unsymmetrized) centres, own base rule ===")
run_inversion_test(raw_centres, "raw")

print("\n=== inversion-symmetrized centres ===")
run_inversion_test(symmetrized_centres, "symmetrized")
