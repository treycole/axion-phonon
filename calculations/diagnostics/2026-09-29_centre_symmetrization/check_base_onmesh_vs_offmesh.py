#!/usr/bin/env python3
"""Is base's inversion-symmetry residual an aliasing/interpolation artifact?

Hypothesis: the finite-mesh representative-choice ambiguity in Fourier-
transforming A(R)/H(R) provably cannot affect H(k)/A(k) AT the original
ab-initio mesh points (see the (e^{ik_j N}-1)=0 argument worked through in
conversation) -- it can only show up when interpolating to k NOT on that
mesh. Base's own inversion symmetry (base_inversion) is exact (no O(Q)
breaking, unlike the mode pairs), so any residual here is pure floor, not
physics.

Test: evaluate the SAME base_inversion covariance check (as in
check_magnetic_symmetry.ipynb) at two different kinds of k-points:
  (a) "on-mesh": actual ab-initio mp_grid=(8,8,8) k-points. A Gamma-centered
      MP mesh is always closed under k -> -k, so both k and -k are exact
      mesh points here -- if the residual is an interpolation artifact, it
      should vanish (to floating-point precision) at these points.
  (b) "off-mesh": the SAME generic points check_magnetic_symmetry.ipynb uses
      (rng seed 20260825), which are deliberately NOT on the ab-initio mesh.

Run for BOTH A(R) formulas (transl_inv_full_from_chk and right_centre_from_chk,
own centres -- base is its own reference so own/shared coincide here).

Run in the axion conda env: python check_base_onmesh_vs_offmesh.py
"""

import sys
from pathlib import Path

import numpy as np
from pythtb import W90

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from modules import curvature as cv
from wannier.position_matrix import (
    _position_matrix_for_model,
    apply_orbital_dependent_R_mapping,
    position_matrix_mmn_formula,
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


def build_model(formula, label):
    """Base model, own centres, H(R)+A(R) both built by the given A(R) formula."""
    w = W90(str(BASE_DIRECTORY), PREFIX)
    checkpoint = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
    centres = np.asarray(checkpoint["wannier_centres"], dtype=float)
    apply_orbital_dependent_R_mapping(w, MP_GRID, centers_cartesian=centres, verbose=False)
    embedding = w.lattice.orb_vecs.copy()
    model = w.model(orb_vecs=embedding, **OPTIONS)

    if formula == "transl_inv_full_from_chk":
        position = position_matrix_transl_inv_full(
            checkpoint, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centres, verbose=False
        )
        keys = sorted(position)
        data = {
            "AA": np.stack([position.pop(R).transpose(1, 2, 0) for R in keys]),
            "iRvec": np.asarray(keys, dtype=int),
            "wcc_cart": centres,
        }
    else:
        checkpoint.pop("v_matrix", None)
        checkpoint2 = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
        data = position_matrix_mmn_formula(
            checkpoint2, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centres
        )
    model._pos_r, info = _position_matrix_for_model(model, data)
    print(f"  [{label}] {formula}: {len(data['iRvec'])} R vectors")
    return model


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


def inversion_residual(model, pos, kpoints):
    """base_inversion covariance test, matching metrics.csv's exact metric."""
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)
    k_cartesian = kpoints @ reciprocal
    transformed_cartesian = k_cartesian @ INVERSION.T  # antiunitary=False for pure inversion
    transformed = transformed_cartesian @ reciprocal_inverse

    original = trace_curvature_vector(model, pos, kpoints)
    evaluated = trace_curvature_vector(model, pos, transformed)
    out = {}
    for part in ("internal", "cross", "external", "total"):
        predicted = np.linalg.det(INVERSION) * original[part] @ INVERSION.T
        residual = evaluated[part] - predicted
        scale = np.abs(np.concatenate((evaluated[part], predicted)))
        # per-point relative residual (each of len(kpoints)*3 components), not just the max
        per_component = np.abs(residual).ravel() / np.maximum(scale.max(), 1e-30)
        out[part] = {
            "max": float(per_component.max()),
            "median": float(np.median(per_component)),
            "rms": float(np.sqrt(np.mean(per_component ** 2))),
            "n": len(per_component),
        }
    return out


# ============================================================== generic (off-mesh) points
# 60 points, same range as the notebook (Cartesian, A^-1), larger sample size
rng = np.random.default_rng(20260825)
generic_cartesian = rng.uniform(-0.23, 0.23, size=(60, 3))

# ============================================================== on-mesh points
# 60 distinct ab-initio mesh points (Gamma-centered MP mesh: closed under k -> -k exactly).
rng2 = np.random.default_rng(20260930)
all_indices = np.array(np.meshgrid(*[np.arange(n) for n in MP_GRID], indexing="ij")).reshape(3, -1).T
chosen = rng2.choice(len(all_indices), size=60, replace=False)
mesh_reduced = all_indices[chosen] / np.asarray(MP_GRID)

for formula in ("transl_inv_full_from_chk", "right_centre_from_chk"):
    print(f"\n=== {formula} ===")
    model = build_model(formula, "base")
    pos = cv.position_terms(model, dtype=np.complex128)
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)

    generic_reduced = generic_cartesian @ reciprocal_inverse
    off_mesh = inversion_residual(model, pos, generic_reduced)
    on_mesh = inversion_residual(model, pos, mesh_reduced)

    print(f"  {'part':10s} {'on max':>9s} {'on med':>9s} {'on rms':>9s}"
          f"  {'off max':>9s} {'off med':>9s} {'off rms':>9s}")
    for part in ("internal", "cross", "external", "total"):
        o, f = on_mesh[part], off_mesh[part]
        print(f"  {part:10s} {o['max']:9.5f} {o['median']:9.5f} {o['rms']:9.5f}"
              f"  {f['max']:9.5f} {f['median']:9.5f} {f['rms']:9.5f}")
