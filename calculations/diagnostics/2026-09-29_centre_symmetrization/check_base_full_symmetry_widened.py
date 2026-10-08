#!/usr/bin/env python3
"""Reproduce check_magnetic_symmetry.ipynb's base table, on the production chk source,
at the now-widened WS-distance tolerance (1e-3 A, wannier/position_matrix.py default).

The existing validation/results/magnetic_symmetry_jaemo_qe76/metrics.csv used
POSITION_SOURCE="jaemo" (Jae-Mo's own fork-written _r.dat, read as written -- it never
calls _map_R_by_orbital_positions, so the Python-side tolerance never touched it and
those numbers cannot change). This script uses right_centre_from_chk (this session's
production-equivalent source, used for every dtheta rerun), which does go through the
tolerance-affected orbital-dependent R mapping.

Same three base tests as the notebook, same methodology (generic_cartesian: rng seed
20260825, N=7; z-axis line: N_POINTS=17, q in linspace(-0.25,0.25,17)), so the numbers
are directly comparable to the existing table.

Run in the axion conda env: python check_base_full_symmetry_widened.py
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
S4Z = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])

OPERATIONS = {
    "base_inversion": (INVERSION, False),
    "base_S4z_control": (S4Z, False),
}


def build_model(formula):
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
        checkpoint2 = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
        data = position_matrix_mmn_formula(
            checkpoint2, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centres
        )
    model._pos_r, info = _position_matrix_for_model(model, data)
    print(f"[{formula}] {len(data['iRvec'])} R vectors")
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


def transformed_kpoints(kpoints, reciprocal_lattice, rotation, antiunitary):
    k_cartesian = np.asarray(kpoints) @ reciprocal_lattice
    sign = -1.0 if antiunitary else 1.0
    transformed_cartesian = sign * k_cartesian @ rotation.T
    return transformed_cartesian @ np.linalg.inv(reciprocal_lattice)


for formula in ("right_centre_from_chk", "transl_inv_full_from_chk"):
    print(f"\n=== {formula} ===")
    model = build_model(formula)
    pos = cv.position_terms(model, dtype=np.complex128)
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)

    rng = np.random.default_rng(20260825)
    generic_cartesian = rng.uniform(-0.23, 0.23, size=(7, 3))
    generic = generic_cartesian @ reciprocal_inverse

    original_base = trace_curvature_vector(model, pos, generic)

    print(f"  {'test':22s} {'part':10s} {'rel.resid':>10s}  {'max_resid':>10s}  {'scale':>10s}")
    for test, (rotation, antiunitary) in OPERATIONS.items():
        transformed = transformed_kpoints(generic, reciprocal, rotation, antiunitary)
        evaluated = trace_curvature_vector(model, pos, transformed)
        factor = (-1.0 if antiunitary else 1.0) * np.linalg.det(rotation)
        for part, original_part in original_base.items():
            predicted = factor * original_part @ rotation.T
            residual = evaluated[part] - predicted
            max_residual = float(np.max(np.abs(residual)))
            scale = float(np.max(np.abs(np.concatenate((evaluated[part], predicted)))))
            rel = max_residual / max(scale, 1e-30)
            print(f"  {test:22s} {part:10s} {rel:10.4%}  {max_residual:10.6f}  {scale:10.6f}")

    # mirror_intersection: z-axis invariant line, omega should vanish identically
    q = np.linspace(-0.25, 0.25, 17)
    direction = np.array([0.0, 0.0, 1.0])
    line_k = (q[:, None] * direction) @ reciprocal_inverse
    base_z = trace_curvature_vector(model, pos, line_k)
    for part in ("internal", "cross", "external", "total"):
        value = base_z[part]
        scale = original_base[part]
        max_residual = float(np.max(np.abs(value)))
        rel = max_residual / max(float(np.max(np.abs(scale))), 1e-30)
        print(f"  {'base_mirror_intersection':22s} {part:10s} {rel:10.4%}  {max_residual:10.6f}"
              f"  {float(np.max(np.abs(scale))):10.6f}")
