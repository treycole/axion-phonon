#!/usr/bin/env python3
"""Does widening the WS-distance tie tolerance reduce base's on-mesh inversion residual?

check_ws_tie_consistency.py found that _map_R_by_orbital_positions's closest-image
selection (tolerance 1e-5 A) is inconsistent under base's own inversion 14% of the
time specifically when the gap to the 2nd-closest candidate is < 0.001 A -- i.e.
near-degenerate ties that the tight tolerance fails to detect as ties, so one image
gets 100% of the weight instead of a symmetric 50/50 split.

This reruns check_base_onmesh_vs_offmesh.py's exact inversion-residual test (same
metric, same on-mesh/off-mesh point sets) at three tolerances: the production 1e-5,
and two widened values (1e-3, 1e-2) that should catch the gap<0.001 A cases (and,
at 1e-2, bleed into the next bucket which was already 0% inconsistent -- a control
for "did widening simply destroy real, well-separated distinctions instead").

H(R) mapping is called directly via _map_R_by_orbital_positions (bypassing
apply_orbital_dependent_R_mapping, which does not forward a tolerance) so this stays
a read-only diagnostic with no production-code changes.

Run in the axion conda env: python check_tie_tolerance_widened.py
"""

import sys
from pathlib import Path

import numpy as np
from pythtb import W90

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from modules import curvature as cv
from wannier.position_matrix import (
    _map_R_by_orbital_positions,
    _position_matrix_for_model,
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
TOLERANCES = (1e-5, 1e-3, 1e-2)


def build_model(formula, tolerance):
    w = W90(str(BASE_DIRECTORY), PREFIX)
    checkpoint = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk", include_m_matrix=True)
    centres = np.asarray(checkpoint["wannier_centres"], dtype=float)

    effective = {R: block["h"] / float(block["deg"]) for R, block in w.ham_r.items()}
    mapped = _map_R_by_orbital_positions(effective, w.lat, centres, MP_GRID, tolerance=tolerance)
    w.ham_r = {R: {"h": value, "deg": 1} for R, value in mapped.items()}

    embedding = w.lattice.orb_vecs.copy()
    model = w.model(orb_vecs=embedding, **OPTIONS)

    if formula == "transl_inv_full_from_chk":
        position = position_matrix_transl_inv_full(
            checkpoint, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centres,
            R_distance_tolerance=tolerance, verbose=False,
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
            checkpoint2, BASE_DIRECTORY / f"{PREFIX}.nnkp", ws_centres_cartesian=centres,
            R_distance_tolerance=tolerance,
        )
    model._pos_r, info = _position_matrix_for_model(model, data)
    return model, len(data["iRvec"])


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
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)
    k_cartesian = kpoints @ reciprocal
    transformed_cartesian = k_cartesian @ INVERSION.T
    transformed = transformed_cartesian @ reciprocal_inverse

    original = trace_curvature_vector(model, pos, kpoints)
    evaluated = trace_curvature_vector(model, pos, transformed)
    out = {}
    for part in ("internal", "cross", "external", "total"):
        predicted = np.linalg.det(INVERSION) * original[part] @ INVERSION.T
        residual = evaluated[part] - predicted
        scale = np.abs(np.concatenate((evaluated[part], predicted)))
        per_component = np.abs(residual).ravel() / np.maximum(scale.max(), 1e-30)
        out[part] = {
            "max": float(per_component.max()),
            "median": float(np.median(per_component)),
            "rms": float(np.sqrt(np.mean(per_component ** 2))),
        }
    return out


rng = np.random.default_rng(20260825)
generic_cartesian = rng.uniform(-0.23, 0.23, size=(60, 3))
rng2 = np.random.default_rng(20260930)
all_indices = np.array(np.meshgrid(*[np.arange(n) for n in MP_GRID], indexing="ij")).reshape(3, -1).T
chosen = rng2.choice(len(all_indices), size=60, replace=False)
mesh_reduced = all_indices[chosen] / np.asarray(MP_GRID)

for formula in ("transl_inv_full_from_chk", "right_centre_from_chk"):
    print(f"\n=== {formula} ===")
    for tol in TOLERANCES:
        model, nR = build_model(formula, tol)
        pos = cv.position_terms(model, dtype=np.complex128)
        reciprocal_inverse = np.linalg.inv(np.asarray(model.recip_lat_vecs, dtype=float))
        generic_reduced = generic_cartesian @ reciprocal_inverse

        off_mesh = inversion_residual(model, pos, generic_reduced)
        on_mesh = inversion_residual(model, pos, mesh_reduced)

        print(f"  tolerance={tol:.0e}  ({nR} R vectors)")
        print(f"    {'part':10s} {'on max':>9s} {'on med':>9s} {'on rms':>9s}"
              f"  {'off max':>9s} {'off med':>9s} {'off rms':>9s}")
        for part in ("internal", "cross", "external", "total"):
            o, f = on_mesh[part], off_mesh[part]
            print(f"    {part:10s} {o['max']:9.5f} {o['median']:9.5f} {o['rms']:9.5f}"
                  f"  {f['max']:9.5f} {f['median']:9.5f} {f['rms']:9.5f}")
