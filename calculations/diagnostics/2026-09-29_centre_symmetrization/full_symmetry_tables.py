#!/usr/bin/env python3
"""Full magnetic-point-group covariance tables for YIO base and modes 1/2/3,
on the production chk source (right_centre_from_chk) at the now-widened
WS-distance tolerance (1e-3 A, wannier/position_matrix.py default).

Operation lists + priming (unprimed/primed) come from
logs/point_groups.npz (copied from the job tmp dir), derived from the
actual DFT geometry (atom positions + converged Ir AIAO moments + each mode's
own displacement pattern) in derive_point_groups.py -- not from literature
recall. base: verified m-3m' (48). modes 1/2/3: verified identical -4'3m' (24)
(expected: all three are normal modes of the same GM2'=A2u irrep, so they
share one isotropy subgroup regardless of which sublattice carries them).

Run in the axion conda env: python full_symmetry_tables.py
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
DATA = REPO / "data"

from modules import curvature as cv
from wannier.wannier_io import load_wannier_pair

PREFIX = "Y2Ir2O7"
BASE = DATA / "Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
MODES = {
    "mode1": DATA / "Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/192345/trial_02_Y_p_Ir_d_O_sp/no-displacement/192908bc889",
    "mode2": DATA / "Y2Ir2O7/phonon/GM2-/mode2/Q_0.01A/output/192346/trial_02_Y_p_Ir_d_O_sp/no-displacement/193196bc889",
    "mode3": DATA / "Y2Ir2O7/phonon/GM2-/mode3/Q_0.01A/output/192347/trial_02_Y_p_Ir_d_O_sp/no-displacement/193215bc889",
}
N_OCCUPIED = 156
BATCH_SIZE = 4
OPTIONS = dict(zero_energy=0, min_hopping_norm=1e-5, max_distance=125, ignorable_imaginary_part=None)

GROUPS = np.load(Path(__file__).with_name("logs") / "point_groups.npz", allow_pickle=True)


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


def op_label(rotation, idx):
    det = np.linalg.det(rotation)
    tr = np.trace(rotation)
    kind = "proper" if det > 0 else "improper"
    return f"#{idx:02d} {kind} tr={tr:+.1f}"


def run_structure(label, model, pos, operations, primings):
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    reciprocal_inverse = np.linalg.inv(reciprocal)

    rng = np.random.default_rng(20260825)
    generic_cartesian = rng.uniform(-0.23, 0.23, size=(7, 3))
    generic = generic_cartesian @ reciprocal_inverse
    original = trace_curvature_vector(model, pos, generic)

    print(f"\n=== {label}: {len(operations)} operations ===")
    print(f"  {'op':24s} {'prime':8s} {'internal':>9s} {'cross':>9s} {'external':>9s} {'total':>9s}")
    rows = []
    for idx, (rotation, priming) in enumerate(zip(operations, primings)):
        antiunitary = (priming == "primed")
        transformed = transformed_kpoints(generic, reciprocal, rotation, antiunitary)
        evaluated = trace_curvature_vector(model, pos, transformed)
        factor = (-1.0 if antiunitary else 1.0) * np.linalg.det(rotation)
        rel = {}
        for part, original_part in original.items():
            predicted = factor * original_part @ rotation.T
            residual = evaluated[part] - predicted
            max_residual = float(np.max(np.abs(residual)))
            scale = float(np.max(np.abs(np.concatenate((evaluated[part], predicted)))))
            rel[part] = max_residual / max(scale, 1e-30)
        label_str = op_label(rotation, idx)
        print(f"  {label_str:24s} {priming:8s} {rel['internal']:8.3%} {rel['cross']:8.3%} "
              f"{rel['external']:8.3%} {rel['total']:8.3%}")
        rows.append((label_str, priming, rel))

    totals = np.array([r[2]["total"] for r in rows])
    print(f"  --- total residual across group: min={totals.min():.3%} median={np.median(totals):.3%} "
          f"max={totals.max():.3%} ---")
    return rows


print("Loading base...")
pair_m1 = load_wannier_pair(
    BASE, MODES["mode1"], prefix=PREFIX, cache_directory=DATA / "Y2Ir2O7/phonon/_AA_cache",
    position_source="right_centre_from_chk", shared_ws_centres=True, model_options=OPTIONS, verbose=False,
)
base_model = pair_m1.base_model
base_pos = cv.position_terms(base_model, dtype=np.complex128)
run_structure("base (m-3m')", base_model, base_pos, GROUPS["base_ops"], GROUPS["base_priming"])

mode1_pos = cv.position_terms(pair_m1.mode_model, dtype=np.complex128)
run_structure("mode1 (-4'3m')", pair_m1.mode_model, mode1_pos, GROUPS["mode1_ops"], GROUPS["mode1_priming"])
del pair_m1, mode1_pos

for mode in ("mode2", "mode3"):
    print(f"\nLoading {mode}...")
    pair = load_wannier_pair(
        BASE, MODES[mode], prefix=PREFIX, cache_directory=DATA / "Y2Ir2O7/phonon/_AA_cache",
        position_source="right_centre_from_chk", shared_ws_centres=True, model_options=OPTIONS, verbose=False,
    )
    pos = cv.position_terms(pair.mode_model, dtype=np.complex128)
    run_structure(f"{mode} (-4'3m')", pair.mode_model, pos, GROUPS[f"{mode}_ops"], GROUPS[f"{mode}_priming"])
    del pair, pos
