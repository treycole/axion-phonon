#!/usr/bin/env python3
"""YIO base full m-3m' covariance table (48 ops) at several WS tie tolerances, passed explicitly to
both the H(R) fold (_map_R_by_orbital_positions) and the A(R) fold (position_matrix_mmn_formula,
i.e. right_centre_from_chk). Same k-points and metric as full_symmetry_tables.py.

For base, every tolerance in [2.05e-3, 2.12e-2) A gives an identical fold (no image's excess
distance falls in that band), so 5e-3 stands for that whole physically allowed window.
"""

import sys
import time
from pathlib import Path

import numpy as np
from pythtb import W90

REPO = Path("/Users/treycole/Repos/axion-phonon")
sys.path.insert(0, str(REPO))

from modules import curvature as cv
from wannier.position_matrix import (
    _map_R_by_orbital_positions,
    _position_matrix_for_model,
    position_matrix_mmn_formula,
)
from wannier.wannier_io import read_wannier_checkpoint

TOLERANCES = (1e-5, 5e-3)
PREFIX = "Y2Ir2O7"
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
MP_GRID = (8, 8, 8)
N_OCCUPIED = 156
BATCH_SIZE = 4
OPTIONS = dict(zero_energy=0, min_hopping_norm=1e-5, max_distance=125, ignorable_imaginary_part=None)
GROUPS = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz",
                 allow_pickle=True)
OUT = Path(__file__).with_name("base_full_group_tolerances.npz")

import pythtb
print(f"pythtb: {pythtb.__file__}", flush=True)
checkpoint = read_wannier_checkpoint(BASE / f"{PREFIX}.chk", include_m_matrix=True)
centres = np.asarray(checkpoint["wannier_centres"], dtype=float)


def build(tol):
    w = W90(str(BASE), PREFIX)
    effective = {R: block["h"] / float(block["deg"]) for R, block in w.ham_r.items()}
    mapped = _map_R_by_orbital_positions(effective, w.lat, centres, MP_GRID, tolerance=tol)
    w.ham_r = {R: {"h": value, "deg": 1} for R, value in mapped.items()}
    model = w.model(orb_vecs=w.lattice.orb_vecs.copy(), **OPTIONS)
    data = position_matrix_mmn_formula(checkpoint, BASE / f"{PREFIX}.nnkp",
                                       ws_centres_cartesian=centres, R_distance_tolerance=tol)
    model._pos_r, _ = _position_matrix_for_model(model, data)
    return model, len(mapped), len(data["iRvec"])


def trace_curvature_vector(model, pos, kpoints):
    rinv = np.linalg.inv(np.asarray(model.recip_lat_vecs, dtype=float))
    pieces = {name: [] for name in ("internal", "cross", "external", "total")}
    for start in range(0, len(kpoints), BATCH_SIZE):
        c = cv.connection(cv.fields(model, pos, kpoints[start:start + BATCH_SIZE]), N_OCCUPIED)
        omega = {"internal": cv.omega_internal(c), "cross": cv.omega_cross(c), "external": cv.omega_external(c)}
        omega["total"] = omega["internal"] + omega["cross"] + omega["external"]
        for name in pieces:
            cart = np.einsum("ia,jb,ab...->ij...", rinv, rinv, omega[name], optimize=True)
            tr = np.trace(cart, axis1=-2, axis2=-1)
            pieces[name].append(np.stack((tr[1, 2], tr[2, 0], tr[0, 1]), axis=-1))
    return {name: np.concatenate(v, axis=0) for name, v in pieces.items()}


ops, primings = GROUPS["base_ops"], GROUPS["base_priming"]
parts = ("internal", "cross", "external", "total")
results = {}
for tol in TOLERANCES:
    t0 = time.time()
    model, nH, nA = build(tol)
    pos = cv.position_terms(model, dtype=np.complex128)
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    rinv = np.linalg.inv(reciprocal)
    generic = np.random.default_rng(20260825).uniform(-0.23, 0.23, size=(7, 3)) @ rinv
    original = trace_curvature_vector(model, pos, generic)
    table = np.zeros((len(ops), 4))
    print(f"\n=== tolerance {tol:.0e} A: H(R) {nH} R, A(R) {nA} R (built in {time.time() - t0:.0f} s) ===")
    print(f"  {'op':24s} {'prime':8s} {'internal':>9s} {'cross':>9s} {'external':>9s} {'total':>9s}", flush=True)
    for idx, (rotation, priming) in enumerate(zip(ops, primings)):
        anti = priming == "primed"
        kt = ((-1.0 if anti else 1.0) * (generic @ reciprocal) @ rotation.T) @ rinv
        evaluated = trace_curvature_vector(model, pos, kt)
        factor = (-1.0 if anti else 1.0) * np.linalg.det(rotation)
        for j, part in enumerate(parts):
            predicted = factor * original[part] @ rotation.T
            scale = float(np.max(np.abs(np.concatenate((evaluated[part], predicted)))))
            table[idx, j] = float(np.max(np.abs(evaluated[part] - predicted))) / max(scale, 1e-30)
        kind = "proper" if np.linalg.det(rotation) > 0 else "improper"
        label = f"#{idx:02d} {kind} tr={np.trace(rotation):+.1f}"
        print(f"  {label:24s} {priming:8s} " + " ".join(f"{v:8.3%}" for v in table[idx]), flush=True)
    for j, part in enumerate(parts):
        v = table[1:, j]
        print(f"  {part:9s} (excl. identity): min={v.min():.3%} median={np.median(v):.3%} max={v.max():.3%}")
    print(f"  done in {time.time() - t0:.0f} s", flush=True)
    results[f"tol_{tol:.0e}"] = table
    np.savez(OUT, **results)
    del model, pos
