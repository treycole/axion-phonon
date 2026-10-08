#!/usr/bin/env python3
"""Does a denser ab-initio mesh (10x10x10 vs 8x8x8) reduce the YIO base symmetry floor?

Full m-3m' covariance table (48 ops) of the Cartesian trace curvature, split into
internal/cross/external/total, same k-points and metric as
2026-09-29_centre_symmetrization/full_symmetry_tables.py (whose base table is the 8x8x8
reference for the chk route: median total 1.563%).

Routes:
  chk   -- H(R) from _hr.dat folded on the .chk centres, A(R) = position_matrix_mmn_formula
           on the .chk m_matrix (= right_centre_from_chk), default WS tolerance.
  jaemo -- _hr.dat/_r.dat from the write_ndegen_applied run, read as written
           (build_model_with_position_matrix, source="rdat_ndegen_applied").

Run in the axion conda env: python base_symmetry_vs_mesh.py
"""

import gc
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
    build_model_with_position_matrix,
    position_matrix_mmn_formula,
)
from wannier.wannier_io import read_wannier_checkpoint

# ---- settings ----------------------------------------------------------------------------
PREFIX = "Y2Ir2O7"
NK8 = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
NK10 = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/mesh_10x10x10/202972bc889"
# (label, route, directory) in run order; jaemo directories hold the write_ndegen_applied files
RUNS = [
    ("nk10", "chk", NK10),
    ("nk10", "jaemo", NK10),
    ("nk8", "jaemo", NK8 / "jaemo"),
]
N_OCCUPIED = 156
BATCH_SIZE = 4
OPTIONS = dict(zero_energy=0, min_hopping_norm=1e-5, max_distance=125, ignorable_imaginary_part=None)
GROUPS = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz",
                 allow_pickle=True)
OUT = Path(__file__).with_name("base_symmetry_vs_mesh.npz")
# -------------------------------------------------------------------------------------------


def build_chk(directory):
    checkpoint = read_wannier_checkpoint(directory / f"{PREFIX}.chk", include_m_matrix=True)
    centres = np.asarray(checkpoint["wannier_centres"], dtype=float)
    mp_grid = tuple(int(n) for n in checkpoint["mp_grid"])
    w = W90(str(directory), PREFIX)
    effective = {R: block["h"] / float(block["deg"]) for R, block in w.ham_r.items()}
    mapped = _map_R_by_orbital_positions(effective, w.lat, centres, mp_grid)
    w.ham_r = {R: {"h": value, "deg": 1} for R, value in mapped.items()}
    model = w.model(orb_vecs=w.lattice.orb_vecs.copy(), **OPTIONS)
    data = position_matrix_mmn_formula(checkpoint, directory / f"{PREFIX}.nnkp", ws_centres_cartesian=centres)
    del checkpoint, w
    gc.collect()
    model._pos_r, _ = _position_matrix_for_model(model, data)
    return model, mp_grid


def build_jaemo(directory):
    model, result = build_model_with_position_matrix(
        W90(str(directory), PREFIX), directory=directory, prefix=PREFIX,
        source="rdat_ndegen_applied", model_options=OPTIONS)
    return model, result.n_R


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
        del c, omega
    return {name: np.concatenate(v, axis=0) for name, v in pieces.items()}


ops, primings = GROUPS["base_ops"], GROUPS["base_priming"]
parts = ("internal", "cross", "external", "total")
results = dict(np.load(OUT)) if OUT.exists() else {}
for label, route, directory in RUNS:
    t0 = time.time()
    model, info = build_chk(directory) if route == "chk" else build_jaemo(directory)
    pos = cv.position_terms(model, dtype=np.complex128)
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    rinv = np.linalg.inv(reciprocal)
    generic = np.random.default_rng(20260825).uniform(-0.23, 0.23, size=(7, 3)) @ rinv
    original = trace_curvature_vector(model, pos, generic)
    table = np.zeros((len(ops), 4))
    print(f"\n=== {label} {route}: {directory}  ({info}; H(R) {len(model._ham_r)} R, "
          f"built in {time.time() - t0:.0f} s) ===")
    print(f"  {'op':24s} {'prime':8s} {'internal':>9s} {'cross':>9s} {'external':>9s} {'total':>9s}", flush=True)
    bands = np.zeros(len(ops))
    e0 = np.linalg.eigvalsh(model.hamiltonian(generic, flatten_spin_axis=True))
    for idx, (rotation, priming) in enumerate(zip(ops, primings)):
        anti = priming == "primed"
        kt = ((-1.0 if anti else 1.0) * (generic @ reciprocal) @ rotation.T) @ rinv
        evaluated = trace_curvature_vector(model, pos, kt)
        factor = (-1.0 if anti else 1.0) * np.linalg.det(rotation)
        for j, part in enumerate(parts):
            predicted = factor * original[part] @ rotation.T
            scale = float(np.max(np.abs(np.concatenate((evaluated[part], predicted)))))
            table[idx, j] = float(np.max(np.abs(evaluated[part] - predicted))) / max(scale, 1e-30)
        bands[idx] = np.max(np.abs(np.linalg.eigvalsh(model.hamiltonian(kt, flatten_spin_axis=True)) - e0))
        kind = "proper" if np.linalg.det(rotation) > 0 else "improper"
        tag = f"#{idx:02d} {kind} tr={np.trace(rotation):+.1f}"
        print(f"  {tag:24s} {priming:8s} " + " ".join(f"{v:8.3%}" for v in table[idx])
              + f"   bands {bands[idx] * 1e3:6.2f} meV", flush=True)
    for j, part in enumerate(parts):
        v = table[1:, j]
        print(f"  {part:9s} (excl. identity): min={v.min():.3%} median={np.median(v):.3%} max={v.max():.3%}")
    print(f"  bands     (excl. identity): median={np.median(bands[1:]) * 1e3:.2f} meV "
          f"max={bands[1:].max() * 1e3:.2f} meV")
    print(f"  done in {time.time() - t0:.0f} s", flush=True)
    results[f"{label}_{route}"] = table
    results[f"{label}_{route}_bands_eV"] = bands
    np.savez(OUT, **results)
    del model, pos
    gc.collect()
