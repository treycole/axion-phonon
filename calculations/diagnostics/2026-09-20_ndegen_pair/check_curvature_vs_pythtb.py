#!/usr/bin/env python3
"""modules/curvature.py against PythTB's own curvature, on a real structure.

PythTB's ``berry_curvature(non_abelian=True)`` is an independent implementation of the
same three k planes (validated against postw90).  On generic k-points:

    omega_internal(c)                  vs  berry_curvature(include_external=False)
    omega_cross(c) + omega_external(c) vs  berry_curvature(include_external=True)
                                           - berry_curvature(include_external=False)

PythTB returns cross and external summed (B^X + B^E), so that pair is compared as a sum.
Each code fixes its own eigenvector phases, and the matrices are in the occupied-band
basis, so only gauge-invariant quantities are compared: the band trace and the Frobenius
norm of every plane.

    python check_curvature_vs_pythtb.py [RUN_DIR] [N_OCC]
"""

import sys
from pathlib import Path

import numpy as np
from pythtb import W90

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))

from modules.curvature import (
    connection, fields, omega_cross, omega_external, omega_internal, position_terms,
)

DATA = Path(__file__).resolve().parents[3] / "data"

RUN = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA / "Y2Ir2O7/base/soc/u_3.0/output/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181538bc889" / "jaemo"
N_OCC = int(sys.argv[2]) if len(sys.argv) > 2 else 156
PLANES = ((0, 1), (1, 2), (2, 0))

w90 = W90(str(RUN), "Y2Ir2O7")
tau = w90.lattice.orb_vecs.copy()
model = w90.model(
    orb_vecs=tau, zero_energy=0, min_hopping_norm=1e-5, max_distance=125,
    ignorable_imaginary_part=None,
)
del w90

k = np.random.default_rng(0).uniform(0.0, 1.0, size=(6, 3))
occupied = np.arange(N_OCC)

# ours
c = connection(fields(model, position_terms(model, dtype=np.complex128), k), N_OCC)
ours = {"internal": omega_internal(c), "cross + external": omega_cross(c) + omega_external(c)}

# PythTB's
internal = model.berry_curvature(k, occ_idxs=occupied, non_abelian=True, include_external=False)
total = model.berry_curvature(k, occ_idxs=occupied, non_abelian=True, include_external=True)
theirs = {"internal": internal, "cross + external": total - internal}


def band_trace(x):
    return np.trace(x, axis1=-2, axis2=-1)


worst = 0.0
for name in ours:
    for mu, nu in PLANES:
        a, b = ours[name][mu, nu], theirs[name][mu, nu]  # (nk, Nv, Nv)
        trace_error = np.abs(band_trace(a) - band_trace(b)).max()
        norm_error = np.abs(np.linalg.norm(a, axis=(-2, -1)) - np.linalg.norm(b, axis=(-2, -1))).max()
        scale = max(np.abs(band_trace(b)).max(), np.linalg.norm(b, axis=(-2, -1)).max())
        relative = max(trace_error, norm_error) / scale
        worst = max(worst, relative)
        print(f"{name:17s} plane ({mu},{nu}):  max|trace diff| {trace_error:.2e}   "
              f"max|Frobenius diff| {norm_error:.2e}   relative to scale {scale:.2e}: {relative:.1e}")

print(f"\nworst relative disagreement over all planes: {worst:.1e}")
raise SystemExit(0 if worst < 1e-6 else 1)
