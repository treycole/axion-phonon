"""Y2Ir2O7 production models obey their magnetic point groups (paper tab:symmetry-yio; port of
analysis/validation/symmetry/check_magnetic_symmetry.ipynb and of the script behind the table,
archive/diagnostics/2026-10-26_yio_symmetry_conv1e-10/symmetry_table_pinned.py). Base: m-3m' (47 operations besides
the identity); modes: -4'3m' (23). Operations and priming: yio_point_groups.npz (from QE's symmetry output).

- DFT eigenvalues on the mesh (.eig, lowest EIG_BANDS bands): e_n(g k) = e_n(k), k -> -g k for primed operations;
- interpolated bands at random k;
- the Cartesian band-traced curvature omega = (Tr Omega_yz, Tr Omega_zx, Tr Omega_xy) (AIAO breaks T, so it is not
  zero): omega(g k) = det(R) R omega(k), or omega(-g k) = -det(R) R omega(k) for primed operations, for the internal,
  cross, external parts and the total;
- control: a fourfold operation with its priming flipped (not a symmetry) must fail (paper: ~200%).

The base is the local production run (windows/f12.422_w22.644_k8/base, wannier90 220455); the mode runs' .chk are on
the cluster only (skipped here).
"""
from pathlib import Path

import numpy as np
import pytest

from axion_phonon import curvature

# ---- settings
GROUPS = Path(__file__).resolve().parent / "yio_point_groups.npz"
EIG_BANDS = 260               # lowest bands in the DFT test (paper)
N_GENERIC = 7
SEED = 20260825
BATCH = 4
THRESHOLDS = dict(            # maxima over the operations other than the identity (paper 10-26: DFT 1.4 ueV;
    dft_ueV=2.0,              # internal 0.000%, cross 0.011%, external 0.011%, total 0.000%)
    bands_eV=1e-6,
    internal=1e-5,
    cross=1e-3,
    external=1e-3,
    total=1e-4,
    control_min=0.5,
)
PARTS = ("internal", "cross", "external", "total")
# --------------
pytestmark = pytest.mark.data


def trace_vector(model, pos, k, reciprocal, n_occ):
    """omega (nk, 3), Cartesian, per part."""
    rinv = np.linalg.inv(reciprocal)
    out = {p: [] for p in PARTS}
    for s in range(0, len(k), BATCH):
        c = curvature.connection(curvature.fields(model, pos, k[s:s + BATCH]), n_occ)
        om = {"internal": curvature.omega_internal(c), "cross": curvature.omega_cross(c),
              "external": curvature.omega_external(c)}
        om["total"] = om["internal"] + om["cross"] + om["external"]
        for p in PARTS:
            tr = np.trace(np.einsum("ia,jb,ab...->ij...", rinv, rinv, om[p], optimize=True), axis1=-2, axis2=-1).real
            out[p].append(np.stack((tr[1, 2], tr[2, 0], tr[0, 1]), axis=-1))
    return {p: np.concatenate(v) for p, v in out.items()}


def transformed(K, reciprocal, op, anti):
    return ((-1.0 if anti else 1.0) * K @ op.T) @ np.linalg.inv(reciprocal)


def control_operation(ops, primed):
    """A fourfold operation (|trace| = 1, not +-E) with its priming flipped, checked not to be in the group."""
    i = next(i for i, R in enumerate(ops) if abs(abs(np.trace(R)) - 1) < 1e-9)
    flipped = not primed[i]
    assert not any(np.allclose(R, ops[i]) and p == flipped for R, p in zip(ops, primed)), "control is in the group"
    return ops[i], flipped


def symmetry_failures(run, model, n_occ, group):
    g = np.load(GROUPS, allow_pickle=True)
    ops, primed = g[f"{group}_ops"], np.array([str(p) == "primed" for p in g[f"{group}_priming"]])
    assert np.allclose(ops[0], np.eye(3)) and not primed[0]
    bad = []
    # DFT eigenvalues on the mesh
    grid = np.array(run["mp_grid"])
    kred, eig, reciprocal = run["kpts"], run["energies"][:, :EIG_BANDS], run["reciprocal"]
    index = {tuple(np.rint(k * grid).astype(int) % grid): i for i, k in enumerate(kred)}
    dft = 0.0
    for R, anti in zip(ops[1:], primed[1:]):
        n = np.rint(transformed(kred @ reciprocal, reciprocal, R, anti) * grid).astype(int)
        dft = max(dft, np.abs(eig[[index[tuple(v % grid)] for v in n]] - eig).max() * 1e6)
    if dft > THRESHOLDS["dft_ueV"]:
        bad.append(f"DFT eigenvalues {dft:.2f} ueV")
    # interpolated bands and traced curvature at random k
    K = np.random.default_rng(SEED).uniform(-0.23, 0.23, size=(N_GENERIC, 3))
    k0 = K @ np.linalg.inv(reciprocal)
    pos = curvature.position_terms(model, dtype=np.complex128)
    e0 = np.linalg.eigvalsh(model.hamiltonian(k0, flatten_spin_axis=True))
    w0 = trace_vector(model, pos, k0, reciprocal, n_occ)
    worst = dict.fromkeys(PARTS + ("bands_eV",), 0.0)
    cases = [(R, anti, False) for R, anti in zip(ops[1:], primed[1:])] + [(*control_operation(ops, primed), True)]
    for R, anti, is_control in cases:
        kt = transformed(K, reciprocal, R, anti)
        w = trace_vector(model, pos, kt, reciprocal, n_occ)
        factor = (-1.0 if anti else 1.0) * np.linalg.det(R)
        res = {}
        for p in PARTS:
            pred = factor * w0[p] @ R.T
            res[p] = np.abs(w[p] - pred).max() / max(np.abs(np.concatenate((w[p], pred))).max(), 1e-30)
        res["bands_eV"] = np.abs(np.linalg.eigvalsh(model.hamiltonian(kt, flatten_spin_axis=True)) - e0).max()
        if is_control:
            print(f"{group}: control (fourfold, priming flipped) total {res['total']:.2e}")
            if res["total"] < THRESHOLDS["control_min"]:
                bad.append(f"control does not fail: total {res['total']:.2e}")
        else:
            worst = {p: max(worst[p], res[p]) for p in worst}
    print(f"{group}: DFT {dft:.2f} ueV; max over {len(ops) - 1} operations: "
          + ", ".join(f"{p} {v:.2e}" for p, v in worst.items()))
    for p, v in worst.items():
        if v > THRESHOLDS[p]:
            bad.append(f"{p} {v:.2e} > {THRESHOLDS[p]}")
    return bad


def test_base(yio_base):
    bad = symmetry_failures(yio_base["run"], yio_base["model"], yio_base["n_occ"], "base")
    assert not bad, "Y2Ir2O7 base (m-3m'):\n  " + "\n  ".join(bad)
