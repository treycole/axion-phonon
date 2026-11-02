r"""Symmetry tests of a Wannier model against the magnetic point group of its structure (the paper's symmetry tables).

The group is read from the structure's pw.x output (``io.qe.read_symmetry``: QE's Cartesian matrices and 'Time
Reversal' flags), never typed in. ``anti_translation`` adds every operation again with the opposite priming, for a
structure whose anti-translation {1'|t} QE does not list (MnBi2Te4 base: fractional translations are disabled for the
magnetic supercell).

An operation g (Cartesian R; primed = combined with time reversal) maps k -> R k, or k -> -R k when primed. Tested:

- DFT eigenvalues on the Wannier mesh: e_n(g k) = e_n(k)  (``dft_residuals``);
- interpolated bands: e_n(g k) = e_n(k) at off-mesh k  (``table``, column ``bands_meV``);
- the occupied-band curvature, Cartesian, as the dual vector F_a = eps_abc Omega_bc / 2 (an n_occ x n_occ matrix per
  component), through one of two gauge-invariant quantities (``table``, columns per part of eq (24)):
  * ``"vector"``: omega_a = Re Tr F_a (the Berry curvature, a pseudovector): omega(g k) = det(R) R omega(k),
    and omega(-R k) = -det(R) R omega(k) for a primed g (time reversal flips it);
  * ``"tensor"``: T_ab = Re Tr[F_a F_b], even in F: T(g k) = R T(k) R^T for every g. Used where PT forces Tr F = 0
    (MnBi2Te4), so the vector test would compare zero with zero.

Residuals are max |evaluated - predicted| over the k points and components, relative to the larger of the two
maxima, per part (internal, cross, external, total).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from . import curvature
from .io import qe

PARTS = ("internal", "cross", "external", "total")


def group(scf_out: str | Path, anti_translation: bool = False):
    """(operations (n, 3, 3) Cartesian, primed (n,) bool, names, listed (n,) bool: True if QE lists it, False if added
    as the product with the anti-translation)."""
    s = qe.read_symmetry(scf_out)
    ops, primed = s["cart"], s["time_reversal"]
    names = [n + (" T" if p else "") for n, p in zip(s["names"], primed)]
    if not (np.allclose(ops[0], np.eye(3)) and not primed[0]):
        raise ValueError(f"{scf_out}: the first operation is not the identity")
    listed = np.ones(len(ops), bool)
    if anti_translation:
        ops, primed = np.concatenate([ops, ops]), np.concatenate([primed, ~primed])
        names = names + [n + " {1'|t}" for n in names]
        listed = np.r_[listed, np.zeros(len(listed), bool)]
    keys = [(tuple(np.round(R, 6).ravel()), bool(p)) for R, p in zip(ops, primed)]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{scf_out}: repeated operations (anti_translation={anti_translation})")
    return ops, primed, names, listed


def controls(ops, primed, candidates):
    """The (name, R, primed) of ``candidates`` that are not in the group: operations that must fail."""
    return [c for c in candidates if not any(np.allclose(R, c[1]) and p == c[2] for R, p in zip(ops, primed))]


def to_reduced(K, reciprocal, R, primed):
    """Reduced k of R K, or -R K for a primed operation; K Cartesian rows (1/A), reciprocal rows b_i."""
    return ((-1.0 if primed else 1.0) * K @ R.T) @ np.linalg.inv(reciprocal)


def dft_residuals(run: dict, ops, primed, n_bands: int):
    """max |e_n(g k) - e_n(k)| (micro-eV) over the Wannier mesh for each operation: (n_ops, 2), the lowest n_bands
    bands and all bands. ``run`` from ``io.wannier90.read_run`` (kpts, energies, reciprocal, mp_grid)."""
    kred, eig, grid = run["kpts"], run["energies"], np.array(run["mp_grid"])
    index = {tuple(np.rint(k * grid).astype(int) % grid): i for i, k in enumerate(kred)}
    out = []
    for R, p in zip(ops, primed):
        n = to_reduced(kred @ run["reciprocal"], run["reciprocal"], R, p) * grid
        if not np.allclose(n, np.rint(n), atol=1e-6):
            raise ValueError("an operation does not map the Wannier mesh onto itself")
        d = np.abs(eig[[index[tuple(v % grid)] for v in np.rint(n).astype(int)]] - eig)
        out.append((d[:, :n_bands].max() * 1e6, d.max() * 1e6))
    return np.array(out)


def curvature_dual(model, pos, k, reciprocal, n_occ: int, batch: int = 4) -> dict:
    """F_a = (Omega_yz, Omega_zx, Omega_xy) of the occupied bands, Cartesian (A^2): (3, nk, n_occ, n_occ) per part
    of eq (24)."""
    rinv = np.linalg.inv(reciprocal)
    out = {p: [] for p in PARTS}
    for s in range(0, len(k), batch):
        c = curvature.connection(curvature.fields(model, pos, k[s:s + batch]), n_occ)
        om = {"internal": curvature.omega_internal(c), "cross": curvature.omega_cross(c),
              "external": curvature.omega_external(c)}
        om["total"] = om["internal"] + om["cross"] + om["external"]
        for p in PARTS:
            cart = np.einsum("ia,jb,ab...->ij...", rinv, rinv, om[p], optimize=True)    # reduced -> Cartesian
            out[p].append(np.stack((cart[1, 2], cart[2, 0], cart[0, 1])))
    return {p: np.concatenate(v, axis=1) for p, v in out.items()}


def invariant(F: np.ndarray, kind: str) -> np.ndarray:
    """``"vector"``: Re Tr F_a, (nk, 3); ``"tensor"``: Re Tr[F_a F_b], (nk, 3, 3)."""
    if kind == "vector":
        return np.einsum("akmm->ka", F).real
    if kind == "tensor":
        return np.einsum("akmn,bknm->kab", F, F, optimize=True).real
    raise ValueError(kind)


def predict(X0: np.ndarray, R: np.ndarray, primed: bool, kind: str) -> np.ndarray:
    """What the invariant at g k must be, given its value X0 at k."""
    if kind == "vector":
        return (-1.0 if primed else 1.0) * np.linalg.det(R) * X0 @ R.T
    return np.einsum("ab,kbc,dc->kad", R, X0, R)


def table(model, k, reciprocal, cases, n_occ: int, kind: str, batch: int = 4):
    """Rows (one per case (name, R, primed)) of [max band change (meV), residual internal, cross, external, total], and
    the invariant at k per part. The first case should be the identity (a row of zeros: the check of the bookkeeping)."""
    pos = curvature.position_terms(model, dtype=np.complex128)
    K = k @ reciprocal
    e0 = np.linalg.eigvalsh(model.hamiltonian(k, flatten_spin_axis=True))
    X0 = {p: invariant(F, kind) for p, F in curvature_dual(model, pos, k, reciprocal, n_occ, batch).items()}
    rows = []
    for name, R, p in cases:
        kt = to_reduced(K, reciprocal, R, p)
        de = np.abs(np.linalg.eigvalsh(model.hamiltonian(kt, flatten_spin_axis=True)) - e0).max() * 1e3
        X = {q: invariant(F, kind) for q, F in curvature_dual(model, pos, kt, reciprocal, n_occ, batch).items()}
        res = []
        for q in PARTS:
            pred = predict(X0[q], R, p, kind)
            res.append(np.abs(X[q] - pred).max() / max(np.abs(pred).max(), np.abs(X[q]).max(), 1e-300))
        rows.append([de] + res)
    return np.array(rows), X0


def random_k(reciprocal, n: int, half_width: float, rng):
    """n reduced k from Cartesian components uniform in [-half_width, half_width] (1/A); rng a numpy Generator."""
    return rng.uniform(-half_width, half_width, size=(n, 3)) @ np.linalg.inv(reciprocal)


def near_axis_k(reciprocal, n: int, radius: float, rng):
    """n reduced k at distance 0.3-1 x radius (1/A) from the axis b1 + b2 + b3 (MnBi2Te4: the k_z line, along c), at
    a random height between Gamma and its half (Z). Draws from rng after random_k, as the 10-22 run did."""
    axis = reciprocal.sum(0)
    zhat = axis / np.linalg.norm(axis)
    if np.abs(zhat[:2]).max() > 1e-9:
        raise ValueError("b1 + b2 + b3 is not along z")
    phi = rng.uniform(0, 2 * np.pi, n)
    perp = np.stack([np.cos(phi), np.sin(phi), np.zeros(n)], 1) * rng.uniform(0.3, 1, (n, 1)) * radius
    return (perp + rng.uniform(0, 1, (n, 1)) * axis / 2) @ np.linalg.inv(reciprocal)
