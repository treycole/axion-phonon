#!/usr/bin/env python3
r"""Axion response d(theta)/d(beta): finite difference in beta, c2, BZ integral.

Stage 3 of three (see ``modules/curvature.py`` for the map).  Sections 8 and 9 of
``notes/berry_curvature_derivation.md`` (the "note"); equation numbers are its.

What it does
------------

The note's section 8.3 lists the three new ingredients of the phonon direction:
``dH/dbeta`` and ``dX/dbeta`` are **finite differences of two structures**
(two ``_hr.dat`` and two ``_r.dat``), and ``Lambda(R)`` is optional.  So:

1. ``fields`` (from ``curvature.py``) for the base and the displaced structure,
   at the *same* k-points and with the *same* nominal ``tau`` (assumption A3);
2. ``beta_terms``: ``dH/dbeta = (H_mode - H_base)/dbeta``, ``dA/dbeta`` likewise,
   and the mixed curl ``Omega_{l,beta} = -dA_l/dbeta`` in the frozen Wannier
   gauge (8.5), or ``dA_beta,l/dk - dA_l/dbeta`` when ``Lambda(R)`` is given (35);
3. ``curvature.connection`` and ``curvature.omega`` build the four-dimensional
   curvature ``B_{mu nu}`` (kx, ky, kz, beta);
4. ``c2_density``: ``c2 = (1/16 pi) eps^{mu nu rho sigma} tr[B_mu nu B_rho sigma]``
   (38); ``integrate_c2``: ``d(theta)/d(beta) = integral_BZ c2 d^3k``.

Two estimates, not an error bar
-------------------------------

The result has two entries, ``dtheta[0]`` and ``dtheta[1]``.  Both use the
**same** secant ``dH/dbeta`` and ``dA/dbeta``; they differ only in whether the
spatial fields (``H``, ``dH/dk``, ``A``, curl) are evaluated at beta = 0 (the
base) or at beta = dbeta (the mode).  There is no "right" endpoint: they agree
only if ``c2`` is linear in beta, and their spread measures how much ``c2``
varies across the step.  It is not a validity criterion.  The secant is the
derivative at the midpoint to O(dbeta^2), so the **mean** of the two is the
natural O(dbeta^2) estimate and half their difference is the leading curvature
term d2theta/dbeta2.  Only a sign flip between them is diagnostic: dtheta/dbeta
cannot reverse within one small step.

Use it yourself
---------------

::

    from pythtb import W90
    from modules.axion import dtheta

    base, mode = W90("base_dir", "seed"), W90("mode_dir", "seed")
    tau = base.lattice.orb_vecs                       # one fixed tau for both (A3)
    opts = dict(min_hopping_norm=1e-5)
    model_base = base.model(orb_vecs=tau, **opts)
    model_mode = mode.model(orb_vecs=tau, **opts)

    result = dtheta(model_base, model_mode, nk=12, dbeta=0.01, n_occ=156)
    result["dtheta"]          # (2,)  base-end and mode-end estimates, rad / Angstrom
    result["c2_density"]      # (nk, nk, nk, 2)
"""

from __future__ import annotations

import csv
import itertools
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from modules.curvature import (
    BetaTerms,
    Fields,
    PositionTerms,
    connection,
    fields,
    omega,
    position_terms,
)

__all__ = [
    "beta_terms",
    "c2_density",
    "dtheta",
    "integrate_c2",
    "parametric_connection",
    "prepare_lambda",
    "save_dtheta_results",
]


def _levi_civita_4() -> np.ndarray:
    """The rank-4 Levi-Civita symbol with eps[0, 1, 2, 3] = +1 (eps^{xyz beta} = +1)."""
    eps = np.zeros((4, 4, 4, 4))
    for perm in itertools.permutations(range(4)):
        inversions = sum(
            perm[i] > perm[j] for i in range(4) for j in range(i + 1, 4)
        )
        eps[perm] = -1.0 if inversions % 2 else 1.0
    return eps


_EPS4 = _levi_civita_4()


# --------------------------------------------------------------------------- #
# The fourth direction                                                        #
# --------------------------------------------------------------------------- #
def beta_terms(
    f_base: Fields,
    f_mode: Fields,
    dbeta: float,
    A_beta: np.ndarray | None = None,
    dA_beta: np.ndarray | None = None,
) -> BetaTerms:
    r"""``dH/dbeta``, ``dA/dbeta`` and the mixed curl from two structures (8.3).

    ``dH/dbeta = (H_mode - H_base)/dbeta`` and ``dA/dbeta`` likewise are the
    finite differences of two ``_hr.dat`` and two ``_r.dat``.

    Frozen Wannier gauge (default, section 8.5, ``A_beta = 0``): the mixed curl is
    ``Omega_{l,beta} = -dA_l/dbeta``.  Given the parametric connection ``A_beta``
    (33) and its k-derivative ``dA_beta`` (both ``(3, nk, J, J)`` for the
    derivative), the curl is ``d_l A_beta - d_beta A_l`` (35).
    """

    dA = (f_mode.A - f_base.A) / dbeta
    dH = (f_mode.H - f_base.H) / dbeta
    if A_beta is None:
        return BetaTerms(dH=dH, curl=-dA, A=None)
    return BetaTerms(dH=dH, curl=dA_beta - dA, A=A_beta)


def prepare_lambda(Lambda_R, dtype):
    """Impose Lambda(R) = -Lambda(-R)^dag (34) and flatten for the Bloch sum.

    ``Lambda_R`` is ``{R: Lambda(R)}``.  Pass the result to ``parametric_connection``.
    """
    raw = {}
    for vector, block in Lambda_R.items():
        raw[tuple(int(value) for value in vector)] = np.asarray(block, dtype=dtype)
    if not raw:
        raise ValueError("Lambda_R is empty; pass None for the frozen gauge.")
    for R in list(raw):
        raw.setdefault(tuple(-value for value in R), -raw[R].conj().T)
    ordered = sorted(raw)
    stack = np.stack(
        [
            0.5 * (raw[R] - raw[tuple(-value for value in R)].conj().T)
            for R in ordered
        ]
    )
    return np.asarray(ordered, dtype=float), stack.reshape(len(ordered), -1).astype(dtype)


def parametric_connection(pos: PositionTerms, prepared, k: np.ndarray):
    r"""``A_beta(k) = i sum_R e^{i k.b} Lambda(R)`` (33) and its k-derivative.

    ``prepared`` is the output of ``prepare_lambda``.  Returns ``(A_beta, dA_beta)``
    with shapes ``(nk, J, J)`` and ``(3, nk, J, J)``, ready for ``beta_terms``.
    """
    R, Lambda = prepared
    dtype = pos.dtype
    nk, J = len(k), pos.tau.shape[0]
    cell = np.exp(2j * np.pi * (k @ R.T)).astype(dtype)
    orbital = np.exp(2j * np.pi * (k @ pos.tau.T))
    bond = (np.conj(orbital)[:, :, None] * orbital[:, None, :]).astype(dtype)
    Lambda_k = (cell @ Lambda).reshape(nk, J, J)
    A_beta = 1j * Lambda_k * bond
    dA_beta = np.empty((3, nk, J, J), dtype=dtype)
    for axis in range(3):
        R_term = ((cell * R[None, :, axis]) @ Lambda).reshape(nk, J, J)
        tau_difference = pos.tau[None, :, axis] - pos.tau[:, None, axis]  # tau_t - tau_s
        dA_beta[axis] = -2.0 * np.pi * (R_term + tau_difference * Lambda_k) * bond
    return A_beta.astype(dtype, copy=False), dA_beta


# --------------------------------------------------------------------------- #
# c2 and its integral                                                         #
# --------------------------------------------------------------------------- #
def c2_density(Omega: np.ndarray) -> np.ndarray:
    r"""``c2 = (1/16 pi) eps^{mu nu rho sigma} tr[B_{mu nu} B_{rho sigma}]``, eq (38).

    ``Omega`` is ``(4, 4, nk, Nv, Nv)`` from ``curvature.omega``.  The full
    ``M x M`` matrix is kept until after the product: ``tr(BB) != tr(B) tr(B)``.
    Every surviving term has exactly one beta leg, so this equals the reduced
    form (39), ``(1/2 pi) sum_l tr[B~_l B_{l beta}]`` with ``B~_l = eps^{lij} B_ij / 2``.
    """
    return np.einsum("ijkl,ij...mn,kl...nm->...", _EPS4, Omega, Omega) / (16.0 * np.pi)


def integrate_c2(c2: np.ndarray, include_endpoint: bool = False):
    r"""``d(theta)/d(beta) = integral_BZ c2 d^3k`` on a uniform mesh.

    ``c2`` is ``(nk, nk, nk, ...)``.  Returns ``(sum rule, Simpson, d3k)``: the
    plain mesh sum times ``d3k`` (exact for a periodic integrand up to the
    aliasing of the mesh), and Simpson with the periodic endpoint appended.  The
    BZ is the unit cube of reduced coordinates, so the volume is 1.
    """
    nk = c2.shape[0]
    if include_endpoint:  # the mesh already contains the k = 1 face
        integration_grid = c2[:-1, :-1, :-1]
        d3k = 1.0 / (nk - 1) ** 3
        simpson_grid, dk = c2, 1.0 / (nk - 1)
    else:
        integration_grid = c2
        d3k = 1.0 / nk**3
        simpson_grid = np.concatenate([c2, c2[:1]], axis=0)
        simpson_grid = np.concatenate([simpson_grid, simpson_grid[:, :1]], axis=1)
        simpson_grid = np.concatenate([simpson_grid, simpson_grid[:, :, :1]], axis=2)
        dk = 1.0 / nk
    total = np.sum(integration_grid, axis=(0, 1, 2)) * d3k
    simpson_total = simpson(simpson_grid, dx=dk, axis=0)
    simpson_total = simpson(simpson_total, dx=dk, axis=0)
    simpson_total = simpson(simpson_total, dx=dk, axis=0)
    return total, simpson_total, d3k


# --------------------------------------------------------------------------- #
# The whole calculation                                                       #
# --------------------------------------------------------------------------- #
def dtheta(
    model_base,
    model_mode,
    nk: int,
    dbeta: float,
    n_occ: int,
    *,
    Lambda_R=None,
    positions: tuple[PositionTerms, PositionTerms] | None = None,
    include_endpoint: bool = False,
    batch: int = 100,
    progress=None,
) -> dict:
    r"""``d(theta)/d(beta)`` from a base and a displaced structure.

    Parameters
    ----------
    model_base, model_mode :
        PythTB models built with ``W90.model(orb_vecs=tau)`` from the base and
        displaced Wannier90 runs.  **``tau`` must be the same array** (A3); this
        is checked.
    nk :
        The k-mesh is ``nk x nk x nk`` in reduced coordinates.
    dbeta :
        The phonon amplitude between the two structures (Angstrom).
    n_occ :
        Number of occupied bands (bands ``0 .. n_occ-1``).
    Lambda_R :
        Optional ``{R: Lambda(R)}`` (section 8.4).  ``None`` selects the frozen
        Wannier gauge ``A_beta = 0`` -- an explicit choice; know which one a run
        is using.
    positions :
        Optional precomputed ``(position_terms(base), position_terms(mode))``,
        to reuse across a sweep in ``nk`` (each holds about 1.6 GB for J = 176).
    progress :
        Optional callable taking a message string.

    Returns
    -------
    dict with ``dtheta`` and ``dtheta_simpson`` (each ``(2,)``: base-end and
    mode-end estimates sharing one secant, see the module docstring),
    ``c2_density (nk,nk,nk,2)``, ``minimum_gap (2,)``, ``minimum_gap_kpoint
    (2,3)``, ``d3k`` and ``nk``.
    """

    log = progress if progress is not None else (lambda message: None)
    pos_base, pos_mode = (
        positions
        if positions is not None
        else (position_terms(model_base), position_terms(model_mode))
    )
    if pos_base.tau.shape != pos_mode.tau.shape or not np.allclose(
        pos_base.tau, pos_mode.tau, atol=1e-12
    ):
        raise ValueError(
            "Base and mode use different nominal tau.  Assumption A3 needs one "
            "fixed tau: build both models with W90.model(orb_vecs=tau)."
        )
    prepared = None if Lambda_R is None else prepare_lambda(Lambda_R, pos_base.dtype)

    k_mesh = model_base.k_uniform_mesh([nk, nk, nk], include_endpoints=include_endpoint)
    n_k = len(k_mesh)
    c2 = np.zeros((n_k, 2), dtype=np.complex64)
    minimum_gap = np.full(2, np.inf)
    minimum_gap_kpoint = np.full((2, 3), np.nan)

    for start in range(0, n_k, batch):
        stop = min(start + batch, n_k)
        k = k_mesh[start:stop]
        log(f"  k points {start} to {stop} of {n_k}")

        f_base = fields(model_base, pos_base, k)
        f_mode = fields(model_mode, pos_mode, k)
        if prepared is None:
            beta = beta_terms(f_base, f_mode, dbeta)
        else:
            A_beta, dA_beta = parametric_connection(pos_base, prepared, k)
            beta = beta_terms(f_base, f_mode, dbeta, A_beta, dA_beta)

        for endpoint, f in enumerate((f_base, f_mode)):
            c = connection(f, n_occ, beta)  # same secant at both endpoints
            c2[start:stop, endpoint] = c2_density(omega(c))
            gap = c.E[:, n_occ] - c.E[:, n_occ - 1]
            local = int(np.argmin(gap))
            if gap[local] < minimum_gap[endpoint]:
                minimum_gap[endpoint] = float(gap[local])
                minimum_gap_kpoint[endpoint] = k[local]

    c2_grid = c2.reshape((nk, nk, nk, 2))
    total, simpson_total, d3k = integrate_c2(c2_grid, include_endpoint)
    return {
        "nk": nk,
        "d3k": d3k,
        "dtheta": total,
        "dtheta_simpson": simpson_total,
        "minimum_gap": minimum_gap,
        "minimum_gap_kpoint": minimum_gap_kpoint,
        "c2_density": c2_grid,
    }


# --------------------------------------------------------------------------- #
# Saving                                                                      #
# --------------------------------------------------------------------------- #


def save_dtheta_results(
    output_directory: str | Path,
    results: Sequence[Mapping[str, object]],
    metadata: Mapping[str, object],
) -> tuple[Path, Path]:
    """Save completed mesh results and rebuild the NPZ/CSV sweep summaries."""

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    ordered = sorted(results, key=lambda result: int(result["nk"]))
    if not ordered:
        raise ValueError("No dtheta results to save.")

    for result in ordered:
        nk = int(result["nk"])
        np.savez_compressed(
            output_directory / f"nk_{nk}.npz",
            nk=np.int32(nk),
            d3k=np.float64(result["d3k"]),
            dtheta=result["dtheta"],
            dtheta_simpson=result["dtheta_simpson"],
            minimum_gap=result["minimum_gap"],
            minimum_gap_kpoint=result["minimum_gap_kpoint"],
            c2_density=result["c2_density"],
            # File-level aliases keep existing analysis notebooks readable.
            min_gap=result["minimum_gap"],
            min_gap_k=result["minimum_gap_kpoint"],
            chern2_density=result["c2_density"],
            **metadata,
        )

    nks = np.array([int(result["nk"]) for result in ordered], dtype=np.int32)
    dtheta = np.array([result["dtheta"] for result in ordered], dtype=np.complex64)
    dtheta_simpson = np.array(
        [result["dtheta_simpson"] for result in ordered], dtype=np.complex64
    )
    minimum_gap = np.array(
        [result["minimum_gap"] for result in ordered], dtype=np.float64
    )
    minimum_gap_kpoint = np.array(
        [result["minimum_gap_kpoint"] for result in ordered], dtype=np.float64
    )
    d3k = np.array([result["d3k"] for result in ordered], dtype=np.float64)

    summary_npz = output_directory / "summary.npz"
    np.savez_compressed(
        summary_npz,
        nks=nks,
        dtheta=dtheta,
        dtheta_simpson=dtheta_simpson,
        minimum_gap=minimum_gap,
        minimum_gap_kpoint=minimum_gap_kpoint,
        d3k=d3k,
        min_gap=minimum_gap,
        min_gap_k=minimum_gap_kpoint,
        **metadata,
    )

    summary_csv = output_directory / "summary.csv"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "nk",
                "dtheta_base_real",
                "dtheta_base_imag",
                "dtheta_mode_real",
                "dtheta_mode_imag",
                "dtheta_simpson_base_real",
                "dtheta_simpson_mode_real",
                "minimum_gap_base",
                "minimum_gap_mode",
                "minimum_gap_kpoint_base_x",
                "minimum_gap_kpoint_base_y",
                "minimum_gap_kpoint_base_z",
                "minimum_gap_kpoint_mode_x",
                "minimum_gap_kpoint_mode_y",
                "minimum_gap_kpoint_mode_z",
                "d3k",
            ]
        )
        for result in ordered:
            dtheta_endpoint = np.asarray(result["dtheta"])
            dtheta_simp = np.asarray(result["dtheta_simpson"])
            gaps = np.asarray(result["minimum_gap"])
            gap_k = np.asarray(result["minimum_gap_kpoint"])
            writer.writerow(
                [
                    int(result["nk"]),
                    dtheta_endpoint[0].real,
                    dtheta_endpoint[0].imag,
                    dtheta_endpoint[1].real,
                    dtheta_endpoint[1].imag,
                    dtheta_simp[0].real,
                    dtheta_simp[1].real,
                    gaps[0],
                    gaps[1],
                    *gap_k[0].tolist(),
                    *gap_k[1].tolist(),
                    result["d3k"],
                ]
            )

    return summary_npz, summary_csv
