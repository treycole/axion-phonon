#!/usr/bin/env python3
r"""Everything that produces the Wannier position matrix, and one door into it.

STAGE 1, in-house route.  Builds ``A(R)`` in Python from ``.chk`` and ``.mmn`` and
holds ``ws_pair_consistency``.  The current path does not need it: Wannier90 writes
``_r.dat`` directly (see README.md).  Kept for the ``.mmn`` route and the
validation notebooks.

``A(R)_{ij,a} = <0i| r_a |Rj>`` is the *only* input the external Berry-curvature
terms need beyond ``H(R)``.  Getting it wrong does not raise -- it produces a
plausible number -- so this module exists to keep every way of obtaining it in
one place, behind one entry point, with the trade-offs written down.

Use :func:`build_model_with_position_matrix`.  Do not call the per-source
builders directly unless you are adding a new source.

    from wannier.position_matrix import build_model_with_position_matrix

    model, result = build_model_with_position_matrix(
        w90, directory=run_dir, prefix="SrTiO3", source="mmn",
    )
    # model now carries H(R) and A(R) on one shared R-set;
    # result.diagnostics holds the numbers worth printing.

Then hand ``model`` to ``position_terms`` in ``modules/curvature.py``.

The sources
-----------

Every source below reproduces the ab-initio data *exactly* on the source
k-mesh.  They differ only *between* mesh points, because a finite ``N1 x N2 x
N3`` mesh fixes ``A(R)`` only up to its aliasing class.  That makes the choice
a choice of **interpolating function**, not a gauge -- see
``notes/log/2026-08-20_handoff_wannier_conventions.md``.  Two consequences run through this
whole module:

1. ``H(R)`` and ``A(R)`` must use the **same** R-set and the same
   representative rule.  The curvature formula assumes ``A(k)`` is the
   connection of the eigenstates of ``H(k)``; mixing two rules breaks that
   everywhere except on the source mesh.  Every builder here asserts it.
2. The representative rule must itself commute with the space group, or the
   interpolant silently breaks crystal symmetry.  ``min |R|`` does not;
   the pair-dependent ``min |R + tau_j - tau_i|`` does.  That is why the
   ``*_pair`` routes exist.

=================  =========================  ==================================
source             reads                      character
=================  =========================  ==================================
``mmn``            ``.chk`` + ``.mmn``        Production.  Sums the neighbour
                                              links into one finite difference,
                                              *then* maps to Wigner-Seitz.
                                              Right-orbital link phase
                                              ``exp(i b.tau_j)``.
``rfull``          ``_r_full.dat``            postw90 ``transl_inv_full=T``,
                                              ``use_ws_distance=T``.  Symmetric
                                              link phase, each link kept
                                              separate *through* the Fourier
                                              transform, mapped afterwards.
                                              A genuinely different
                                              discretisation, not a
                                              reimplementation of ``mmn``.
``rdat_pair``      ``_r.dat``                 Cheapest: strips Wannier90's
                                              scalar WS weights and applies the
                                              pair rule in house.  Needs no
                                              ``.mmn`` and no extra run.
``rdat_direct``    ``_r.dat``                 Legacy control only.  Wannier90's
                                              own ``min |R|`` R-set, no
                                              reassignment.  Cannot be paired
                                              with an orbital-dependent
                                              ``H(R)``; kept to reproduce old
                                              numbers.
``rdat_ndegen_``   ``_r.dat`` + ``_hr.dat``   Requires a run with the
``applied``        + ``_wsvec.dat``           ``write_ndegen_applied`` .win
                                              flag (jaemolihm/wannier90,
                                              ``plan5-write-ndegen-applied``):
                                              every WS weight is divided out
                                              at write time, so this route is
                                              ``rdat_direct`` with no Python
                                              remap needed -- and, unlike it,
                                              refuses a file that does not
                                              carry the flag.  H(R) verified
                                              to match the pair rule at
                                              machine precision; A(R) to
                                              relative_l2 ~1e-2 against
                                              ``mmn`` (2026-09-15).
=================  =========================  ==================================

``mmn`` and ``rfull`` do not commute at finite mesh density: summing the links
before mapping is not the same as mapping before summing.  Both are defensible;
``notes/log/2026-08-27_rfull_vs_mmn_position_matrix.md`` records why production keeps ``mmn``
(closer to DFPT, and its cache is 4.7x smaller than the ``_r_full.dat``
export).

Adding a source
---------------

Write a builder with the signature of :func:`_build_mmn` and add one entry to
:data:`POSITION_SOURCES`.  Nothing else in the codebase needs to change --
``validation/symmetry/check_convention_symmetry.ipynb`` will pick it up, and so
will anything else that goes through the entry point.

Caching
-------

The ``.mmn`` files are tens of GB per structure.  ``build_position_matrix_from_mmn``
writes a provenance-stamped ``.npz`` archive (~25x smaller) that reproduces
``A(R)`` bit-for-bit, and will use it after the ``.mmn`` has been deleted.  The
stamp covers the inputs *and* the options, so changing a convention invalidates
the cache rather than silently returning the old matrix.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from collections.abc import Mapping, Sequence

import numpy as np
from scipy.io import FortranFile

# One-way import: this module reads files through wannier_io, never the reverse.
# wannier_io.load_wannier_pair imports *us* from inside the function body.
from wannier.wannier_io import (
    _mmn_header,
    _mmn_records,
    _nnkp_block,
    _plain_or_gzip_path,
    check_mmn_compatibility,
    read_wannier_checkpoint,
)

__all__ = [
    "POSITION_SOURCES",
    "PositionMatrixResult",
    "apply_orbital_dependent_R_mapping",
    "apply_rdat_pair_convention",
    "build_model_with_position_matrix",
    "build_position_matrix_from_mmn",
    "connection_from_mmn",
    "finite_difference_weights",
    "install_postw90_position_matrix",
    "position_matrix_mmn_formula",
    "position_matrix_transl_inv_full",
    "read_postw90_position_matrix",
    "remap_A_R_to_centres",
    "ws_pair_consistency",
]
def _postw90_header_logical(header: str, keyword: str) -> bool:
    """Read one required ``keyword=T/F`` field from an exporter header."""

    match = re.search(
        rf"(?:^|;)\s*{re.escape(keyword)}\s*=\s*([TF])(?:\s*;|\s*$)",
        header,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise ValueError(
            f"postw90 position-matrix header has no {keyword}=T/F field"
        )
    return match.group(1).upper() == "T"


def _wsvec_write_ndegen_applied(directory: Path, prefix: str) -> bool | None:
    """Read the ``write_ndegen_applied=`` token Wannier90 stamps into
    ``seedname_wsvec.dat``'s header (added alongside the ``.win`` flag of the
    same name; see notes/log/2026-09-13_jaemo_transl_inv_full_ws_distance.md).

    Returns ``None`` -- "unknown", not "false" -- when the file is absent or
    predates the flag (a plain ``... with use_ws_distance=.true.`` header,
    no such token).
    """

    path = _plain_or_gzip_path(directory / f"{prefix}_wsvec.dat")
    if not path.is_file():
        return None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", errors="replace") as handle:
        header = handle.readline()
    match = re.search(
        r"write_ndegen_applied\s*=\s*\.?(true|false)\.?", header, re.IGNORECASE
    )
    if match is None:
        return None
    return match.group(1).lower() == "true"


def read_postw90_position_matrix(path: str | Path) -> dict:
    r"""Read ``seedname_r_full.dat`` without assuming ``_hr.dat`` weights.

    With ``use_ws_distance=T``, postw90 has already distributed each orbital
    pair among the representatives minimizing
    ``|R + tau_right - tau_left|``.  Those weights cannot be reconstructed
    from the scalar degeneracy list in ``_hr.dat``; the exporter therefore
    writes effective blocks and marks them with ``effective_ws_weights=T``.
    """

    path = Path(path)
    with path.open("r", encoding="ascii") as handle:
        header = handle.readline().strip()
        try:
            num_wann = int(handle.readline())
            num_R = int(handle.readline())
        except ValueError as error:
            raise ValueError(f"{path}: malformed num_wann/num_R header") from error
        if num_wann <= 0 or num_R <= 0:
            raise ValueError(f"{path}: num_wann and num_R must be positive")

        blocks: dict[tuple[int, int, int], np.ndarray] = {}
        rows_per_R = num_wann * num_wann
        for ir in range(num_R):
            lines = list(itertools.islice(handle, rows_per_R))
            if len(lines) != rows_per_R:
                raise ValueError(
                    f"{path}: R block {ir + 1}/{num_R} has "
                    f"{len(lines)} of {rows_per_R} rows"
                )
            values = np.fromstring("".join(lines), sep=" ")
            if values.size != rows_per_R * 11:
                raise ValueError(
                    f"{path}: R block {ir + 1}/{num_R} has "
                    f"{values.size} numeric fields; expected {rows_per_R * 11}"
                )
            values = values.reshape(rows_per_R, 11)
            integer_fields = np.rint(values[:, :5]).astype(int)
            if not np.array_equal(values[:, :5], integer_fields):
                raise ValueError(f"{path}: non-integer R/orbital index in block {ir + 1}")
            if not np.all(integer_fields[:, :3] == integer_fields[0, :3]):
                raise ValueError(f"{path}: mixed R vectors in block {ir + 1}")
            left = integer_fields[:, 3] - 1
            right = integer_fields[:, 4] - 1
            if (
                np.any(left < 0)
                or np.any(left >= num_wann)
                or np.any(right < 0)
                or np.any(right >= num_wann)
                or len(set(zip(left.tolist(), right.tolist()))) != rows_per_R
            ):
                raise ValueError(f"{path}: invalid or duplicate orbital indices")
            vector = tuple(int(value) for value in integer_fields[0, :3])
            if vector in blocks:
                raise ValueError(f"{path}: duplicate R vector {vector}")
            position = values[:, 5::2] + 1j * values[:, 6::2]
            block = np.empty((3, num_wann, num_wann), dtype=complex)
            block[:, left, right] = position.T
            blocks[vector] = block

        if any(line.strip() for line in handle):
            raise ValueError(f"{path}: unexpected data after {num_R} R blocks")

    transl_inv_full = _postw90_header_logical(header, "transl_inv_full")
    use_ws_distance = _postw90_header_logical(header, "use_ws_distance")
    effective_ws_weights = _postw90_header_logical(
        header, "effective_ws_weights"
    )
    if effective_ws_weights != use_ws_distance:
        raise ValueError(
            f"{path}: inconsistent use_ws_distance={use_ws_distance} and "
            f"effective_ws_weights={effective_ws_weights}"
        )
    return {
        "header": header,
        "num_wann": num_wann,
        "num_R": num_R,
        "transl_inv_full": transl_inv_full,
        "use_ws_distance": use_ws_distance,
        "effective_ws_weights": effective_ws_weights,
        "position_R": blocks,
    }


def install_postw90_position_matrix(
    model,
    path: str | Path,
    *,
    centers_cartesian: np.ndarray | None = None,
    tolerance: float = 1e-8,
    centres_tolerance: float = 1e-6,
) -> dict:
    """Install a pair-weighted postw90 export on an already pair-mapped model.

    ``tolerance`` bounds the dimensionless ``A(R) = A(-R)^dag`` relative defect.
    ``centres_tolerance`` bounds the ``R = 0`` diagonal against the checkpoint
    Wannier centres and carries units of Angstrom, so it cannot share a value
    with the former.  It matches the 1e-6 A used by the ``.mmn`` route in
    ``_position_matrix_for_model``: postw90 accumulates the centres from the
    ``transl_inv_full`` sum over neighbour shells while ``wannier90.x`` computes
    them from its own Marzari--Vanderbilt expression, so the two agree only to
    the summation-order noise of the k-sum (~1e-8 A here), not to machine
    precision.
    """

    data = read_postw90_position_matrix(path)
    if not data["transl_inv_full"]:
        raise ValueError(f"{path}: transl_inv_full must be T")
    if not data["use_ws_distance"] or not data["effective_ws_weights"]:
        raise ValueError(
            f"{path}: the production r_full path requires use_ws_distance=T "
            "and effective_ws_weights=T"
        )
    if int(data["num_wann"]) != model.norb:
        raise ValueError(
            f"{path}: num_wann={data['num_wann']} but model has {model.norb} orbitals"
        )
    if any(abs(float(value) - 1.0) > 1e-14 for value in model._ham_deg.values()):
        raise ValueError(
            "H(R) must receive orbital-dependent R mapping first, with all "
            "degeneracies equal to one"
        )

    position_R = data["position_R"]
    h_R = set(model._ham_r)
    a_R = set(position_R)
    h_only = h_R - a_R
    if h_only:
        raise ValueError(f"{path}: {len(h_only)} H(R) vectors have no A(R) block")
    sample_h = next(iter(model._ham_r.values()))
    for R in a_R - h_R:
        model._ham_r[R] = np.zeros_like(sample_h)
        model._ham_deg[R] = 1.0

    error_squared = 0.0
    norm_squared = 0.0
    maximum = 0.0
    for R, value in position_R.items():
        negative = tuple(-component for component in R)
        if negative not in position_R:
            raise ValueError(f"{path}: A(R) is not closed under R -> -R")
        defect = value - np.conj(np.transpose(position_R[negative], (0, 2, 1)))
        error_squared += float(np.vdot(defect, defect).real)
        norm_squared += float(np.vdot(value, value).real)
        maximum = max(maximum, float(np.abs(defect).max()))
    relative = np.sqrt(error_squared / norm_squared) if norm_squared else 0.0
    if not np.isfinite(relative) or relative > tolerance:
        raise ValueError(
            f"{path}: A(R)=A(-R)^dag relative defect {relative:.3e} "
            f"exceeds {tolerance:.3e}"
        )

    centre_error = float("nan")
    if centers_cartesian is not None:
        centers_cartesian = np.asarray(centers_cartesian, dtype=float)
        origin = position_R.get(_ZERO_R)
        if origin is None:
            raise ValueError(f"{path}: no R=(0,0,0) block")
        diagonal = np.arange(model.norb)
        exported = np.real(origin[:, diagonal, diagonal]).T
        centre_error = float(np.abs(exported - centers_cartesian).max())
        if centre_error > centres_tolerance:
            raise ValueError(
                f"{path}: R=0 diagonal differs from the Wannier centres by "
                f"{centre_error:.3e} Angstrom, above {centres_tolerance:.3e}"
            )

    model._pos_r = position_R
    return {
        **{key: value for key, value in data.items() if key != "position_R"},
        "n_R_pos": len(position_R),
        "zero_H_blocks_added": len(a_R - h_R),
        "pair_hermiticity_max": maximum,
        "pair_hermiticity_relative": float(relative),
        "centres_max_error": centre_error,
    }


def finite_difference_weights(b_cart: np.ndarray, tolerance: float = 1e-6):
    r"""Return ``w_b`` satisfying ``sum_b w_b b_a b_c = delta_ac``."""

    b_cart = np.asarray(b_cart, dtype=float)
    norms = np.linalg.norm(b_cart, axis=1)
    shells: list[float] = []
    labels = np.full(len(norms), -1, dtype=int)
    for index, norm in enumerate(norms):
        for shell, reference in enumerate(shells):
            if abs(norm - reference) < 1e-5:
                labels[index] = shell
                break
        else:
            labels[index] = len(shells)
            shells.append(norm)

    pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    equations = np.zeros((6, len(shells)))
    for shell in range(len(shells)):
        vectors = b_cart[labels == shell]
        for row, (axis_a, axis_b) in enumerate(pairs):
            equations[row, shell] = (
                vectors[:, axis_a] * vectors[:, axis_b]
            ).sum()
    target = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    shell_weights, *_ = np.linalg.lstsq(equations, target, rcond=None)
    if np.abs(equations @ shell_weights - target).max() > tolerance:
        raise ValueError("finite-difference b-shell completeness is not satisfied")
    return shell_weights[labels]


def connection_from_mmn(
    checkpoint: dict,
    path: str | Path,
    *,
    translational_invariant_MV: bool = False,
    reciprocal_lattice: np.ndarray | None = None,
) -> dict:
    r"""Build ``A_a(k)`` from MMN links with ``exp(i b.tau_right)``."""

    path = Path(path)
    stream = _mmn_records(path)
    num_bands, num_kpoints, num_neighbors = next(stream)
    for key, value in (
        ("num_bands", num_bands),
        ("num_kpts", num_kpoints),
        ("nntot", num_neighbors),
    ):
        if int(checkpoint[key]) != value:
            raise ValueError(
                f"{path}: {key}={value} disagrees with checkpoint {checkpoint[key]}"
            )

    kpoints = np.asarray(checkpoint["kpt_red"], dtype=float)
    reciprocal = np.asarray(
        checkpoint["recip_lattice"]
        if reciprocal_lattice is None
        else reciprocal_lattice,
        dtype=float,
    )
    V = np.asarray(checkpoint["v_matrix"])
    centers = np.asarray(checkpoint["wannier_centres"], dtype=float)
    num_wann = int(checkpoint["num_wann"])
    A = np.zeros((num_kpoints, num_wann, num_wann, 3), dtype=complex)
    diagonal = np.arange(num_wann)
    neighbors = np.empty((num_kpoints, num_neighbors), dtype=int)
    shifts = np.empty((num_kpoints, num_neighbors, 3), dtype=int)
    b_cart = np.empty((num_kpoints, num_neighbors, 3))
    weights = np.empty((num_kpoints, num_neighbors))
    completeness_error = 0.0

    shell = []
    for index, record in enumerate(stream):
        ik, ib = divmod(index, num_neighbors)
        neighbors[ik, ib], shifts[ik, ib] = _mmn_header(
            path, record, ik, num_kpoints
        )
        shell.append(record)
        if ib < num_neighbors - 1:
            continue

        b_reduced = kpoints[neighbors[ik]] + shifts[ik] - kpoints[ik]
        b_cart[ik] = b_reduced @ reciprocal
        weights[ik] = finite_difference_weights(b_cart[ik])
        completeness = np.einsum(
            "b,ba,bc->ac", weights[ik], b_cart[ik], b_cart[ik]
        )
        completeness_error = max(
            completeness_error,
            float(np.abs(completeness - np.eye(3)).max()),
        )

        V_left = V[ik].conj().T
        for neighbor_index, block in enumerate(shell):
            body = block[5:].reshape((num_bands, num_bands, 2))
            overlap = (body[..., 0] + 1j * body[..., 1]).T
            overlap_wannier = (
                V_left @ overlap @ V[neighbors[ik, neighbor_index]]
            )
            weight = weights[ik, neighbor_index]
            b_vector = b_cart[ik, neighbor_index]
            term = (
                1j
                * overlap_wannier[:, :, None]
                * weight
                * b_vector[None, None]
            )
            if translational_invariant_MV:
                term[diagonal, diagonal] = (
                    -np.log(overlap_wannier.diagonal()).imag[:, None]
                    * weight
                    * b_vector
                )
            term *= np.exp(1j * centers.dot(b_vector))[None, :, None]
            A[ik] += term
        shell = []

    return {
        "A": A,
        "neighbors": neighbors,
        "reciprocal_shifts": shifts,
        "b_cart": b_cart,
        "weights": weights,
        "completeness_error": completeness_error,
    }


_ZERO_R = (0, 0, 0)
_CACHE_SCHEMA = 3
_DEFAULT_WS_DISTANCE_TOLERANCE = 1e-3
_A_R_ALGORITHM = "inhouse-mmn-connection-v1"


def _map_R_by_orbital_positions(
    matrices: Mapping[tuple[int, int, int], np.ndarray],
    lattice_vectors: np.ndarray,
    centers_cartesian: np.ndarray,
    mp_grid: Sequence[int],
    tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
) -> dict[tuple[int, int, int], np.ndarray]:
    r"""Choose the representative minimizing ``|R+tau_right-tau_left|``.

    The matrix orbital indices must be the final two axes.  Equal-distance
    representatives are retained with equal weight.
    """

    lattice_vectors = np.asarray(lattice_vectors, dtype=float)
    centers_cartesian = np.asarray(centers_cartesian, dtype=float)
    num_orbitals = centers_cartesian.shape[0]
    input_R = list(matrices)
    R_cartesian = np.asarray(input_R, dtype=float) @ lattice_vectors
    images = np.array(
        [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)],
        dtype=float,
    )
    shifts_reduced = images * np.asarray(mp_grid, dtype=float)
    shifts_cartesian = shifts_reduced @ lattice_vectors

    R_index = {R: index for index, R in enumerate(input_R)}
    extra_R: list[tuple[int, int, int]] = []
    targets = np.empty((len(input_R), len(images)), dtype=int)
    for row, R in enumerate(input_R):
        for column, shift in enumerate(shifts_reduced):
            target_R = tuple(int(round(R[axis] + shift[axis])) for axis in range(3))
            if target_R not in R_index:
                R_index[target_R] = len(input_R) + len(extra_R)
                extra_R.append(target_R)
            targets[row, column] = R_index[target_R]

    output_R = input_R + extra_R
    sample = matrices[input_R[0]]
    leading_shape = sample.shape[:-2]
    output = np.zeros(
        (len(output_R),) + leading_shape + (num_orbitals, num_orbitals),
        dtype=complex,
    )
    stack = np.asarray([matrices[R] for R in input_R])
    for left in range(num_orbitals):
        for right in range(num_orbitals):
            displacement = centers_cartesian[right] - centers_cartesian[left]
            distances = np.linalg.norm(
                R_cartesian[:, None] + shifts_cartesian[None] + displacement,
                axis=2,
            )
            selected = distances <= distances.min(axis=1, keepdims=True) + tolerance
            multiplicity = selected.sum(axis=1)
            values = stack[..., left, right] / multiplicity.reshape(
                (-1,) + (1,) * len(leading_shape)
            )
            rows, columns = np.nonzero(selected)
            np.add.at(
                output,
                (targets[rows, columns], Ellipsis, left, right),
                values[rows],
            )

    keep = np.abs(output).reshape(len(output_R), -1).max(axis=1) > 0.0
    return {R: output[index] for index, R in enumerate(output_R) if keep[index]}


def apply_rdat_pair_convention(
    w90,
    mp_grid: Sequence[int],
    *,
    centers_cartesian: np.ndarray | None = None,
    R_distance_tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
    centres_tolerance: float = 1e-6,
    verbose: bool = True,
):
    r"""Put ``_r.dat``'s ``A(R)`` and ``_hr.dat``'s ``H(R)`` in one pair-dependent set.

    This is the ``use_ws_distance`` treatment applied to the position matrix
    Wannier90 already writes, with no extra files and no ``.mmn``.

    ``wannier90.x`` always writes ``_r.dat`` on the centre-*independent*
    ``min |R|`` Wigner--Seitz set; ``use_ws_distance = .true.`` does not change
    that file, it emits the pair-dependent remapping separately in
    ``seedname_wsvec.dat``.  The rule there is ``min |R + tau_j - tau_i|``,
    which is exactly ``_map_R_by_orbital_positions``, so the remapping can be
    done in house from ``_r.dat`` alone.

    Call this **before** ``W90.model`` and instead of
    ``apply_orbital_dependent_R_mapping``; it does the ``H(R)`` mapping too, so
    the two operators cannot end up in different conventions.  Returns
    ``(position_R, diagnostics)``; assign ``model._pos_r = position_R`` after
    building the model.
    """

    if w90.pos_r is None:
        raise ValueError(
            "This Wannier90 run has no position matrix: _r.dat / _tb.dat is "
            "absent, so there is nothing to remap."
        )
    degeneracies = {R: float(block["deg"]) for R, block in w90.ham_r.items()}
    if degeneracies and all(value == 1.0 for value in degeneracies.values()):
        raise ValueError(
            "Every H(R) degeneracy is already 1. Either (a) this W90 object "
            "has already been remapped -- call apply_rdat_pair_convention on "
            "a freshly read W90 object, remapping twice silently corrupts "
            "A(R) -- or (b) this run used the write_ndegen_applied .win "
            "flag, whose R-list Wannier90 already redistributed at write "
            "time in a structure this function does not expect; use "
            "source='rdat_ndegen_applied' instead, not this route."
        )
    if set(w90.pos_r) != set(w90.ham_r):
        raise ValueError("_r.dat and _hr.dat cover different R-vector sets.")

    if centers_cartesian is None:
        centers_cartesian = (
            np.asarray(w90.lattice.orb_vecs, dtype=float)
            @ np.asarray(w90.lat, dtype=float)
        )
    centers_cartesian = np.asarray(centers_cartesian, dtype=float)

    # Undo Wannier90's scalar Wigner--Seitz weights before re-distributing.
    raw = {
        R: np.asarray(block, dtype=complex) / degeneracies[R]
        for R, block in w90.pos_r.items()
    }
    centres_before = float(
        np.abs(
            np.real(np.einsum("ass->sa", raw[_ZERO_R])) - centers_cartesian
        ).max()
    )
    if centres_before > centres_tolerance:
        raise ValueError(
            f"_r.dat R=0 diagonal differs from the Wannier centres by "
            f"{centres_before:.3e} Angstrom, above {centres_tolerance:.3e}"
        )

    mapped = _map_R_by_orbital_positions(
        raw, w90.lat, centers_cartesian, mp_grid,
        tolerance=R_distance_tolerance,
    )

    # The remapping moves each (i, j) block independently, so A(R) and
    # A(-R)^dag land on partner vectors only if both are present.  Project onto
    # the Hermitian part, exactly as the .mmn pipeline does.
    sample = next(iter(mapped.values()))
    zero = np.zeros_like(sample)
    all_R = set(mapped)
    all_R.update(tuple(-component for component in R) for R in mapped)
    position_R = {
        R: 0.5
        * (
            mapped.get(R, zero)
            + np.conj(
                np.transpose(
                    mapped.get(tuple(-component for component in R), zero),
                    (0, 2, 1),
                )
            )
        )
        for R in sorted(all_R)
    }

    defect_squared = norm_squared = 0.0
    maximum = 0.0
    for R, value in mapped.items():
        partner = mapped.get(tuple(-component for component in R))
        if partner is None:
            maximum = float("inf")
            defect_squared = float("inf")
            break
        defect = value - np.conj(np.transpose(partner, (0, 2, 1)))
        defect_squared += float(np.vdot(defect, defect).real)
        norm_squared += float(np.vdot(value, value).real)
        maximum = max(maximum, float(np.abs(defect).max()))
    relative = (
        np.sqrt(defect_squared / norm_squared)
        if np.isfinite(defect_squared) and norm_squared
        else float("inf")
    )

    apply_orbital_dependent_R_mapping(
        w90, mp_grid, centers_cartesian, verbose=verbose
    )

    centres_after = float(
        np.abs(
            np.real(np.einsum("ass->sa", position_R[_ZERO_R])) - centers_cartesian
        ).max()
    )
    diagnostics = {
        "n_R_rdat": len(raw),
        "n_R_pos": len(position_R),
        "centres_vs_rdat_before": centres_before,
        "centres_vs_rdat_after": centres_after,
        "rdat_pair_hermiticity_max": maximum,
        "rdat_pair_hermiticity_relative": float(relative),
    }
    if verbose:
        print(
            f"    A(R) from _r.dat, pair-dependent remap: {len(raw)} -> "
            f"{len(position_R)} R-vectors; pre-projection Hermiticity rel "
            f"{relative:.2e}; centres error {centres_after:.2e} A"
        )
    return position_R, diagnostics


def apply_orbital_dependent_R_mapping(
    w90,
    mp_grid: Sequence[int],
    centers_cartesian: np.ndarray | None = None,
    *,
    verbose: bool = True,
):
    r"""Map ``H(R)`` in place using ``min |R+tau_right-tau_left|``.

    Call this before ``W90.model``.  The returned Hamiltonian blocks already
    include the old Wigner--Seitz weights, so all resulting degeneracies are 1.
    """

    if centers_cartesian is None:
        centers_cartesian = (
            np.asarray(w90.lattice.orb_vecs, dtype=float)
            @ np.asarray(w90.lat, dtype=float)
        )
    effective = {
        R: block["h"] / float(block["deg"]) for R, block in w90.ham_r.items()
    }
    mapped = _map_R_by_orbital_positions(
        effective, w90.lat, centers_cartesian, mp_grid
    )
    outside_weight = sum(
        np.linalg.norm(value) ** 2 for R, value in mapped.items() if R not in effective
    )
    total_weight = sum(np.linalg.norm(value) ** 2 for value in mapped.values())
    w90.ham_r = {R: {"h": value, "deg": 1} for R, value in mapped.items()}
    if verbose:
        moved = np.sqrt(outside_weight / total_weight) if total_weight else 0.0
        print(
            "    H(R) orbital-dependent mapping "
            f"(min |R + tau_right - tau_left|): {len(effective)} -> "
            f"{len(mapped)} R-vectors, {moved:.3%} of the weight moved "
            "outside the input set"
        )
    return w90


def ws_pair_consistency(
    base_ham_R: Mapping[tuple[int, int, int], np.ndarray],
    mode_ham_R: Mapping[tuple[int, int, int], np.ndarray],
    mp_grid: Sequence[int],
    *,
    significance: float = 1e-4,
    support: float = 1e-6,
) -> dict[str, float | int]:
    r"""Did two structures fold ``H(R)`` onto the same Wigner-Seitz images?

    A finite mesh fixes ``H(R)`` only up to the aliasing class
    ``[R] = {R + sum_i n_i N_i a_i}``; how a hopping is spread over the images
    in a class is a choice of interpolant.  With ``write_ndegen_applied``
    Wannier90 makes that choice itself, per run, from that run's own Wannier
    centres, and shares weight equally between images tied to within 1e-5 A.
    A symmetric base has exact ties.  A displaced structure has moved centres
    and generally breaks them, so a bond spread as (1/2, 1/2) over two images
    in the base can sit entirely on one of them in the mode.  Both are valid
    interpolants of the same mesh data, but ``(H_mode - H_base)/dbeta`` then
    contains the difference between two interpolants as well as the response.

    **``H(R)`` only.**  The images of an ``H(R)`` class are copies of one
    another, so class sums are the mesh-determined invariant.  The images of a
    ``transl_inv_full`` ``A(R)`` class are *not*: ``exp(-i b.R/2)`` is applied
    per b at the final lattice vector, so images differ by ``(-1)**n_i`` and a
    class sum is not meaningful (notes/log/2026-09-20_jaemo_ndegen_yio_mode1.md).
    Do not pass ``A(R)`` here, and do not redistribute ``A(R)`` from class sums.

    Key-set equality of the two ``R`` lists is *not* the right test: the sets
    legitimately differ (a hopping the displacement creates, or the
    ``min_hopping_norm`` cut).  Instead, for every class and orbital pair whose
    class-summed hopping exceeds ``significance`` in **both** structures --
    which excludes a hopping the displacement merely creates -- compare which
    images carry weight (``|H| > support``; the files are printed to six
    decimals).  Any difference there is a change of interpolant.

    Parameters
    ----------
    base_ham_R, mode_ham_R :
        ``{R: (n, n) array}`` with Wannier90's scalar degeneracy already
        divided out (``h / deg``; ``deg`` is 1 for ``write_ndegen_applied``).
    mp_grid :
        The ab-initio k-mesh, which defines the aliasing classes.
    significance :
        Class-summed ``|H|`` below which an entry is not compared (eV).
    support :
        ``|H|`` above which an image counts as carrying weight (eV).

    Returns
    -------
    dict with
    ``entries_compared``, ``entries_mismatched``, ``mismatch_fraction`` :
        How many (class, m, n) entries were compared and how many changed
        support.
    ``delta_H_contamination`` :
        ``||dH on mismatched entries||_F / ||dH||_F`` with
        ``dH = H_mode - H_base`` per image: the share of the finite-difference
        numerator that is a change of interpolant, not of the physics.
    ``base_weight_fraction`` :
        Share of ``||H_base||_F^2`` sitting on the mismatched entries.
    ``n_R_base``, ``n_R_mode``, ``n_R_common`` :
        Sizes of the two R lists, for context only.
    """

    grid = np.asarray(mp_grid, dtype=int)
    base_R, mode_R = set(base_ham_R), set(mode_ham_R)
    if not base_R or not mode_R:
        raise ValueError("Both structures need a non-empty H(R) to compare.")
    sample = np.asarray(next(iter(base_ham_R.values())))
    zero = np.zeros(sample.shape, dtype=complex)

    images_by_class: dict[tuple[int, ...], list[tuple[int, int, int]]] = {}
    for R in base_R | mode_R:
        key = tuple(int(x) for x in np.mod(np.asarray(R, dtype=int), grid))
        images_by_class.setdefault(key, []).append(R)

    compared = mismatched = 0
    delta_sq = delta_mismatched_sq = 0.0
    base_sq = base_mismatched_sq = 0.0
    for images in images_by_class.values():
        base = np.array([np.asarray(base_ham_R.get(R, zero)) for R in images])
        mode = np.array([np.asarray(mode_ham_R.get(R, zero)) for R in images])
        significant = (np.abs(base.sum(axis=0)) > significance) & (
            np.abs(mode.sum(axis=0)) > significance
        )
        changed = ((np.abs(base) > support) != (np.abs(mode) > support)).any(axis=0)
        bad = significant & changed
        delta = np.abs(mode - base) ** 2
        compared += int(significant.sum())
        mismatched += int(bad.sum())
        delta_sq += float(delta.sum())
        delta_mismatched_sq += float(delta[:, bad].sum())
        base_sq += float((np.abs(base) ** 2).sum())
        base_mismatched_sq += float((np.abs(base[:, bad]) ** 2).sum())

    return {
        "entries_compared": compared,
        "entries_mismatched": mismatched,
        "mismatch_fraction": mismatched / compared if compared else 0.0,
        "delta_H_contamination": (
            float(np.sqrt(delta_mismatched_sq / delta_sq)) if delta_sq else 0.0
        ),
        "base_weight_fraction": base_mismatched_sq / base_sq if base_sq else 0.0,
        "n_R_base": len(base_R),
        "n_R_mode": len(mode_R),
        "n_R_common": len(base_R & mode_R),
    }


def _without_input_paths(metadata: dict) -> dict:
    """Cache metadata with the recorded input paths dropped, for comparison.

    Each input is identified by its role (``chk``, ``mmn``, ...), size and
    mtime; moving the run directory (a rename keeps both) must not make a warm
    cache read as stale, since its ``.mmn`` is usually deleted and cannot be
    re-read to rebuild it.
    """
    inputs = {
        key: {field: value for field, value in stamp.items() if field != "path"}
        for key, stamp in metadata.get("inputs", {}).items()
    }
    return {**metadata, "inputs": inputs}


def _A_R_cache_metadata(
    checkpoint_directory: Path,
    mmn_directory: Path,
    prefix: str,
    options: dict,
    recorded: dict | None = None,
) -> dict:
    """Provenance stamp for one A(R) cache entry.

    ``recorded`` is the metadata a existing cache already carries.  When it is
    given, an input that no longer exists on disk keeps the stamp the cache
    recorded for it instead of raising.  That is what lets a warm cache be used
    after the multi-GB ``.mmn`` has been deleted: every input that *is* still
    present is compared byte-stamp for byte-stamp exactly as before, so a stale
    cache is still caught -- only the deleted file's own staleness becomes
    unverifiable, and ``build_position_matrix_from_mmn`` says so out loud.
    """

    def stamp(path: Path) -> dict:
        path = path.resolve()
        stat = path.stat()
        return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

    inputs = {
        "chk": _plain_or_gzip_path(checkpoint_directory / f"{prefix}.chk"),
        "mmn": _plain_or_gzip_path(mmn_directory / f"{prefix}.mmn"),
    }
    recorded_inputs = (recorded or {}).get("inputs", {})
    for label, path in (
        ("chk_nnkp", checkpoint_directory / f"{prefix}.nnkp"),
        ("mmn_nnkp", mmn_directory / f"{prefix}.nnkp"),
        ("mmn_eig", mmn_directory / f"{prefix}.eig"),
    ):
        # These are optional, but once a cache has recorded one, dropping it
        # silently would change the metadata and read as a stale cache.  Keep
        # the recorded stamp so deleting them is allowed too.
        if path.exists() or label in recorded_inputs:
            inputs[label] = path

    stamps = {}
    for key, path in sorted(inputs.items()):
        if path.exists():
            stamps[key] = stamp(path)
        elif key in recorded_inputs:
            stamps[key] = recorded_inputs[key]
        else:
            raise FileNotFoundError(
                f"Required Wannier90 input does not exist: {path}"
            )
    return {
        "schema": _CACHE_SCHEMA,
        "inputs": stamps,
        "options": options,
        "algorithm": _A_R_ALGORITHM,
    }


def _connection_to_A_R(
    checkpoint: dict,
    A_k: np.ndarray,
    distance_tolerance: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Fourier transform ``A(k)`` and map each orbital pair to its physical R."""

    mp_grid = tuple(int(value) for value in checkpoint["mp_grid"])
    num_kpoints = int(checkpoint["num_kpts"])
    if num_kpoints != int(np.prod(mp_grid)):
        raise ValueError(
            f"Checkpoint has {num_kpoints} k-points but mp_grid={mp_grid}"
        )
    grid_coordinates_float = (
        np.asarray(checkpoint["kpt_red"]) * np.asarray(mp_grid)[None]
    )
    grid_coordinates = np.rint(grid_coordinates_float).astype(int)
    if not np.allclose(grid_coordinates, grid_coordinates_float, atol=1e-7):
        raise ValueError("Checkpoint k-points are not on its Monkhorst-Pack grid")
    grid_coordinates %= np.asarray(mp_grid)
    if len({tuple(kpoint) for kpoint in grid_coordinates}) != num_kpoints:
        raise ValueError("Checkpoint k-point grid contains duplicates")
    if not np.any(np.all(grid_coordinates == 0, axis=1)):
        raise ValueError("The A(R) construction requires a Gamma-centered grid")

    num_wann = int(checkpoint["num_wann"])
    grid = np.zeros(mp_grid + (num_wann, num_wann, 3), dtype=complex)
    for index, kpoint in enumerate(grid_coordinates):
        grid[tuple(kpoint)] = A_k[index]
    del A_k
    grid_R = np.fft.fftn(grid, axes=(0, 1, 2)) / num_kpoints
    del grid

    # The mapping expects orbital indices last, whereas A(R) stores Cartesian last.
    raw = {
        R: grid_R[R].transpose(2, 0, 1) for R in np.ndindex(mp_grid)
    }
    mapped_first = _map_R_by_orbital_positions(
        raw,
        checkpoint["real_lattice"],
        checkpoint["wannier_centres"],
        mp_grid,
        tolerance=distance_tolerance,
    )
    del raw, grid_R
    return _hermitian_A_R(mapped_first, num_wann)


def _hermitian_A_R(
    mapped_first: Mapping[tuple[int, int, int], np.ndarray], num_wann: int
) -> tuple[np.ndarray, np.ndarray]:
    """Stack ``{R: (3, n, n)}`` as ``A(R)`` (Cartesian last), closed under
    ``R -> -R`` and paired as ``[A(R) + A(-R)^dag] / 2``."""

    mapped = {
        R: value.transpose(1, 2, 0) for R, value in mapped_first.items()
    }

    zero = np.zeros((num_wann, num_wann, 3), dtype=complex)
    all_R = set(mapped)
    all_R.update(tuple(-value for value in R) for R in mapped)
    R_vectors = np.array(sorted(all_R), dtype=int)
    A_R = np.empty((len(R_vectors), num_wann, num_wann, 3), dtype=complex)
    for index, vector in enumerate(R_vectors):
        R = tuple(int(value) for value in vector)
        negative_R = tuple(-value for value in R)
        A_R[index] = 0.5 * (
            mapped.get(R, zero)
            + mapped.get(negative_R, zero).swapaxes(0, 1).conj()
        )
    return A_R, R_vectors


def _nnkp_neighbours(nnkp: Path, num_kpoints: int) -> tuple[np.ndarray, np.ndarray]:
    """``(neighbour index, reciprocal shift G)`` per ``(k, b)`` from ``.nnkp``."""

    from wannier.wannier_io import _nnkp_block

    values = np.fromstring(_nnkp_block(nnkp, "nnkpts"), sep=" ", dtype=float)
    nntot = int(values[0])
    rows = values[1:].astype(int).reshape(num_kpoints, nntot, 5)
    if not np.array_equal(rows[:, :, 0], np.repeat(np.arange(1, num_kpoints + 1), nntot).reshape(num_kpoints, nntot)):
        raise ValueError(f"{nnkp}: nnkpts rows are not ordered by k-point")
    return rows[:, :, 1] - 1, rows[:, :, 2:]


def position_matrix_transl_inv_full(
    checkpoint: dict,
    nnkp: str | Path,
    *,
    ws_centres_cartesian: np.ndarray | None = None,
    R_distance_tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
    verbose: bool = True,
) -> dict[tuple[int, int, int], np.ndarray]:
    r"""``<0i|r|Rj>`` by Wannier90's ``transl_inv_full`` + ``write_ndegen_applied``
    formula, from the checkpoint's ``m_matrix``: ``{R: (3, J, J)}``, weights applied.

    A transcription of ``hamiltonian_get_rmn`` (Wannier90 ``src/hamiltonian.F90``,
    Jae-Mo Lihm's ``write_ndegen_applied`` branch), one b vector at a time:

    1. k-space half of the phase: ``M'_ij(k, b) = M_ij(k, b) exp(i b.(r_i + r_j)/2)``
       with this run's own centres ``r``;
    2. class sums ``C_b[R] = (1/N) sum_k M'(k, b) exp(-i k.R)`` on the mesh;
    3. each class onto its images ``min |R + r_j - r_i|`` (ties split equally);
    4. real-space half at the *final* image, ``i w_b b exp(-i b.R/2)``, then sum b;
    5. the ``R = 0`` diagonal set to the centres.

    Step 3 uses ``ws_centres_cartesian`` when given (the shared base centres of a
    finite difference, assumption A3'); otherwise the run's own centres, which is
    what Wannier90 does and what reproduces its ``_r.dat``.  Step 1 always uses the
    run's own centres: that dependence is continuous in the displacement, not a tie.
    Needs ``read_wannier_checkpoint(..., include_m_matrix=True)``.
    """

    if "m_matrix" not in checkpoint:
        raise ValueError("Read the checkpoint with include_m_matrix=True.")
    mp_grid = tuple(int(value) for value in checkpoint["mp_grid"])
    num_kpoints = int(checkpoint["num_kpts"])
    num_wann = int(checkpoint["num_wann"])
    kpoints = np.asarray(checkpoint["kpt_red"], dtype=float)
    lattice = np.asarray(checkpoint["real_lattice"], dtype=float)
    reciprocal = np.asarray(checkpoint["recip_lattice"], dtype=float)
    centres = np.asarray(checkpoint["wannier_centres"], dtype=float)
    ws_centres = centres if ws_centres_cartesian is None else np.asarray(
        ws_centres_cartesian, dtype=float
    )
    M = checkpoint["m_matrix"]

    neighbours, shifts = _nnkp_neighbours(Path(nnkp), num_kpoints)
    b_cart = (kpoints[neighbours] + shifts - kpoints[:, None, :]) @ reciprocal
    reference = b_cart[0]
    weights = finite_difference_weights(reference)
    # For each k, the column of m_matrix that carries reference b vector g
    # (Wannier90's nnord).
    order = np.empty((num_kpoints, len(reference)), dtype=int)
    for ik in range(num_kpoints):
        match = np.abs(b_cart[ik][None, :, :] - reference[:, None, :]).max(axis=2) < 1e-6
        if not np.all(match.sum(axis=1) == 1):
            raise ValueError(f"k-point {ik}: b vectors do not match those of k-point 0")
        order[ik] = match.argmax(axis=1)

    grid_float = kpoints * np.asarray(mp_grid)
    grid_index = np.rint(grid_float).astype(int) % np.asarray(mp_grid)
    if not np.allclose(grid_float, np.rint(grid_float), atol=1e-7):
        raise ValueError("Checkpoint k-points are not on its Monkhorst-Pack grid")

    total: dict[tuple[int, int, int], np.ndarray] = {}
    for g, b in enumerate(reference):
        half = np.exp(0.5j * centres @ b)
        grid = np.zeros(mp_grid + (num_wann, num_wann), dtype=complex)
        for ik in range(num_kpoints):
            grid[tuple(grid_index[ik])] = half[:, None] * M[ik, order[ik, g]] * half[None, :]
        grid = np.fft.fftn(grid, axes=(0, 1, 2)) / num_kpoints
        classes = {R: grid[R] for R in np.ndindex(mp_grid)}
        mapped = _map_R_by_orbital_positions(
            classes, lattice, ws_centres, mp_grid, tolerance=R_distance_tolerance
        )
        del classes, grid
        prefactor = 1j * weights[g] * b  # (3,)
        for R, value in mapped.items():
            phase = np.exp(-0.5j * b @ (np.asarray(R, dtype=float) @ lattice))
            term = (prefactor * phase)[:, None, None] * value[None]
            if R in total:
                total[R] += term
            else:
                total[R] = term
        del mapped
        if verbose:
            print(f"    b {g + 1}/{len(reference)} done ({len(total)} R vectors)")

    origin = total.setdefault(_ZERO_R, np.zeros((3, num_wann, num_wann), dtype=complex))
    diagonal = np.arange(num_wann)
    origin[:, diagonal, diagonal] = centres.T
    return total


def position_matrix_mmn_formula(
    checkpoint: dict,
    nnkp: str | Path,
    *,
    ws_centres_cartesian: np.ndarray | None = None,
    R_distance_tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
    verbose: bool = False,
) -> dict:
    r"""The production ``.mmn`` A(R), built from the checkpoint's ``m_matrix``.

    ``A(k) = sum_b i w_b b M(k, b) exp(i b.r_j)`` (right-hand centre), Fourier
    transform to class sums, fold, Hermitian pairing.  ``connection_from_mmn``
    forms ``M = V(k)^dag M_raw(k, b) V(k+b)`` from the ``.mmn``; the checkpoint
    stores exactly that rotated overlap as ``m_matrix``, so no ``.mmn`` is needed.
    Returns the cache layout (``AA``, ``iRvec``, ``wcc_cart``, ``real_lattice``)
    for ``_position_matrix_for_model``.  ``ws_centres_cartesian`` folds directly
    onto shared centres (exact: the images of a class are copies).

    One b vector at a time (mirroring ``position_matrix_transl_inv_full``),
    unlike the old implementation which summed all b's into a dense ``(k, i,
    j, 3)`` array before the Wigner-Seitz fold: that made the fold's ``(1097 ->
    ~13824 candidate images)`` expansion 3x bigger per call than it needs to be
    (an extra Cartesian axis carried through the expensive step), and, under
    ``ws_centres_cartesian``, ran it *twice* per structure (this function's own
    fold, then a separate ``remap_A_R_to_centres`` re-fold) -- 77 GB peak for a
    176-orbital YIO pair.  Folding per b before assembling the Cartesian sum
    avoids both: one ~6 GB call per b, 8 in total, freed between, and the given
    centres are used directly so no second fold is needed.  Checked against
    ``compare_mmn_formula_chk_vs_cache.py`` (own centres): reproduces the old
    implementation's agreement with the ``.mmn`` cache exactly (relative L2
    9.3e-7, matching the pre-refactor number to machine precision).
    """

    if "m_matrix" not in checkpoint:
        raise ValueError("Read the checkpoint with include_m_matrix=True.")
    mp_grid = tuple(int(value) for value in checkpoint["mp_grid"])
    num_kpoints = int(checkpoint["num_kpts"])
    num_wann = int(checkpoint["num_wann"])
    kpoints = np.asarray(checkpoint["kpt_red"], dtype=float)
    lattice = np.asarray(checkpoint["real_lattice"], dtype=float)
    reciprocal = np.asarray(checkpoint["recip_lattice"], dtype=float)
    centres = np.asarray(checkpoint["wannier_centres"], dtype=float)
    ws_centres = centres if ws_centres_cartesian is None else np.asarray(
        ws_centres_cartesian, dtype=float
    )
    M = checkpoint["m_matrix"]

    neighbours, shifts = _nnkp_neighbours(Path(nnkp), num_kpoints)
    b_cart = (kpoints[neighbours] + shifts - kpoints[:, None, :]) @ reciprocal
    reference = b_cart[0]
    weights = finite_difference_weights(reference)
    # For each k, the column of m_matrix that carries reference b vector g
    # (same matching as position_matrix_transl_inv_full; the star of b vectors
    # is the same at every k, just not necessarily in the same order).
    order = np.empty((num_kpoints, len(reference)), dtype=int)
    for ik in range(num_kpoints):
        match = np.abs(b_cart[ik][None, :, :] - reference[:, None, :]).max(axis=2) < 1e-6
        if not np.all(match.sum(axis=1) == 1):
            raise ValueError(f"k-point {ik}: b vectors do not match those of k-point 0")
        order[ik] = match.argmax(axis=1)

    grid_float = kpoints * np.asarray(mp_grid)
    grid_index = np.rint(grid_float).astype(int) % np.asarray(mp_grid)
    if not np.allclose(grid_float, np.rint(grid_float), atol=1e-7):
        raise ValueError("Checkpoint k-points are not on its Monkhorst-Pack grid")

    total: dict[tuple[int, int, int], np.ndarray] = {}
    for g, b in enumerate(reference):
        right_phase = np.exp(1j * centres @ b)  # exp(i b.r_j), right-centre link phase
        grid = np.zeros(mp_grid + (num_wann, num_wann), dtype=complex)
        for ik in range(num_kpoints):
            grid[tuple(grid_index[ik])] = M[ik, order[ik, g]] * right_phase[None, :]
        grid = np.fft.fftn(grid, axes=(0, 1, 2)) / num_kpoints
        classes = {R: grid[R] for R in np.ndindex(mp_grid)}
        mapped = _map_R_by_orbital_positions(
            classes, lattice, ws_centres, mp_grid, tolerance=R_distance_tolerance
        )
        del classes, grid
        prefactor = 1j * weights[g] * b  # (3,), no per-image phase in this formula
        for R, value in mapped.items():
            term = prefactor[:, None, None] * value[None]
            if R in total:
                total[R] = total[R] + term
            else:
                total[R] = term
        del mapped
        if verbose:
            print(f"    b {g + 1}/{len(reference)} done ({len(total)} R vectors)")

    A_R, R_vectors = _hermitian_A_R(total, num_wann)
    return {
        "AA": A_R,
        "iRvec": R_vectors,
        "wcc_cart": centres,
        "real_lattice": np.asarray(checkpoint["real_lattice"], dtype=float),
    }


def remap_A_R_to_centres(
    data: dict,
    centers_cartesian: np.ndarray,
    mp_grid: Sequence[int],
    tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
) -> dict:
    r"""Re-fold an ``.mmn`` ``A(R)`` onto the images that *other* centres choose.

    ``build_position_matrix_from_mmn`` folds each aliasing class
    ``[R] = {R + sum_i n_i N_i a_i}`` onto ``min |R + tau_t - tau_s|`` using the
    run's own Wannier centres, splitting ties equally.  In this route the link
    phase ``exp(i b.tau_right)`` is applied before the Fourier transform, so it
    does not depend on R and the images of a class are exact copies.  Summing a
    class and redistributing it with a different set of centres is therefore
    exact; ``_map_R_by_orbital_positions`` does both in one pass, because every
    input image of a class lands on the same new targets.

    Used to put the base and the displaced structure on **one** rule (assumption
    A3', notes on Wigner-Seitz ties): pass the same ``centers_cartesian`` for both.
    Not valid for a ``transl_inv_full`` ``A(R)`` (``rfull``, Jae-Mo's
    ``_r.dat``), whose images differ by ``(-1)**n_i``.  The ``R = 0`` diagonal is
    untouched (its bond is 0, so it is never tied); ``_position_matrix_for_model``
    still sets it to the run's own centres.
    """

    A_R = np.asarray(data["AA"])
    num_wann = A_R.shape[1]
    centers_cartesian = np.asarray(centers_cartesian, dtype=float)
    if centers_cartesian.shape != (num_wann, 3):
        raise ValueError(
            f"centers_cartesian has shape {centers_cartesian.shape}; "
            f"A(R) needs ({num_wann}, 3)."
        )
    matrices = {
        tuple(int(value) for value in R): A_R[index].transpose(2, 0, 1)
        for index, R in enumerate(data["iRvec"])
    }
    mapped_first = _map_R_by_orbital_positions(
        matrices, data["real_lattice"], centers_cartesian, mp_grid, tolerance
    )
    del matrices
    A_R_new, R_vectors = _hermitian_A_R(mapped_first, num_wann)
    return {**data, "AA": A_R_new, "iRvec": R_vectors}


def _position_matrix_for_model(model, data: dict, tolerance: float = 1e-6):
    """Convert cached A(R) to PythTB's weighted position-matrix convention."""

    A_R = np.asarray(data["AA"])
    R_vectors = [tuple(int(value) for value in R) for R in data["iRvec"]]
    centers = np.asarray(data["wcc_cart"])
    if A_R.shape[1] != model.norb:
        raise ValueError(
            f"num_wann mismatch: A(R) has {A_R.shape[1]}, model has {model.norb}."
        )
    if model._pos_r is None or model._ham_deg is None:
        raise ValueError("The model must carry both position and Hamiltonian blocks.")

    hamiltonian_R = set(model._ham_deg)
    position_support = set(R_vectors)
    missing = hamiltonian_R - position_support
    if missing:
        raise ValueError(
            f"{len(missing)} H(R) vectors are absent from the MMN-derived A(R)."
        )
    unpaired = {
        R for R in hamiltonian_R if tuple(-value for value in R) not in hamiltonian_R
    }
    if unpaired:
        raise ValueError(
            f"H(R) is not closed under R -> -R ({len(unpaired)} unpaired vectors)."
        )

    index_by_R = {R: index for index, R in enumerate(R_vectors)}
    diagonal = np.arange(model.norb)
    residual = np.real(
        np.einsum("ssa->sa", A_R[index_by_R[_ZERO_R]])
    )
    outside = [index for R, index in index_by_R.items() if R not in hamiltonian_R]
    total_norm_squared = float(np.vdot(A_R, A_R).real)
    outside_norm_squared = (
        float(np.vdot(A_R[outside], A_R[outside]).real) if outside else 0.0
    )
    outside_weight = (
        np.sqrt(outside_norm_squared / total_norm_squared)
        if total_norm_squared
        else 0.0
    )

    position_R = {}
    for R in sorted(position_support):
        block = A_R[index_by_R[R]].transpose(2, 0, 1).astype(complex)
        if R == _ZERO_R:
            block[:, diagonal, diagonal] = centers.T
        position_R[R] = block * float(model._ham_deg.get(R, 1.0))

    old_R0 = model._pos_r[_ZERO_R] / float(model._ham_deg.get(_ZERO_R, 1.0))
    centers_difference = np.abs(
        np.real(np.einsum("ass->sa", old_R0)) - centers
    ).max()
    if centers_difference > tolerance:
        raise ValueError(
            "Wannier centers from _r.dat and .chk differ by "
            f"{centers_difference:.2e} Angstrom."
        )

    def hermiticity(matrices: dict) -> tuple[float, float]:
        error_squared = 0.0
        norm_squared = 0.0
        maximum = 0.0
        for R, value in matrices.items():
            partner = matrices.get(tuple(-component for component in R))
            if partner is None:
                return float("inf"), float("inf")
            defect = value - np.conj(np.transpose(partner, (0, 2, 1)))
            error_squared += float(np.vdot(defect, defect).real)
            norm_squared += float(np.vdot(value, value).real)
            maximum = max(maximum, float(np.abs(defect).max()))
        relative = np.sqrt(error_squared / norm_squared) if norm_squared else 0.0
        return maximum, float(relative)

    rdat_max, rdat_relative = hermiticity(model._pos_r)
    mmn_max, mmn_relative = hermiticity(position_R)
    if not np.isfinite(mmn_relative) or mmn_relative > 1e-10:
        raise ValueError(
            "MMN-derived A(R) failed A(R)=A(-R)^dag: "
            f"relative defect {mmn_relative:.2e}, max {mmn_max:.2e} Angstrom."
        )

    # A finite Fourier grid can contain position-matrix blocks at R vectors
    # where the separately written H(R) file has no block.  Those A(R) terms
    # still contribute to the connection and must not be truncated to H's
    # support.  Add explicit zero Hamiltonian blocks so both the MMN and the
    # postw90 r_full installers represent that independent support identically.
    sample_h = next(iter(model._ham_r.values()))
    position_only_R = position_support - hamiltonian_R
    for R in position_only_R:
        model._ham_r[R] = np.zeros_like(sample_h)
        model._ham_deg[R] = 1.0

    diagnostics = {
        "n_R_pos": len(position_R),
        "n_R_mmn": len(R_vectors),
        "zero_H_blocks_added": len(position_only_R),
        "weight_outside_ham_R": float(outside_weight),
        "residual_R0_max": float(np.abs(residual).max()),
        "residual_R0_sum": residual.sum(axis=0),
        "centres_vs_rdat": float(centers_difference),
        "rdat_pair_hermiticity_max": rdat_max,
        "rdat_pair_hermiticity_relative": rdat_relative,
        "mmn_pair_hermiticity_max": mmn_max,
        "mmn_pair_hermiticity_relative": mmn_relative,
    }
    return position_R, diagnostics


def build_position_matrix_from_mmn(
    model,
    checkpoint_directory: str | Path,
    prefix: str,
    *,
    mmn_directory: str | Path | None = None,
    cache: str | Path | None = None,
    translational_invariant_JM: bool = False,
    translational_invariant_MV: bool = False,
    R_distance_tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
    ws_centres_cartesian: np.ndarray | None = None,
    mp_grid: Sequence[int] | None = None,
) -> dict:
    r"""Build orbital-aware ``A(R)`` from MMN and install it on ``model``.

    The MMN links receive ``exp(i b.tau_right)`` before the finite-difference
    sum.  The Fourier coefficients are then mapped with the same
    ``min |R+tau_right-tau_left|`` rule used for H(R).

    The cache always holds the fold onto the run's own centres.  Pass
    ``ws_centres_cartesian`` (and ``mp_grid``) to re-fold it onto other centres
    before installing -- see ``remap_A_R_to_centres``; H(R) must then be mapped
    with the same array.
    """
    if ws_centres_cartesian is not None and mp_grid is None:
        raise ValueError("ws_centres_cartesian needs mp_grid.")

    checkpoint_directory = Path(checkpoint_directory)
    mmn_directory = (
        Path(mmn_directory) if mmn_directory is not None else checkpoint_directory
    )
    if translational_invariant_JM and translational_invariant_MV:
        raise ValueError("The JM and MV schemes are mutually exclusive.")
    if translational_invariant_JM:
        raise NotImplementedError(
            "The JM scheme requires an R-dependent half-link phase."
        )
    if not np.isfinite(R_distance_tolerance) or R_distance_tolerance <= 0:
        raise ValueError("R_distance_tolerance must be positive and finite.")
    cache_path = Path(cache) if cache is not None else None
    absent = [
        path
        for path in (
            checkpoint_directory / f"{prefix}.chk",
            mmn_directory / f"{prefix}.mmn",
        )
        if not _plain_or_gzip_path(path).is_file()
    ]
    # A warm cache already holds the finished A(R); the raw inputs are then only
    # needed for provenance.  Allow them to have been deleted -- the .mmn is
    # tens of GB per structure and the cache is ~25x smaller.
    recorded = None
    if absent:
        if cache_path is None or not cache_path.exists():
            raise FileNotFoundError(
                "Required Wannier90 input does not exist: "
                + ", ".join(str(path) for path in absent)
            )
        try:
            with np.load(cache_path, allow_pickle=False) as archive:
                recorded = json.loads(str(archive["cache_meta"].item()))
        except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
            raise FileNotFoundError(
                f"{', '.join(str(path) for path in absent)} is missing and the "
                f"cache {cache_path} carries no usable provenance ({error})."
            ) from error
        print(
            "    "
            + ", ".join(path.name for path in absent)
            + f" absent; reusing {cache_path.name} on its recorded provenance "
            "(staleness of the deleted input cannot be re-verified)"
        )
    else:
        check_mmn_compatibility(checkpoint_directory, mmn_directory, prefix)

    # Keep these cache field names stable: the large production caches were
    # built with this provenance contract and should remain reusable.
    options = {
        "transl_inv_JM": bool(translational_invariant_JM),
        "transl_inv_MV": bool(translational_invariant_MV),
        "ws_dist_tol": float(R_distance_tolerance),
    }
    metadata = _A_R_cache_metadata(
        checkpoint_directory, mmn_directory, prefix, options, recorded
    )
    metadata_json = json.dumps(metadata, sort_keys=True, separators=(",", ":"))
    data = None
    if cache_path is not None and cache_path.exists():
        try:
            with np.load(cache_path, allow_pickle=False) as archive:
                if "cache_meta" not in archive.files:
                    reason = "legacy cache has no provenance metadata"
                elif _without_input_paths(
                    json.loads(str(archive["cache_meta"].item()))
                ) != _without_input_paths(metadata):
                    reason = "input, option, or algorithm provenance changed"
                else:
                    data = {
                        key: archive[key]
                        for key in archive.files
                        if key != "cache_meta"
                    }
                    reason = None
        except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
            reason = f"cache is unreadable ({error})"
        if data is None:
            print(f"    ignoring stale A(R) cache {cache_path.name}: {reason}")

    if data is None:
        checkpoint = read_wannier_checkpoint(
            checkpoint_directory / f"{prefix}.chk"
        )
        nnkp = checkpoint_directory / f"{prefix}.nnkp"
        reciprocal_lattice = None
        if nnkp.is_file():
            values = np.fromstring(_nnkp_block(nnkp, "recip_lattice"), sep=" ")
            if values.size != 9:
                raise ValueError(f"Expected a 3x3 recip_lattice block in {nnkp}")
            reciprocal_lattice = values.reshape(3, 3)
        connection = connection_from_mmn(
            checkpoint,
            mmn_directory / f"{prefix}.mmn",
            translational_invariant_MV=translational_invariant_MV,
            reciprocal_lattice=reciprocal_lattice,
        )
        A_R, R_vectors = _connection_to_A_R(
            checkpoint,
            connection.pop("A"),
            R_distance_tolerance,
        )
        if translational_invariant_MV:
            origin = np.flatnonzero(np.all(R_vectors == 0, axis=1))[0]
            diagonal = np.arange(int(checkpoint["num_wann"]))
            A_R[origin, diagonal, diagonal] = 0.0
        data = {
            "AA": A_R,
            "iRvec": R_vectors,
            "wcc_cart": np.asarray(checkpoint["wannier_centres"]),
            "real_lattice": np.asarray(checkpoint["real_lattice"]),
            "num_wann": np.int32(checkpoint["num_wann"]),
            "transl_inv_JM": np.int8(translational_invariant_JM),
            "transl_inv_MV": np.int8(translational_invariant_MV),
            "ws_dist_tol": np.float64(R_distance_tolerance),
            "bshell_completeness_error": np.float64(
                connection["completeness_error"]
            ),
        }
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{cache_path.name}.",
                suffix=".npz",
                dir=cache_path.parent,
            )
            os.close(descriptor)
            temporary_path = Path(temporary_name)
            try:
                np.savez_compressed(
                    temporary_path,
                    **data,
                    cache_meta=np.array(metadata_json),
                )
                os.replace(temporary_path, cache_path)
            except OSError as error:
                print(f"    WARNING: could not write A(R) cache: {error}")
            finally:
                if temporary_path.exists():
                    temporary_path.unlink()

    if ws_centres_cartesian is not None:
        data = remap_A_R_to_centres(
            data, ws_centres_cartesian, mp_grid, tolerance=R_distance_tolerance
        )
    position_R, diagnostics = _position_matrix_for_model(model, data)
    model._pos_r = position_R
    return {**data, **diagnostics}


# =========================================================================== #
# The entry point
# =========================================================================== #


@dataclass(frozen=True)
class PositionMatrixResult:
    """What a source did, and how well.

    ``diagnostics`` keys vary by source -- each builder documents its own -- so
    read it with ``.get()``.  The four fields above it are always present.
    """

    #: Registry key that produced this, e.g. ``"mmn"``.
    source: str
    #: Files actually read, for provenance in logs and result archives.
    reads: str
    #: The real-space representative rule, in words.  Two models can only be
    #: finite-differenced against each other if this string matches.
    ws_rule: str
    #: Number of R vectors shared by H(R) and A(R) after installation.
    n_R: int
    diagnostics: dict = field(default_factory=dict)


def _build_mmn(w90, model_options, *, directory, prefix, mmn_directory,
               cache, R_distance_tolerance, centres_cartesian, **_):
    """Production route: rebuild A(R) from the .chk/.mmn link overlaps.

    The R-mapping is applied to ``w90`` *before* the model is built, because
    ``w90.model()`` freezes H(R) at whatever R-set the W90 object holds.
    """
    mmn_directory = Path(mmn_directory) if mmn_directory else directory
    if mmn_directory != directory:
        # Borrowing a .mmn from another trial is legitimate -- it depends only
        # on the nscf, not on the projections -- but it must be checked.
        check_mmn_compatibility(directory, mmn_directory, prefix)
    mp_grid = read_wannier_checkpoint(directory / f"{prefix}.chk")["mp_grid"]
    apply_orbital_dependent_R_mapping(w90, mp_grid, verbose=False)
    model = w90.model(**model_options)
    info = build_position_matrix_from_mmn(
        model, directory, prefix,
        mmn_directory=mmn_directory,
        cache=cache,
        R_distance_tolerance=R_distance_tolerance,
    )
    return model, ".chk/.mmn", "min|R+tau_j-tau_i| (links summed first)", info


def _build_rfull(w90, model_options, *, directory, prefix, centres_cartesian, **_):
    """postw90 ``_r_full.dat`` written with transl_inv_full=T, use_ws_distance=T."""
    mp_grid = read_wannier_checkpoint(directory / f"{prefix}.chk")["mp_grid"]
    apply_orbital_dependent_R_mapping(w90, mp_grid, verbose=False)
    model = w90.model(**model_options)
    info = install_postw90_position_matrix(
        model, directory / f"{prefix}_r_full.dat",
        centers_cartesian=centres_cartesian(w90),
    )
    return model, "_r_full.dat", "postw90 min|R+tau_j-tau_i| (links kept separate)", info


def _build_rdat_pair(w90, model_options, *, directory, prefix,
                     R_distance_tolerance, centres_cartesian, **_):
    """Cheap route: ordinary _r.dat with the use_ws_distance rule applied here.

    Note the position matrix is assigned to ``model._pos_r`` *after* the model
    is built -- unlike the mmn/rfull routes, the reassignment is computed from
    ``w90`` directly rather than by remapping the W90 object first.
    """
    mp_grid = read_wannier_checkpoint(directory / f"{prefix}.chk")["mp_grid"]
    # The R-mapping always uses physical centres -- that is what makes it
    # point-group covariant.  The *embedding* (tau in H, v and the external
    # fields) is a separate choice, set through model_options["orb_vecs"].
    _position_attr(w90, "pos_r")  # fail with a readable message, not AttributeError
    position_R, info = apply_rdat_pair_convention(
        w90, mp_grid,
        centers_cartesian=centres_cartesian(w90),
        R_distance_tolerance=R_distance_tolerance,
        verbose=False,
    )
    model = w90.model(**model_options)
    model._pos_r = position_R
    return model, "_r.dat", "in-house min|R+tau_j-tau_i| (use_ws_distance rule)", info


def _build_rdat_direct(w90, model_options, *, directory, prefix, **_):
    """Legacy control: Wannier90's own min|R| R-set, no reassignment at all."""
    model = w90.model(**model_options)
    if _position_attr(model, "_pos_r") is None:
        raise ValueError(f"{directory} has no {prefix}_r.dat position matrix")
    return model, "_r.dat", "direct Wannier90 min|R|", {}


#: Angstrom.  ``_r.dat`` and ``_centres.xyz`` are both printed to six decimals,
#: so the honest floor is ~1e-6; a wrong b-vector ordering is off by >>0.01.
_NDEGEN_CENTRES_TOLERANCE = 1e-5


def _build_rdat_ndegen_applied(
    w90, model_options, *, directory, prefix, centres_cartesian, **_
):
    """``_r.dat``/``_hr.dat`` written with the ``write_ndegen_applied`` .win flag.

    Upstream (jaemolihm/wannier90, branch ``plan5-write-ndegen-applied``; see
    notes/log/2026-09-13_jaemo_transl_inv_full_ws_distance.md), this flag makes wannier90.x divide every
    Wigner-Seitz weight out *before* writing: ``_hr.dat``'s and ``_tb.dat``'s
    ndegen block becomes all 1s, and ``_r.dat`` -- which never had a ndegen
    block of its own -- becomes self-contained. ``O_mn(k) = sum_R exp(i
    k.R) O_mn(R)`` with no ndegen, no ``_wsvec.dat``. So unlike
    ``rdat_direct``, this route's ``w90.model()`` H(R)/A(R) are taken exactly
    as parsed, with no Python-side remap, and it *requires* the flag rather
    than silently assuming it (a stock run would silently give the old wrong
    ``rdat_direct`` numbers here).

    Validated on SrTiO3 Ti_Q_2A trial04 with ``use_ws_distance``,
    ``transl_inv_full`` and ``write_ndegen_applied`` all set
    (notes/log/2026-09-16_jaemo_rdat_production_closure.md): C4z/Mx covariance of
    Tr Omega at the ``mmn`` level, A(R) Hermiticity ~1e-12, and agreement with
    an independently built postw90 ``rfull`` export to ~1e-8.  It differs from
    ``mmn`` by the known rfull-vs-mmn interpolant gap (~6%/38%), which is a
    property of the two families, not a defect.  Raw H(R) matches the in-house
    pair remap to machine precision (2026-09-15).  Without ``transl_inv_full``
    the Hermiticity defect is ~7e-3 -- do not use such a run.

    The ``R = 0`` diagonal of A(R) is checked against the Wannier centres:
    ``<0n|r|0n>`` *is* the centre, so a wrong b-vector ordering (for example a
    ``.chk`` written by one Wannier90 version and re-read by another under
    ``restart = plot``) shows up there as an error of order Angstrom.
    ``info["centres_max_error"]`` reports the measured value.
    """
    flag = _wsvec_write_ndegen_applied(Path(directory), prefix)
    if flag is not True:
        found = (
            "no write_ndegen_applied= token"
            if flag is None
            else "write_ndegen_applied=.false."
        )
        raise ValueError(
            f"{directory}/{prefix}_wsvec.dat does not confirm "
            f"write_ndegen_applied=.true. (found {found}). This source "
            "needs a Wannier90 build with that .win flag "
            "(jaemolihm/wannier90, branch plan5-write-ndegen-applied); use "
            "source='rdat_pair' or 'mmn' for an ordinary run."
        )
    degeneracies = {R: float(block["deg"]) for R, block in w90.ham_r.items()}
    if degeneracies and not all(value == 1.0 for value in degeneracies.values()):
        raise ValueError(
            f"{directory}/{prefix}_hr.dat has non-trivial degeneracies even "
            "though write_ndegen_applied=.true.; inputs are inconsistent, "
            "refusing to guess which is right."
        )
    model = w90.model(**model_options)
    if _position_attr(model, "_pos_r") is None:
        raise ValueError(f"{directory} has no {prefix}_r.dat position matrix")

    pos_r = model._pos_r
    if _ZERO_R not in pos_r:
        raise ValueError(f"{directory}/{prefix}_r.dat has no R = 0 block")
    centres = np.asarray(centres_cartesian(w90), dtype=float)
    diagonal = np.real(np.einsum("ass->sa", np.asarray(pos_r[_ZERO_R])))
    centres_max_error = float(np.abs(diagonal - centres).max())
    if centres_max_error > _NDEGEN_CENTRES_TOLERANCE:
        raise ValueError(
            f"{directory}/{prefix}_r.dat R=0 diagonal differs from the "
            f"Wannier centres by {centres_max_error:.3e} Angstrom, above "
            f"{_NDEGEN_CENTRES_TOLERANCE:.1e}. <0n|r|0n> is the centre, so "
            "the position matrix was built with the wrong b-vector ordering "
            "or the wrong .chk."
        )

    defect_sq = norm_sq = 0.0
    maximum = 0.0
    for R, block in pos_r.items():
        partner = pos_r.get(tuple(-component for component in R))
        if partner is None:
            continue
        defect = block - np.conj(np.transpose(partner, (0, 2, 1)))
        defect_sq += float(np.vdot(defect, defect).real)
        norm_sq += float(np.vdot(block, block).real)
        maximum = max(maximum, float(np.abs(defect).max()))
    info = {
        "ndegen_applied_hermiticity_max": maximum,
        "ndegen_applied_hermiticity_relative": (
            float(np.sqrt(defect_sq / norm_sq)) if norm_sq else float("nan")
        ),
        "centres_max_error": centres_max_error,
    }
    return (
        model,
        "_r.dat/_hr.dat (write_ndegen_applied)",
        "min|R+tau_j-tau_i| (baked in at write time, no Python remap)",
        info,
    )


#: The registry.  Add a source by writing a builder above and one entry here.
#: ``requires`` is checked before anything is read, so a missing export fails
#: with a filename rather than deep inside a Fortran reader.
POSITION_SOURCES: dict[str, dict] = {
    "mmn": {
        "builder": _build_mmn,
        # Needed only on a cold cache.  A warm _AA_cache/*.npz reproduces A(R)
        # bit-for-bit with the (tens of GB) .mmn deleted, which is the whole
        # point of the cache -- so these are waived when one is supplied.
        "requires": (),
        "requires_cold": ("{prefix}.chk", "{prefix}.mmn", "{prefix}.nnkp"),
        "summary": "production; links summed, then mapped to Wigner-Seitz",
    },
    "rfull": {
        "builder": _build_rfull,
        "requires": ("{prefix}.chk", "{prefix}_r_full.dat"),
        "summary": "postw90 transl_inv_full; links kept separate through the FT",
    },
    "rdat_pair": {
        "builder": _build_rdat_pair,
        "requires": ("{prefix}.chk", "{prefix}_r.dat"),
        "summary": "cheap; _r.dat with the pair rule applied in house",
    },
    "rdat_direct": {
        "builder": _build_rdat_direct,
        "requires": ("{prefix}_r.dat",),
        "summary": "legacy control; Wannier90 min|R|, no reassignment",
    },
    "rdat_ndegen_applied": {
        "builder": _build_rdat_ndegen_applied,
        "requires": ("{prefix}_r.dat", "{prefix}_hr.dat", "{prefix}_wsvec.dat"),
        "summary": (
            "_r.dat/_hr.dat from a write_ndegen_applied=T run; no .mmn, "
            "no Python remap needed -- see notes/log/2026-09-13_jaemo_transl_inv_full_ws_distance.md"
        ),
    },
}


_PYTHTB_HINT = (
    "This build of pythtb has no Wannier position-matrix API. "
    "A(R) support lives on the pythtb 'external' branch (it adds W90.pos_r "
    "and TBModel._pos_r by reading prefix_r.dat); 'main' and the topic "
    "branches do not have it. Check `git -C <pythtb> branch --show-current`."
)


def _position_attr(obj, attr: str):
    """Read ``obj.<attr>``, or explain that pythtb is missing the feature.

    Without this the failure surfaces as a bare ``AttributeError`` on a
    private name, which reads like a bug in this repository rather than a
    checked-out pythtb branch that predates the position-matrix support.
    """
    try:
        return getattr(obj, attr)
    except AttributeError as error:
        raise AttributeError(
            f"{type(obj).__name__} has no {attr!r}. {_PYTHTB_HINT}"
        ) from error


def _default_centres(w90):
    """Cartesian Wannier centres of a W90 object, in Angstrom."""
    return (
        np.asarray(w90.lattice.orb_vecs, dtype=float)
        @ np.asarray(w90.lat, dtype=float)
    )


def build_model_with_position_matrix(
    w90,
    *,
    directory: str | Path,
    prefix: str,
    source: str = "mmn",
    model_options: Mapping[str, object] | None = None,
    mmn_directory: str | Path | None = None,
    cache: str | Path | None = None,
    R_distance_tolerance: float = _DEFAULT_WS_DISTANCE_TOLERANCE,
) -> tuple[object, PositionMatrixResult]:
    r"""Build a PythTB model carrying a consistent ``H(R)`` and ``A(R)``.

    This is the single entry point for the position matrix.  Every route ends
    with H(R) and A(R) on one shared R-set, which is asserted, not assumed.

    Parameters
    ----------
    w90 :
        A ``pythtb.W90`` for one run directory.  It is **mutated** by the
        ``mmn`` and ``rfull`` routes (their R-mapping has to happen before
        ``w90.model()`` is called), so pass a fresh one per source if you are
        comparing several.
    directory, prefix :
        The run directory and Wannier90 seedname.
    source :
        A key of :data:`POSITION_SOURCES`.  See the module docstring for what
        each one does and when it is the right choice.
    model_options :
        Passed to ``w90.model()``.  Pass ``orb_vecs=`` here to pin the
        convention-I embedding -- required when two structures are going to be
        finite-differenced, so that ``(H_mode - H_base)`` carries no gauge
        rotation.  Never set tau on ``w90.lattice``: ``TBModel.lattice``
        returns a copy, so that assignment is a silent no-op and would leave
        H and v at the Wannier centres while only the external terms moved.
    mmn_directory :
        Borrow the ``.mmn`` from another trial built on the same nscf.
        Compatibility is checked, not assumed.  ``mmn`` source only.
    cache :
        Provenance-stamped ``.npz`` archive of ``A(R)``.  Written if absent,
        reused if warm -- including after the ``.mmn`` has been deleted.
        ``mmn`` source only.
    R_distance_tolerance :
        Ties in ``|R + tau_j - tau_i|`` closer than this share the weight.
        Too tight and symmetry-equivalent images are split; the default is the
        production value.

    Returns
    -------
    (model, result)
        ``model`` is ready for ``curvature.position_terms``; ``result`` is a
        :class:`PositionMatrixResult`.

    Raises
    ------
    ValueError
        Unknown ``source``, or H(R) and A(R) ended on different R-sets.
    FileNotFoundError
        A file the source needs is absent.
    """
    if source not in POSITION_SOURCES:
        raise ValueError(
            f"Unknown position-matrix source {source!r}. "
            f"Available: {', '.join(sorted(POSITION_SOURCES))}."
        )
    directory = Path(directory)
    spec = POSITION_SOURCES[source]

    needed = list(spec["requires"])
    cache_is_warm = cache is not None and Path(cache).exists()
    if not cache_is_warm:
        needed += list(spec.get("requires_cold", ()))
    missing = [
        name for name in (p.format(prefix=prefix) for p in needed)
        if not _plain_or_gzip_path(directory / name).is_file()
    ]
    if missing:
        hint = "" if cache_is_warm else (
            "  (a warm cache= archive would waive the .chk/.mmn requirement)"
            if spec.get("requires_cold") else ""
        )
        raise FileNotFoundError(
            f"source={source!r} needs {', '.join(missing)} in {directory}{hint}"
        )

    model, reads, ws_rule, info = spec["builder"](
        w90,
        dict(model_options or {}),
        directory=directory,
        prefix=prefix,
        mmn_directory=mmn_directory,
        cache=cache,
        R_distance_tolerance=R_distance_tolerance,
        centres_cartesian=_default_centres,
    )

    # The invariant the whole external-curvature formula rests on.  A(k) must
    # be the connection of the eigenstates of H(k); if the two matrices live on
    # different R-sets that only holds on the source mesh.
    if _position_attr(model, "_pos_r") is None:
        raise ValueError(f"source={source!r} produced no position matrix")
    h_R, a_R = set(model._ham_r), set(model._pos_r)

    # H(R) without A(R) is a real failure: the position blocks are missing,
    # not merely small, and the external terms would silently drop them.
    if h_R - a_R:
        raise ValueError(
            f"source={source!r}: {len(h_R - a_R)} R vectors have H(R) but no "
            "A(R). The position blocks are missing, not merely small."
        )

    # The reverse is benign and expected.  Pair-dependent reassignment prunes
    # real-space blocks whose entire H(R) matrix is exactly zero, but an
    # occupied-only model can still carry nonzero A(R) there.  Restoring
    # explicit zero Hamiltonian blocks makes the two operators use literally
    # the same R set without changing H(k) by anything.
    zero_H_blocks_added = len(a_R - h_R)
    if zero_H_blocks_added:
        sample = next(iter(model._ham_r.values()))
        for R in a_R - h_R:
            model._ham_r[R] = np.zeros_like(sample)
            model._ham_deg[R] = 1.0
        if set(model._ham_r) != a_R:
            raise RuntimeError(f"source={source!r}: failed to align R sets")

    diagnostics = dict(info)
    diagnostics.setdefault("zero_H_blocks_added", zero_H_blocks_added)
    return model, PositionMatrixResult(
        source=source, reads=reads, ws_rule=ws_rule,
        n_R=len(a_R), diagnostics=diagnostics,
    )
