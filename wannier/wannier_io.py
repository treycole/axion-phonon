#!/usr/bin/env python3
"""I/O boundary for the Wannier-interpolated axion calculation.

STAGE 1, in-house route.  Python readers for Wannier90's ``.chk``, ``.mmn`` and
``.nnkp``, the ``Lambda(R)`` loader, and the older base/mode pair loader.  The
current path does not use them: Wannier90 writes ``_hr.dat`` and ``_r.dat``
(see README.md) and PythTB reads those.  Kept for the ``.mmn`` route and the
validation notebooks.

The driver should not know how Wannier90 files, ``A(R)`` caches, or result
archives are laid out.  It loads one ``WannierPair`` and saves a list of mesh
results through the two public functions in this module.

This module owns the **file readers** (``.chk``, ``.mmn``, ``.nnkp``,
``Lambda(R)``) and the base/mode pair loader.  Everything that *constructs*
the position matrix ``A(R)`` -- the four sources, the Wigner-Seitz
reassignment rules, and the ``A(R)`` cache -- lives in
``wannier/position_matrix.py``, behind one entry point.  Go there to add or
compare a position-matrix source.

The result writer ``save_dtheta_results`` lives in ``modules/axion.py``.
"""

from __future__ import annotations

import gzip
import hashlib
import itertools
import json
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Mapping, Sequence

import numpy as np
from pythtb import W90
from scipy.io import FortranFile

__all__ = [
    "WannierPair",
    "check_mmn_compatibility",
    "load_Lambda_R",
    "load_wannier_pair",
    "read_mmn_overlaps",
    "read_wannier_checkpoint",
]

# Wannier90 file readers -----------------------------------------------------


def _plain_or_gzip_path(path: str | Path) -> Path:
    """Return ``path`` or its ``.gz`` counterpart when only that exists."""

    path = Path(path)
    if path.is_file():
        return path
    compressed = Path(f"{path}.gz")
    return compressed if compressed.is_file() else path


class _SequentialFortranReader:
    """Minimal little-endian sequential-record reader for gzipped checkpoints."""

    def __init__(self, handle, path: Path):
        self._handle = handle
        self._path = path

    def read_record(self, dtype):
        # Records over 2 GiB (e.g. m_matrix on a 10x10x10 YIO mesh) are split into
        # subrecords; a negative leading marker means another subrecord follows.
        chunks = []
        while True:
            leading = self._handle.read(4)
            if len(leading) != 4:
                raise EOFError(f"{self._path}: truncated Fortran record header")
            marker = int(np.frombuffer(leading, dtype="<i4")[0])
            size = abs(marker)
            payload = self._handle.read(size)
            trailing = self._handle.read(4)
            if len(payload) != size or len(trailing) != 4:
                raise EOFError(f"{self._path}: truncated Fortran record")
            if abs(int(np.frombuffer(trailing, dtype="<i4")[0])) != size:
                raise ValueError(f"{self._path}: inconsistent Fortran record markers")
            chunks.append(payload)
            if marker >= 0:
                break
        if len(chunks) == 1:
            return np.frombuffer(chunks[0], dtype=np.dtype(dtype)).copy()
        buffer = np.empty(sum(len(chunk) for chunk in chunks), dtype=np.uint8)
        offset = 0
        while chunks:
            chunk = chunks.pop(0)
            buffer[offset:offset + len(chunk)] = np.frombuffer(chunk, dtype=np.uint8)
            offset += len(chunk)
        return buffer.view(np.dtype(dtype))

    def close(self):
        self._handle.close()


def read_wannier_checkpoint(path: str | Path, *, include_m_matrix: bool = False) -> dict:
    """Read the Wannier90 unformatted checkpoint needed by the MMN pipeline.

    ``include_m_matrix=True`` also returns ``m_matrix[k, b, i, j]``, the overlaps
    ``<u_ik|u_j,k+b>`` rotated into the Wannier gauge (``num_wann`` x ``num_wann``,
    ``b`` in the ``.nnkp`` neighbour order).  It is ~2 GB for the Y2Ir2O7 runs.
    """

    def read_string(handle) -> str:
        return b"".join(handle.read_record("c")).decode("ascii").strip()

    path = _plain_or_gzip_path(path)
    handle = (
        _SequentialFortranReader(gzip.open(path, "rb"), path)
        if path.suffix == ".gz"
        else FortranFile(str(path), "r")
    )
    read_int = lambda: handle.read_record("i4")
    read_float = lambda: handle.read_record("f8")

    def read_complex():
        values = read_float()
        return values[::2] + 1j * values[1::2]

    header = read_string(handle)
    num_bands = int(read_int()[0])
    _num_exclude = int(read_int()[0])
    exclude_bands = read_int()
    real_lattice = read_float().reshape((3, 3), order="F")
    reciprocal_lattice = read_float().reshape((3, 3), order="F")
    num_kpoints = int(read_int()[0])
    mp_grid = tuple(int(value) for value in read_int())
    kpoints = read_float().reshape((num_kpoints, 3))
    num_neighbors = int(read_int()[0])
    num_wann = int(read_int()[0])
    read_string(handle)  # checkpoint label
    disentangled = bool(read_int()[0])

    if disentangled:
        read_float()  # omega invariant
        lwindow = read_int().reshape((num_kpoints, num_bands)).astype(bool)
        ndimwin = read_int()
        U_dis = read_complex().reshape(
            (num_kpoints, num_wann, num_bands)
        ).swapaxes(1, 2)
    U = read_complex().reshape(
        (num_kpoints, num_wann, num_wann)
    ).swapaxes(1, 2)
    m_record = handle.read_record("f8")
    m_matrix = None
    if include_m_matrix:
        # Fortran m_matrix(i, j, nn, k), written as consecutive (re, im) pairs.
        # (re, im) pairs are already complex128's memory layout: view, don't copy (~4 GB at 10x10x10).
        m_matrix = m_record.view(np.complex128).reshape(
            (num_kpoints, num_neighbors, num_wann, num_wann)
        ).swapaxes(2, 3)
    del m_record
    centers = read_float().reshape((num_wann, 3))
    spreads = read_float().reshape((num_wann,))
    handle.close()

    if disentangled:
        V = np.zeros((num_kpoints, num_bands, num_wann), dtype=complex)
        for ik in range(num_kpoints):
            dimension = ndimwin[ik]
            V[ik, lwindow[ik]] = U_dis[ik][:dimension] @ U[ik]
    else:
        V = U

    lattice_product = real_lattice @ reciprocal_lattice.T / (2 * np.pi)
    if np.linalg.norm(lattice_product - np.eye(3)) > 1e-10:
        raise ValueError(f"{path}: real and reciprocal lattices are inconsistent")

    return {
        "header": header,
        "num_bands": num_bands,
        "num_wann": num_wann,
        "num_kpts": num_kpoints,
        "nntot": num_neighbors,
        "mp_grid": mp_grid,
        "exclude_bands": exclude_bands,
        "kpt_red": kpoints,
        "real_lattice": real_lattice,
        "recip_lattice": reciprocal_lattice,
        "v_matrix": V,
        "wannier_centres": centers,
        "wannier_spreads": spreads,
        **({"m_matrix": m_matrix} if include_m_matrix else {}),
    }


def _mmn_records(path: str | Path, chunk_bytes: int = 1 << 26):
    """Stream formatted MMN records without materializing the full file."""

    path = _plain_or_gzip_path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        handle.readline()
        try:
            num_bands, num_kpoints, num_neighbors = (
                int(value) for value in handle.readline().split()
            )
        except ValueError as error:
            raise ValueError(f"{path}: malformed .mmn dimensions") from error
        yield num_bands, num_kpoints, num_neighbors

        record_size = 5 + 2 * num_bands * num_bands
        pending = np.empty(0, dtype=np.float64)
        trailing_bytes = b""
        records_seen = 0

        def complete_records(values):
            nonlocal pending, records_seen
            pending = values if pending.size == 0 else np.concatenate((pending, values))
            count, remainder = divmod(pending.size, record_size)
            for index in range(count):
                yield pending[index * record_size : (index + 1) * record_size].copy()
            records_seen += count
            pending = pending[pending.size - remainder :].copy()

        while True:
            chunk = handle.read(chunk_bytes)
            if not chunk:
                break
            newline = chunk.rfind(b"\n")
            if newline < 0:
                trailing_bytes += chunk
                continue
            values = np.fromstring(trailing_bytes + chunk[:newline], sep=" ")
            trailing_bytes = chunk[newline + 1 :]
            yield from complete_records(values)
        if trailing_bytes.strip():
            yield from complete_records(np.fromstring(trailing_bytes, sep=" "))

    expected = num_kpoints * num_neighbors
    if pending.size or records_seen != expected:
        raise ValueError(
            f"{path}: expected {expected} MMN records, got {records_seen} "
            f"with {pending.size} values left over"
        )


def _mmn_header(
    path: str | Path, record: np.ndarray, ik: int, num_kpoints: int
) -> tuple[int, np.ndarray]:
    header = record[:5]
    rounded = np.rint(header)
    if not np.array_equal(header, rounded):
        raise ValueError(f"{path}: an MMN block header contains a non-integer")
    header = rounded.astype(int)
    if header[0] - 1 != ik:
        raise ValueError(f"{path}: MMN blocks are not in Wannier90 (ik,ib) order")
    neighbor = int(header[1]) - 1
    if not 0 <= neighbor < num_kpoints:
        raise ValueError(f"{path}: MMN neighbor index is out of range")
    return neighbor, header[2:5]


def read_mmn_overlaps(path: str | Path) -> dict:
    """Read a small MMN file, preserving link order and matrix orientation."""

    stream = _mmn_records(path)
    num_bands, num_kpoints, num_neighbors = next(stream)
    overlaps = np.empty(
        (num_kpoints, num_neighbors, num_bands, num_bands), dtype=complex
    )
    neighbors = np.empty((num_kpoints, num_neighbors), dtype=int)
    reciprocal_shifts = np.empty((num_kpoints, num_neighbors, 3), dtype=int)
    for index, record in enumerate(stream):
        ik, ib = divmod(index, num_neighbors)
        neighbors[ik, ib], reciprocal_shifts[ik, ib] = _mmn_header(
            path, record, ik, num_kpoints
        )
        body = record[5:].reshape((num_bands, num_bands, 2))
        overlaps[ik, ib] = (body[..., 0] + 1j * body[..., 1]).T
    return {
        "num_bands": num_bands,
        "num_kpts": num_kpoints,
        "nntot": num_neighbors,
        "M": overlaps,
        "neigh": neighbors,
        "G": reciprocal_shifts,
    }


# Orbital-position-aware real-space representation -------------------------


def load_Lambda_R(path: str | Path) -> dict[tuple[int, int, int], np.ndarray]:
    r"""Load ``Lambda(R)=<0s|d/dbeta|Rt>`` from an NPZ archive."""

    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        if "R" not in data or "Lambda" not in data:
            raise ValueError(f"{path} must contain arrays named 'R' and 'Lambda'.")
        R_vectors = np.asarray(data["R"])
        blocks = np.asarray(data["Lambda"])
    if (
        R_vectors.ndim != 2
        or R_vectors.shape[1] != 3
        or not np.allclose(R_vectors, np.rint(R_vectors))
    ):
        raise ValueError("Lambda R vectors must be integers with shape (n_R, 3).")
    if (
        blocks.ndim != 3
        or blocks.shape[0] != R_vectors.shape[0]
        or blocks.shape[1] != blocks.shape[2]
    ):
        raise ValueError("Lambda must have shape (n_R, J, J).")
    result = {}
    for vector, block in zip(np.rint(R_vectors).astype(int), blocks):
        R = tuple(int(value) for value in vector)
        if R in result:
            raise ValueError(f"Duplicate Lambda block for R={R} in {path}.")
        result[R] = block
    return result


def _nnkp_block(path: Path, name: str) -> str:
    text = path.read_text(errors="ignore")
    match = re.search(rf"begin\s+{name}(.*?)end\s+{name}", text, re.I | re.S)
    if match is None:
        raise ValueError(f"No '{name}' block in {path}")
    return "\n".join(
        line.strip() for line in match.group(1).splitlines() if line.strip()
    )


def check_mmn_compatibility(
    checkpoint_directory: str | Path,
    mmn_directory: str | Path,
    prefix: str,
) -> None:
    """Verify that a borrowed MMN and the checkpoint share one NSCF run."""

    checkpoint_directory = Path(checkpoint_directory)
    mmn_directory = Path(mmn_directory)
    if checkpoint_directory.resolve() == mmn_directory.resolve():
        return

    # An MMN kept at the NSCF root rather than inside one trial's directory is
    # its natural home: the overlaps depend on the NSCF alone, and every trial
    # underneath shares it.  Such a directory carries no .nnkp/.eig of its own
    # to compare against, but the containment *is* the provenance -- the
    # borrower is a descendant of the run that produced the file.  num_bands,
    # num_kpts and nntot are still checked against the checkpoint downstream in
    # connection_from_mmn, so a mismatched file cannot pass silently.
    mmn_nnkp = mmn_directory / f"{prefix}.nnkp"
    if not mmn_nnkp.is_file():
        if mmn_directory.resolve() in checkpoint_directory.resolve().parents:
            return
        raise ValueError(
            f"Cannot verify borrowing {prefix}.mmn from {mmn_directory}: it "
            f"has no {prefix}.nnkp to compare against, and it is not a parent "
            f"of {checkpoint_directory} (which would establish the shared "
            "NSCF by containment)."
        )

    for block in ("kpoints", "nnkpts"):
        checkpoint_text = _nnkp_block(
            checkpoint_directory / f"{prefix}.nnkp", block
        )
        mmn_text = _nnkp_block(mmn_nnkp, block)
        if checkpoint_text != mmn_text:
            raise ValueError(
                f"Cannot borrow {prefix}.mmn from {mmn_directory}: its "
                f"'{block}' block differs from {checkpoint_directory}."
            )
    checkpoint_eig = checkpoint_directory / f"{prefix}.eig"
    mmn_eig = mmn_directory / f"{prefix}.eig"
    if checkpoint_eig.exists() and mmn_eig.exists():
        hashes = [
            hashlib.md5(path.read_bytes()).hexdigest()
            for path in (checkpoint_eig, mmn_eig)
        ]
        if hashes[0] != hashes[1]:
            raise ValueError(
                f"Cannot borrow {prefix}.mmn from {mmn_directory}: the two "
                f"{prefix}.eig files differ."
            )


@dataclass(frozen=True)
class WannierPair:
    """Base/mode models expressed in one convention-I embedding."""

    base_model: object
    mode_model: object
    embedding: np.ndarray
    Lambda_R: Mapping[tuple[int, int, int], np.ndarray] | None
    diagnostics: dict[str, dict[str, float | int]]


def load_wannier_pair(
    base_directory: str | Path,
    mode_directory: str | Path,
    *,
    prefix: str,
    cache_directory: str | Path,
    A_R_cache_files: Mapping[str, str | Path] | None = None,
    position_source: str = "mmn",
    orbital_dependent_R: bool = True,
    include_position_terms: bool = True,
    Lambda_file: str | Path | None = None,
    max_ws_pair_contamination: float | None = 0.10,
    shared_ws_centres: bool = False,
    model_options: Mapping[str, object] | None = None,
    verbose: bool = True,
) -> WannierPair:
    """Load base and mode models with consistent ``H(R)``, ``A(R)``, and tau.

    ``position_source='mmn'`` rebuilds the connection from the large overlap
    file. ``position_source='rfull'`` loads the compact postw90 export made
    with ``transl_inv_full=T`` and ``use_ws_distance=T``. Both require the
    orbital-dependent real-space mapping. ``position_source='rdat'`` is kept
    only as an explicit legacy comparison.

    ``position_source='rdat_ndegen_applied'`` reads ``_r.dat``/``_hr.dat`` from
    a Jae-Mo ``write_ndegen_applied`` run (all three of ``use_ws_distance``,
    ``transl_inv_full`` and ``write_ndegen_applied`` set; no ``.mmn`` or
    ``.chk`` is read).  Wannier90 has already chosen the Wigner-Seitz images,
    so ``orbital_dependent_R`` must be ``False``.  Because it chooses them per
    run from that run's own centres, the base and mode can spread a
    tied hopping differently; :func:`wannier.position_matrix.ws_pair_consistency`
    measures how much of ``H_mode - H_base`` that accounts for, reports it in
    ``diagnostics["ws_pair"]``, and raises above ``max_ws_pair_contamination``
    (``None`` reports without raising).

    The 0.10 default is a tripwire for a gross mismatch of rules, **not** an
    accuracy bound.  On Y2Ir2O7 mode 1 (Q = 0.01 A) the measured values were 0
    for a shared rule, 5.5% for two structures each folded from their own
    centres (what this route and the ``.mmn`` route both do), and 15.7% for
    mixed rules.  How much a given contamination moves ``dtheta`` was not
    calibrated -- see notes/progress/2026-09-20_jaemo_ndegen_yio_mode1.md.

    ``shared_ws_centres=True`` (``.mmn`` only) folds H(R) and A(R) of **both**
    structures onto the images chosen by the **base** Wannier centres -- the
    same array as the embedding tau -- so the two share one interpolant and the
    contamination above is exactly zero (assumption A3': tau *and* the image
    weights fixed across the finite difference).  The default, ``False``, folds
    each structure onto its own centres, as every run before 2026-09-23 did.
    Valid when the two structures share the lattice, the Wannier functions
    correspond one to one, and the centres move little compared with bond
    lengths; see notes/progress/2026-09-23_ws_ties_in_beta_derivative.md.
    """

    # Imported here, not at module scope: position_matrix imports this module
    # for its file readers, so a top-level import would be circular.
    from wannier.position_matrix import (
        apply_orbital_dependent_R_mapping,
        build_model_with_position_matrix,
        build_position_matrix_from_mmn,
        install_postw90_position_matrix,
        ws_pair_consistency,
    )

    base_directory = Path(base_directory)
    mode_directory = Path(mode_directory)
    cache_directory = Path(cache_directory)
    Lambda_file = Path(Lambda_file) if Lambda_file is not None else None

    if position_source not in (
        "mmn", "rfull", "rdat", "rdat_ndegen_applied", "transl_inv_full_from_chk", "right_centre_from_chk"
    ):
        raise ValueError(
            "position_source must be 'mmn', 'rfull', 'rdat', "
            "'rdat_ndegen_applied', 'transl_inv_full_from_chk' or 'right_centre_from_chk'."
        )
    if shared_ws_centres:
        if include_position_terms and position_source not in ("mmn", "transl_inv_full_from_chk", "right_centre_from_chk"):
            raise ValueError(
                "shared_ws_centres re-folds A(R) from class sums, which is exact "
                "only for the .mmn A(R) (images of a class are copies). A "
                "transl_inv_full A(R) read from a file ('rfull', "
                "'rdat_ndegen_applied') has images that differ by (-1)**n_i; use "
                "'transl_inv_full_from_chk', which applies that phase after the fold."
            )
        if not orbital_dependent_R:
            raise ValueError(
                "shared_ws_centres replaces the centres of the in-house "
                "orbital-dependent mapping; it needs orbital_dependent_R=True."
            )
    if position_source == "rdat_ndegen_applied":
        if not include_position_terms:
            raise ValueError(
                "position_source='rdat_ndegen_applied' exists to supply A(R); "
                "it needs include_position_terms=True."
            )
        if orbital_dependent_R:
            raise ValueError(
                "rdat_ndegen_applied H(R)/A(R) already carry the Wigner-Seitz "
                "weights Wannier90 applied at write time; a second "
                "orbital-dependent mapping would apply them twice. Set "
                "orbital_dependent_R=False."
            )
    if (
        include_position_terms
        and position_source in ("mmn", "rfull", "transl_inv_full_from_chk", "right_centre_from_chk")
        and not orbital_dependent_R
    ):
        raise ValueError(
            f"The {position_source}-derived A(R) uses orbital-dependent "
            "representatives; H(R) must use the same mapping. Set "
            "orbital_dependent_R=True."
        )
    if include_position_terms and position_source == "rdat" and orbital_dependent_R:
        raise ValueError(
            "The legacy _r.dat has the old R-set and cannot be paired with an "
            "orbital-dependent H(R). Use position_source='mmn'/'rfull' or "
            "disable the mapping."
        )

    def A_R_cache(label: str, directory: Path) -> Path:
        if A_R_cache_files is not None and label in A_R_cache_files:
            return Path(A_R_cache_files[label])
        return cache_directory / f"A_R_{label}_{directory.name}.npz"

    required = []
    for label, directory in (("base", base_directory), ("mode", mode_directory)):
        required.extend(
            [directory / f"{prefix}.win", directory / f"{prefix}_hr.dat"]
        )
        if include_position_terms and position_source in ("mmn", "rdat", "transl_inv_full_from_chk", "right_centre_from_chk"):
            required.append(directory / f"{prefix}_r.dat")
        if include_position_terms and position_source in ("transl_inv_full_from_chk", "right_centre_from_chk"):
            required.extend(directory / f"{prefix}.{ext}" for ext in ("chk", "nnkp"))
        if include_position_terms and position_source == "rfull":
            required.append(directory / f"{prefix}_r_full.dat")
        if position_source == "rdat_ndegen_applied":
            required.extend(
                [
                    directory / f"{prefix}_r.dat",
                    _plain_or_gzip_path(directory / f"{prefix}_wsvec.dat"),
                ]
            )
        # A warm A(R) cache stands in for the deleted .chk/.mmn; the builder
        # checks its provenance.
        if (
            include_position_terms
            and position_source == "mmn"
            and not A_R_cache(label, directory).is_file()
        ):
            required.extend(
                directory / f"{prefix}.{extension}"
                for extension in ("chk", "mmn", "nnkp")
            )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing Wannier inputs:\n" + "\n".join(missing))

    def log(message: str) -> None:
        if verbose:
            print(message)

    log(f"Loading Wannier90 data:\n  base: {base_directory}\n  mode: {mode_directory}")
    w90_base = W90(str(base_directory), prefix)
    w90_mode = W90(str(mode_directory), prefix)

    if orbital_dependent_R or position_source == "rdat_ndegen_applied":
        grids = []
        for directory in (base_directory, mode_directory):
            text = (directory / f"{prefix}.win").read_text(errors="ignore")
            match = re.search(
                r"^\s*mp_grid\s*[=:]?\s*(\d+)\s+(\d+)\s+(\d+)",
                text,
                re.I | re.M,
            )
            if match is None:
                raise ValueError(f"No mp_grid in {directory / f'{prefix}.win'}")
            grids.append(tuple(int(value) for value in match.groups()))
        if grids[0] != grids[1]:
            raise ValueError(
                f"Base and mode mp_grid differ ({grids[0]} vs {grids[1]})."
            )
        mp_grid = grids[0]

    # One array of centres for all four foldings (H and A, base and mode) when
    # shared; None lets each fold use its own run's centres.
    # The centres that choose Wigner-Seitz images.  H(R) and A(R) of a structure
    # must be folded on the SAME array (CLAUDE.md pitfall 1): the .mmn and .chk
    # A(R) builders fold on the checkpoint's centres, so H(R) uses those too.
    # PythTB's centres (read from text files) differ by ~5e-9 A, which is enough
    # to flip near-ties at the 1e-5 A tolerance edge (20 H(R) entries on YIO mode 1).
    def fold_centres(label: str, directory: Path, w90) -> np.ndarray:
        if include_position_terms and position_source == "mmn":
            cache = A_R_cache(label, directory)
            if cache.is_file():
                with np.load(cache, allow_pickle=False) as archive:
                    return np.asarray(archive["wcc_cart"], float)
        if include_position_terms and position_source in ("mmn", "transl_inv_full_from_chk", "right_centre_from_chk"):
            return np.asarray(
                read_wannier_checkpoint(directory / f"{prefix}.chk")["wannier_centres"], float
            )
        return np.asarray(w90.lattice.orb_vecs, float) @ np.asarray(w90.lat, float)

    own_centres = {
        "base": fold_centres("base", base_directory, w90_base),
        "mode": fold_centres("mode", mode_directory, w90_mode),
    } if orbital_dependent_R else {"base": None, "mode": None}

    ws_centres = None
    if shared_ws_centres:
        if not np.allclose(w90_base.lat, w90_mode.lat, rtol=0.0, atol=1e-8):
            raise ValueError(
                "shared_ws_centres needs one lattice: base and mode cells differ, "
                "so their aliasing classes are not the same sets."
            )
        if len(w90_base.lattice.orb_vecs) != len(w90_mode.lattice.orb_vecs):
            raise ValueError("Base and mode have different numbers of Wannier functions.")
        ws_centres = own_centres["base"]
        shift = np.linalg.norm(own_centres["mode"] - ws_centres, axis=1)
        diagnostics_ws = {
            "shared": 1,
            "max_centre_shift": float(shift.max()),
            "n_centres_moved_beyond_tolerance": int((shift > 1e-5).sum()),
        }
        log(
            "Shared Wigner-Seitz rule: both structures folded on the base centres "
            f"(mode centres differ by up to {shift.max():.2e} A; "
            f"{diagnostics_ws['n_centres_moved_beyond_tolerance']} beyond the 1e-5 A tie tolerance)"
        )
    if orbital_dependent_R:
        log(f"Applying orbital-dependent R mapping on mp_grid={mp_grid}")
        for label, w90 in (("base", w90_base), ("mode", w90_mode)):
            apply_orbital_dependent_R_mapping(
                w90,
                mp_grid,
                centers_cartesian=ws_centres if shared_ws_centres else own_centres[label],
                verbose=verbose,
            )

    options = dict(
        zero_energy=0,
        min_hopping_norm=1e-5,
        max_distance=125,
        ignorable_imaginary_part=None,
    )
    if model_options is not None:
        options.update(model_options)

    diagnostics: dict[str, dict[str, float | int]] = {}
    if shared_ws_centres:
        diagnostics["ws_centres"] = diagnostics_ws

    # Take tau from the base endpoint once, then build both models in that one
    # convention-I frame so (H_mode-H_base)/dbeta contains no gauge rotation.
    if position_source == "rdat_ndegen_applied":
        # Wannier90 applied the Wigner-Seitz weights when it wrote the files, so
        # H(R) and A(R) are taken as parsed; the entry point checks the R=0
        # diagonal of A(R) against the centres and H(R)/A(R) on one R-set.
        base_model, base_result = build_model_with_position_matrix(
            w90_base,
            directory=base_directory,
            prefix=prefix,
            source=position_source,
            model_options=options,
        )
        embedding = base_model.orb_vecs.copy()
        mode_model, mode_result = build_model_with_position_matrix(
            w90_mode,
            directory=mode_directory,
            prefix=prefix,
            source=position_source,
            model_options={**options, "orb_vecs": embedding},
        )
        for label, result in (("base", base_result), ("mode", mode_result)):
            diagnostics[label] = {
                key: value
                for key, value in result.diagnostics.items()
                if isinstance(value, (int, float, np.integer, np.floating))
            }
            log(
                f"  {label} A(R): {result.n_R} R vectors; Hermiticity "
                f"rel={result.diagnostics['ndegen_applied_hermiticity_relative']:.2e}"
                f"; R=0 diagonal vs centres "
                f"{result.diagnostics['centres_max_error']:.2e} A"
            )
        # ``deg`` is 1 for every R (the builder refuses anything else), so the
        # parsed blocks are the effective H(R) and need no copy.
        pair = ws_pair_consistency(
            {R: block["h"] for R, block in w90_base.ham_r.items()},
            {R: block["h"] for R, block in w90_mode.ham_r.items()},
            mp_grid,
        )
        diagnostics["ws_pair"] = pair
        log(
            f"  base/mode Wigner-Seitz images: "
            f"{pair['entries_mismatched']}/{pair['entries_compared']} entries "
            f"changed support ({pair['mismatch_fraction']:.2%}); interpolant "
            f"change is {pair['delta_H_contamination']:.2%} of |H_mode-H_base|; "
            f"R lists {pair['n_R_base']}/{pair['n_R_mode']} "
            f"({pair['n_R_common']} common)"
        )
        if (
            max_ws_pair_contamination is not None
            and pair["delta_H_contamination"] > max_ws_pair_contamination
        ):
            raise ValueError(
                "Base and mode spread Wigner-Seitz-tied hoppings over different "
                f"images: {pair['delta_H_contamination']:.2%} of "
                "|H_mode - H_base| is a change of interpolant, not of the "
                f"physics (limit {max_ws_pair_contamination:.2%}). The finite "
                "difference would mix the two. Pass "
                "max_ws_pair_contamination=None to inspect "
                "diagnostics['ws_pair'] anyway."
            )
    else:
        embedding = w90_base.model(**options).orb_vecs.copy()
        base_model = w90_base.model(orb_vecs=embedding, **options)
        mode_model = w90_mode.model(orb_vecs=embedding, **options)
        if orbital_dependent_R:
            # The same image check as the ndegen route, on the in-house fold
            # (deg is 1 after the mapping).  Exactly 0 with shared centres.
            pair = ws_pair_consistency(
                {R: block["h"] for R, block in w90_base.ham_r.items()},
                {R: block["h"] for R, block in w90_mode.ham_r.items()},
                mp_grid,
            )
            diagnostics["ws_pair"] = pair
            log(
                f"  base/mode Wigner-Seitz images: {pair['entries_mismatched']}/"
                f"{pair['entries_compared']} entries changed support; interpolant "
                f"change is {pair['delta_H_contamination']:.2%} of |H_mode-H_base|"
            )
            if shared_ws_centres and pair["entries_mismatched"] != 0:
                raise RuntimeError(
                    "shared_ws_centres was set but base and mode H(R) still use "
                    f"different images ({pair['entries_mismatched']} entries)."
                )
            if (
                max_ws_pair_contamination is not None
                and pair["delta_H_contamination"] > max_ws_pair_contamination
            ):
                raise ValueError(
                    f"{pair['delta_H_contamination']:.2%} of |H_mode - H_base| is a "
                    "change of Wigner-Seitz images (limit "
                    f"{max_ws_pair_contamination:.2%}); consider shared_ws_centres=True."
                )
    for label, model in (("base", base_model), ("mode", mode_model)):
        if not np.allclose(model.orb_vecs, embedding, atol=1e-12):
            raise RuntimeError(f"{label} model did not accept the shared embedding.")
    log(
        "Built both models in one embedding "
        f"(max |tau|={np.abs(embedding).max():.6f} reduced)"
    )

    if include_position_terms and position_source != "rfull":
        for label, model in (("base", base_model), ("mode", mode_model)):
            if model._pos_r is None:
                raise ValueError(f"{label} model has no position matrix.")

    if include_position_terms and position_source == "rfull":
        for label, model, w90, directory in (
            ("base", base_model, w90_base, base_directory),
            ("mode", mode_model, w90_mode, mode_directory),
        ):
            centers_cartesian = (
                np.asarray(w90.lattice.orb_vecs, dtype=float)
                @ np.asarray(w90.lat, dtype=float)
            )
            info = install_postw90_position_matrix(
                model,
                directory / f"{prefix}_r_full.dat",
                centers_cartesian=centers_cartesian,
            )
            diagnostics[label] = {
                "n_R_pos": int(info["n_R_pos"]),
                "zero_H_blocks_added": int(info["zero_H_blocks_added"]),
                "pair_hermiticity_max": float(info["pair_hermiticity_max"]),
                "pair_hermiticity_relative": float(
                    info["pair_hermiticity_relative"]
                ),
                "centres_max_error": float(info["centres_max_error"]),
            }
            log(
                f"  {label} A(R): loaded {info['n_R_pos']} pair-weighted "
                f"blocks from {prefix}_r_full.dat; Hermiticity "
                f"rel={info['pair_hermiticity_relative']:.2e}"
            )

    if include_position_terms and position_source == "mmn":
        cache_directory.mkdir(parents=True, exist_ok=True)
        for label, model, directory in (
            ("base", base_model, base_directory),
            ("mode", mode_model, mode_directory),
        ):
            info = build_position_matrix_from_mmn(
                model,
                directory,
                prefix,
                cache=A_R_cache(label, directory),
                ws_centres_cartesian=ws_centres,
                mp_grid=mp_grid if ws_centres is not None else None,
            )
            diagnostics[label] = {
                key: int(info[key]) if key.startswith("n_R") else float(info[key])
                for key in (
                    "n_R_pos",
                    "n_R_mmn",
                    "bshell_completeness_error",
                    "weight_outside_ham_R",
                    "rdat_pair_hermiticity_relative",
                    "mmn_pair_hermiticity_relative",
                    "centres_vs_rdat",
                    "residual_R0_max",
                )
            }
            diag = diagnostics[label]
            log(
                f"  {label} A(R): {diag['n_R_pos']} shared R vectors; "
                f"outside-H weight={diag['weight_outside_ham_R']:.2e}; "
                f"curl shell error={diag['bshell_completeness_error']:.2e}"
            )
            # ``info`` also carries the multi-gigabyte uncompressed AA array.
            # The model now owns the converted position matrix, so release it.
            del info

    if include_position_terms and position_source in ("transl_inv_full_from_chk", "right_centre_from_chk"):
        # Both formulas rebuilt from the checkpoint's m_matrix, folded on
        # ws_centres (shared) or each run's own .chk centres.  transl_inv_full_from_chk:
        # Wannier90's transl_inv_full + write_ndegen_applied formula, per-b phase
        # exp(-i b.R/2) at the final image.  right_centre_from_chk: the production .mmn formula.
        from wannier.position_matrix import (
            _position_matrix_for_model,
            position_matrix_mmn_formula,
            position_matrix_transl_inv_full,
        )

        for label, model, directory in (
            ("base", base_model, base_directory),
            ("mode", mode_model, mode_directory),
        ):
            checkpoint = read_wannier_checkpoint(
                directory / f"{prefix}.chk", include_m_matrix=True
            )
            checkpoint.pop("v_matrix", None)
            if position_source == "right_centre_from_chk":
                data = position_matrix_mmn_formula(
                    checkpoint, directory / f"{prefix}.nnkp", ws_centres_cartesian=ws_centres
                )
                keys = data["iRvec"]
            else:
                position = position_matrix_transl_inv_full(
                    checkpoint,
                    directory / f"{prefix}.nnkp",
                    ws_centres_cartesian=ws_centres,
                    verbose=False,
                )
                keys = sorted(position)
                data = {
                    "AA": np.stack([position.pop(R).transpose(1, 2, 0) for R in keys]),
                    "iRvec": np.asarray(keys, dtype=int),
                    "wcc_cart": np.asarray(checkpoint["wannier_centres"]),
                }
                del position
            del checkpoint
            model._pos_r, info = _position_matrix_for_model(model, data)
            del data
            diagnostics[label] = {
                key: float(value) for key, value in info.items()
                if isinstance(value, (int, float, np.integer, np.floating))
            }
            log(f"  {label} A(R): {position_source} from .chk m_matrix, {len(keys)} R vectors")

    if include_position_terms:
        for label, model in (("base", base_model), ("mode", mode_model)):
            if model._pos_r is None:
                raise ValueError(f"{label} model has no position matrix.")
            if set(model._ham_r) != set(model._pos_r):
                raise ValueError(f"{label} model has inconsistent H(R)/A(R) sets.")

    Lambda_R = None
    if include_position_terms and Lambda_file is not None and Lambda_file.is_file():
        Lambda_R = load_Lambda_R(Lambda_file)
        log(f"Loaded Lambda(R) from {Lambda_file}")
    elif include_position_terms:
        log("No Lambda(R): using the frozen Wannier gauge A_beta=0")

    return WannierPair(
        base_model=base_model,
        mode_model=mode_model,
        embedding=embedding,
        Lambda_R=Lambda_R,
        diagnostics=diagnostics,
    )
