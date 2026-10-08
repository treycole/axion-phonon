"""Validate a Wannier90 position matrix against the MMN construction.

This is a run-level regression test for the large SrTiO3 artifacts, which are
not suitable for the ordinary unit-test suite.  It compares every complex
component of ``<0i|r|Rj>`` after putting the two inputs in the same convention:

* Wigner--Seitz degeneracies are removed from ``seedname_r.dat``;
* orbital-dependent R representatives are selected with
  ``min |R + tau_j - tau_i|``; and
* the cached MMN residual has its R=0 diagonal replaced by the Wannier centres.

The optional old run is a control: it detects whether adding
``transl_inv_full`` changed the position matrix written by ``wannier90.x``.
Wannier90 documents that keyword as a ``postw90.x`` option, so this control is
important when interpreting a failed equivalence check.

Pure library: no argparse, no ``__main__``. Used directly by
``validate_transl_inv_full.ipynb`` and imported by
``validation/unit/test_position_matrix_compare.py`` (``compare_matrices``).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO, Mapping

import numpy as np
from pythtb import W90

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import (
    STO_AA_CACHE,
    STO_BASE_RUN,
)
from validation._paths import DATA
from wannier.wannier_io import read_wannier_checkpoint
from wannier.position_matrix import _map_R_by_orbital_positions

PREFIX = "SrTiO3"
TRIAL = STO_BASE_RUN / "trial_03_Sr_sp_Ti_pd_O_sp"
DEFAULT_OUTPUT = TRIAL / "188725bc889"  # a transl_inv_full rerun of set c; no longer in data/, only its A(R) cache
DEFAULT_REFERENCE = DATA / "SrTiO3/base/soc/output/181642bc889/trial_03_Sr_sp_Ti_pd_O_sp/185116bc889"
DEFAULT_CACHE = (
    DATA / "SrTiO3/_AA_cache/AA_trial_03_Sr_sp_Ti_pd_O_sp_base_185116bc889_wbnone_ws1em05.npz"
)
ZERO_R = (0, 0, 0)


@dataclass(frozen=True)
class MatrixComparison:
    """Numerical comparison of two R-indexed matrix dictionaries."""

    equivalent: bool
    max_abs_error_angstrom: float
    rms_error_angstrom: float
    relative_frobenius_error: float
    candidate_frobenius_norm_angstrom: float
    reference_frobenius_norm_angstrom: float
    n_R_candidate: int
    n_R_reference: int
    n_R_union: int
    worst_R: tuple[int, int, int]
    worst_component: int
    worst_left_orbital: int
    worst_right_orbital: int
    worst_candidate_real: float
    worst_candidate_imag: float
    worst_reference_real: float
    worst_reference_imag: float


def effective_rdat(w90: W90) -> dict[tuple[int, int, int], np.ndarray]:
    """Return position blocks with Wannier90 degeneracy weights removed."""

    if w90.pos_r is None:
        raise ValueError("The Wannier90 run has no position matrix (_r.dat/_tb.dat).")
    if set(w90.pos_r) != set(w90.ham_r):
        raise ValueError("The position and Hamiltonian R-vector sets differ.")
    return {
        R: np.asarray(block, dtype=complex) / float(w90.ham_r[R]["deg"])
        for R, block in w90.pos_r.items()
    }


def _relocated(recorded: str) -> Path:
    """A cache's recorded input path, followed to ``data/`` if the runs were moved there."""
    path = Path(recorded)
    moved = Path(recorded.replace("/calculations/", "/data/", 1))
    return moved if not path.exists() and moved.exists() else path


def position_from_cache(
    path: Path,
) -> tuple[dict[tuple[int, int, int], np.ndarray], np.ndarray, dict]:
    """Load the MMN-derived residual and restore its physical R=0 diagonal."""

    with np.load(path, allow_pickle=False) as archive:
        required = {"AA", "iRvec", "wcc_cart", "cache_meta"}
        missing = required - set(archive.files)
        if missing:
            raise ValueError(f"{path} is missing cache fields: {sorted(missing)}")
        A_R = np.asarray(archive["AA"])
        R_vectors = np.asarray(archive["iRvec"])
        centres = np.asarray(archive["wcc_cart"], dtype=float)
        metadata = json.loads(str(archive["cache_meta"].item()))

    if A_R.shape != (len(R_vectors), len(centres), len(centres), 3):
        raise ValueError(
            f"Unexpected AA shape {A_R.shape}; expected "
            f"({len(R_vectors)}, {len(centres)}, {len(centres)}, 3)."
        )
    matrices = {
        tuple(int(value) for value in R): block.transpose(2, 0, 1).copy()
        for R, block in zip(R_vectors, A_R)
    }
    if ZERO_R not in matrices:
        raise ValueError("The MMN-derived cache has no R=(0,0,0) block.")
    diagonal = np.arange(len(centres))
    matrices[ZERO_R][:, diagonal, diagonal] = centres.T
    return matrices, centres, metadata


def hermitianize(
    matrices: Mapping[tuple[int, int, int], np.ndarray],
) -> dict[tuple[int, int, int], np.ndarray]:
    """Apply the same R/-R Hermitian projection as the MMN pipeline."""

    sample = next(iter(matrices.values()))
    zero = np.zeros_like(sample)
    all_R = set(matrices)
    all_R.update(tuple(-component for component in R) for R in matrices)
    return {
        R: 0.5
        * (
            matrices.get(R, zero)
            + np.conj(
                np.transpose(
                    matrices.get(tuple(-component for component in R), zero),
                    (0, 2, 1),
                )
            )
        )
        for R in sorted(all_R)
    }


def align_rdat(
    matrices: Mapping[tuple[int, int, int], np.ndarray],
    lattice: np.ndarray,
    centres: np.ndarray,
    mp_grid: tuple[int, int, int],
    distance_tolerance: float,
) -> dict[tuple[int, int, int], np.ndarray]:
    """Use the MMN pipeline's orbital-dependent R representation."""

    mapped = _map_R_by_orbital_positions(
        matrices,
        lattice,
        centres,
        mp_grid,
        tolerance=distance_tolerance,
    )
    return hermitianize(mapped)


def compare_matrices(
    candidate: Mapping[tuple[int, int, int], np.ndarray],
    reference: Mapping[tuple[int, int, int], np.ndarray],
    *,
    max_abs_tolerance: float,
    relative_tolerance: float,
) -> MatrixComparison:
    """Compare R-indexed matrices without stacking the large production data."""

    if not candidate or not reference:
        raise ValueError("Both position-matrix dictionaries must be non-empty.")
    sample = np.asarray(next(iter(reference.values())))
    if sample.ndim != 3 or sample.shape[0] != 3:
        raise ValueError("Position blocks must have shape (3, num_wann, num_wann).")
    zero = np.zeros_like(sample)
    candidate_norm_squared = 0.0
    reference_norm_squared = 0.0
    error_norm_squared = 0.0
    count = 0
    maximum = -1.0
    worst = (ZERO_R, 0, 0, 0, 0j, 0j)

    all_R = sorted(set(candidate) | set(reference))
    for R in all_R:
        left = np.asarray(candidate.get(R, zero))
        right = np.asarray(reference.get(R, zero))
        if left.shape != sample.shape or right.shape != sample.shape:
            raise ValueError(f"Inconsistent position-block shape at R={R}.")
        difference = left - right
        absolute = np.abs(difference)
        local_index = np.unravel_index(int(np.argmax(absolute)), absolute.shape)
        local_maximum = float(absolute[local_index])
        if local_maximum > maximum:
            maximum = local_maximum
            component, orbital_i, orbital_j = local_index
            worst = (
                R,
                int(component),
                int(orbital_i),
                int(orbital_j),
                complex(left[local_index]),
                complex(right[local_index]),
            )
        candidate_norm_squared += float(np.vdot(left, left).real)
        reference_norm_squared += float(np.vdot(right, right).real)
        error_norm_squared += float(np.vdot(difference, difference).real)
        count += difference.size

    candidate_norm = float(np.sqrt(candidate_norm_squared))
    reference_norm = float(np.sqrt(reference_norm_squared))
    error_norm = float(np.sqrt(error_norm_squared))
    relative = error_norm / reference_norm if reference_norm else error_norm
    rms = float(np.sqrt(error_norm_squared / count))
    equivalent = maximum <= max_abs_tolerance and relative <= relative_tolerance
    R, component, orbital_i, orbital_j, left_value, right_value = worst
    return MatrixComparison(
        equivalent=equivalent,
        max_abs_error_angstrom=maximum,
        rms_error_angstrom=rms,
        relative_frobenius_error=relative,
        candidate_frobenius_norm_angstrom=candidate_norm,
        reference_frobenius_norm_angstrom=reference_norm,
        n_R_candidate=len(candidate),
        n_R_reference=len(reference),
        n_R_union=len(all_R),
        worst_R=R,
        worst_component=component,
        worst_left_orbital=orbital_i,
        worst_right_orbital=orbital_j,
        worst_candidate_real=float(left_value.real),
        worst_candidate_imag=float(left_value.imag),
        worst_reference_real=float(right_value.real),
        worst_reference_imag=float(right_value.imag),
    )


def binary_reader(path: Path) -> BinaryIO:
    return gzip.open(path, "rb") if path.suffix == ".gz" else path.open("rb")


def content_digest(path: Path) -> str:
    """Hash uncompressed content, making .mmn and .mmn.gz comparable."""

    digest = hashlib.md5()  # nosec B324: file identity, not cryptography
    with binary_reader(path) as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def active_command_mentions(path: Path, executable: str) -> bool:
    if not path.is_file():
        return False
    pattern = re.compile(rf"(^|[ /]){re.escape(executable)}([ ]|$)")
    return any(
        pattern.search(line) is not None
        for line in path.read_text(errors="ignore").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )


def print_comparison(label: str, result: MatrixComparison) -> None:
    status = "PASS" if result.equivalent else "FAIL"
    print(f"{label}: {status}")
    print(f"  relative Frobenius error : {result.relative_frobenius_error:.9e}")
    print(f"  maximum absolute error   : {result.max_abs_error_angstrom:.9e} Angstrom")
    print(f"  RMS element error        : {result.rms_error_angstrom:.9e} Angstrom")
    print(
        "  R-vector counts         : "
        f"candidate={result.n_R_candidate}, reference={result.n_R_reference}, "
        f"union={result.n_R_union}"
    )
    print(
        "  worst element            : "
        f"R={result.worst_R}, alpha={result.worst_component}, "
        f"i={result.worst_left_orbital}, j={result.worst_right_orbital} "
        "(zero-based; Wannier90 i/j="
        f"{result.worst_left_orbital + 1}/{result.worst_right_orbital + 1})"
    )
    print(
        "  candidate/reference      : "
        f"({result.worst_candidate_real:.9e}{result.worst_candidate_imag:+.9e}j) / "
        f"({result.worst_reference_real:.9e}{result.worst_reference_imag:+.9e}j)"
    )


def run_equivalence_check(
    output_dir: Path,
    cache_path: Path,
    reference_dir: Path | None,
    prefix: str = PREFIX,
    max_abs: float = 1e-6,
    max_relative: float = 1e-6,
    r_distance_tolerance: float = 1e-5,
) -> dict:
    """Run the full equivalence check and return the JSON-serializable report."""

    output_dir = output_dir.resolve()
    cache_path = cache_path.resolve()
    checkpoint_path = output_dir / f"{prefix}.chk"
    win_path = output_dir / f"{prefix}.win"
    job_path = output_dir / "send_job.sh"
    mmn_candidates = (
        output_dir / f"{prefix}.mmn",
        output_dir / f"{prefix}.mmn.gz",
    )
    output_mmn = next((path for path in mmn_candidates if path.is_file()), None)
    for path in (output_dir, cache_path, checkpoint_path, win_path):
        if not path.exists():
            raise FileNotFoundError(path)

    cached_position, cached_centres, cache_metadata = position_from_cache(cache_path)
    checkpoint = read_wannier_checkpoint(checkpoint_path)
    w90 = W90(str(output_dir), prefix)
    output_position = effective_rdat(w90)
    aligned_output = align_rdat(
        output_position,
        np.asarray(checkpoint["real_lattice"]),
        cached_centres,
        tuple(int(value) for value in checkpoint["mp_grid"]),
        r_distance_tolerance,
    )
    equivalence = compare_matrices(
        aligned_output,
        cached_position,
        max_abs_tolerance=max_abs,
        relative_tolerance=max_relative,
    )

    cache_inputs = cache_metadata.get("inputs", {})
    source_chk = _relocated(cache_inputs.get("chk", {}).get("path", ""))
    source_mmn = _relocated(cache_inputs.get("mmn", {}).get("path", ""))
    gauge_max_error = float("nan")
    centre_max_error = float("nan")
    if source_chk.is_file():
        cached_checkpoint = read_wannier_checkpoint(source_chk)
        gauge_max_error = float(
            np.max(
                np.abs(
                    np.asarray(checkpoint["v_matrix"])
                    - np.asarray(cached_checkpoint["v_matrix"])
                )
            )
        )
        centre_max_error = float(
            np.max(
                np.abs(
                    np.asarray(checkpoint["wannier_centres"])
                    - np.asarray(cached_checkpoint["wannier_centres"])
                )
            )
        )

    mmn_content_equal: bool | None = None
    if output_mmn is not None and source_mmn.is_file():
        mmn_content_equal = content_digest(output_mmn) == content_digest(source_mmn)

    reference_comparison: MatrixComparison | None = None
    if reference_dir is not None and reference_dir.is_dir():
        old_w90 = W90(str(reference_dir.resolve()), prefix)
        reference_comparison = compare_matrices(
            output_position,
            effective_rdat(old_w90),
            max_abs_tolerance=max_abs,
            relative_tolerance=max_relative,
        )

    win_text = win_path.read_text(errors="ignore")
    keyword_enabled = re.search(
        r"^\s*transl_inv_full\s*=\s*(true|t|\.true\.)\s*(?:!.*)?$",
        win_text,
        flags=re.IGNORECASE | re.MULTILINE,
    ) is not None
    ran_wannier90 = active_command_mentions(job_path, "wannier90.x")
    ran_postw90 = active_command_mentions(job_path, "postw90.x")

    print("Position-matrix equivalence test")
    print(f"  output _r.dat           : {output_dir}")
    print(f"  MMN-derived reference   : {cache_path}")
    print(
        "  tolerances              : "
        f"max={max_abs:.1e} Angstrom, relative={max_relative:.1e}"
    )
    print(f"  transl_inv_full enabled : {keyword_enabled}")
    print(f"  active wannier90.x call : {ran_wannier90}")
    print(f"  active postw90.x call   : {ran_postw90}")
    print(f"  MMN content matches cache input : {mmn_content_equal}")
    print(f"  checkpoint gauge max error      : {gauge_max_error:.9e}")
    print(f"  checkpoint centre max error     : {centre_max_error:.9e} Angstrom")
    print_comparison("\n_r.dat vs MMN-derived position matrix", equivalence)
    if reference_comparison is not None:
        print_comparison(
            "\n_r.dat vs old no-transl_inv_full control", reference_comparison
        )

    report = {
        "output_directory": str(output_dir),
        "cache": str(cache_path),
        "cache_metadata": cache_metadata,
        "keyword_enabled": keyword_enabled,
        "active_wannier90_call": ran_wannier90,
        "active_postw90_call": ran_postw90,
        "mmn_content_equal_to_cache_input": mmn_content_equal,
        "checkpoint_gauge_max_error": gauge_max_error,
        "checkpoint_centre_max_error_angstrom": centre_max_error,
        "equivalence_tolerances": {
            "max_abs_angstrom": max_abs,
            "relative_frobenius": max_relative,
        },
        "rdat_vs_mmn": asdict(equivalence),
        "rdat_vs_old_control": (
            asdict(reference_comparison) if reference_comparison is not None else None
        ),
    }
    return report


def save_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nWrote report: {path.resolve()}")


def load_report(path: Path) -> dict:
    return json.loads(path.read_text())
