"""Checks added around the ``rdat_ndegen_applied`` position source.

* ``ws_pair_consistency`` -- do base and mode spread a Wigner-Seitz-tied
  hopping over the same images?
* the ``R = 0`` diagonal check of ``A(R)`` against the Wannier centres.
* the argument guards in ``load_wannier_pair``.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from wannier.position_matrix import (
    _NDEGEN_CENTRES_TOLERANCE,
    _build_rdat_ndegen_applied,
    ws_pair_consistency,
)
from wannier.wannier_io import load_wannier_pair

GRID = (4, 4, 4)
# R_A and R_B differ by 4 a_1: one aliasing class on a 4x4x4 mesh.
R_A = (2, 0, 0)
R_B = (-2, 0, 0)
R_ZERO = (0, 0, 0)


def _hopping(value: complex) -> np.ndarray:
    """A 2-orbital block with a Hermitian pair of off-diagonal entries."""
    return np.array([[0.0, value], [np.conj(value), 0.0]], dtype=complex)


def _onsite() -> np.ndarray:
    return np.diag([1.0, 2.0]).astype(complex)


def _model(weight_a: float, weight_b: float, hopping: float = 0.3) -> dict:
    """H(R) with one bond of size ``hopping`` split (weight_a, weight_b)."""
    return {
        R_ZERO: _onsite(),
        R_A: _hopping(weight_a * hopping),
        R_B: _hopping(weight_b * hopping),
    }


class WsPairConsistencyTests(unittest.TestCase):
    def test_identical_structures_are_consistent(self):
        base = _model(0.5, 0.5)
        result = ws_pair_consistency(base, base, GRID)
        self.assertEqual(result["entries_mismatched"], 0)
        self.assertEqual(result["delta_H_contamination"], 0.0)
        self.assertEqual(result["base_weight_fraction"], 0.0)

    def test_same_tie_pattern_with_a_physical_shift_is_consistent(self):
        base = _model(0.5, 0.5, hopping=0.30)
        mode = _model(0.5, 0.5, hopping=0.306)
        result = ws_pair_consistency(base, mode, GRID)
        self.assertEqual(result["entries_mismatched"], 0)
        self.assertEqual(result["delta_H_contamination"], 0.0)

    def test_broken_tie_is_flagged_and_dominates_the_difference(self):
        # Base: (1/2, 1/2).  Mode: the tie is broken, all weight on R_A, plus
        # the same 2% physical change as above.
        base = _model(0.5, 0.5, hopping=0.30)
        mode = _model(1.0, 0.0, hopping=0.306)
        result = ws_pair_consistency(base, mode, GRID)
        # The (0, 1) and (1, 0) entries of the bond's class.
        self.assertEqual(result["entries_mismatched"], 2)
        self.assertGreater(result["delta_H_contamination"], 0.99)
        self.assertGreater(result["base_weight_fraction"], 0.0)

    def test_a_hopping_the_displacement_creates_is_not_a_mismatch(self):
        # Symmetry forbids the bond in the base; the mode makes it 0.02 eV.
        # That is the response, not a change of interpolant.
        base = {R_ZERO: _onsite()}
        mode = {R_ZERO: _onsite(), (1, 1, 0): _hopping(0.02)}
        result = ws_pair_consistency(base, mode, GRID)
        self.assertEqual(result["entries_mismatched"], 0)
        self.assertEqual(result["n_R_common"], 1)
        self.assertEqual(result["n_R_mode"], 2)

    def test_a_tie_change_below_significance_is_ignored(self):
        base = _model(0.5, 0.5, hopping=5e-5)
        mode = _model(1.0, 0.0, hopping=5e-5)
        result = ws_pair_consistency(base, mode, GRID)
        self.assertEqual(result["entries_mismatched"], 0)

    def test_empty_input_is_rejected(self):
        with self.assertRaises(ValueError):
            ws_pair_consistency({}, _model(0.5, 0.5), GRID)


class _FakeW90:
    """Just enough of a ``pythtb.W90`` for the ndegen builder."""

    def __init__(self, centres: np.ndarray, pos_diagonal: np.ndarray, *, with_R0=True):
        self.ham_r = {R_ZERO: {"h": _onsite(), "deg": 1}}
        self._pos_diagonal = pos_diagonal
        self._with_R0 = with_R0

    def model(self, **_options):
        n = self._pos_diagonal.shape[0]
        block = np.zeros((3, n, n), dtype=complex)
        for s in range(n):
            block[:, s, s] = self._pos_diagonal[s]
        pos_r = {R_ZERO: block} if self._with_R0 else {(1, 0, 0): block}
        return SimpleNamespace(_pos_r=pos_r)


class NdegenCentresCheckTests(unittest.TestCase):
    CENTRES = np.array([[0.1, 0.2, 0.3], [1.5, -0.7, 2.25]])

    def _run(self, pos_diagonal, *, with_R0=True):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / "seed_wsvec.dat").write_text(
                "## written on 20Sep2026 with use_ws_distance=.true.  "
                "write_ndegen_applied=.true.\n"
            )
            w90 = _FakeW90(self.CENTRES, pos_diagonal, with_R0=with_R0)
            return _build_rdat_ndegen_applied(
                w90,
                {},
                directory=directory,
                prefix="seed",
                centres_cartesian=lambda _w90: self.CENTRES,
            )

    def test_matching_centres_pass_and_the_error_is_reported(self):
        _model_, _reads, _rule, info = self._run(self.CENTRES + 3e-7)
        self.assertLess(info["centres_max_error"], _NDEGEN_CENTRES_TOLERANCE)
        self.assertAlmostEqual(info["centres_max_error"], 3e-7, places=12)

    def test_shifted_centre_is_rejected(self):
        wrong = self.CENTRES.copy()
        wrong[1, 2] += 0.4  # a wrong b-vector ordering moves it by ~Angstrom
        with self.assertRaisesRegex(ValueError, "R=0 diagonal differs"):
            self._run(wrong)

    def test_missing_R0_block_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "no R = 0 block"):
            self._run(self.CENTRES, with_R0=False)


class LoaderGuardTests(unittest.TestCase):
    def _call(self, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            return load_wannier_pair(
                tmp,
                tmp,
                prefix="seed",
                cache_directory=tmp,
                position_source="rdat_ndegen_applied",
                **kwargs,
            )

    def test_orbital_dependent_mapping_is_refused(self):
        with self.assertRaisesRegex(ValueError, "twice"):
            self._call(orbital_dependent_R=True)

    def test_position_terms_are_required(self):
        with self.assertRaisesRegex(ValueError, "include_position_terms"):
            self._call(orbital_dependent_R=False, include_position_terms=False)

    def test_missing_files_are_named(self):
        with self.assertRaisesRegex(FileNotFoundError, "seed_wsvec.dat"):
            self._call(orbital_dependent_R=False)


if __name__ == "__main__":
    unittest.main()
