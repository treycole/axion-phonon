"""``shared_ws_centres``: fold base and mode onto one Wigner-Seitz rule.

* ``remap_A_R_to_centres`` re-folds an ``.mmn`` ``A(R)`` (images of a class are
  copies) exactly: folding on centres P and re-folding on Q equals folding on Q.
* After a shared fold, a displaced structure that broke a base tie no longer
  changes images (``ws_pair_consistency`` finds nothing).
* ``load_wannier_pair`` refuses the option where it is not valid.
"""

from __future__ import annotations

import tempfile
import unittest

import numpy as np

from wannier.position_matrix import (
    _hermitian_A_R,
    _map_R_by_orbital_positions,
    remap_A_R_to_centres,
    ws_pair_consistency,
)
from wannier.wannier_io import load_wannier_pair

GRID = (4, 4, 4)
LATTICE = np.eye(3) * 3.0
# The class period along x is 4 cells = 12 A, so R = (2,0,0) and (-2,0,0) tie
# exactly for every orbital pair whose centres share x (both here; the y offset
# keeps the pair distinct).  Moving orbital 1 along x in the mode breaks it.
BASE_CENTRES = np.array([[0.0, 0.0, 0.0], [0.0, 1.5, 0.0]])
MODE_CENTRES = BASE_CENTRES + np.array([[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]])


def _class_sums(seed: int, leading=(3,)):
    """Random ``{R: (..., 2, 2)}`` on the raw 4x4x4 grid (one entry per class)."""
    rng = np.random.default_rng(seed)
    return {
        R: rng.normal(size=leading + (2, 2)) + 1j * rng.normal(size=leading + (2, 2))
        for R in np.ndindex(GRID)
    }


def _fold(raw, centres):
    return _map_R_by_orbital_positions(raw, LATTICE, centres, GRID)


class RemapTests(unittest.TestCase):
    def test_refold_equals_direct_fold(self):
        raw = _class_sums(0)
        via_mode = _hermitian_A_R(_fold(raw, MODE_CENTRES), 2)
        direct = _hermitian_A_R(_fold(raw, BASE_CENTRES), 2)
        data = {"AA": via_mode[0], "iRvec": via_mode[1], "real_lattice": LATTICE}
        refolded = remap_A_R_to_centres(data, BASE_CENTRES, GRID)
        got = {tuple(R): A for R, A in zip(refolded["iRvec"], refolded["AA"])}
        want = {tuple(R): A for R, A in zip(direct[1], direct[0])}
        for R in set(got) | set(want):
            np.testing.assert_allclose(
                got.get(R, 0.0), want.get(R, 0.0), atol=1e-12, err_msg=str(R)
            )

    def test_the_two_folds_really_differ(self):
        # Guard against a vacuous test: own-centre folds must break a tie.
        raw = _class_sums(1, leading=())
        base, mode = _fold(raw, BASE_CENTRES), _fold(raw, MODE_CENTRES)
        own = ws_pair_consistency(base, mode, GRID, significance=0.0)
        self.assertGreater(own["entries_mismatched"], 0)

    def test_shared_fold_removes_the_image_change(self):
        raw_base = _class_sums(2, leading=())
        # A physical change of every class sum, small compared with the sums.
        raw_mode = {R: H + 1e-2 * np.ones_like(H) for R, H in raw_base.items()}
        shared = ws_pair_consistency(
            _fold(raw_base, BASE_CENTRES),
            _fold(raw_mode, BASE_CENTRES),
            GRID,
            significance=0.0,
        )
        self.assertEqual(shared["entries_mismatched"], 0)
        self.assertEqual(shared["delta_H_contamination"], 0.0)

    def test_wrong_shape_is_rejected(self):
        data = {"AA": np.zeros((1, 2, 2, 3)), "iRvec": np.zeros((1, 3), int),
                "real_lattice": LATTICE}
        with self.assertRaisesRegex(ValueError, "shape"):
            remap_A_R_to_centres(data, np.zeros((3, 3)), GRID)


class LoaderGuardTests(unittest.TestCase):
    def _call(self, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            return load_wannier_pair(
                tmp, tmp, prefix="seed", cache_directory=tmp,
                shared_ws_centres=True, **kwargs,
            )

    def test_transl_inv_full_sources_are_refused(self):
        for source in ("rfull", "rdat_ndegen_applied"):
            with self.subTest(source=source), self.assertRaisesRegex(
                ValueError, "copies"
            ):
                self._call(position_source=source, orbital_dependent_R=True)

    def test_needs_the_in_house_mapping(self):
        with self.assertRaisesRegex(ValueError, "orbital_dependent_R=True"):
            self._call(position_source="mmn", orbital_dependent_R=False)


if __name__ == "__main__":
    unittest.main()
