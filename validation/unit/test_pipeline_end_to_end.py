"""File boundaries: what the pair loader refuses, and what the result writer writes.

The curvature and dtheta algebra now live in ``test_curvature.py``.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from modules.axion import save_dtheta_results
from wannier.wannier_io import load_wannier_pair


class TestWannierIOBoundary(unittest.TestCase):
    def test_rejects_mixed_R_conventions_before_reading_files(self):
        with self.assertRaisesRegex(ValueError, "legacy _r.dat"):
            load_wannier_pair(
                "missing-base",
                "missing-mode",
                prefix="test",
                cache_directory="missing-cache",
                position_source="rdat",
                orbital_dependent_R=True,
            )

    def test_rejects_rfull_without_pair_dependent_H_mapping(self):
        with self.assertRaisesRegex(ValueError, r"rfull-derived A\(R\)"):
            load_wannier_pair(
                "missing-base",
                "missing-mode",
                prefix="test",
                cache_directory="missing-cache",
                position_source="rfull",
                orbital_dependent_R=False,
            )

    def test_result_writer_uses_readable_names_and_legacy_aliases(self):
        result = {
            "nk": 2,
            "d3k": 0.125,
            "dtheta": np.array([1.0 + 2.0j, 3.0 + 4.0j]),
            "dtheta_simpson": np.array([1.1 + 0j, 3.1 + 0j]),
            "minimum_gap": np.array([0.5, 0.6]),
            "minimum_gap_kpoint": np.zeros((2, 3)),
            "c2_density": np.zeros((2, 2, 2, 2)),
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            summary_npz, summary_csv = save_dtheta_results(
                output, [result], {"position_source": np.array(["mmn"])}
            )
            self.assertTrue(summary_csv.is_file())
            with np.load(summary_npz, allow_pickle=False) as data:
                self.assertIn("minimum_gap", data.files)
                self.assertIn("min_gap", data.files)
                np.testing.assert_array_equal(data["minimum_gap"], data["min_gap"])


if __name__ == "__main__":
    unittest.main()
