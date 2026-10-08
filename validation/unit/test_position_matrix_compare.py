from __future__ import annotations

import unittest

import numpy as np

from validation.inputs.validate_transl_inv_full_lib import compare_matrices


class TestPositionMatrixComparison(unittest.TestCase):
    def test_identical_complex_position_matrices_pass(self):
        block = np.zeros((3, 2, 2), dtype=complex)
        block[1, 0, 1] = 0.2 - 0.3j
        matrices = {(0, 0, 0): block, (1, 0, 0): 0.5 * block}

        result = compare_matrices(
            matrices,
            matrices,
            max_abs_tolerance=1e-10,
            relative_tolerance=1e-10,
        )

        self.assertTrue(result.equivalent)
        self.assertEqual(result.max_abs_error_angstrom, 0.0)
        self.assertEqual(result.relative_frobenius_error, 0.0)

    def test_reports_the_largest_complex_matrix_element_mismatch(self):
        reference = {(0, 0, 0): np.zeros((3, 2, 2), dtype=complex)}
        candidate = {(0, 0, 0): reference[(0, 0, 0)].copy()}
        candidate[(0, 0, 0)][2, 1, 0] = 0.02 - 0.01j

        result = compare_matrices(
            candidate,
            reference,
            max_abs_tolerance=1e-6,
            relative_tolerance=1e-6,
        )

        self.assertFalse(result.equivalent)
        self.assertAlmostEqual(result.max_abs_error_angstrom, np.sqrt(5e-4))
        self.assertEqual(result.worst_R, (0, 0, 0))
        self.assertEqual(result.worst_component, 2)
        self.assertEqual(result.worst_left_orbital, 1)
        self.assertEqual(result.worst_right_orbital, 0)


if __name__ == "__main__":
    unittest.main()
