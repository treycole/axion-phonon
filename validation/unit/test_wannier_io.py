from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from wannier.wannier_io import read_mmn_overlaps
from wannier.position_matrix import (
    _A_R_cache_metadata,
    _position_matrix_for_model,
    connection_from_mmn,
    finite_difference_weights,
    install_postw90_position_matrix,
    read_postw90_position_matrix,
)


def _write_mmn(path: Path, matrices, neighbours, reciprocal_shifts) -> None:
    matrices = np.asarray(matrices)
    num_kpoints, num_neighbors, num_bands, num_bands_2 = matrices.shape
    assert num_bands == num_bands_2
    lines = [
        "synthetic test mmn",
        f"{num_bands} {num_kpoints} {num_neighbors}",
    ]
    for ik in range(num_kpoints):
        for neighbor in range(num_neighbors):
            shift = reciprocal_shifts[ik, neighbor]
            lines.append(
                f"{ik + 1} {neighbours[ik, neighbor] + 1} "
                f"{shift[0]} {shift[1]} {shift[2]}"
            )
            for value in matrices[ik, neighbor].T.flat:
                lines.append(f"{value.real:.17g} {value.imag:.17g}")
    path.write_text("\n".join(lines) + "\n")


def _write_r_full(
    path: Path,
    matrices,
    *,
    transl_inv_full: bool = True,
    use_ws_distance: bool = True,
    effective_ws_weights: bool = True,
) -> None:
    num_wann = next(iter(matrices.values())).shape[1]
    logical = lambda value: "T" if value else "F"
    lines = [
        "synthetic; "
        f"transl_inv_full={logical(transl_inv_full)}; "
        f"use_ws_distance={logical(use_ws_distance)}; "
        f"effective_ws_weights={logical(effective_ws_weights)}",
        str(num_wann),
        str(len(matrices)),
    ]
    for R, block in matrices.items():
        for right in range(num_wann):
            for left in range(num_wann):
                value = block[:, left, right]
                fields = [
                    *R,
                    left + 1,
                    right + 1,
                    value[0].real,
                    value[0].imag,
                    value[1].real,
                    value[1].imag,
                    value[2].real,
                    value[2].imag,
                ]
                lines.append(" ".join(str(item) for item in fields))
    path.write_text("\n".join(lines) + "\n")


class TestMMNInput(unittest.TestCase):
    def test_preserves_link_order_and_matrix_orientation(self):
        matrices = np.empty((2, 2, 2, 2), dtype=complex)
        for ik in range(2):
            for neighbor in range(2):
                base = 100 * ik + 10 * neighbor
                matrices[ik, neighbor] = np.array(
                    [
                        [base + 1 + 1j, base + 2 + 2j],
                        [base + 3 + 3j, base + 4 + 4j],
                    ]
                )
        neighbours = np.array([[1, 0], [1, 0]])
        reciprocal_shifts = np.array(
            [[[0, 0, 0], [1, 0, 0]], [[-1, 0, 0], [0, 0, 0]]]
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.mmn"
            _write_mmn(path, matrices, neighbours, reciprocal_shifts)
            result = read_mmn_overlaps(path)

        np.testing.assert_array_equal(result["neigh"], neighbours)
        np.testing.assert_array_equal(result["G"], reciprocal_shifts)
        np.testing.assert_allclose(result["M"], matrices, atol=0, rtol=0)

    def test_reads_gzip_mmn_through_plain_name_fallback(self):
        matrices = np.array([[[[1 + 2j]]]])
        neighbours = np.zeros((1, 1), dtype=int)
        reciprocal_shifts = np.zeros((1, 1, 3), dtype=int)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.mmn"
            _write_mmn(path, matrices, neighbours, reciprocal_shifts)
            compressed = Path(f"{path}.gz")
            compressed.write_bytes(gzip.compress(path.read_bytes()))
            path.unlink()
            result = read_mmn_overlaps(path)

        np.testing.assert_array_equal(result["neigh"], neighbours)
        np.testing.assert_array_equal(result["G"], reciprocal_shifts)
        np.testing.assert_allclose(result["M"], matrices, atol=0, rtol=0)

    def test_connection_pairs_each_overlap_with_its_b_vector_and_phase(self):
        b_vectors = np.array(
            [
                [1, 0, 0],
                [-1, 0, 0],
                [0, 1, 0],
                [0, -1, 0],
                [0, 0, 1],
                [0, 0, -1],
            ]
        )
        reciprocal_shifts = b_vectors[None]
        neighbours = np.zeros((1, 6), dtype=int)
        overlaps = np.array(
            [1 + 0.1j, 2 - 0.2j, 3 + 0.3j, 4 - 0.4j, 5 + 0.5j, 6 - 0.6j]
        )
        matrices = overlaps.reshape(1, 6, 1, 1)
        center = np.array([[0.2, -0.1, 0.3]])
        checkpoint = {
            "num_bands": 1,
            "num_kpts": 1,
            "nntot": 6,
            "num_wann": 1,
            "kpt_red": np.zeros((1, 3)),
            "recip_lattice": np.eye(3),
            "v_matrix": np.ones((1, 1, 1), dtype=complex),
            "wannier_centres": center,
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.mmn"
            _write_mmn(path, matrices, neighbours, reciprocal_shifts)
            result = connection_from_mmn(checkpoint, path)

        weights = np.full(6, 0.5)
        expected = 1j * np.einsum(
            "b,b,ba,b->a",
            overlaps,
            weights,
            b_vectors,
            np.exp(1j * b_vectors.dot(center[0])),
        )
        np.testing.assert_allclose(result["A"][0, 0, 0], expected, atol=1e-14)
        self.assertEqual(result["completeness_error"], 0.0)

    def test_finite_difference_weights_satisfy_completeness(self):
        b_vectors = np.array(
            [
                [1, 0, 0],
                [-1, 0, 0],
                [0, 2, 0],
                [0, -2, 0],
                [0, 0, 3],
                [0, 0, -3],
            ],
            dtype=float,
        )
        weights = finite_difference_weights(b_vectors)
        np.testing.assert_allclose(
            np.einsum("b,ba,bc->ac", weights, b_vectors, b_vectors),
            np.eye(3),
            atol=1e-14,
        )


class TestPostw90PositionInput(unittest.TestCase):
    @staticmethod
    def _matrices():
        origin = np.zeros((3, 2, 2), dtype=complex)
        origin[:, 0, 0] = [0.1, 0.2, 0.3]
        origin[:, 1, 1] = [0.4, 0.5, 0.6]
        positive = np.zeros_like(origin)
        positive[:, 0, 1] = [0.02 + 0.03j, -0.04j, 0.05]
        negative = np.conj(np.transpose(positive, (0, 2, 1)))
        return {
            (0, 0, 0): origin,
            (1, 0, 0): positive,
            (-1, 0, 0): negative,
        }

    def test_reads_pair_weighted_blocks_and_orientation(self):
        matrices = self._matrices()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test_r_full.dat"
            _write_r_full(path, matrices)
            result = read_postw90_position_matrix(path)

        self.assertTrue(result["transl_inv_full"])
        self.assertTrue(result["use_ws_distance"])
        self.assertTrue(result["effective_ws_weights"])
        self.assertEqual(result["num_R"], 3)
        for R, expected in matrices.items():
            np.testing.assert_array_equal(result["position_R"][R], expected)

    def test_installs_effective_blocks_with_unit_degeneracy(self):
        matrices = self._matrices()

        class Model:
            norb = 2
            _ham_r = {(0, 0, 0): np.eye(2, dtype=complex)}
            _ham_deg = {(0, 0, 0): 1.0}
            _pos_r = None

        model = Model()
        centers = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test_r_full.dat"
            _write_r_full(path, matrices)
            result = install_postw90_position_matrix(
                model, path, centers_cartesian=centers
            )

        self.assertEqual(result["zero_H_blocks_added"], 2)
        self.assertEqual(set(model._ham_r), set(matrices))
        self.assertEqual(set(model._pos_r), set(matrices))
        self.assertTrue(all(value == 1.0 for value in model._ham_deg.values()))
        self.assertEqual(result["pair_hermiticity_relative"], 0.0)

    def test_centre_check_uses_its_own_angstrom_tolerance(self):
        """The centre check is in Angstrom and must not share the Hermiticity bound.

        postw90 accumulates the R=0 diagonal from the transl_inv_full neighbour
        sum while wannier90.x computes the centres from its own expression, so
        the two agree only to k-sum rounding.  Real exports land near 1e-8 A.
        """
        matrices = self._matrices()

        def model():
            class Model:
                norb = 2
                _ham_r = {(0, 0, 0): np.eye(2, dtype=complex)}
                _ham_deg = {(0, 0, 0): 1.0}
                _pos_r = None

            return Model()

        origin = matrices[(0, 0, 0)]
        exact = np.real(origin[:, np.arange(2), np.arange(2)]).T
        offset = exact + 1e-7

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test_r_full.dat"
            _write_r_full(path, matrices)

            # Inside the default 1e-6 A budget, outside a 1e-8 one.
            result = install_postw90_position_matrix(
                model(), path, centers_cartesian=offset
            )
            self.assertLess(result["centres_max_error"], 1e-6)
            with self.assertRaisesRegex(ValueError, "R=0 diagonal differs"):
                install_postw90_position_matrix(
                    model(), path, centers_cartesian=offset,
                    centres_tolerance=1e-8,
                )

            # A tight Hermiticity bound must not reject an acceptable centre.
            install_postw90_position_matrix(
                model(), path, centers_cartesian=offset, tolerance=1e-14
            )

    def test_rejects_ordinary_scalar_degeneracy_export(self):
        matrices = self._matrices()

        class Model:
            norb = 2
            _ham_r = {(0, 0, 0): np.eye(2, dtype=complex)}
            _ham_deg = {(0, 0, 0): 1.0}

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test_r_full.dat"
            _write_r_full(
                path,
                matrices,
                use_ws_distance=False,
                effective_ws_weights=False,
            )
            with self.assertRaisesRegex(ValueError, "requires use_ws_distance=T"):
                install_postw90_position_matrix(Model(), path)

    def test_rejects_inconsistent_weight_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test_r_full.dat"
            _write_r_full(
                path,
                self._matrices(),
                use_ws_distance=True,
                effective_ws_weights=False,
            )
            with self.assertRaisesRegex(ValueError, "inconsistent use_ws_distance"):
                read_postw90_position_matrix(path)


class TestMMNPositionInstallation(unittest.TestCase):
    def test_preserves_position_support_outside_hamiltonian_support(self):
        origin = np.zeros((3, 2, 2), dtype=complex)
        origin[:, 0, 0] = [0.1, 0.2, 0.3]
        origin[:, 1, 1] = [0.4, 0.5, 0.6]
        positive = np.zeros_like(origin)
        positive[:, 0, 1] = [0.02 + 0.03j, -0.04j, 0.05]
        negative = np.conj(np.transpose(positive, (0, 2, 1)))
        matrices = {
            (-1, 0, 0): negative,
            (0, 0, 0): origin,
            (1, 0, 0): positive,
        }
        R_vectors = np.array(sorted(matrices), dtype=int)
        A_R = np.array([matrices[tuple(R)].transpose(1, 2, 0) for R in R_vectors])
        centers = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])

        class Model:
            norb = 2

            def __init__(self):
                self._ham_r = {(0, 0, 0): np.eye(2, dtype=complex)}
                self._ham_deg = {(0, 0, 0): 1.0}
                self._pos_r = {(0, 0, 0): origin.copy()}

        model = Model()
        position_R, diagnostics = _position_matrix_for_model(
            model,
            {"AA": A_R, "iRvec": R_vectors, "wcc_cart": centers},
        )

        self.assertEqual(set(position_R), set(matrices))
        self.assertEqual(set(model._ham_r), set(matrices))
        self.assertEqual(diagnostics["zero_H_blocks_added"], 2)
        np.testing.assert_array_equal(model._ham_r[(-1, 0, 0)], np.zeros((2, 2)))
        np.testing.assert_array_equal(model._ham_r[(1, 0, 0)], np.zeros((2, 2)))
        self.assertGreater(diagnostics["weight_outside_ham_R"], 0.0)


if __name__ == "__main__":
    unittest.main()


class TestCacheProvenanceWithDeletedInputs(unittest.TestCase):
    """A warm A(R) cache must survive deleting the multi-GB .mmn.

    The cache already holds the finished matrix; the raw inputs only carry
    provenance.  Reusing the recorded stamp for a deleted file must not weaken
    staleness detection on the files that are still there.
    """

    OPTIONS = {"transl_inv_JM": False, "transl_inv_MV": False, "ws_dist_tol": 1e-5}

    def _directory(self, root: Path) -> Path:
        root.mkdir(parents=True, exist_ok=True)
        (root / "test.chk").write_bytes(b"checkpoint")
        (root / "test.mmn").write_bytes(b"overlaps")
        return root

    def test_recorded_stamp_stands_in_for_a_deleted_mmn(self):
        with tempfile.TemporaryDirectory() as temporary:
            d = self._directory(Path(temporary) / "run")
            full = _A_R_cache_metadata(d, d, "test", self.OPTIONS)
            self.assertIn("mmn", full["inputs"])

            (d / "test.mmn").unlink()
            with self.assertRaises(FileNotFoundError):
                _A_R_cache_metadata(d, d, "test", self.OPTIONS)

            relaxed = _A_R_cache_metadata(d, d, "test", self.OPTIONS, recorded=full)
            self.assertEqual(relaxed, full)

    def test_a_surviving_input_is_still_checked_for_staleness(self):
        with tempfile.TemporaryDirectory() as temporary:
            d = self._directory(Path(temporary) / "run")
            full = _A_R_cache_metadata(d, d, "test", self.OPTIONS)
            (d / "test.mmn").unlink()

            # The .chk changing must still invalidate, even though the .mmn is
            # being taken on trust.  Otherwise the relaxation would hide a
            # genuinely stale cache.
            (d / "test.chk").write_bytes(b"a different checkpoint entirely")
            relaxed = _A_R_cache_metadata(d, d, "test", self.OPTIONS, recorded=full)
            self.assertNotEqual(relaxed, full)
            self.assertEqual(relaxed["inputs"]["mmn"], full["inputs"]["mmn"])

    def test_options_still_invalidate_when_an_input_is_absent(self):
        with tempfile.TemporaryDirectory() as temporary:
            d = self._directory(Path(temporary) / "run")
            full = _A_R_cache_metadata(d, d, "test", self.OPTIONS)
            (d / "test.mmn").unlink()
            changed = _A_R_cache_metadata(
                d, d, "test", {**self.OPTIONS, "transl_inv_MV": True}, recorded=full
            )
            self.assertNotEqual(changed, full)
