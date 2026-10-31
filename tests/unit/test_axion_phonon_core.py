"""axion_phonon pieces that need no data: the tube quadrature, the Wigner-Seitz fold, the b-shell weights, .amn I/O."""
from __future__ import annotations

import itertools
import tempfile
import unittest
from pathlib import Path

import numpy as np

from axion_phonon import tube, wannier
from axion_phonon.io import wannier90


class TestTube(unittest.TestCase):
    def test_constant_integrand(self):
        """c2 = 1 integrates to 1 (the BZ volume) at every in-plane level and plane count."""
        kl = tube.klist(nk=8, m=4, half=1, n3=8)
        levels, n3s, table = tube.integrate(np.ones((len(kl["k"]), 2)), kl)
        self.assertEqual(levels, [1, 2, 4])
        self.assertEqual(n3s, [8, 4, 2])
        np.testing.assert_allclose(table, 1.0, atol=1e-12)

    def test_tube_is_kz_line(self):
        """The box points nearest the tube include k = t (1, 1, 1) in every plane."""
        kl = tube.klist(nk=8, m=4, half=1, n3=8)
        for p in range(8):
            k = kl["k"][kl["plane"] == p]
            t = p / 8
            self.assertTrue(np.any(np.all(np.abs(k - t) < 1e-12, axis=1)))


class TestFold(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.lattice = np.array([[1.0, 0.1, 0.0], [0.0, 1.2, 0.2], [0.3, 0.0, 2.0]])
        self.grid = (3, 3, 2)
        self.positions = rng.uniform(0, 1, (3, 3)) @ self.lattice
        self.reps = wannier.class_representatives(self.grid, self.lattice)
        self.H = {R: rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for R in self.reps}

    def test_one_representative_per_class(self):
        classes = {tuple(np.mod(R, self.grid)) for R in self.reps}
        self.assertEqual(len(classes), np.prod(self.grid))

    def test_class_sums_conserved(self):
        """A fold only redistributes each class over its images."""
        folded = wannier.fold(self.H, self.lattice, self.positions, self.grid)
        sums = {}
        for R, block in folded.items():
            c = tuple(np.mod(R, self.grid))
            sums[c] = sums.get(c, 0) + block
        for R, block in self.H.items():
            np.testing.assert_allclose(sums[tuple(np.mod(R, self.grid))], block, atol=1e-13)

    def test_entry_on_nearest_image(self):
        """Without ties, each entry lands on the image minimising |R + tau_j - tau_i|."""
        folded = wannier.fold(self.H, self.lattice, self.positions, self.grid, tolerance=0.0)
        L = np.array(list(itertools.product((-1, 0, 1), repeat=3))) * np.array(self.grid)
        for R0 in self.reps[:5]:
            for i, j in itertools.product(range(3), repeat=2):
                cand = np.asarray(R0) + L
                d = np.linalg.norm(cand @ self.lattice + self.positions[j] - self.positions[i], axis=1)
                if np.sort(d)[1] - d.min() < 1e-9:      # a tie: shared between images
                    continue
                best = tuple(int(x) for x in cand[d.argmin()])
                self.assertAlmostEqual(folded[best][i, j], self.H[R0][i, j], places=13)


class TestShells(unittest.TestCase):
    def test_cubic_completeness(self):
        b = np.array([v for v in itertools.product((-1, 0, 1), repeat=3) if sum(map(abs, v)) == 1], float) * 0.5
        w = wannier.finite_difference_weights(b)
        np.testing.assert_allclose(np.einsum("b,ba,bc->ac", w, b, b), np.eye(3), atol=1e-12)


class TestAmn(unittest.TestCase):
    def test_write_layout(self):
        A = np.arange(2 * 3 * 2).reshape(2, 3, 2) * (1 + 0.5j)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "x.amn"
            wannier90.write_amn(path, A, "test")
            rows = np.loadtxt(path, skiprows=2)
        for m, n, k, re, im in rows:
            self.assertEqual(A[int(k) - 1, int(m) - 1, int(n) - 1], re + 1j * im)


if __name__ == "__main__":
    unittest.main()
