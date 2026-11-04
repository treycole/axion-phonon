"""axion_phonon/paw/ pieces that need no QE data: QE's real Y_lm, simpson, G matching.

The data check (agreement with the code that made the corrected .amn, S-orthonormality of QE's wavefunctions) is
archive/validation/restructure/accept_paw.py.
"""
from __future__ import annotations

import unittest

import numpy as np

from axion_phonon.io.qe import Pseudo
from axion_phonon.paw.cross_structure import angular_gradient_tables, gradient_matrix, match_g
from axion_phonon.paw.qe_port import channels, simpson, ylm_real


class TestYlm(unittest.TestCase):
    def test_orthonormal_on_sphere(self):
        lmax = 3
        x, w = np.polynomial.legendre.leggauss(24)                      # cos(theta)
        phi = np.linspace(0, 2 * np.pi, 48, endpoint=False)
        ct, ph = np.meshgrid(x, phi, indexing="ij")
        st = np.sqrt(1 - ct ** 2)
        g = np.stack([st * np.cos(ph), st * np.sin(ph), ct], axis=-1).reshape(-1, 3) * 1.7
        weights = (w[:, None] * np.full(len(phi), 2 * np.pi / len(phi))[None, :]).ravel()
        y = ylm_real(g, lmax)
        gram = (y * weights[:, None]).T @ y
        np.testing.assert_allclose(gram, np.eye((lmax + 1) ** 2), atol=1e-12)

    def test_qe_order_l1(self):
        # QE ylmr2: lm = 2 -> z, 3 -> -x, 4 -> -y (times sqrt(3/4pi))
        g = np.array([[0.3, -0.5, 0.8]])
        y = ylm_real(g, 1)[0]
        c = np.sqrt(3 / (4 * np.pi)) / np.linalg.norm(g)
        np.testing.assert_allclose(y[1:], c * np.array([0.8, -0.3, 0.5]), atol=1e-14)


class TestGradientMatrix(unittest.TestCase):
    def test_against_3d_quadrature(self):
        # <R_i Y_i| d_a |R_j Y_j> from gradient_matrix vs a direct 3D integral (central differences in x, y, z)
        x = np.arange(-9.0, 2.2, 0.0125)                               # QE-like log mesh r = e^x
        r = np.exp(x)
        rab = r * 0.0125
        lll = np.array([0, 1, 2])
        u = [r ** (l + 1) * np.exp(-a * r ** 2) for l, a in zip(lll, (1.0, 0.7, 0.9))]
        n = len(r)
        ps = Pseudo(r=r, rab=rab, beta=[np.zeros(n)] * 3, lll=lll, jjj=None, kkbeta=n, qqq=np.zeros((3, 3)),
                    qfunc0=np.zeros((3, 3, n)), aewfc=u, pswfc=[np.zeros(n)] * 3)
        ch = channels(ps, lspinorb=False)
        G, H = angular_gradient_tables(2)
        D = gradient_matrix(ps, ch, G, H)
        # reference: radial Gauss-Legendre on [0, 6] x angular Gauss-Legendre/uniform
        rq, rw = np.polynomial.legendre.leggauss(80)
        rq, rw = 3.0 * (rq + 1), 3.0 * rw
        ct, cw = np.polynomial.legendre.leggauss(12)
        ph = np.linspace(0, 2 * np.pi, 24, endpoint=False)
        CT, PH = np.meshgrid(ct, ph, indexing="ij")
        ST = np.sqrt(1 - CT ** 2)
        rhat = np.stack([ST * np.cos(PH), ST * np.sin(PH), CT], -1).reshape(-1, 3)
        aw = (cw[:, None] * np.full(len(ph), 2 * np.pi / len(ph))[None]).ravel()
        pts = (rq[:, None, None] * rhat[None]).reshape(-1, 3)
        wts = (rw[:, None] * rq[:, None] ** 2 * aw[None]).ravel()
        radial = [lambda s, l=l, a=a: s ** l * np.exp(-a * s ** 2) for l, a in zip(lll, (1.0, 0.7, 0.9))]

        def field(p, ih):                                               # R(|p|) Y(p) for channel ih
            s = np.linalg.norm(p, axis=1)
            return radial[ch.indv[ih]](s) * ylm_real(p, 2)[:, ch.nhtolm[ih]]

        h = 1e-5
        for ih in range(len(ch.indv)):
            fi = field(pts, ih)
            for jh in range(len(ch.indv)):
                if abs(ch.nhtol[ih] - ch.nhtol[jh]) != 1:
                    continue
                for a in range(3):
                    e = np.zeros(3)
                    e[a] = h
                    dfj = (field(pts + e, jh) - field(pts - e, jh)) / (2 * h)
                    ref = np.sum(wts * fi * dfj)
                    self.assertAlmostEqual(D[a, ih, jh], ref, delta=2e-6)
        np.testing.assert_allclose(D, -D.transpose(0, 2, 1), atol=2e-6)   # d is anti-Hermitian


class TestNumerics(unittest.TestCase):
    def test_simpson_odd_and_even(self):
        for n in (101, 100):
            x = np.linspace(0, 2, n)
            rab = np.full(n, x[1] - x[0])
            self.assertAlmostEqual(simpson(x ** 3, rab), 4.0, places=6)

    def test_match_g(self):
        rng = np.random.default_rng(0)
        m1 = rng.integers(-20, 20, size=(500, 3))
        m1 = np.unique(m1, axis=0)
        perm = rng.permutation(len(m1))
        i1, i2 = match_g(m1, m1[perm])
        np.testing.assert_array_equal(m1[i1], m1[perm][i2])
        self.assertEqual(len(i1), len(m1))


if __name__ == "__main__":
    unittest.main()
