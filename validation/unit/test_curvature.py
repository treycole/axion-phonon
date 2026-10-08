"""Tests for ``modules/curvature.py`` and ``modules/axion.py``.

Three groups, none needing first-principles data:

* analytic covariance checks: each piece of Omega, and the total, against exactly
  solvable moving frames (``TestAnalyticCovariance``);
* small algebraic properties: Hermiticity, antisymmetry, the Kubo equation, eq 38 ==
  eq 39 (``TestAlgebra``);
* the k-space inputs: Bloch sums against a brute-force sum over R, and the
  parametric connection ``A_beta`` against its own central difference
  (``TestBlochSums``, ``TestParametricConnection``).

Agreement with PythTB and WannierBerri on real Wannier90 output is checked in
``validation/external_curvature``.
"""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from modules.axion import (
    beta_terms,
    c2_density,
    dtheta,
    integrate_c2,
    parametric_connection,
    prepare_lambda,
)
from modules.curvature import (
    BetaTerms,
    Fields,
    connection,
    fields,
    omega,
    omega_cross,
    omega_external,
    omega_internal,
    position_terms,
)
TIGHT = dict(atol=1e-12, rtol=1e-12)


def _analytic_two_level_data(kpoints: np.ndarray):
    r"""Return a gapped model whose occupied spinor is known analytically."""

    kpoints = np.asarray(kpoints, dtype=float)
    x = 2.0 * np.pi * kpoints[:, 0]
    y = 2.0 * np.pi * kpoints[:, 1]
    theta_0 = 0.47
    amplitude = 0.23
    theta = theta_0 + amplitude * np.sin(x)
    dtheta_dx = 2.0 * np.pi * amplitude * np.cos(x)
    dphase_dy = 2.0 * np.pi

    occupied_spinor = np.stack(
        [
            np.cos(theta),
            np.exp(1j * y) * np.sin(theta),
        ],
        axis=1,
    )
    dspinor_dx = np.stack(
        [
            -np.sin(theta) * dtheta_dx,
            np.exp(1j * y) * np.cos(theta) * dtheta_dx,
        ],
        axis=1,
    )
    dspinor_dy = np.stack(
        [
            np.zeros_like(theta),
            1j * dphase_dy * np.exp(1j * y) * np.sin(theta),
        ],
        axis=1,
    )

    projector = np.einsum(
        "ki,kj->kij", occupied_spinor, occupied_spinor.conj()
    )
    hamiltonian = np.eye(2)[None] - 2.0 * projector
    dH = np.zeros((3, len(kpoints), 2, 2), dtype=np.complex128)
    for axis, derivative in enumerate((dspinor_dx, dspinor_dy)):
        dH[axis] = -2.0 * (
            np.einsum("ki,kj->kij", derivative, occupied_spinor.conj())
            + np.einsum(
                "ki,kj->kij", occupied_spinor, derivative.conj()
            )
        )

    omega_xy = (
        -2.0
        * dphase_dy
        * np.sin(theta)
        * np.cos(theta)
        * dtheta_dx
    )
    expected = np.zeros(
        (3, 3, len(kpoints), 1, 1), dtype=np.complex128
    )
    expected[0, 1, :, 0, 0] = omega_xy
    expected[1, 0, :, 0, 0] = -omega_xy
    return hamiltonian, dH, theta, expected


def _pauli_rotation_and_derivative(
    generator: np.ndarray, angle: float, derivative: float
):
    r"""Exponentiate a Pauli block with eigenvalues ``(-1, 0, +1)``."""

    identity = np.eye(len(generator), dtype=np.complex128)
    active_projector = generator @ generator
    rotation = (
        identity
        + (np.cos(angle / 2.0) - 1.0) * active_projector
        - 1j * np.sin(angle / 2.0) * generator
    )
    drotation = derivative * (
        -0.5 * np.sin(angle / 2.0) * active_projector
        - 0.5j * np.cos(angle / 2.0) * generator
    )
    return rotation, drotation


def _rotate_frame(
    kpoints: np.ndarray,
    hamiltonian: np.ndarray,
    dH: np.ndarray,
    generator_x: np.ndarray,
    generator_y: np.ndarray,
):
    r"""Apply an analytic moving basis, returning exact H, dH, A, and curl."""

    rotated_hamiltonian = np.empty_like(hamiltonian)
    rotated_dH = np.zeros_like(dH)
    dimension = hamiltonian.shape[-1]
    connection = np.zeros(
        (len(kpoints), 3, dimension, dimension), dtype=np.complex128
    )

    for ik, (kx, ky) in enumerate(kpoints[:, :2]):
        x = 2.0 * np.pi * kx
        y = 2.0 * np.pi * ky
        alpha = 0.71 * np.sin(x)
        beta = 0.53 * np.sin(y)
        dalpha_dx = 0.71 * 2.0 * np.pi * np.cos(x)
        dbeta_dy = 0.53 * 2.0 * np.pi * np.cos(y)

        rotation_x, drotation_x = _pauli_rotation_and_derivative(
            generator_x, alpha, dalpha_dx
        )
        rotation_y, drotation_y = _pauli_rotation_and_derivative(
            generator_y, beta, dbeta_dy
        )
        frame = rotation_x @ rotation_y
        dframe_dx = drotation_x @ rotation_y
        dframe_dy = rotation_x @ drotation_y

        rotated_hamiltonian[ik] = (
            frame.conj().T @ hamiltonian[ik] @ frame
        )
        for axis, dframe in enumerate((dframe_dx, dframe_dy)):
            rotated_dH[axis, ik] = (
                dframe.conj().T @ hamiltonian[ik] @ frame
                + frame.conj().T @ dH[axis, ik] @ frame
                + frame.conj().T @ hamiltonian[ik] @ dframe
            )
            connection[ik, axis] = 1j * frame.conj().T @ dframe

    # A = i G^dagger dG is a pure-gauge connection in the complete frame,
    # so dA = i[A_x, A_y] exactly.
    curl = np.zeros_like(connection)
    curl[:, 2] = 1j * (
        connection[:, 0] @ connection[:, 1]
        - connection[:, 1] @ connection[:, 0]
    )
    return rotated_hamiltonian, rotated_dH, connection, curl


def _rotate_two_level_frame(
    kpoints: np.ndarray,
    hamiltonian: np.ndarray,
    dH: np.ndarray,
):
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_y = np.array([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)
    return _rotate_frame(
        kpoints, hamiltonian, dH, sigma_x, sigma_y
    )


def _fields(H, dH, A=None, curl=None):
    """Old-style arrays (A and curl with k first) -> the new direction-first Fields."""
    nk, J = H.shape[:2]
    zeros = np.zeros((3, nk, J, J), dtype=complex)
    return Fields(
        H=H,
        dH=dH[:3],
        A=zeros if A is None else A.swapaxes(0, 1),
        curl=zeros if curl is None else curl.swapaxes(0, 1),
    )


def _beta(dH, mixed_curl=None, A_beta=None):
    nk, J = dH.shape[:2]
    curl = (
        np.zeros((3, nk, J, J), dtype=complex)
        if mixed_curl is None
        else mixed_curl.swapaxes(0, 1)
    )
    return BetaTerms(dH=dH, curl=curl, A=A_beta)


def _pieces(c):
    return omega_internal(c), omega_cross(c), omega_external(c), omega(c)


class TestAnalyticCovariance(unittest.TestCase):
    def setUp(self):
        self.kpoints = np.array(
            [[0.13, 0.19, 0.0], [0.31, 0.37, 0.0], [0.67, 0.73, 0.0]]
        )
        (self.H, self.dH, self.theta, self.expected) = _analytic_two_level_data(
            self.kpoints
        )

    def test_curvature_moves_from_internal_to_external_when_J_equals_M(self):
        internal, cross, external, total = _pieces(
            connection(_fields(self.H, self.dH), n_occ=1)
        )
        nk = len(self.kpoints)
        A = np.zeros((nk, 3, 1, 1), dtype=complex)
        A[:, 1, 0, 0] = -2.0 * np.pi * np.sin(self.theta) ** 2
        curl = np.zeros_like(A)
        curl[:, 2, 0, 0] = self.expected[0, 1, :, 0, 0]
        occ_internal, occ_cross, occ_external, occ_total = _pieces(
            connection(
                _fields(-np.ones((nk, 1, 1), complex), np.zeros((3, nk, 1, 1), complex), A, curl),
                n_occ=1,
            )
        )
        np.testing.assert_allclose(internal, self.expected, **TIGHT)
        np.testing.assert_array_equal(cross, 0.0)
        np.testing.assert_array_equal(external, 0.0)
        np.testing.assert_array_equal(occ_internal, 0.0)
        np.testing.assert_array_equal(occ_cross, 0.0)
        np.testing.assert_allclose(occ_external, self.expected, **TIGHT)
        np.testing.assert_allclose(total, occ_total, **TIGHT)

    def test_total_is_invariant_when_all_three_terms_are_nonzero(self):
        H, dH, A, curl = _rotate_two_level_frame(self.kpoints, self.H, self.dH)
        internal, cross, external, total = _pieces(connection(_fields(H, dH, A, curl), 1))
        for piece in (internal, cross, external):
            self.assertGreater(np.max(np.abs(piece)), 1e-2)
        np.testing.assert_allclose(total, self.expected, **TIGHT)

    def test_nonabelian_covariance_including_occupied_commutator(self):
        nk = len(self.kpoints)
        H = np.zeros((nk, 3, 3), dtype=complex)
        H[:, 0, 0] = -2.0
        H[:, 1:, 1:] = self.H
        dH = np.zeros((3, nk, 3, 3), dtype=complex)
        dH[:, :, 1:, 1:] = self.dH
        gx = np.zeros((3, 3), dtype=complex)
        gy = np.zeros((3, 3), dtype=complex)
        gx[:2, :2] = [[0.0, 1.0], [1.0, 0.0]]
        gy[:2, :2] = [[0.0, -1j], [1j, 0.0]]
        Hr, dHr, A, curl = _rotate_frame(self.kpoints, H, dH, gx, gy)

        internal, cross, external, total = _pieces(connection(_fields(Hr, dHr, A, curl), 2))
        expected = np.zeros_like(total)
        expected[0, 1, :, 1, 1] = self.expected[0, 1, :, 0, 0]
        expected[1, 0, :, 1, 1] = self.expected[1, 0, :, 0, 0]
        for piece in (internal, cross, external):
            self.assertGreater(np.max(np.abs(piece)), 1e-2)
        np.testing.assert_allclose(total, expected, **TIGHT)

    def test_mixed_curvature_moves_from_internal_to_external(self):
        """The same geometry may live in dH_beta or in A_beta."""
        nk = len(self.kpoints)
        mixed_dH = np.zeros((4, nk, 2, 2), dtype=complex)
        mixed_dH[0], mixed_dH[3] = self.dH[0], self.dH[1]
        expected = np.zeros((4, 4, nk, 1, 1), dtype=complex)
        expected[0, 3], expected[3, 0] = self.expected[0, 1], self.expected[1, 0]

        internal, cross, external, total = _pieces(
            connection(_fields(self.H, mixed_dH), 1, _beta(mixed_dH[3]))
        )

        mixed_curl = np.zeros((nk, 3, 1, 1), dtype=complex)
        mixed_curl[:, 0, 0, 0] = self.expected[0, 1, :, 0, 0]
        A_beta = np.zeros((nk, 1, 1), dtype=complex)
        A_beta[:, 0, 0] = -2.0 * np.pi * np.sin(self.theta) ** 2
        one_band = _fields(-np.ones((nk, 1, 1), complex), np.zeros((4, nk, 1, 1), complex))
        occ_internal, occ_cross, occ_external, occ_total = _pieces(
            connection(one_band, 1, _beta(np.zeros((nk, 1, 1), complex), mixed_curl, A_beta))
        )
        np.testing.assert_allclose(internal, expected, **TIGHT)
        np.testing.assert_array_equal(cross, 0.0)
        np.testing.assert_array_equal(external, 0.0)
        np.testing.assert_array_equal(occ_internal, 0.0)
        np.testing.assert_array_equal(occ_cross, 0.0)
        np.testing.assert_allclose(occ_external, expected, **TIGHT)
        np.testing.assert_allclose(total, occ_total, **TIGHT)

    def test_mixed_total_is_covariant_with_all_terms_nonzero(self):
        """Reindex the exact moving-frame xy model as an x-beta model."""
        H, dHr, A, curl = _rotate_two_level_frame(self.kpoints, self.H, self.dH)
        nk = len(self.kpoints)
        mixed_dH = np.zeros((4, nk, 2, 2), dtype=complex)
        mixed_dH[0], mixed_dH[3] = dHr[0], dHr[1]
        spatial_A = np.zeros_like(A)
        spatial_A[:, 0] = A[:, 0]
        mixed_curl = np.zeros_like(A)
        mixed_curl[:, 0] = curl[:, 2]

        internal, cross, external, total = _pieces(
            connection(
                _fields(H, mixed_dH, spatial_A),
                1,
                _beta(mixed_dH[3], mixed_curl, A_beta=A[:, 1]),
            )
        )
        expected = np.zeros_like(total)
        expected[0, 3], expected[3, 0] = self.expected[0, 1], self.expected[1, 0]
        for piece in (internal, cross, external):
            self.assertGreater(np.max(np.abs(piece)), 1e-2)
        np.testing.assert_allclose(total, expected, **TIGHT)


class TestAlgebra(unittest.TestCase):
    def _random_connection(self, seed=3, J=6, n_occ=3, nk=4, with_beta=True):
        rng = np.random.default_rng(seed)

        def hermitian(*shape):
            m = rng.normal(size=shape) + 1j * rng.normal(size=shape)
            return m + m.conj().swapaxes(-1, -2)

        H = hermitian(nk, J, J)
        dH = hermitian(4, nk, J, J)
        A = hermitian(nk, 3, J, J)
        curl = hermitian(nk, 3, J, J)
        mixed_curl = hermitian(nk, 3, J, J)
        A_beta = hermitian(nk, J, J)
        return H, dH, A, curl, mixed_curl, A_beta

    def test_pieces_are_hermitian_and_antisymmetric_and_sum_to_omega(self):
        H, dH, A, curl, mixed_curl, A_beta = self._random_connection()
        c = connection(_fields(H, dH, A, curl), 3, _beta(dH[3], mixed_curl, A_beta))
        internal, cross, external, total = _pieces(c)
        for piece in (internal, cross, external, total):
            np.testing.assert_allclose(piece, -piece.swapaxes(0, 1), atol=1e-10)
            np.testing.assert_allclose(piece, piece.conj().swapaxes(-1, -2), atol=1e-10)
            self.assertGreater(np.max(np.abs(piece)), 0.0)
        np.testing.assert_allclose(total, internal + cross + external, atol=0.0, rtol=0.0)

    def test_without_external_terms_only_the_internal_piece_survives(self):
        H, dH, *_ = self._random_connection()
        internal, cross, external, total = _pieces(connection(_fields(H, dH), 3))
        np.testing.assert_array_equal(cross, 0.0)
        np.testing.assert_array_equal(external, 0.0)
        np.testing.assert_allclose(total, internal)

    def test_internal_matches_the_kubo_equation(self):
        energies = np.array([-1.0, 0.5, 2.0])
        H = np.diag(energies)[None].astype(complex)
        rng = np.random.default_rng(7)
        dH = np.empty((4, 1, 3, 3), dtype=complex)
        for mu in range(4):
            m = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
            dH[mu, 0] = m + m.conj().T
        c = connection(_fields(H, dH), 1, _beta(dH[3]))
        expected = np.zeros_like(omega_internal(c))
        for mu in range(4):
            for nu in range(4):
                for empty in (1, 2):
                    expected[mu, nu, 0, 0, 0] += 1j * (
                        dH[mu, 0, 0, empty] * dH[nu, 0, empty, 0]
                        - dH[nu, 0, 0, empty] * dH[mu, 0, empty, 0]
                    ) / (energies[0] - energies[empty]) ** 2
        np.testing.assert_allclose(omega_internal(c), expected, atol=1e-12)
        self.assertTrue(c.frozen_gauge)

    def test_frozen_gauge_is_recorded_only_when_beta_has_no_A(self):
        H, dH, A, curl, mixed_curl, A_beta = self._random_connection()
        f = _fields(H, dH, A, curl)
        self.assertFalse(connection(f, 3).frozen_gauge)
        self.assertTrue(connection(f, 3, _beta(dH[3], mixed_curl)).frozen_gauge)
        self.assertFalse(connection(f, 3, _beta(dH[3], mixed_curl, A_beta)).frozen_gauge)

    def test_eq38_equals_eq39(self):
        H, dH, A, curl, mixed_curl, A_beta = self._random_connection(seed=11)
        Om = omega(connection(_fields(H, dH, A, curl), 3, _beta(dH[3], mixed_curl, A_beta)))
        # (39): (1/2 pi) sum_l tr[ B~_l B_{l beta} ],  B~_l = eps^{lij} B_ij / 2
        B_tilde = np.stack(
            [Om[1, 2], Om[2, 0], Om[0, 1]]  # B~_x = B_yz, B~_y = B_zx, B~_z = B_xy
        )
        c2_39 = sum(
            np.einsum("kmn,knm->k", B_tilde[l], Om[l, 3]) for l in range(3)
        ) / (2.0 * np.pi)
        np.testing.assert_allclose(c2_density(Om), c2_39, atol=1e-10)

    def test_levi_civita_matches_pythtb(self):
        from pythtb.utils import levi_civita

        from modules.axion import _EPS4

        np.testing.assert_array_equal(_EPS4, levi_civita(4, 4))

    def test_integral_of_a_constant_is_that_constant(self):
        nk = 5
        grid = np.full((nk, nk, nk, 2), 3.0 + 0j)
        total, simpson_total, d3k = integrate_c2(grid)
        np.testing.assert_allclose(total, 3.0)
        np.testing.assert_allclose(simpson_total, 3.0)
        self.assertAlmostEqual(d3k, 1.0 / nk**3)


def _fake_model(seed=5, J=4, upper_energy=None, positions=True):
    """A stand-in for a PythTB model: just what ``position_terms`` and ``fields`` read."""
    rng = np.random.default_rng(seed)
    lat = np.array([[3.0, 0.2, 0.0], [0.4, 3.5, 0.1], [0.0, 0.3, 4.1]])
    recip = 2.0 * np.pi * np.linalg.inv(lat).T
    tau = rng.uniform(-0.3, 0.3, size=(J, 3))
    pos_r, ham_deg = {}, {}
    for R in [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, -1), (-1, -1, 1)]:
        ham_deg[R] = 1
        block = rng.normal(size=(3, J, J)) + 1j * rng.normal(size=(3, J, J))
        pos_r[R] = block if positions else np.zeros_like(block)
    model = SimpleNamespace(
        norb=J,
        lattice=SimpleNamespace(lat_vecs=lat),
        recip_lat_vecs=recip,
        orb_vecs=tau,
        _pos_r=pos_r,
        _ham_deg=ham_deg,
    )
    if upper_energy is not None:  # k-independent two-level model
        def hamiltonian(k, flatten_spin_axis=True):
            H = np.zeros((len(k), J, J), dtype=complex)
            H[:, 0, 0], H[:, 1, 1] = -1.0, upper_energy
            return H

        model.hamiltonian = hamiltonian
        model.velocity = lambda k, flatten_spin_axis=True: np.zeros((3, len(k), J, J), complex)
        model.k_uniform_mesh = lambda mesh, include_endpoints=False: np.array(
            np.meshgrid(*[np.arange(n) / n for n in mesh], indexing="ij")
        ).reshape(3, -1).T
    else:
        model.hamiltonian = lambda k, flatten_spin_axis=True: np.zeros((len(k), J, J), complex)
        model.velocity = lambda k, flatten_spin_axis=True: np.zeros((3, len(k), J, J), complex)
    return model


class TestBlochSums(unittest.TestCase):
    def test_position_terms_subtract_tau_from_the_R0_diagonal(self):
        model = _fake_model()
        pos = position_terms(model, dtype=np.complex128)
        origin = [i for i, R in enumerate(map(tuple, pos.R)) if R == (0, 0, 0)][0]
        raw = model._pos_r[(0, 0, 0)]
        hermitian = 0.5 * (raw + raw.conj().transpose(0, 2, 1))
        centre = np.einsum("ass->sa", hermitian).real  # <0s|r|0s>
        tau_cartesian = model.orb_vecs @ model.lattice.lat_vecs
        np.testing.assert_allclose(
            np.einsum("ass->sa", pos.X[origin]).real, centre - tau_cartesian, atol=1e-12
        )

    def test_A_and_curl_match_a_brute_force_sum_over_R(self):
        """Eqs (8) and (9) written as plain loops, from the raw stored blocks."""
        model = _fake_model()
        J = model.norb
        lat, tau = model.lattice.lat_vecs, model.orb_vecs
        k = np.array([[0.11, 0.23, 0.37], [0.5, 0.0, 0.25], [0.81, 0.4, 0.66]])
        f = fields(model, position_terms(model, dtype=np.complex128), k)

        tau_cartesian = tau @ lat
        A = np.zeros((len(k), 3, J, J), dtype=complex)     # Cartesian components
        curl = np.zeros((len(k), 3, 3, J, J), dtype=complex)  # Cartesian, (mu, nu)
        for R, raw in model._pos_r.items():
            partner = model._pos_r[tuple(-x for x in R)]
            X = 0.5 * (raw + partner.conj().transpose(0, 2, 1))  # Hermitian part
            if R == (0, 0, 0):
                for s in range(J):
                    X[:, s, s] -= tau_cartesian[s]               # X(0) = <0s|r - tau_s|0s>
            for s in range(J):
                for t in range(J):
                    b_red = np.array(R) + tau[t] - tau[s]        # b_st(R), reduced
                    b = b_red @ lat                              # Cartesian
                    phase = np.exp(2j * np.pi * (k @ b_red))
                    A[:, :, s, t] += phase[:, None] * X[:, s, t]
                    for mu in range(3):
                        for nu in range(3):
                            weight = 1j * (b[mu] * X[nu, s, t] - b[nu] * X[mu, s, t])
                            curl[:, mu, nu, s, t] += phase * weight

        # Reduced components follow k = sum_i xi_i b_i, so a covector is A_i = b_i . A
        # and the dual two-form is curl_l = Omega_{(l+1)(l+2)} = b_{l+1}^mu b_{l+2}^nu Omega_{mu nu}
        # (b_i are the reciprocal vectors, 2 pi included).
        recip = model.recip_lat_vecs
        A_red = np.einsum("im,kmst->ikst", recip, A)
        curl_red = np.stack([
            np.einsum("m,n,kmnst->kst", recip[(l + 1) % 3], recip[(l + 2) % 3], curl)
            for l in range(3)
        ])
        np.testing.assert_allclose(f.A, A_red, atol=1e-11)
        np.testing.assert_allclose(f.curl, curl_red, atol=1e-11)


class TestParametricConnection(unittest.TestCase):
    def _setup(self):
        model = _fake_model(seed=9)
        J = model.norb
        rng = np.random.default_rng(2)
        Lambda_R = {
            R: rng.normal(size=(J, J)) + 1j * rng.normal(size=(J, J))
            for R in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, -1)]
        }
        pos = position_terms(model, dtype=np.complex128)
        return model, pos, prepare_lambda(Lambda_R, np.complex128)

    def test_A_beta_is_hermitian(self):
        """Lambda(R) = -Lambda(-R)^dag (34) makes A_beta(k) a Hermitian connection."""
        _, pos, prepared = self._setup()
        k = np.array([[0.11, 0.23, 0.37], [0.5, 0.0, 0.25]])
        A_beta, _ = parametric_connection(pos, prepared, k)
        np.testing.assert_allclose(A_beta, A_beta.conj().transpose(0, 2, 1), atol=1e-12)

    def test_dA_beta_is_the_k_derivative_of_A_beta(self):
        _, pos, prepared = self._setup()
        k = np.array([[0.11, 0.23, 0.37], [0.5, 0.0, 0.25]])
        _, dA_beta = parametric_connection(pos, prepared, k)
        h = 1e-5
        for axis in range(3):
            step = np.zeros(3)
            step[axis] = h
            plus, _ = parametric_connection(pos, prepared, k + step)
            minus, _ = parametric_connection(pos, prepared, k - step)
            np.testing.assert_allclose(dA_beta[axis], (plus - minus) / (2 * h), atol=1e-7)

    def test_beta_terms_frozen_gauge_and_with_lambda(self):
        model, pos, prepared = self._setup()
        J = model.norb
        k = np.array([[0.11, 0.23, 0.37], [0.5, 0.0, 0.25]])
        rng = np.random.default_rng(3)
        zeros = np.zeros((3, len(k), J, J), complex)
        f_base = Fields(H=np.zeros((len(k), J, J), complex), dH=zeros, A=zeros, curl=zeros)
        dA = rng.normal(size=(3, len(k), J, J)) + 0j
        f_mode = Fields(H=f_base.H, dH=zeros, A=dA, curl=zeros)

        frozen = beta_terms(f_base, f_mode, 1.0)        # A_beta = 0: curl = -dA/dbeta (8.5)
        np.testing.assert_array_equal(frozen.curl, -dA)
        self.assertIsNone(frozen.A)

        A_beta, dA_beta = parametric_connection(pos, prepared, k)
        with_lambda = beta_terms(f_base, f_mode, 1.0, A_beta, dA_beta)   # (35)
        np.testing.assert_allclose(with_lambda.curl, dA_beta - dA, atol=1e-13)
        np.testing.assert_array_equal(with_lambda.A, A_beta)


class TestDtheta(unittest.TestCase):
    def test_k_independent_model_gives_zero(self):
        base = _fake_model(J=2, upper_energy=1.0, positions=False)
        mode = _fake_model(J=2, upper_energy=1.1, positions=False)
        result = dtheta(base, mode, nk=2, dbeta=0.01, n_occ=1)
        np.testing.assert_allclose(result["c2_density"], 0.0, atol=0.0)
        np.testing.assert_allclose(result["dtheta"], 0.0, atol=0.0)
        np.testing.assert_allclose(result["minimum_gap"], [2.0, 2.1])

    def test_different_nominal_tau_is_refused(self):
        base = _fake_model(J=2, upper_energy=1.0, positions=False, seed=1)
        mode = _fake_model(J=2, upper_energy=1.1, positions=False, seed=2)  # different tau
        with self.assertRaisesRegex(ValueError, "nominal tau"):
            dtheta(base, mode, nk=2, dbeta=0.01, n_occ=1)


if __name__ == "__main__":
    unittest.main()
