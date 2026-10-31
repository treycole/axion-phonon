r"""The Chern-Simons axion angle of a Wannier model, as an independent check of dtheta/dQ (log 2026-10-14).

    theta = -1/(4 pi) int_{[0,1]^3} d^3k eps^{mnr} Tr[ 1/2 A~_m Omega~_nr + i/3 A~_m A~_n A~_r ]      (reduced k)

theta depends on the gauge only by multiples of 2 pi, but the k sum converges only in a smooth gauge. Here: the
projection gauge |psi~(k)> = P(k) g(k) O^{-1/2}, O = g^+ P g, with k-independent trials G in the home-cell orbital
basis, whose convention-I coefficients are g(k) = exp(-2 pi i k.tau) G. In that gauge

    A~ = i O^{-1/2} X O^{-1/2} + i O^{1/2} d(O^{-1/2})  +  W^+ A_vv W,        Omega~ = W^+ Omega W,
    X_m = S^+ (<v|d_m H|c>/(E_v - E_c) S_c + U_v^+ d_m g),   S = U_v^+ g,  S_c = U_c^+ g,  W = S O^{-1/2},

with A_vv (the Wannier position matrix in the band basis) and the full occupied Omega (internal + cross + external)
from ``curvature`` -- the same quantities whose Tr[Omega Omega] gives c2. Everything is analytic in k (no finite
differences), so the integrand is smooth wherever det O != 0.

Use: dtheta/dQ = [theta(structure) - theta(base)] / Q with ONE trial set G (from the base) and one tau for both.
Unlike the c2 route this needs no frozen gauge (A_beta = 0). It works for theta = 0 materials (YIO: agrees with c2
to <= 0.25%); for theta = pi the trials must break inversion and the k sum needs nk >~ 16 (log 2026-10-14 §1).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import qr

from axion_phonon.curvature import PositionTerms, connection, fields, omega

EPS3 = np.zeros((3, 3, 3))
for _a, _b, _c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    EPS3[_a, _b, _c], EPS3[_a, _c, _b] = 1, -1


def mesh(nk: int, shift: float = 0.0) -> np.ndarray:
    """The nk^3 grid in reduced coordinates, shifted by ``shift`` in each direction."""
    return (np.stack(np.meshgrid(*[np.arange(nk) / nk] * 3, indexing="ij"), -1).reshape(-1, 3) + shift) % 1.0


def occupied_home_cell(model, k: np.ndarray, n_occ: int) -> np.ndarray:
    """Occupied eigenvectors in the home-cell orbital basis (convention-I coefficients times exp(+2 pi i k.tau)):
    (n_k, n_orb, n_occ). P(k) = U U^+ there."""
    _, U = np.linalg.eigh(model.hamiltonian(k, flatten_spin_axis=True))
    tau = np.asarray(model.orb_vecs, float)
    return np.exp(2j * np.pi * k @ tau.T)[:, :, None] * U[:, :, :n_occ]


def scdm_trials(model, nk: int, n_occ: int) -> np.ndarray:
    """Trials G (n_orb, n_occ): the SCDM columns of the occupied projector averaged over the nk^3 mesh, orthonormalized."""
    U = occupied_home_cell(model, mesh(nk), n_occ)
    P = np.einsum("kin,kjn->ij", U, U.conj()) / len(U)
    _, _, piv = qr(P, pivoting=True, mode="economic")
    C = P[:, piv[:n_occ]]
    w, V = np.linalg.eigh(C.conj().T @ C)
    return C @ V @ np.diag(w ** -0.5) @ V.conj().T


def smallest_sigma(model, G: np.ndarray, k: np.ndarray, n_occ: int) -> float:
    """min over k of the smallest singular value of <u_occ|g>: the gauge is smooth only if this stays away from 0."""
    U = occupied_home_cell(model, k, n_occ)
    return np.linalg.svd(U.conj().transpose(0, 2, 1) @ G, compute_uv=False)[:, -1].min()


def _inv_sqrt_and_derivative(O: np.ndarray, dO: np.ndarray):
    """O^{-1/2} and its derivative along each dO[m] (Daleckii-Krein, in the eigenbasis of O)."""
    lam, V = np.linalg.eigh(O)
    f = lam ** -0.5
    gap = lam[:, None] - lam[None, :]
    close = np.abs(gap) < 1e-10 * lam.max()
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.where(close, -0.5 * lam[:, None] ** -1.5, (f[:, None] - f[None, :]) / np.where(close, 1.0, gap))
    Vh = V.conj().T
    return (V * f) @ Vh, np.array([V @ (F * (Vh @ dO[m] @ V)) @ Vh for m in range(len(dO))])


def projection_gauge(model, pos: PositionTerms, G: np.ndarray, k: np.ndarray, n_occ: int):
    """A~ (3, n_k, n_occ, n_occ), Omega~ (3, 3, n_k, n_occ, n_occ) in the projection gauge of G, and min sigma."""
    f = fields(model, pos, k)
    c = connection(f, n_occ)
    Om = omega(c)                                            # (3, 3, n_k, v, v), band gauge
    _, U = np.linalg.eigh(np.asarray(f.H, dtype=np.result_type(f.H.dtype, np.complex64)))
    tau = pos.tau
    A_out = np.empty((3, len(k), n_occ, n_occ), complex)
    Om_out = np.empty((3, 3, len(k), n_occ, n_occ), complex)
    smin = np.inf
    for i in range(len(k)):
        Uv, Uc = U[i, :, :n_occ], U[i, :, n_occ:]
        g = np.exp(-2j * np.pi * tau @ k[i])[:, None] * G
        dg = np.array([-2j * np.pi * tau[:, m][:, None] * g for m in range(3)])
        S, Sc = Uv.conj().T @ g, Uc.conj().T @ g
        smin = min(smin, np.linalg.svd(S, compute_uv=False)[-1])
        Vvc = c.d[:, i].conj().transpose(0, 2, 1)            # <v|d_m H|c> / (E_v - E_c)
        X = np.array([S.conj().T @ (Vvc[m] @ Sc + Uv.conj().T @ dg[m]) for m in range(3)])
        Oi, dOi = _inv_sqrt_and_derivative(S.conj().T @ S, X + X.conj().transpose(0, 2, 1))
        Oh = np.linalg.inv(Oi)
        W = S @ Oi
        Wh = W.conj().T
        for m in range(3):
            A_out[m, i] = 1j * (Oi @ X[m] @ Oi + Oh @ dOi[m]) + Wh @ c.A_vv[m, i] @ W
            for n in range(3):
                Om_out[m, n, i] = Wh @ Om[m, n, i] @ W
    return A_out, Om_out, smin


def theta_sum(model, pos: PositionTerms, G: np.ndarray, k: np.ndarray, n_occ: int, batch: int = 32):
    """sum over the points k of eps Tr[A~ Omega~ / 2 + i A~ A~ A~ / 3], and min sigma. The projection gauge is local
    in k, so a k list can be split into slices and the sums added."""
    total, smin = 0.0, np.inf
    for i in range(0, len(k), batch):
        A, Om, s = projection_gauge(model, pos, G, k[i:i + batch], n_occ)
        smin = min(smin, s)
        AOm = np.einsum("mkab,nrkba->mnr", A, Om, optimize=True)
        AAA = np.einsum("mkab,nkbc,rkca->mnr", A, A, A, optimize=True)
        total += np.einsum("mnr,mnr->", EPS3, 0.5 * AOm + (1j / 3) * AAA)
    return total, smin


def theta(model, pos: PositionTerms, G: np.ndarray, nk: int, n_occ: int, shift: float = 0.0, batch: int = 32):
    """theta_CS on the nk^3 mesh (complex: the imaginary part measures the error), and min sigma."""
    k = mesh(nk, shift)
    total, smin = theta_sum(model, pos, G, k, n_occ, batch)
    return -(total / len(k)) / (4 * np.pi), smin
