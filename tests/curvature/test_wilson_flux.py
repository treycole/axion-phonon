"""The interpolated curvature against the DFT overlaps it came from (port of
analysis/validation/external_curvature/check_vs_wilson_flux.ipynb, which used the .mmn of SrTiO3; production keeps
the overlaps in the .chk). On the Y2Ir2O7 production base (its band-traced curvature is not zero; MnBi2Te4's is, by
PT, so there the test would be empty):

for a sample of plaquettes of the Wannier mesh in each reduced plane (mu, nu), the integral of Tr Omega_{mu nu}
(axion_phonon.curvature, all three parts) over the plaquette, by Gauss-Legendre quadrature, against the phase of the
occupied Wilson loop around it, built from the .chk overlaps M(k, b) in the Wannier gauge and the occupied
eigenvectors of V^dag E V at the four corners: arg det W = -integral Tr Omega, up to the finite-link (1/8 of the zone)
error. The quadrature must be converged (order 2 vs 3), and the residual small against the flux.
"""
import numpy as np
import pytest

from axion_phonon import curvature

# ---- settings
N_PLAQUETTES = 48             # per plane, random (seed 5)
ORDERS = (2, 3)               # Gauss-Legendre points per side
BATCH = 64
QUADRATURE = 5e-3             # max |I(order 3) - I(order 2)| / max |I|; well below RESIDUAL_RMS
RESIDUAL_RMS = 0.05           # rms |wrap(-I - phase)| / rms |phase|, per plane (10-31: 0.025-0.027)
RESIDUAL_MAX = 0.15           # max |wrap(-I - phase)| / max |phase|, per plane (10-31: 0.026; SrTiO3 .mmn: 0.02-0.10)
PLANES = ((0, 1), (1, 2), (2, 0))
# --------------
pytestmark = pytest.mark.data


def occupied_vectors(run, n_occ):
    HW = np.einsum("kbi,kb,kbj->kij", run["V"].conj(), run["energies"], run["V"], optimize=True)
    return np.linalg.eigh(HW)[1][:, :, :n_occ]


def wilson_phases(run, n_occ, starts):
    """arg det of the occupied Wilson loop around the plaquette (k, k + e_mu, k + e_mu + e_nu, k + e_nu), per plane."""
    grid = np.array(run["mp_grid"])
    U = occupied_vectors(run, n_occ)
    steps = np.rint((run["kpts"][run["neighbours"]] + run["G"] - run["kpts"][:, None]) * grid).astype(int)
    index = {tuple(np.rint(k * grid).astype(int) % grid): i for i, k in enumerate(run["kpts"])}

    def link(i, step):
        j = int(np.flatnonzero(np.all(steps[i] == step, axis=1))[0])          # the neighbour along this step
        return U[i].conj().T @ run["M"][i, j] @ U[run["neighbours"][i, j]]

    def move(i, step):
        return index[tuple((np.rint(run["kpts"][i] * grid).astype(int) + step) % grid)]

    out = np.empty((len(PLANES), len(starts)))
    for p, (mu, nu) in enumerate(PLANES):
        e_mu, e_nu = np.eye(3, dtype=int)[mu], np.eye(3, dtype=int)[nu]
        for s, i in enumerate(starts):
            i1 = move(i, e_mu)
            i2 = move(i1, e_nu)
            i3 = move(i, e_nu)
            W = link(i, e_mu) @ link(i1, e_nu) @ link(i2, -e_mu) @ link(i3, -e_nu)
            out[p, s] = np.angle(np.linalg.det(W))
    return out


def plaquette_integrals(model, n_occ, kpts, starts, grid, order):
    """integral of Tr Omega_{mu nu} (reduced) over each plaquette, per plane."""
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x + 1) / 2, w / 2
    pos = curvature.position_terms(model, dtype=np.complex128)
    out = np.empty((len(PLANES), len(starts)))
    for p, (mu, nu) in enumerate(PLANES):
        pts, wts, owner = [], [], []
        for s, i in enumerate(starts):
            for a, wa in zip(x, w):
                for b, wb in zip(x, w):
                    k = kpts[i].copy()
                    k[mu] += a / grid[mu]
                    k[nu] += b / grid[nu]
                    pts.append(k)
                    wts.append(wa * wb / (grid[mu] * grid[nu]))
                    owner.append(s)
        pts, wts, owner = np.array(pts), np.array(wts), np.array(owner)
        tr = np.concatenate([np.trace(curvature.omega(curvature.connection(
            curvature.fields(model, pos, pts[b:b + BATCH]), n_occ))[mu, nu], axis1=-2, axis2=-1).real
            for b in range(0, len(pts), BATCH)])
        out[p] = np.bincount(owner, weights=tr * wts, minlength=len(starts))
    return out


def test_wilson_flux(yio_base):
    run, model, n_occ = yio_base["run"], yio_base["model"], yio_base["n_occ"]
    grid = np.array(run["mp_grid"])
    starts = np.random.default_rng(5).choice(len(run["kpts"]), N_PLAQUETTES, replace=False)
    phase = wilson_phases(run, n_occ, starts)
    I = {o: plaquette_integrals(model, n_occ, run["kpts"], starts, grid, o) for o in ORDERS}
    lo, hi = I[ORDERS[0]], I[ORDERS[-1]]
    assert np.abs(hi - lo).max() < QUADRATURE * np.abs(hi).max(), "quadrature not converged"
    residual = np.angle(np.exp(1j * (-hi - phase)))
    for p, plane in enumerate(PLANES):
        rms = np.sqrt((residual[p] ** 2).mean()) / np.sqrt((phase[p] ** 2).mean())
        mx = np.abs(residual[p]).max() / np.abs(phase[p]).max()
        print(f"plane {plane}: flux rms {np.sqrt((phase[p] ** 2).mean()):.3e} rad, residual rms {rms:.3f}, max {mx:.3f}")
        assert rms < RESIDUAL_RMS and mx < RESIDUAL_MAX, f"plane {plane}: residual rms {rms:.3f}, max {mx:.3f}"
