r"""The overlap of two structures' wavefunctions with the PAW metric, and the displaced structure's .amn (bug #23).

    A(k) = O(k) V_0(k),    O_mn(k) = <psi_Q,mk | psi_0,nk>

V_0 = U_dis U of the base Wannier run (its .chk), so the trial functions of the displaced structure are the base
Wannier functions themselves. Three parts, by where they come from (``notes/theory/paw_overlap.md``):
- the plane-wave dot product on the shared G vectors (``match_g``, ``pw_overlap``) and A = O V_0 (``mode_amn``):
  the same as QEWavefunctions.jl's ``compute_wfc_overlap`` and ``generate_amn.jl --no-displacement``;
- the augmentation term in S of one structure: ``qe_port.py``;
- new here (in neither QE nor QEWavefunctions.jl): between two structures each side's projectors sit on its own
  atoms, and the one-centre factor q_ij becomes q_ij + d . D_ij with d = tau_Q - tau_0 and
  D_ij = <phi_i|grad|phi_j> - <phi~_i|grad|phi~_j> (``gradient_matrix``, ``displacement_term``). This makes
  ``overlap`` exact to first order in d. At d = 0 it is QE's S.
Units: bohr, bohr^-1.
"""
from __future__ import annotations

from functools import cache

import numpy as np

from axion_phonon.io.qe import Pseudo
from axion_phonon.paw.qe_port import Channels, PAWSystem, simpson, ylm_real


# ---------------------------------------------------------------- new: partial-wave gradients
def angular_gradient_tables(lmax: int) -> tuple[np.ndarray, np.ndarray]:
    """``G[a, i, j] = int Y_i rhat_a Y_j dOmega`` and ``H[a, i, j] = int Y_i (r d_a Y_j) dOmega`` (QE real Y_lm).

    With these, <R_i Y_i| d_a |R_j Y_j> = G_a int r^2 R_i R_j' dr + H_a int r R_i R_j dr.
    Exact quadrature (Gauss-Legendre in cos(theta) x uniform phi); r d_a Y by central differences of the
    degree-0 extension of Y.
    """
    x, w = np.polynomial.legendre.leggauss(2 * lmax + 6)
    phi = np.linspace(0, 2 * np.pi, 4 * lmax + 8, endpoint=False)
    ct, ph = np.meshgrid(x, phi, indexing="ij")
    st = np.sqrt(1 - ct ** 2)
    rhat = np.stack([st * np.cos(ph), st * np.sin(ph), ct], axis=-1).reshape(-1, 3)
    wt = (w[:, None] * np.full(len(phi), 2 * np.pi / len(phi))[None, :]).ravel()
    y = ylm_real(rhat, lmax)
    h = 1e-5
    G = np.einsum("p,pi,pa,pj->aij", wt, y, rhat, y)
    H = np.empty_like(G)
    for a in range(3):
        e = np.zeros(3)
        e[a] = h
        dy = (ylm_real(rhat + e, lmax) - ylm_real(rhat - e, lmax)) / (2 * h)
        H[a] = np.einsum("p,pi,pj->ij", wt, y, dy)
    return G, H


def _radial_derivative(u: np.ndarray, rab: np.ndarray) -> np.ndarray:
    """du/dr on a radial mesh (du/di / rab), 4th-order central differences, 2nd order at the two ends."""
    d = np.gradient(u)
    d[2:-2] = (-u[4:] + 8 * u[3:-1] - 8 * u[1:-3] + u[:-4]) / 12.0
    return d / rab


def gradient_matrix(ps: Pseudo, ch: Channels, G: np.ndarray, H: np.ndarray) -> np.ndarray:
    """``D_a(ih, jh) = <phi_i|d_a|phi_j> - <phi~_i|d_a|phi~_j>`` (bohr^-1), one-centre, (3, nh, nh) real.

    The O(d) term of the cross-structure augmentation: <phi_i(r - d)|phi_j(r)> = q_ij + d . D_ij + O(d^2).
    Large components only (PP_AEWFC); the small-component (PP_AEWFC_REL) part of q_ij is not differentiated.
    """
    if not ps.augmented:                                                # norm-conserving: no partial waves
        return np.zeros((3, len(ch.indv), len(ch.indv)))
    k = ps.kkbeta
    r, rab = ps.r[:k], ps.rab[:k]
    nb = len(ps.lll)
    i1 = np.zeros((nb, nb))                                             # int u_i (u_j' - u_j / r) dr, AE - PS
    i2 = np.zeros((nb, nb))                                             # int u_i u_j / r dr, AE - PS
    for a in range(nb):
        for b in range(nb):
            if abs(ps.lll[a] - ps.lll[b]) != 1:
                continue
            f1 = f2 = 0.0
            for u, sign in ((ps.aewfc, 1.0), (ps.pswfc, -1.0)):
                ua, ub = u[a][:k], u[b][:k]
                dub = _radial_derivative(ub, rab)
                f1 = f1 + sign * ua * (dub - ub / r)
                f2 = f2 + sign * ua * ub / r
            i1[a, b], i2[a, b] = simpson(f1, rab), simpson(f2, rab)
    lm = ch.nhtolm
    return (G[:, lm][:, :, lm] * i1[ch.indv][:, ch.indv][None] + H[:, lm][:, :, lm] * i2[ch.indv][:, ch.indv][None])


@cache
def gradient_so(system: PAWSystem) -> dict[str, np.ndarray]:
    """D_ij of every species in QE's spinor channel basis (``Channels.to_so``): {species: (3, 4, nh, nh)}; built once."""
    G, H = angular_gradient_tables(system.lmax)
    return {name: np.stack([system.chan[name].to_so(d) for d in gradient_matrix(ps, system.chan[name], G, H)])
            for name, ps in system.pseudo.items()}


def displacement_term(system: PAWSystem, becp1: np.ndarray, becp2: np.ndarray, disp: np.ndarray) -> np.ndarray:
    """sum_a sum_ij,ss' conj(becp1[m,s,i]) (d_a . D_so^a)[ss',ij] becp2[n,s',j]  ->  (nbnd1, nbnd2).

    ``disp`` (nat, 3; bohr) is each atom's displacement of side 1 relative to side 2: the O(d) change of the
    one-centre overlap of two partial waves whose spheres sit d apart. Atoms with d = 0 contribute nothing.
    """
    dso = gradient_so(system)
    out = np.zeros((becp1.shape[0], becp2.shape[0]), complex)
    for a, sp, off in system.blocks:
        if not np.any(disp[a]):
            continue
        dq = np.tensordot(disp[a], dso[sp], axes=(0, 0))
        sl = slice(off, off + dq.shape[1])
        b1, b2 = becp1[:, :, sl], becp2[:, :, sl]
        for ijs, (s1, s2) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            out += b1[:, s1].conj() @ dq[ijs] @ b2[:, s2].T
    return out


# ---------------------------------------------------------------- the overlap of two structures
def match_g(mill1: np.ndarray, mill2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Indices (i1, i2) with mill1[i1] == mill2[i2]."""
    key = lambda m: (m[:, 0].astype(np.int64) * 4096 + m[:, 1]) * 4096 + m[:, 2]
    k1, k2 = key(mill1 + 2048), key(mill2 + 2048)
    common, i1, i2 = np.intersect1d(k1, k2, assume_unique=True, return_indices=True)
    return i1, i2


def pw_overlap(w1: dict, w2: dict) -> np.ndarray:
    """Plane-wave part <psi1_m|psi2_n> on the shared G vectors (both spin components)."""
    i1, i2 = match_g(w1["mill"], w2["mill"])
    g1, g2 = w1["igwx"], w2["igwx"]
    a = np.concatenate([w1["evc"][:, i1], w1["evc"][:, g1 + i1]], axis=1)
    b = np.concatenate([w2["evc"][:, i2], w2["evc"][:, g2 + i2]], axis=1)
    return a.conj() @ b.T


def overlap(sys1: PAWSystem, w1: dict, sys2: PAWSystem, w2: dict) -> np.ndarray:
    """O_mn = <psi1_m|psi2_n> with the PAW augmentation, exact to first order in d = tau1 - tau2; same k in both.

    ``w1``, ``w2`` from ``io.qe.read_wfc``. The projectors of each side sit on its own atoms; species, channels and
    tables are the same for both structures, so ``sys1``'s serve both. With sys1 = sys2 it is QE's <psi1|S|psi2>.
    """
    if not np.allclose(w1["xk"], w2["xk"], atol=1e-8):
        raise ValueError("k points differ")
    tau1, tau2 = sys1.struct["tau"], sys2.struct["tau"]
    p1 = sys1.vkb(w1["xk"], w1["mill"], w1["b"], tau1)
    p2 = sys1.vkb(w2["xk"], w2["mill"], w2["b"], tau2)
    b1, b2 = sys1.becp(p1, w1["evc"]), sys1.becp(p2, w2["evc"])
    return pw_overlap(w1, w2) + sys1.augmentation(b1, b2) + displacement_term(sys1, b1, b2, tau1 - tau2)


def mode_amn(O: np.ndarray, V_base: np.ndarray) -> np.ndarray:
    """A(k) = O(k) V_0(k): (nk, num_bands, num_wann), the displaced structure's .amn."""
    return O @ V_base
