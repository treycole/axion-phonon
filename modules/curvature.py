#!/usr/bin/env python3
r"""Non-Abelian Berry curvature of a Wannier-interpolated insulator.

Where this sits
---------------

::

    wannier/            STAGE 1  Wannier90 writes  seedname_hr.dat, seedname_r.dat
        |
        v               PythTB reads them (W90.model) and gives H(k), dH/dk
    modules/curvature.py   STAGE 2  THIS FILE.  Builds each piece of Omega.
        |
        v
    modules/axion.py       STAGE 3  finite difference in beta, c2, BZ integral

Everything here transcribes ``notes/berry_curvature_derivation.md`` (the
"note"); equation numbers in the docstrings are the note's.  The note splits the
curvature of the M occupied bands, as an M x M matrix, into three pieces

    B = B_int + B_cross + B_ext                                      (24)

that are built from two ingredients in the band (Hamiltonian) gauge:

    d_mu  = internal connection, from H and dH only                  (31)
    a_mu  = external connection, from the position matrix X(R)       (28), 7.2

The seven functions
-------------------

=================  ==========================================  ==============
function           what it returns                             note
=================  ==========================================  ==============
``position_terms``  X(R) = <0s|r - tau_s|Rt>, bond vectors b(R)  (8), (9)
``fields``          Wannier gauge, k space: H, dH/dk, A, curl A  (7), (8), (9)
``connection``      band gauge: d (internal); a, A_vv, Omega_vv  (10), (28), (31)
                    (external)
``omega_internal``  i[d_mu^dag d_nu - d_nu^dag d_mu]             (29)
``omega_cross``     (d_mu^dag a_nu + h.c.) - (mu <-> nu)         (30)
``omega_external``  Omega_vv - i[A_vv_mu, A_vv_nu]               (27)
``omega``           internal + cross + external                  (24)
=================  ==========================================  ==============

Use it yourself
---------------

::

    from pythtb import W90
    from modules.curvature import (position_terms, fields, connection,
                                   omega, omega_internal, omega_cross, omega_external)

    w90   = W90("run_dir", "seed")            # reads seed_hr.dat, seed_r.dat, seed_wsvec.dat
    tau   = w90.lattice.orb_vecs              # nominal orbital positions, FIXED (assumption A3)
    model = w90.model(orb_vecs=tau, min_hopping_norm=1e-5)

    pos   = position_terms(model)             # once per structure
    f     = fields(model, pos, k)             # k: (nk, 3) reduced coordinates
    c     = connection(f, n_occ)              # n_occ = number of occupied bands
    Om    = omega(c)                          # (3, 3, nk, n_occ, n_occ)
    B_xy  = np.trace(Om[0, 1], axis1=-2, axis2=-1)   # Berry curvature, band-summed

Conventions
-----------

* The **direction index comes first**: ``dH`` is ``(3, nk, J, J)``, ``A`` is
  ``(3, nk, J, J)``, ``Omega`` is ``(D, D, nk, Nv, Nv)``.  ``D = 3`` (kx, ky, kz)
  or ``D = 4`` (kx, ky, kz, beta).  J = Wannier functions, Nv = occupied bands,
  Nc = J - Nv.
* Directions are the **reduced** components (derivatives with respect to the
  reduced k), the convention of ``TBModel.velocity``.
* ``tau`` is the *nominal* orbital position of convention I.  It must be the
  same array for every structure that will be differenced (assumption A3); the
  physical motion of the Wannier centres enters through X(R) instead.
* Insulator: bands ``0 .. n_occ-1`` are occupied and gapped from the rest (A1).
* Arrays are complex64 for the large ones (as in production) unless the input
  carries complex128, which is preserved.

If you want a different way of getting X(R), that decision belongs in stage 1
(what Wannier90 writes), not here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "BetaTerms",
    "Connection",
    "Fields",
    "PositionTerms",
    "connection",
    "fields",
    "omega",
    "omega_cross",
    "omega_external",
    "omega_internal",
    "position_terms",
]

# The three (k_i, k_j) planes; the dual index l of a curl component labels
# the plane (l+1, l+2), i.e. curl[l] = Omega_{l+1, l+2}.
_KK_PLANES = ((1, 2), (2, 0), (0, 1))
_ZERO_R = (0, 0, 0)


def _adjoint(matrix: np.ndarray) -> np.ndarray:
    """Conjugate transpose of the last two axes."""
    return matrix.conj().swapaxes(-1, -2)


# --------------------------------------------------------------------------- #
# Containers                                                                  #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class PositionTerms:
    r"""Real-space position data of one structure (output of ``position_terms``).

    ``X[R, a, s, t] = <0s| r_a - tau_s |Rt>`` in Cartesian components (Angstrom),
    ``curl_weight`` is the R-space weight of the plain curl (9), and ``tau`` the
    nominal positions the Bloch phases will use.
    """

    R: np.ndarray            # (nR, 3) integer lattice vectors
    X: np.ndarray            # (nR, 3, J, J)
    curl_weight: np.ndarray  # (nR, 3, J, J)  i (b_mu X_nu - b_nu X_mu), dual index
    tau: np.ndarray          # (J, 3) reduced nominal positions (fixed embedding)
    one_form: np.ndarray     # (3, 3) Cartesian -> reduced, for a vector (A)
    two_form: np.ndarray     # (3, 3) Cartesian -> reduced, for a 2-form (curl A)
    dtype: np.dtype


@dataclass(frozen=True)
class Fields:
    r"""Wannier-gauge, k-space ingredients for one batch of k-points.

    ``H (nk,J,J)``; ``dH (3,nk,J,J)`` = dH/dk; ``A (3,nk,J,J)`` the external
    connection (8); ``curl (3,nk,J,J)`` its plain curl (9), ``curl[l]`` being
    the (l+1, l+2) plane.  None of these is gauge covariant on its own.
    """

    H: np.ndarray
    dH: np.ndarray
    A: np.ndarray
    curl: np.ndarray


@dataclass(frozen=True)
class BetaTerms:
    r"""What the fourth direction beta adds (section 8); built by ``axion.py``.

    ``dH (nk,J,J)`` = dH/dbeta and ``curl (3,nk,J,J)`` = Omega^W_{l,beta}, both
    from the finite difference of two structures.  ``A`` is the parametric
    connection A_beta (33), or ``None`` for the **frozen Wannier gauge**
    (section 8.5), in which every term containing A_beta is absent.
    """

    dH: np.ndarray
    curl: np.ndarray
    A: np.ndarray | None = None


@dataclass(frozen=True)
class Connection:
    r"""Band-gauge connection, internal and external (output of ``connection``).

    ``d``       (D, nk, Nc, Nv)  internal, d_{mu,cv} of (31)
    ``a``       (D, nk, Nc, Nv)  external, A^{ext,cv}_mu of (28)
    ``A_vv``    (D, nk, Nv, Nv)  external, A^{ext,vv}_mu
    ``Omega_vv``(D, D, nk, Nv, Nv)  plain curl of A^ext in the occupied block
    ``E``       (nk, J)          band energies, ascending
    """

    d: np.ndarray
    a: np.ndarray
    A_vv: np.ndarray
    Omega_vv: np.ndarray
    E: np.ndarray
    n_occ: int
    frozen_gauge: bool


# --------------------------------------------------------------------------- #
# 1. A(R): the position matrix, once per structure                            #
# --------------------------------------------------------------------------- #
def position_terms(model, dtype=np.complex64) -> PositionTerms:
    r"""Prepare ``X(R) = <0s| r - tau_s |Rt>`` and the curl weights, eqs (8), (9).

    ``model`` is a PythTB model built with ``W90.model(orb_vecs=tau)`` from a run
    whose ``seedname_r.dat`` is present.  Reads the position blocks PythTB stored
    from that file, then:

    * takes the Hermitian part, ``X(R) -> [X(R) + X(-R)^dag]/2``, and divides by
      the Wigner-Seitz degeneracy (1 for a ``write_ndegen_applied`` run);
    * subtracts the nominal position from the ``R = 0`` diagonal, so the
      remainder is the *physical* offset of each Wannier centre from ``tau``;
    * forms the bond vectors ``b_st(R) = R + tau_t - tau_s`` and the curl weights
      ``i (b_mu X_nu - b_nu X_mu)`` for the three (k_i, k_j) planes.

    This is the only function that touches PythTB's stored position data
    (``_pos_r``, ``_ham_deg``); it does the same reads as PythTB's own external
    curvature.  Call it once per structure and reuse the result across batches.
    """

    pos_r = getattr(model, "_pos_r", None)
    if pos_r is None:
        raise ValueError(
            "The model has no Wannier position matrix.  Stage 1 must write "
            "seedname_r.dat (write_rmn = true), and PythTB must be on its "
            "'external' branch."
        )
    lattice = np.asarray(model.lattice.lat_vecs, dtype=float)
    tau = np.asarray(model.orb_vecs, dtype=float)
    tau_cartesian = tau @ lattice
    num_orbitals = len(tau)

    R_keys = list(pos_r)
    R = np.asarray(R_keys, dtype=int)
    R_cartesian = R @ lattice

    X = np.empty((len(R_keys), 3, num_orbitals, num_orbitals), dtype=complex)
    for index, key in enumerate(R_keys):
        partner = pos_r.get(tuple(-value for value in key), np.zeros_like(pos_r[key]))
        X[index] = 0.5 * (pos_r[key] + np.conj(np.transpose(partner, (0, 2, 1))))
        X[index] /= float(model._ham_deg.get(key, 1.0))

    origin = [index for index, key in enumerate(R_keys) if key == _ZERO_R]
    if len(origin) != 1:
        raise ValueError("The position matrix needs exactly one R = 0 block.")
    diagonal = np.arange(num_orbitals)
    for axis in range(3):
        X[origin[0], axis, diagonal, diagonal] -= tau_cartesian[:, axis]

    # b_st(R) = R + tau_t - tau_s, indexed [R, s, t, cartesian]
    bond = (
        R_cartesian[:, None, None, :]
        + tau_cartesian[None, None, :, :]
        - tau_cartesian[None, :, None, :]
    )
    curl_weight = np.empty_like(X)
    for dual, (mu, nu) in enumerate(_KK_PLANES):
        curl_weight[:, dual] = 1j * (
            bond[..., mu] * X[:, nu] - bond[..., nu] * X[:, mu]
        )

    # Cartesian -> reduced components: a vector transforms with the reciprocal
    # vectors, a 2-form (the curl) with det(B) B^{-T}.
    reciprocal = np.asarray(model.recip_lat_vecs, dtype=float)
    real_dtype = np.float32 if np.dtype(dtype) == np.complex64 else float
    return PositionTerms(
        R=R,
        X=X.astype(dtype),
        curl_weight=curl_weight.astype(dtype),
        tau=tau,
        one_form=reciprocal.astype(real_dtype),
        two_form=(np.linalg.det(reciprocal) * np.linalg.inv(reciprocal).T).astype(
            real_dtype
        ),
        dtype=np.dtype(dtype),
    )


# --------------------------------------------------------------------------- #
# 2. Wannier gauge, k space                                                   #
# --------------------------------------------------------------------------- #
def fields(model, pos: PositionTerms, k: np.ndarray) -> Fields:
    r"""H(k), dH/dk, A(k) and its plain curl in the Wannier gauge, eqs (7)-(9).

    ``H(k)`` and ``dH/dk`` come straight from PythTB (eq 7 applied to
    ``seedname_hr.dat``).  The external connection is the Bloch sum

        A_{mu,st}(k) = sum_R exp(i k . b_st(R)) X_{mu,st}(R)          (8)

    and its plain curl (9) is the same sum with the curl weights.  Both carry
    the convention-I bond phase ``exp(2 pi i k . (tau_t - tau_s))`` with the
    *fixed* ``tau`` of ``pos``.  ``k`` is ``(nk, 3)`` in reduced coordinates.
    """

    k = np.atleast_2d(np.asarray(k, dtype=float))
    if k.shape[1] != 3:
        raise ValueError("k must have shape (nk, 3).")
    dtype = pos.dtype
    nk = len(k)
    num_orbitals = pos.tau.shape[0]

    H = model.hamiltonian(k, flatten_spin_axis=True).astype(dtype, copy=False)
    dH = model.velocity(k, flatten_spin_axis=True).astype(dtype, copy=False)

    cell_phase = np.exp(2j * np.pi * (k @ pos.R.T)).astype(dtype)  # e^{2 pi i k.R}
    orbital_phase = np.exp(2j * np.pi * (k @ pos.tau.T))
    bond_phase = (
        np.conj(orbital_phase)[:, :, None] * orbital_phase[:, None, :]
    ).astype(dtype)  # e^{2 pi i k.(tau_t - tau_s)}

    shape = (nk, 3, num_orbitals, num_orbitals)
    A_cartesian = (cell_phase @ pos.X.reshape(len(pos.R), -1)).reshape(shape)
    A_cartesian *= bond_phase[:, None]
    curl_cartesian = (
        cell_phase @ pos.curl_weight.reshape(len(pos.R), -1)
    ).reshape(shape)
    curl_cartesian *= bond_phase[:, None]

    return Fields(
        H=H,
        dH=dH,
        A=np.einsum("ua,kaij->ukij", pos.one_form, A_cartesian),
        curl=np.einsum("ua,kaij->ukij", pos.two_form, curl_cartesian),
    )


# --------------------------------------------------------------------------- #
# 3. Band gauge: internal and external connection                             #
# --------------------------------------------------------------------------- #
def connection(f: Fields, n_occ: int, beta: BetaTerms | None = None) -> Connection:
    r"""Rotate to the band gauge and build the internal and external connection.

    Diagonalise ``H = U E U^dag`` (10): ``U`` are the occupied columns and
    ``Ubar`` the empty ones, ``J = Nv + Nc``.

    **Internal** (from H and dH only), first-order perturbation theory (31):

        d_{mu,cv} = [Ubar^dag (dH_mu) U]_{cv} / (E_v - E_c)

    **External** (from the position matrix), the rotations of (8) and (9) that
    section 7.2 says need no further derivation:

        a_mu = Ubar^dag A_mu U      (28)   cross-gap block, c rows, v columns
        A_vv = U^dag A_mu U                 occupied block
        Omega_vv = U^dag (curl A) U

    ``beta=None`` gives ``D = 3`` directions (kx, ky, kz).  Passing ``BetaTerms``
    adds the fourth direction: ``d_beta`` from ``beta.dH``, ``Omega_vv[l, beta]``
    from ``beta.curl``, and ``a_beta``, ``A_vv_beta`` from ``beta.A`` -- or zero
    if ``beta.A is None`` (the frozen Wannier gauge, recorded in
    ``frozen_gauge``).  Requires a gap between bands ``n_occ-1`` and ``n_occ``.
    """

    dtype = np.result_type(f.H.dtype, np.complex64)
    energies, U = np.linalg.eigh(np.asarray(f.H, dtype=dtype))  # columns = eigenvectors
    nk = energies.shape[0]
    Uv, Uc = U[..., :n_occ], U[..., n_occ:]
    Uv_h, Uc_h = _adjoint(Uv), _adjoint(Uc)
    n_empty = Uc.shape[-1]
    n_dir = 3 if beta is None else 4

    # ---- internal: d_{mu,cv}, eq (31) ------------------------------------
    dH = f.dH if beta is None else np.concatenate([f.dH, beta.dH[None]], axis=0)
    gap = energies[:, None, :n_occ] - energies[:, n_occ:, None]  # E_v - E_c, (nk, Nc, Nv)
    d = (Uc_h[None] @ dH @ Uv[None]) / gap[None]

    # ---- external: a, A_vv, Omega_vv, eqs (28), (8), (9) -----------------
    a = np.zeros((n_dir, nk, n_empty, n_occ), dtype=d.dtype)
    A_vv = np.zeros((n_dir, nk, n_occ, n_occ), dtype=d.dtype)
    Omega_vv = np.zeros((n_dir, n_dir, nk, n_occ, n_occ), dtype=d.dtype)

    a[:3] = Uc_h[None] @ f.A @ Uv[None]
    A_vv[:3] = Uv_h[None] @ f.A @ Uv[None]
    curl_vv = Uv_h[None] @ f.curl @ Uv[None]
    for dual, (mu, nu) in enumerate(_KK_PLANES):
        Omega_vv[mu, nu] = curl_vv[dual]
        Omega_vv[nu, mu] = -curl_vv[dual]

    frozen = False
    if beta is not None:
        mixed_vv = Uv_h[None] @ beta.curl @ Uv[None]  # plane (l, beta)
        for axis in range(3):
            Omega_vv[axis, 3] = mixed_vv[axis]
            Omega_vv[3, axis] = -mixed_vv[axis]
        if beta.A is None:
            frozen = True  # A_beta = 0: a[3], A_vv[3] stay zero
        else:
            a[3] = Uc_h @ beta.A @ Uv
            A_vv[3] = Uv_h @ beta.A @ Uv

    return Connection(
        d=d, a=a, A_vv=A_vv, Omega_vv=Omega_vv, E=energies,
        n_occ=n_occ, frozen_gauge=frozen,
    )


# --------------------------------------------------------------------------- #
# 4-7. The three pieces of Omega, and their sum                               #
# --------------------------------------------------------------------------- #
def omega_internal(c: Connection) -> np.ndarray:
    r"""``B_int = i [d_mu^dag d_nu - d_nu^dag d_mu]``, eq (29).

    With no external terms this is the textbook non-Abelian Kubo formula (32).
    Returns ``(D, D, nk, Nv, Nv)``; the same shape for the other pieces.
    """
    q = _adjoint(c.d)[:, None] @ c.d[None]  # q[mu, nu] = d_mu^dag d_nu
    return 1j * (q - _adjoint(q))  # adjoint(q)[mu, nu] = q[nu, mu]


def omega_cross(c: Connection) -> np.ndarray:
    r"""``B_cross = (d_mu^dag a_nu + h.c.) - (mu <-> nu)``, eq (30).

    There is **no factor of i**: the i of (26) cancels against the i hidden in
    ``A^{int,cv} = i d``.  Hermitian by construction.
    """
    t = _adjoint(c.d)[:, None] @ c.a[None]  # t[mu, nu] = d_mu^dag a_nu
    P = t + _adjoint(t)
    return P - P.swapaxes(0, 1)


def omega_external(c: Connection) -> np.ndarray:
    r"""``B_ext = Omega_vv - i [A_vv_mu, A_vv_nu]``, eq (27).

    The commutator is what makes it gauge covariant in the occupied block.
    """
    AA = c.A_vv[:, None] @ c.A_vv[None]  # AA[mu, nu] = A_mu A_nu
    return c.Omega_vv - 1j * (AA - AA.swapaxes(0, 1))


def omega(c: Connection) -> np.ndarray:
    r"""The full occupied-space curvature ``B = B_int + B_cross + B_ext``, eq (24).

    Only the sum is independent of how many Wannier functions are used; the
    individual pieces exchange weight as the Wannier space grows (section 5.4).
    """
    return omega_internal(c) + omega_cross(c) + omega_external(c)
