r"""The PAW metric S of one structure, ported from Quantum ESPRESSO (nothing here is new; ``notes/theory/paw_overlap.md``).

With PAW, QE's stored pseudo-wavefunctions are orthonormal only in the metric

    S = 1 + sum_a sum_ij |p_i^a> q_ij^a <p_j^a|        (spin-orbit: q_ij -> qq_so[ij, s s'])

p_i^a are the projectors (QE's beta functions, ``vkb``) of atom a, one per channel (beta, l, m; j with spin-orbit),
and q_ij = <phi_i|phi_j> - <phi~_i|phi~_j> is the norm the pseudo partial waves miss inside the sphere.
Every function reproduces a QE routine, named in its docstring:
- projectors as ``init_us_2``: (-i)^l Y_lm(k+G) beta_l(|k+G|) e^{-i(k+G).tau} / sqrt(Omega), beta_l(q) from QE's table
  (dq = 0.01 bohr^-1, 4 pi int r beta(r) j_l(qr) dr by QE's ``simpson`` up to kkbeta, cubic Lagrange interpolation),
  Y_lm from ``ylmr2``, ordering species -> atoms of that species -> (beta, m);
- q_ij as ``compute_qqr`` at q = 0 (the L = 0 augmentation integral), rotated to qq_so with ``fcoef``
  (``init_us_1``, ``transform_qq_so``) for spin-orbit species;
- S applied as ``s_psi_nc`` (``PAWSystem.augmentation``).
Check: QE's own wavefunctions come out S-orthonormal to ~1e-14 (``analysis/pipeline/overlap.py`` prints it).
The overlap between two structures, which QE never computes, is ``cross_structure.py``.
Units: bohr, bohr^-1. Inputs come from ``axion_phonon/io/qe.py`` (wavefunctions, structure, UPFs).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import spherical_jn

from axion_phonon.io.qe import Pseudo

DQ = 0.01                                            # QE beta_mod table step (bohr^-1)
LMAXX = 3                                            # QE upf_spinorb lmaxx

# ---------------------------------------------------------------- QE numerics
def simpson(func: np.ndarray, rab: np.ndarray) -> float:
    """QE ``simpson`` (upflib/simpsn.f90), including its even-mesh end formula."""
    mesh = len(func)
    i = np.arange(2, mesh)                                              # Fortran i = 2 .. mesh-1
    s = np.sum(np.where(i % 2 == 0, 4.0, 2.0) * func[1:mesh - 1] * rab[1:mesh - 1])
    if mesh % 2 == 1:
        return (s + func[0] * rab[0] + func[-1] * rab[-1]) / 3.0
    return (s + func[0] * rab[0] - func[-3] * rab[-3] * 0.25 + func[-2] * rab[-2] + func[-1] * rab[-1] * 1.25) / 3.0


def simpson_weights(rab: np.ndarray) -> np.ndarray:
    """Weights w with sum(w * f) == simpson(f, rab)."""
    mesh = len(rab)
    w = np.zeros(mesh)
    w[1:mesh - 1] = np.where(np.arange(2, mesh) % 2 == 0, 4.0, 2.0)
    w[0] = 1.0
    if mesh % 2 == 1:
        w[-1] = 1.0
    else:
        w[-3] -= 0.25
        w[-2] += 1.0
        w[-1] = 1.25
    return w * rab / 3.0


def ylm_real(g: np.ndarray, lmax: int) -> np.ndarray:
    """QE ``ylmr2``: real spherical harmonics, columns lm = l^2 + 1 .. (l+1)^2 in QE order; g is (n, 3)."""
    n = len(g)
    gg = np.sum(g ** 2, axis=1)
    gmod = np.sqrt(gg)
    cost = np.where(gmod < 1e-9, 0.0, g[:, 2] / np.where(gmod < 1e-9, 1.0, gmod))
    sent = np.sqrt(np.maximum(0.0, 1.0 - cost ** 2))
    gx = np.where(np.abs(g[:, 0]) > 1e-9, g[:, 0], 1.0)
    phi = np.where(g[:, 0] > 1e-9, np.arctan(g[:, 1] / gx),
                   np.where(g[:, 0] < -1e-9, np.arctan(g[:, 1] / gx) + np.pi, np.where(g[:, 1] < 0, -np.pi / 2, np.pi / 2)))
    q = np.zeros((n, (lmax + 1) ** 2))                                  # 0-based column = Fortran lm - 1
    q[:, 0] = 1.0
    if lmax == 0:
        return q / np.sqrt(4 * np.pi)
    q[:, 1] = cost
    q[:, 3] = -sent / np.sqrt(2.0)
    for l in range(2, lmax + 1):
        for m in range(0, l - 1):
            lm, lm1, lm2 = l ** 2 + 2 * m, (l - 1) ** 2 + 2 * m, (l - 2) ** 2 + 2 * m
            q[:, lm] = cost * (2 * l - 1) / np.sqrt(l * l - m * m) * q[:, lm1] \
                - np.sqrt((l - 1) ** 2 - m * m) / np.sqrt(l * l - m * m) * q[:, lm2]
        lm, lm1, lm2 = l ** 2 + 2 * l, l ** 2 + 2 * (l - 1), (l - 1) ** 2 + 2 * (l - 1)
        q[:, lm1] = cost * np.sqrt(2 * l - 1) * q[:, lm2]
        q[:, lm] = -np.sqrt(2 * l - 1) / np.sqrt(2 * l) * sent * q[:, lm2]
    y = np.empty_like(q)
    y[:, 0] = q[:, 0] / np.sqrt(4 * np.pi)
    lm = 0
    for l in range(1, lmax + 1):
        c = np.sqrt((2 * l + 1) / (4 * np.pi))
        lm += 1
        y[:, lm] = c * q[:, lm]
        for m in range(1, l + 1):
            lm += 2
            y[:, lm - 1] = c * np.sqrt(2.0) * q[:, lm] * np.cos(m * phi)
            y[:, lm] = c * np.sqrt(2.0) * q[:, lm] * np.sin(m * phi)
    return y


def _spinor(l, j, m, spin):
    """QE ``spinor``: Clebsch-Gordan coefficients of the spin-angle functions."""
    denom = 1.0 / (2 * l + 1)
    if abs(j - l - 0.5) < 1e-8:
        return np.sqrt((l + m + 1.0) * denom) if spin == 1 else np.sqrt((l - m) * denom)
    if m < -l + 1:
        return 0.0
    return np.sqrt((l - m + 1.0) * denom) if spin == 1 else -np.sqrt((l + m) * denom)


def _sph_ind(l, j, m, spin):
    """QE ``sph_ind``."""
    if abs(j - l - 0.5) < 1e-8:
        s = m if spin == 1 else m + 1
    else:
        s = 0 if m < -l + 1 else (m - 1 if spin == 1 else m)
    return 0 if (s < -l or s > l) else s


def _rot_ylm():
    """QE ``init_us_1``: rot_ylm(m', mi) for l = lmaxx (rows m' = -l..l, Fortran index m' + l + 1)."""
    l = LMAXX
    rot = np.zeros((2 * l + 1, 2 * l + 1), complex)
    rot[l, 0] = 1.0
    for n1 in range(2, 2 * l + 2, 2):
        m = n1 // 2
        rot[l - m, n1 - 1] = (-1) ** m / np.sqrt(2)
        rot[l - m, n1] = -1j * (-1) ** m / np.sqrt(2)
        rot[l + m, n1 - 1] = 1 / np.sqrt(2)
        rot[l + m, n1] = 1j / np.sqrt(2)
    return rot


@dataclass
class Channels:
    """Per-species projector channels, as QE's nh, indv, nhtol, nhtolm, nhtoj, qq_nt, qq_so."""
    indv: np.ndarray                                # beta index of each ih
    nhtol: np.ndarray
    nhtolm: np.ndarray                              # 0-based column of ylm_real
    nhtoj: np.ndarray | None
    qq_so: np.ndarray                               # (4, nh, nh): ijs = (1,1), (1,2), (2,1), (2,2)
    fcoef: np.ndarray | None                        # (nh, nh, 2, 2), spin-orbit species only

    def to_so(self, op: np.ndarray) -> np.ndarray:
        """A spin-independent one-centre matrix (nh, nh) in QE's spinor channel basis, as transform_qq_so."""
        out = np.zeros((4,) + op.shape, complex)
        if self.fcoef is None:
            out[0] = out[3] = op
            return out
        for ijs, (s1, s2) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            out[ijs] = sum(self.fcoef[:, :, s1, s] @ op @ self.fcoef[:, :, s, s2] for s in (0, 1))
        return out


def channels(ps: Pseudo, lspinorb: bool = True) -> Channels:
    """QE ``init_us_1``: the channels (nh, indv, nhtol, nhtolm, nhtoj), qq_nt and, with spin-orbit, fcoef and qq_so."""
    indv, nhtol, nhtolm, nhtoj = [], [], [], []
    for nb, l in enumerate(ps.lll):
        for m in range(2 * l + 1):
            indv.append(nb)
            nhtol.append(l)
            nhtolm.append(l * l + m)
            if ps.jjj is not None:
                nhtoj.append(ps.jjj[nb])
    indv, nhtol, nhtolm = np.array(indv), np.array(nhtol), np.array(nhtolm)
    nh = len(indv)
    # qq_nt (compute_qqr at q = 0): L = 0 augmentation integral, QE simpson up to kkbeta, same (l, m) only
    q0 = np.array([[simpson(ps.qfunc0[a, b, :ps.kkbeta], ps.rab[:ps.kkbeta]) for b in range(len(ps.lll))]
                   for a in range(len(ps.lll))])
    qq = np.where(nhtolm[:, None] == nhtolm[None, :], q0[indv][:, indv], 0.0)
    fcoef = None
    if lspinorb and ps.jjj is not None:
        nhtoj = np.array(nhtoj)
        rot = _rot_ylm()
        fcoef = np.zeros((nh, nh, 2, 2), complex)
        for ih in range(nh):
            li, ji, mi = nhtol[ih], nhtoj[ih], nhtolm[ih] - nhtol[ih] ** 2 + 1          # mi Fortran 1-based
            for kh in range(nh):
                lk, jk, mk = nhtol[kh], nhtoj[kh], nhtolm[kh] - nhtol[kh] ** 2 + 1
                if li != lk or abs(ji - jk) > 1e-7:
                    continue
                for s1 in (1, 2):
                    for s2 in (1, 2):
                        c = 0j
                        for m in range(-li - 1, li + 1):
                            m0 = _sph_ind(li, ji, m, s1) + LMAXX
                            m1 = _sph_ind(lk, jk, m, s2) + LMAXX
                            c += rot[m0, mi - 1] * _spinor(li, ji, m, s1) * np.conj(rot[m1, mk - 1]) \
                                * _spinor(lk, jk, m, s2)
                        fcoef[ih, kh, s1 - 1, s2 - 1] = c
        same_beta = indv[:, None] == indv[None, :]
        fcoef = fcoef * same_beta[:, :, None, None]                     # init_us_1: zero for vi /= vj
    ch = Channels(indv=indv, nhtol=nhtol, nhtolm=nhtolm, nhtoj=np.array(nhtoj) if len(nhtoj) else None,
                  qq_so=np.zeros((4, nh, nh), complex), fcoef=fcoef)
    ch.qq_so = ch.to_so(qq)          # transform_qq_so: qq_so(kh,lh,ijs) = sum qq(ih,jh) fcoef(kh,ih,s1,s) fcoef(jh,lh,s,s2)
    return ch


def beta_table(ps: Pseudo, qmax: float) -> np.ndarray:
    """QE ``init_tab_beta``: 4 pi int r beta(r) j_l(q r) dr on q = 0, dq, ... (rows: projectors)."""
    nq = int(qmax / DQ + 4) + 1
    q = np.arange(nq) * DQ
    k = ps.kkbeta
    r, rab = ps.r[:k], ps.rab[:k]
    w = simpson_weights(rab)
    tab = np.empty((len(ps.lll), nq))
    for l in set(ps.lll.tolist()):
        jl = spherical_jn(l, np.outer(q, r))                           # (nq, kkbeta)
        for nb in np.flatnonzero(ps.lll == l):
            tab[nb] = 4 * np.pi * jl @ (ps.beta[nb][:k] * r * w)
    return tab


def interp_beta(tab: np.ndarray, qg: np.ndarray, omega: float) -> np.ndarray:
    """QE ``interp_beta`` (cubic Lagrange on the table), divided by sqrt(Omega); returns (nbeta, len(qg))."""
    px = qg / DQ - np.floor(qg / DQ)
    ux, vx, wx = 1.0 - px, 2.0 - px, 3.0 - px
    i0 = np.floor(qg / DQ).astype(int)
    if i0.max() + 3 >= tab.shape[1]:
        raise ValueError("beta table too short")
    v = (tab[:, i0] * ux * vx * wx / 6.0 + tab[:, i0 + 1] * px * vx * wx / 2.0
         - tab[:, i0 + 2] * px * ux * wx / 2.0 + tab[:, i0 + 3] * px * ux * vx / 6.0)
    return v / np.sqrt(omega)


# ---------------------------------------------------------------- projectors and S
class PAWSystem:
    """Projectors of one structure: species channels and beta tables (built once), vkb per k point."""

    def __init__(self, structure: dict, pseudos: dict[str, Pseudo], qmax: float = 30.0, lspinorb: bool = True):
        """``structure`` from ``io.qe.read_structure``, ``pseudos`` from ``io.qe.read_pseudos``."""
        self.struct = structure
        self.species = [name for name, _ in structure["species"]]
        self.pseudo = pseudos
        self.chan = {name: channels(ps, lspinorb) for name, ps in self.pseudo.items()}
        self.tab = {name: beta_table(ps, qmax) for name, ps in self.pseudo.items()}
        self.lmax = max(int(ps.lll.max()) for ps in self.pseudo.values())
        # QE order of the projector blocks: species, then atoms of that species
        self.blocks = []                                                # (atom index, species, offset)
        off = 0
        for sp in self.species:
            for a, name in enumerate(self.struct["names"]):
                if name == sp:
                    self.blocks.append((a, sp, off))
                    off += len(self.chan[sp].indv)
        self.nkb = off

    def vkb(self, xk: np.ndarray, mill: np.ndarray, b: np.ndarray, tau: np.ndarray | None = None) -> np.ndarray:
        """Projectors (nkb, nG) at k = xk (bohr^-1), G = mill @ b; ``tau`` (bohr) overrides this structure's."""
        tau = self.struct["tau"] if tau is None else tau
        kpg = xk[None, :] + mill @ b
        qg = np.linalg.norm(kpg, axis=1)
        ylm = ylm_real(kpg, self.lmax)
        out = np.empty((self.nkb, len(kpg)), complex)
        cache = {}
        for a, sp, off in self.blocks:
            ch = self.chan[sp]
            if sp not in cache:
                vq = interp_beta(self.tab[sp], qg, self.struct["omega"])
                cache[sp] = ylm[:, ch.nhtolm].T * vq[ch.indv] * ((-1j) ** ch.nhtol)[:, None]
            out[off:off + len(ch.indv)] = cache[sp] * np.exp(-1j * kpg @ tau[a])[None, :]
        return out

    def becp(self, vkb: np.ndarray, evc: np.ndarray) -> np.ndarray:
        """<beta|psi> per spin: (nbnd, 2, nkb) for spinor evc (nbnd, 2 nG)."""
        ng = vkb.shape[1]
        return np.stack([evc[:, :ng] @ vkb.conj().T, evc[:, ng:] @ vkb.conj().T], axis=1)


    def augmentation(self, becp1: np.ndarray, becp2: np.ndarray) -> np.ndarray:
        """QE ``s_psi_nc``'s one-centre term: sum_a sum_ij,ss' conj(becp1[m,s,i]) qq_so^a[ss',ij] becp2[n,s',j].

        Returns (nbnd1, nbnd2). With both becp from one structure, <psi1|S|psi2> = plane-wave overlap + this.
        """
        out = np.zeros((becp1.shape[0], becp2.shape[0]), complex)
        for a, sp, off in self.blocks:
            qq = self.chan[sp].qq_so
            sl = slice(off, off + qq.shape[1])
            b1, b2 = becp1[:, :, sl], becp2[:, :, sl]
            for ijs, (s1, s2) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
                out += b1[:, s1].conj() @ qq[ijs] @ b2[:, s2].T
        return out
