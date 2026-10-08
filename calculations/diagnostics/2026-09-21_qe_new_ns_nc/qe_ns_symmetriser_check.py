#!/usr/bin/env python3
"""Does QE's noncollinear DFT+U occupation-matrix symmetriser (new_ns_nc) preserve a symmetric state?

Quantum ESPRESSO 7.5, PW/src/new_ns.f90::new_ns_nc.  The script re-implements the algorithm with QE's OWN inputs,
read from a pw.x output printed with verbosity = 'high':  the symmetry operations (crystal s, fractional
translation, Cartesian sr, time-reversal flag), irt from checksym's rule, find_u / comp_dspinldau for the spin
matrices, d_matrix for the orbital rotations, and QE's analytic AIAO starting ns (init_ns_nc).

Result on the Y2Ir2O7 all-in-all-out base run (symmetry_with_labels = .true., four Ir species): applied to a
state that is invariant under all 48 operations, new_ns_nc as written does not return it (28 of the 48 terms
are wrong: exactly the operations of order >= 3) and shrinks the site moment by a factor 3.  Using the atom
that isym maps INTO na, nb = irt(invs(isym), na), instead of irt(isym, na), makes it an exact projector.
Section 6 tests the other way to remove the mismatch (keep nb, transpose the contraction) and the version of it in the
open upstream draft MR !2512, which does NOT restore the symmetry on these operations (see the note).

    python qe_ns_symmetriser_check.py PATH/TO/Y2Ir2O7.scf.out

The atomic layout (four Ir species, AIAO angles, fcc lattice) is that of the pyrochlore run and is set below;
adapt AXIS_ANGLES / IR / at for another structure.  See notes/progress/2026-09-21_qe_new_ns_nc_bug.md.
"""
import io, re, sys
import numpy as np
from pathlib import Path

import re
import sys
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else sys.exit(__doc__)
lines = OUT.read_text().splitlines()

# ------------------------------------------------------------------------------------------
# 1. the 48 operations exactly as QE printed them (summary.f90)
# ------------------------------------------------------------------------------------------
row_re = re.compile(r"\(\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*\)(?:\s+(?:f =)?\(\s*(-?[\d.]+)\s*\))?")
ops = []
i = 0
while i < len(lines):
    m = re.match(r"\s+isym =\s+(\d+)\s+(.*?)\s*$", lines[i])
    if m and "Time Reversal" in lines[i + 2]:
        isym, name = int(m.group(1)), m.group(2).strip()
        t_rev = int(lines[i + 2].split()[-1])
        cry, ft = [], []
        j = i + 3
        assert "cryst." in lines[j]
        for k in range(3):
            r = row_re.search(lines[j + k])
            cry.append([int(float(x)) for x in r.groups()[:3]])
            ft.append(float(r.group(4)) if r.group(4) is not None else 0.0)
        j += 4  # 3 rows + blank
        assert "cart." in lines[j], lines[j]
        cart = []
        for k in range(3):
            r = row_re.search(lines[j + k])
            cart.append([float(x) for x in r.groups()[:3]])
        ops.append(dict(isym=isym, name=name, t_rev=t_rev, s=np.array(cry), ft=np.array(ft), sr=np.array(cart)))
        i = j + 3
    else:
        i += 1
assert len(ops) == 48, len(ops)
for o in ops:
    o["sr"] = np.round(o["sr"]).astype(float)  # entries are 0, +-1 (cubic frame)
print(f"parsed {len(ops)} operations, {sum(o['t_rev'] for o in ops)} with time reversal")

# ------------------------------------------------------------------------------------------
# 2. lattice and Ir positions (alat units), atoms 5..8 = Ir1..Ir4 (0-based 4..7)
# ------------------------------------------------------------------------------------------
s2 = 1.0 / np.sqrt(2.0)
at = s2 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)  # columns are a1,a2,a3 (alat units)
tau, label = {}, {}
start = next(k for k, l in enumerate(lines) if 'positions (alat units)' in l)
for ln in lines[start + 1:start + 23]:
    m = re.match(r"\s+(\d+)\s+(\S+)\s+tau\(\s*\d+\) = \(\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*\)", ln)
    if m:
        tau[int(m.group(1))] = np.array(list(map(float, m.groups()[2:])))
        label[int(m.group(1))] = m.group(2)
assert len(tau) == 22
xau = {a: np.linalg.solve(at, t) for a, t in tau.items()}       # crystal coordinates
chem = lambda a: label[a][:2].rstrip("0123456789")               # QE chem_symb(): 'Ir1' -> 'Ir'
IR = [5, 6, 7, 8]


def eqvect(x, y, f, acc=1e-5):
    d = x - y - f
    return np.all(np.abs(d - np.rint(d)) < acc)


# QE checksym(): rau = s^T xau ; irt(isym, na) = nb with rau(na) - xau(nb) - ft in Z^3 and same chemical symbol
def compute_irt(o):
    irt = {}
    for na in tau:
        rau = o["s"].T @ xau[na]
        found = [nb for nb in tau if chem(nb) == chem(na) and eqvect(rau, xau[nb], o["ft"])]
        if len(found) != 1:
            return None
        irt[na] = found[0]
    return irt


bad = 0
for o in ops:
    o["irt"] = compute_irt(o)
    if o["irt"] is None:
        bad += 1
print(f"irt found for {48 - bad}/48 operations using QE's rule (rau = s^T xau, xau(nb) = rau - ft)")
assert bad == 0

# ------------------------------------------------------------------------------------------
# 3. port of find_u / versor / angle_rot (PW/src/divide_class*.f90) and comp_dspinldau
# ------------------------------------------------------------------------------------------
EPS = 1e-7


def is_180(sm):
    return abs(np.trace(sm) + 1.0) < 1e-6


def versor(sm):
    if is_180(sm):
        ax = np.zeros(3)
        for ip in range(3):
            if abs(sm[ip, ip] - 1.0) < EPS:
                ax[ip] = 1.0
        if np.linalg.norm(ax) > EPS:
            return ax
        a1 = np.sqrt(np.abs(np.diag(sm) + 1.0) / 2.0)
        for ip in range(3):
            for jp in range(ip + 1, 3):
                if abs(a1[ip] * a1[jp]) > EPS:
                    a1[ip] = 0.5 * sm[ip, jp] / a1[jp]
    else:
        a1 = np.array([-sm[1, 2] + sm[2, 1], -sm[2, 0] + sm[0, 2], -sm[0, 1] + sm[1, 0]])
    if a1[2] < -EPS:
        a1 = -a1
    elif abs(a1[2]) < EPS and a1[0] < -EPS:
        a1 = -a1
    elif abs(a1[2]) < EPS and abs(a1[0]) < EPS and a1[1] < -EPS:
        a1 = -a1
    return a1 / np.linalg.norm(a1)


def angle_rot(sm):
    if is_180(sm):
        return 180.0
    a1 = np.array([-sm[1, 2] + sm[2, 1], -sm[2, 0] + sm[0, 2], -sm[0, 1] + sm[1, 0]])
    sint = 0.5 * np.linalg.norm(a1)
    sint = min(sint, 1.0)
    ax = a1.copy()
    if ax[2] < -EPS:
        ax = -ax
    elif abs(ax[2]) < EPS and ax[1] < -EPS:
        ax = -ax
    elif abs(ax[2]) < EPS and abs(ax[1]) < EPS and ax[0] < -EPS:
        ax = -ax
    if abs(a1[0]) > EPS:
        sint = np.copysign(sint, a1[0] / ax[0])
    elif abs(a1[1]) > EPS:
        sint = np.copysign(sint, a1[1] / ax[1])
    elif abs(a1[2]) > EPS:
        sint = np.copysign(sint, a1[2] / ax[2])
    ax = a1 / (2.0 * sint)
    if abs(ax[0] ** 2 - 1.0) > EPS:
        cost = (sm[0, 0] - ax[0] ** 2) / (1.0 - ax[0] ** 2)
    elif abs(ax[1] ** 2 - 1.0) > EPS:
        cost = (sm[1, 1] - ax[1] ** 2) / (1.0 - ax[1] ** 2)
    else:
        cost = (sm[2, 2] - ax[2] ** 2) / (1.0 - ax[2] ** 2)
    ang = np.degrees(np.arcsin(sint))
    if ang < 0:
        ang = -ang + 180.0 if cost < 0 else 360.0 + ang
    elif cost < 0:
        ang = -ang + 180.0
    return ang


def find_u(s):
    det = np.linalg.det(s)
    saux = -s if abs(det + 1.0) < 1e-8 else s.copy()
    if np.allclose(saux, np.eye(3), atol=1e-8):
        return np.eye(2, dtype=complex)
    ax, ang = versor(saux), 0.5 * np.radians(angle_rot(saux))
    c, sn = np.cos(ang), np.sin(ang)
    u = np.zeros((2, 2), complex)
    u[0, 0] = complex(c, -ax[2] * sn)
    u[0, 1] = complex(-ax[1] * sn, -ax[0] * sn)
    u[1, 0] = -np.conj(u[0, 1])
    u[1, 1] = np.conj(u[0, 0])
    if c < -1e-8:
        u = -u
    return u


def comp_dspinldau(o):
    u = find_u(o["sr"])
    if o["t_rev"] == 1:                       # d -> i sigma_y d^*
        a, b = np.conj(u[0, 0]), np.conj(u[0, 1])
        u = np.array([[np.conj(u[1, 0]), np.conj(u[1, 1])], [-a, -b]])
    return u


SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]


def so3_from_su2(u):
    """R_ij = 1/2 Tr(sigma_i U sigma_j U^dag): the ACTIVE rotation U implements."""
    return np.array([[0.5 * np.trace(SIG[i] @ u @ SIG[j] @ u.conj().T).real for j in range(3)] for i in range(3)])



# 4. QE's initial ns (init_ns_nc): the reference invariant state.  ns(1..4) = (uu, ud, du, dd)
# ------------------------------------------------------------------------------------------
ang = {5: (125.264389682755, 225.0), 6: (54.735610317245, 135.0), 7: (54.735610317245, 315.0), 8: (125.264389682755, 45.0)}
n_hat = {}
NS0 = {}   # (2x2 spin matrix, QE convention), orbital part is the identity
for a, (a1d, a2d) in ang.items():
    th, ph = np.radians(a1d), np.radians(a2d)
    n_hat[a] = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    cosin, esin = np.cos(th), (np.cos(ph) + 1j * np.sin(ph)) * np.sin(th)
    n, m = 1.4, 0.6
    NS0[a] = np.array([[(n + m * cosin) / 2, m * esin / 2], [m * np.conj(esin) / 2, (n - m * cosin) / 2]])



for o in ops:
    o['irt_inv'] = {o['irt'][a]: a for a in o['irt']}
    o['dS'] = comp_dspinldau(o)
rng = np.random.default_rng(1)
# ---- d_matrix (PW/src/d_matrix.f90): Y_m(S r) = sum_m' D[m,m'] Y_m'(r), 5 random points, real l=2 harmonics ---
def y2(r):
    r = r / np.linalg.norm(r, axis=0)
    x, y, z = r
    c = np.array([np.sqrt(5 / (16 * np.pi)), np.sqrt(15 / (4 * np.pi)), np.sqrt(15 / (4 * np.pi)),
                  np.sqrt(15 / (16 * np.pi)), np.sqrt(15 / (4 * np.pi))])
    return c[:, None] * np.array([3 * z * z - 1, x * z, y * z, x * x - y * y, x * y])


pts = rng.uniform(-0.5, 0.5, size=(3, 5))
Yinv = np.linalg.inv(y2(pts))
for o in ops:
    D = y2(o["sr"] @ pts) @ Yinv
    assert np.allclose(D @ D.T, np.eye(5), atol=1e-9), "D not orthogonal"
    o["D"] = D
print("d_matrix port: all 48 D_S (l=2) orthogonal")

# composite index: (m, s) with m fastest, s slowest, as in QE's proj(offsetU+m+ldim*(is-1)); tensor X[m1,m2,s1,s2]
def Ut(o):
    return np.kron(o["dS"].conj(), o["D"])            # spin outer, orbital inner: index = m + 5*s


def to_mat(X):   # X[m1,m2,s1,s2] -> 10x10 with rows (s1,m1), cols (s2,m2)
    return X.transpose(2, 0, 3, 1).reshape(10, 10)


def from_mat(M):
    return M.reshape(2, 5, 2, 5).transpose(1, 3, 0, 2)


def forward(o, Xa):
    """what the operation does to the ns of one atom (QE's X = rho^T convention)."""
    M = to_mat(Xa)
    U = Ut(o)
    return from_mat(U @ (M.T if o["t_rev"] == 1 else M) @ U.conj().T)


def qe_as_written(X):
    """new_ns_nc literally: nr1(m1,m2,is1,is2,na) += conj(dS(is1,is3)) d(m1,m3) nr[(m3|m4),(is3|is4)](nb) dS(is2,is4) d(m2,m4)/nsym
       with nb = irt(isym,na);  time reversal uses nr(m4,m3,is4,is3)."""
    out = {a: np.zeros((5, 5, 2, 2), complex) for a in IR}
    for a in IR:
        for o in ops:
            nb = o["irt"][a]
            Xb = X[nb].transpose(1, 0, 3, 2) if o["t_rev"] == 1 else X[nb]
            dS, D = o["dS"], o["D"]
            out[a] += np.einsum("ax,by,ik,jl,klxy->ijab", dS.conj(), dS, D, D, Xb) / 48
    return out


def fixed(X):
    """same forms, PRE-image atom: nb = irt^-1(isym, na)."""
    out = {a: np.zeros((5, 5, 2, 2), complex) for a in IR}
    for a in IR:
        for o in ops:
            nb = o["irt_inv"][a]
            Xb = X[nb].transpose(1, 0, 3, 2) if o["t_rev"] == 1 else X[nb]
            dS, D = o["dS"], o["D"]
            out[a] += np.einsum("ax,by,ik,jl,klxy->ijab", dS.conj(), dS, D, D, Xb) / 48
    return out


def diff(A, B):
    return max(np.abs(A[a] - B[a]).max() for a in IR)


# sanity: my einsum reproduces the loop structure  conj(dS[is1,is3]) d[m1,m3] X[m3,m4,is3,is4] dS[is2,is4] d[m2,m4]
# (indices: a=is1,x=is3,b=is2,y=is4,i=m1,k=m3,j=m2,l=m4)

# ---- 1. reference invariant state built independently: group-average of random Hermitian input by the FORWARD action ----
def random_state():
    st = {}
    for a in IR:
        A = rng.normal(size=(10, 10)) + 1j * rng.normal(size=(10, 10))
        st[a] = from_mat(A + A.conj().T)
    return st


# an invariant state: X_a = (1/48) sum_g  Ut_g^dag-type backward image of a random input at atom irt(g,a):
#   invariant  <=>  X_{irt(g,a)} = forward_g(X_a).  Build it as X_a = mean_g backward_g( R_{irt(g,a)} ).
def backward(o, Xb):
    U = Ut(o)
    M = to_mat(Xb)
    Mb = U.conj().T @ M @ U                        # inverse of the unitary part of the forward map
    return from_mat(Mb.T if o["t_rev"] == 1 else Mb)  # forward: X_b = Ut (X_a^T) Ut^dag  ->  X_a = (Ut^dag X_b Ut)^T

R = random_state()
Xtrue = {a: sum(backward(o, R[o["irt"][a]]) for o in ops) / 48 for a in IR}
# hermiticity in the composite index
for a in IR:
    M = to_mat(Xtrue[a]); assert np.allclose(M, M.conj().T, atol=1e-12)

# ---- 2. is Xtrue invariant under the FORWARD action of every one of the 48 operations? (validates the group action)
dev = 0.0
for o in ops:
    for a in IR:
        dev = max(dev, np.abs(forward(o, Xtrue[a]) - Xtrue[o["irt"][a]]).max())
print(f"\ninvariance of the reference state under the forward action of all 48 operations: max deviation {dev:.2e}")

# ---- 3. the same on the analytic AIAO ns (orbital identity) ---------------------------------------------------
NS0f = {a: np.einsum("ij,ab->ijab", np.eye(5), NS0[a]) for a in IR}
dev0 = max(np.abs(forward(o, NS0f[a]) - NS0f[o["irt"][a]]).max() for o in ops for a in IR)
print(f"invariance of QE's initial AIAO ns under the forward action of all 48 operations:  max deviation {dev0:.2e}")

# ---- 4. the symmetrisers on the invariant states --------------------------------------------------------------------
print("\n                                              QE as written (nb = irt)      same forms, nb = irt^-1")
for name, st in (("random invariant orbital x spin state", Xtrue), ("analytic AIAO ns0 (|m| = 3)", NS0f)):
    print(f"  fixed-point test  {name:38s} {diff(qe_as_written(st), st):.3e}                  {diff(fixed(st), st):.3e}")
Y = random_state()
P1, F1 = qe_as_written(Y), fixed(Y)
print(f"  idempotence P(P(X)) - P(X) on a random input   {'':17s} {diff(qe_as_written(P1), P1):.3e}                  {diff(fixed(F1), F1):.3e}")
print(f"  is the output invariant under the forward action?  QE: {max(np.abs(forward(o, P1[a]) - P1[o['irt'][a]]).max() for o in ops for a in IR):.3e}   fixed: {max(np.abs(forward(o, F1[a]) - F1[o['irt'][a]]).max() for o in ops for a in IR):.3e}")

# ---- 5. what it does to the physical moment ---------------------------------------------------------------------
def moment(Xa):
    rho_spin = np.einsum("iiab->ab", Xa).T          # trace over orbital, transpose back to rho
    return np.array([np.trace(s @ rho_spin).real for s in SIG])
print("\nsite moment m = Tr[sigma . rho] after symmetrising the exact AIAO ns0 (input |m| = 3):")
for a in IR:
    print(f"   Ir{a-4}: QE as written |m| = {np.linalg.norm(moment(qe_as_written(NS0f)[a])):.4f}    fixed |m| = {np.linalg.norm(moment(fixed(NS0f)[a])):.4f}")

# ---- 6. the alternative fix: keep nb = irt(isym,na) and transpose the contraction (upstream draft MR !2512) ----------
# MR !2512 (gitlab.com/QEF/q-e, opened 2024-12-11, still a draft with no review on 2026-09-21; raw diff at
# .../merge_requests/2512.diff, lines 50-73 for l = 2) leaves nb alone and changes only the index order and the
# CONJG placement:
#     unitary : CONJG(dS(is3,is1)) d(m3,m1) nr(m3,m4,is3,is4,nb) dS(is4,is2) d(m4,m2)
#     t_rev   :       dS(is3,is1)  d(m3,m1) nr(m4,m3,is4,is3,nb) CONJG(dS(is4,is2)) d(m4,m2)
# "transposed" below is the transposed contraction with the CONJG placement that follows from X_b = Ut X_a Ut^dag,
# Ut = D (x) conj(dS), i.e. the inverse of the forward map (unitary: plain first; t_rev: CONJG first).
E = "xa,yb,ki,lj,klxy->ijab"      # x=is3 a=is1 y=is4 b=is2 k=m3 i=m1 l=m4 j=m2


def _term(kind, o, a, X):
    tr = o["t_rev"] == 1
    nb = o["irt_inv"][a] if kind == "fixed" else o["irt"][a]
    Xb = X[nb].transpose(1, 0, 3, 2) if tr else X[nb]
    dS, D = o["dS"], o["D"]
    if kind in ("qe", "fixed"):
        return np.einsum("ax,by,ik,jl,klxy->ijab", dS.conj(), dS, D, D, Xb)
    conj_first = (not tr) if kind == "mr2512" else tr
    f, s = (dS.conj(), dS) if conj_first else (dS, dS.conj())
    return np.einsum(E, f, s, D, D, Xb)


def sym(kind, X):
    return {a: sum(_term(kind, o, a, X) for o in ops) / 48 for a in IR}


def fwd_dev(P):
    return max(np.abs(forward(o, P[a]) - P[o["irt"][a]]).max() for o in ops for a in IR)


def order(o):
    M, n = np.eye(3), 0
    while n < 12:
        M, n = o["sr"] @ M, n + 1
        if np.allclose(M, np.eye(3), atol=1e-6):
            return n
    return -1


print("\nsymmetrisers compared: deviation from a fixed point on the two invariant states, idempotence, invariance of the")
print("output under the 48 forward actions, and |m| returned from the exact AIAO ns0 (|m| = 3)")
print(f"  {'variant':34s} {'random inv.':>11s} {'AIAO ns0':>9s} {'idempot.':>9s} {'invariant':>9s} {'|m|':>6s}")
for kind, label in (("qe", "QE 7.5 / develop, nb = irt"), ("mr2512", "MR !2512 as drafted, nb = irt"),
                    ("transposed", "transposed + derived CONJG, nb = irt"), ("fixed", "as written, nb = irt(invs)  [patch]")):
    P = sym(kind, Y)
    m = np.mean([np.linalg.norm(moment(sym(kind, NS0f)[a])) for a in IR])
    print(f"  {label:34s} {diff(sym(kind, Xtrue), Xtrue):11.1e} {diff(sym(kind, NS0f), NS0f):9.1e} "
          f"{diff(sym(kind, P), P):9.1e} {fwd_dev(P):9.1e} {m:6.3f}")

print("\nper-operation test on the exact AIAO ns0: how many single terms of the group sum fail to reproduce ns0 (> 1e-8)")
for kind in ("qe", "mr2512", "transposed", "fixed"):
    tally = {}
    for o in ops:
        key = f"order {order(o)}{'T' if o['t_rev'] else ''}"
        bad = max(np.abs(_term(kind, o, a, NS0f) - NS0f[a]).max() for a in IR) > 1e-8
        n, f_ = tally.get(key, (0, 0))
        tally[key] = (n + 1, f_ + bad)
    print(f"  {kind:11s}", "   ".join(f"{k}: {f_}/{n}" for k, (n, f_) in sorted(tally.items())))
