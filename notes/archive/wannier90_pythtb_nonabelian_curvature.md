> ⚠ **Superseded.** This is the earliest draft of the derivation now canonical
> in [`theory/derivation.md`](../theory/derivation.md)
> (2026-08-04). A cleaned-up intermediate revision is
> [`external_derivation.md`](external_derivation.md) (2026-08-03). Kept for
> provenance only — nothing in the codebase should cite this version. See
> [`notes/README.md`](../README.md).

# Wannier90, PythTB, and the gauge-covariant non-Abelian Berry curvature

This note derives the position-matrix corrections needed to recover the
gauge-covariant non-Abelian Berry curvature of an occupied subspace from a
finite Wannier representation. The main goals are to explain:

1. why Wannier position matrix elements are Fourier coefficients of a Berry
   connection;
2. how Wannier90 convention II is converted to PythTB convention I;
3. why the complementary projector appears in the non-Abelian curvature;
4. how that projector separates into states inside and outside the finite
   Wannier subspace;
5. how the internal, cross, and external curvature terms arise; and
6. how the result depends on the number \(J\) of Wannier functions.

Throughout most of the derivation, \(\mathbf{k}\) is a Cartesian wavevector and
\(\partial_i\) means differentiation with respect to a Cartesian momentum
component or another adiabatic coordinate. Reduced-coordinate factors used by
PythTB are discussed separately.

## 1. Wannier representation

Let

$$
|\mathbf{R},s\rangle,
\qquad
s=1,\ldots,J,
$$

be a set of \(J\) orthonormal Wannier functions per cell:

$$
\langle \mathbf{R}',s|\mathbf{R},t\rangle
=
\delta_{\mathbf{R}',\mathbf{R}}\delta_{st}.
$$

Their centers are

$$
\boldsymbol{\tau}_s
=
\langle \mathbf{0},s|\hat{\mathbf r}|\mathbf{0},s\rangle.
$$

Wannier90 writes these centers to `seedname_centres.xyz`. PythTB reads them,
converts them to reduced coordinates, and uses them as the positions of the
tight-binding orbitals.

### 1.1 Wannier90 convention II

The usual Wannier90 Bloch-like states are

$$
|u_{s\mathbf{k}}^{W,\mathrm{II}}\rangle
=
\frac{1}{\sqrt{N}}
\sum_{\mathbf{R}}
e^{-i\mathbf{k}\cdot(\hat{\mathbf r}-\mathbf{R})}
|\mathbf{R},s\rangle.
$$

The corresponding Wannier Hamiltonian is

$$
H_{st}^{W,\mathrm{II}}(\mathbf{k})
=
\langle u_{s\mathbf{k}}^{W,\mathrm{II}}|
\hat H_{\mathbf{k}}
|u_{t\mathbf{k}}^{W,\mathrm{II}}\rangle.
$$

After disentanglement and localization, Wannier90 obtains a \(J\)-dimensional
Wannier subspace. If \(V(\mathbf{k})\) denotes the combined disentanglement and
localization transformation, then

$$
H^{W,\mathrm{II}}(\mathbf{k})
=
V^\dagger(\mathbf{k})
\operatorname{diag}(E_{n\mathbf{k}})
V(\mathbf{k}).
$$

Wannier90 Fourier transforms this Hamiltonian:

$$
H_{st}^{W}(\mathbf{R})
=
\frac{1}{N_k}
\sum_{\mathbf{k}}
e^{-i\mathbf{k}\cdot\mathbf{R}}
H_{st}^{W,\mathrm{II}}(\mathbf{k}).
$$

The real-space matrix elements are

$$
H_{st}^{W}(\mathbf{R})
=
\langle \mathbf{0},s|\hat H|\mathbf{R},t\rangle.
$$

They are written to `seedname_hr.dat`, or to the Hamiltonian section of
`seedname_tb.dat` when `write_tb = true`.

The files store Wigner-Seitz degeneracy weights \(w_{\mathbf R}\). Define the
effective real-space matrix

$$
\widetilde H_{st}(\mathbf R)
=
\frac{H_{st}^{W}(\mathbf R)}{w_{\mathbf R}}.
$$

Then convention II is reconstructed as

$$
H_{st}^{W,\mathrm{II}}(\mathbf k)
=
\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf R}
\widetilde H_{st}(\mathbf R).
$$

An optional energy-zero shift is applied as

$$
H^{W,\mathrm{II}}(\mathbf k)
\longrightarrow
H^{W,\mathrm{II}}(\mathbf k)-E_{\mathrm{zero}}I.
$$

### 1.2 PythTB convention I

PythTB uses convention I:

$$
|u_{s\mathbf{k}}^{W,\mathrm I}\rangle
=
e^{i\mathbf{k}\cdot\boldsymbol{\tau}_s}
|u_{s\mathbf{k}}^{W,\mathrm{II}}\rangle.
$$

Equivalently,

$$
|u_{s\mathbf{k}}^{W,\mathrm I}\rangle
=
\frac{1}{\sqrt N}
\sum_{\mathbf R}
e^{-i\mathbf k\cdot
(\hat{\mathbf r}-\mathbf R-\boldsymbol{\tau}_s)}
|\mathbf R,s\rangle.
$$

Define the bond vector

$$
\mathbf b_{st}(\mathbf R)
=
\mathbf R+\boldsymbol{\tau}_t-\boldsymbol{\tau}_s.
$$

PythTB reconstructs

$$
H_{st}^{W,\mathrm I}(\mathbf k)
=
\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
\widetilde H_{st}(\mathbf R)
-E_{\mathrm{zero}}\delta_{st}.
$$

The convention-I and convention-II Hamiltonians are related by

$$
H^{W,\mathrm I}(\mathbf k)
=
D^\dagger(\mathbf k)
H^{W,\mathrm{II}}(\mathbf k)
D(\mathbf k),
$$

where

$$
D_{st}(\mathbf k)
=
\delta_{st}e^{i\mathbf k\cdot\boldsymbol{\tau}_s}.
$$

They have identical eigenvalues but different coefficient eigenvectors:

$$
U_{sn}^{\mathrm I}(\mathbf k)
=
e^{-i\mathbf k\cdot\boldsymbol{\tau}_s}
U_{sn}^{\mathrm{II}}(\mathbf k).
$$

Consequently, the Hamiltonian matrix reconstructed by PythTB is not literally
the convention-II matrix that Wannier90 Fourier transformed. It is the matrix
of the same projected operator in a different, \(\mathbf{k}\)-dependent basis.

## 2. Why the position matrix is a Berry connection

Wannier90 provides the real-space position matrix

$$
r_{\mu,st}(\mathbf R)
=
\langle \mathbf 0,s|\hat r_\mu|\mathbf R,t\rangle.
$$

Current PythTB reads this from `seedname_tb.dat`, produced by
`write_tb = true`. A legacy `seedname_r.dat` containing the same position-only
data is also supported. As with \(H(\mathbf R)\), define

$$
\widetilde r_{\mu,st}(\mathbf R)
=
\frac{r_{\mu,st}(\mathbf R)}{w_{\mathbf R}}.
$$

The convention-II Wannier-basis Berry connection is

$$
A_{\mu,st}^{W,\mathrm{II}}(\mathbf k)
=
i\langle
u_{s\mathbf k}^{W,\mathrm{II}}
|
\partial_{k_\mu}
u_{t\mathbf k}^{W,\mathrm{II}}
\rangle.
$$

Differentiating the convention-II Bloch sum gives

$$
i\partial_{k_\mu}
|u_{t\mathbf k}^{W,\mathrm{II}}\rangle
=
\frac{1}{\sqrt N}
\sum_{\mathbf R}
(\hat r_\mu-R_\mu)
e^{-i\mathbf k\cdot(\hat{\mathbf r}-\mathbf R)}
|\mathbf R,t\rangle.
$$

Substituting this expression into the connection, using translational
covariance, and setting

$$
\mathbf L=\mathbf R-\mathbf R',
$$

gives

$$
A_{\mu,st}^{W,\mathrm{II}}
=
\frac{1}{N}
\sum_{\mathbf R',\mathbf R}
e^{i\mathbf k\cdot(\mathbf R-\mathbf R')}
\langle\mathbf R',s|
\hat r_\mu-R_\mu
|\mathbf R,t\rangle,
$$

and

$$
\langle\mathbf R',s|
\hat r_\mu-R_\mu
|\mathbf R,t\rangle
=
\langle\mathbf 0,s|
\hat r_\mu-L_\mu
|\mathbf L,t\rangle.
$$

The term proportional to \(L_\mu\) vanishes because

$$
L_\mu
\langle\mathbf 0,s|\mathbf L,t\rangle
=
L_\mu\delta_{\mathbf L,\mathbf 0}\delta_{st}
=0.
$$

The remaining sum over \(\mathbf R'\) cancels the normalization, giving

$$
A_{\mu,st}^{W,\mathrm{II}}(\mathbf k)
=
\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf R}
\widetilde r_{\mu,st}(\mathbf R).
$$

Thus the Wannier position matrix elements are the real-space Fourier
coefficients of the Wannier-gauge Berry connection.

This is also the representation-theoretic statement

$$
P_W\hat r_\mu P_W
\quad\longleftrightarrow\quad
i\partial_{k_\mu}+A_\mu^W(\mathbf k):
$$

position acts as a covariant derivative in momentum space.

## 3. Position connection in convention I

The basis transformation from convention II to convention I gives

$$
A_\mu^{W,\mathrm I}
=
D^\dagger A_\mu^{W,\mathrm{II}}D
+iD^\dagger\partial_{k_\mu}D.
$$

Since

$$
iD^\dagger\partial_{k_\mu}D
=
-\tau_\mu,
$$

the convention-I connection is

$$
A_\mu^{W,\mathrm I}
=
D^\dagger A_\mu^{W,\mathrm{II}}D-\tau_\mu.
$$

Define the residual convention-I position matrix

$$
X_{\mu,st}(\mathbf R)
=
\widetilde r_{\mu,st}(\mathbf R)
-\tau_{s\mu}\delta_{\mathbf R,\mathbf 0}\delta_{st}.
$$

Then

$$
A_{\mu,st}^{W}(\mathbf k)
=
\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
X_{\mu,st}(\mathbf R).
$$

Here and below \(A^W\) denotes the convention-I connection unless stated
otherwise.

The ordinary curl of this connection is

$$
\Omega_{\mu\nu,st}^{W}
=
\partial_{k_\mu}A_{\nu,st}^{W}
-\partial_{k_\nu}A_{\mu,st}^{W}.
$$

Therefore,

$$
\Omega_{\mu\nu,st}^{W}(\mathbf k)
=
i\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
\left[
b_{st,\mu}(\mathbf R)X_{\nu,st}(\mathbf R)
-b_{st,\nu}(\mathbf R)X_{\mu,st}(\mathbf R)
\right].
$$

Under the diagonal-position approximation,

$$
\widetilde r_{\mu,st}(\mathbf R)
\simeq
\tau_{s\mu}\delta_{\mathbf R,\mathbf 0}\delta_{st},
$$

so that

$$
X_{\mu,st}(\mathbf R)\simeq 0
$$

and

$$
A_\mu^W\simeq 0,
\qquad
\Omega_{\mu\nu}^W\simeq 0.
$$

This is why a convention-I point-orbital tight-binding model can compute its
geometry from Hamiltonian eigenvector derivatives alone.

## 4. Occupied and complementary subspaces

Diagonalize the \(J\times J\) convention-I Hamiltonian:

$$
\sum_t
H_{st}^{W}(\mathbf k)U_{tn}(\mathbf k)
=
E_n(\mathbf k)U_{sn}(\mathbf k).
$$

The corresponding eigenstates are

$$
|u_{n\mathbf k}\rangle
=
\sum_{s=1}^{J}
|u_{s\mathbf k}^{W}\rangle U_{sn}(\mathbf k).
$$

For a gapped insulator, let \(M\) of these states be occupied:

$$
n\in O=\{1,\ldots,M\}.
$$

The remaining \(J-M\) states are the complementary states inside the Wannier
subspace:

$$
c\in C=\{M+1,\ldots,J\}.
$$

They are often called "conduction" states in the code, but they are only the
unoccupied eigenstates represented inside the finite \(J\)-dimensional Wannier
space. They are not the complete set of first-principles conduction bands.

Define the projectors

$$
P_O
=
\sum_{n\in O}|u_n\rangle\langle u_n|,
$$

$$
P_C
=
\sum_{c\in C}|u_c\rangle\langle u_c|,
$$

and

$$
P_W
=
\sum_{s=1}^{J}|u_s^W\rangle\langle u_s^W|.
$$

Because the occupied and complementary Hamiltonian-gauge states together span
the Wannier space,

$$
P_W=P_O+P_C.
$$

The projector outside the Wannier space is

$$
Q_W=1-P_W.
$$

Therefore the complement of the occupied space is

$$
Q_O
=
1-P_O
=
P_C+Q_W.
$$

This exact identity is the origin of the two kinds of complementary
contributions:

1. \(P_C\) contains the \(J-M\) unoccupied states represented explicitly by the
   Wannier Hamiltonian.
2. \(Q_W\) contains every Hilbert-space direction omitted from the finite
   Wannier representation.

Equivalently,

$$
\mathcal H
=
\mathcal H_O
\oplus
\mathcal H_C
\oplus
\mathcal H_{\mathrm{out}},
$$

with dimensions

$$
M,
\qquad
J-M,
\qquad
\dim(\mathcal H)-J.
$$

## 5. Gauge-covariant non-Abelian curvature

For occupied states \(m,n\in O\), define the Hermitian non-Abelian connection

$$
\mathcal A_{i,mn}
=
i\langle u_m|\partial_i u_n\rangle.
$$

The index \(i\) may denote \(k_x,k_y,k_z\), or an adiabatic parameter such as a
phonon amplitude \(\lambda\).

The gauge-covariant non-Abelian curvature is

$$
\mathcal F_{ij}
=
\partial_i\mathcal A_j
-\partial_j\mathcal A_i
-i[\mathcal A_i,\mathcal A_j].
$$

Under an occupied-band gauge transformation \(G\),

$$
|u_n\rangle
\longrightarrow
\sum_{m\in O}|u_m\rangle G_{mn},
$$

the connection transforms inhomogeneously,

$$
\mathcal A_i
\longrightarrow
G^\dagger\mathcal A_iG+iG^\dagger\partial_iG,
$$

while the curvature transforms covariantly,

$$
\mathcal F_{ij}
\longrightarrow
G^\dagger\mathcal F_{ij}G.
$$

### 5.1 How the complementary projector appears

First compute the ordinary curl:

$$
\partial_i\mathcal A_{j,mn}
-\partial_j\mathcal A_{i,mn}
=
i\langle\partial_i u_m|\partial_j u_n\rangle
-i\langle\partial_j u_m|\partial_i u_n\rangle.
$$

The second derivatives cancel because

$$
\partial_i\partial_j=\partial_j\partial_i.
$$

For the commutator, use

$$
\langle\partial_i u_m|u_\ell\rangle
=
-\langle u_m|\partial_i u_\ell\rangle.
$$

It follows that

$$
\mathcal A_{i,m\ell}\mathcal A_{j,\ell n}
=
\langle\partial_i u_m|u_\ell\rangle
\langle u_\ell|\partial_j u_n\rangle.
$$

Summing over occupied \(\ell\),

$$
\sum_{\ell\in O}
\mathcal A_{i,m\ell}\mathcal A_{j,\ell n}
=
\langle\partial_i u_m|P_O|\partial_j u_n\rangle.
$$

Substitution into the curvature gives the projector expression

$$
\mathcal F_{ij,mn}
=
i\langle\partial_i u_m|Q_O|\partial_j u_n\rangle
-i\langle\partial_j u_m|Q_O|\partial_i u_n\rangle.
$$

Only the derivative components that leave the occupied space contribute.
Derivative components inside the occupied space merely rotate the occupied
basis and are removed by the commutator.

Using \(Q_O=P_C+Q_W\),

$$
\mathcal F_{ij}
=
\mathcal F_{ij}^{C}
+\mathcal F_{ij}^{Q_W}.
$$

## 6. Internal and external derivatives

Differentiate an occupied eigenstate:

$$
\partial_i|u_n\rangle
=
\sum_s
|u_s^W\rangle\partial_iU_{sn}
+
\sum_s
|\partial_i u_s^W\rangle U_{sn}.
$$

The first term is called internal:

$$
|\partial_i u_n\rangle_{\mathrm{int}}
=
\sum_s|u_s^W\rangle\partial_iU_{sn}.
$$

It changes the coefficient eigenvector inside the fixed Wannier basis.

The second term is called external:

$$
|\partial_i u_n\rangle_{\mathrm{ext}}
=
\sum_s|\partial_i u_s^W\rangle U_{sn}.
$$

It changes the Wannier basis itself.

### 6.1 Internal complementary amplitude

Project the internal derivative onto \(c\in C\):

$$
d_{i,cn}
=
\langle u_c|\partial_i u_n\rangle_{\mathrm{int}}
=
\sum_sU_{sc}^{*}\partial_iU_{sn}.
$$

Differentiating the Hamiltonian eigenvalue equation gives

$$
d_{i,cn}
=
\frac{
\displaystyle
\sum_{s,t}
U_{sc}^{*}
(\partial_iH_{st}^{W})
U_{tn}
}{
E_n-E_c
}.
$$

For a momentum direction,

$$
\partial_{k_\mu}H_{st}^{W}
=
i\sum_{\mathbf R}
b_{st,\mu}(\mathbf R)
e^{i\mathbf k\cdot\mathbf b_{st}(\mathbf R)}
\widetilde H_{st}(\mathbf R).
$$

The array \(d_i^{co}\) has shape

$$
(J-M)\times M.
$$

Here `co` means that its row index is complementary and its column index is
occupied.

### 6.2 External complementary amplitude

Rotate the Wannier-basis connection between complementary and occupied states:

$$
a_{i,cn}
=
\sum_{s,t}
U_{sc}^{*}A_{i,st}^{W}U_{tn}.
$$

Because

$$
\langle u_s^W|\partial_i u_t^W\rangle
=
-iA_{i,st}^{W},
$$

the external derivative amplitude is

$$
\langle u_c|\partial_i u_n\rangle_{\mathrm{ext}}
=
-ia_{i,cn}.
$$

Therefore,

$$
\langle u_c|\partial_i u_n\rangle
=
d_{i,cn}-ia_{i,cn}.
$$

This single equation generates the internal, cross, and external-external
terms.

## 7. Contribution from \(P_C\)

Insert the explicit complementary projector:

$$
\mathcal F_{ij,mn}^{C}
=
i\sum_{c\in C}
\left[
\langle\partial_i u_m|u_c\rangle
\langle u_c|\partial_j u_n\rangle
-
(i\leftrightarrow j)
\right].
$$

Using

$$
\langle u_c|\partial_i u_n\rangle
=
d_{i,cn}-ia_{i,cn},
$$

gives

$$
\mathcal F_{ij,mn}^{C}
=
i\sum_{c\in C}
\left[
(d_{i,cm}-ia_{i,cm})^{*}
(d_{j,cn}-ia_{j,cn})
-
(i\leftrightarrow j)
\right].
$$

Expanding produces three contributions.

### 7.1 Internal term

$$
\mathcal F_{ij,mn}^{\mathrm{int}}
=
i\sum_{c\in C}
\left[
d_{i,cm}^{*}d_{j,cn}
-d_{j,cm}^{*}d_{i,cn}
\right].
$$

This is the Hamiltonian-only non-Abelian Kubo curvature. It depends explicitly
on the \(J-M\) complementary states in the Wannier Hamiltonian.

### 7.2 Internal-external cross term

$$
\mathcal F_{ij,mn}^{\mathrm{cross}}
=
i\sum_{c\in C}
\left[
d_{i,cm}^{*}(-ia_{j,cn})
+
(-ia_{i,cm})^{*}d_{j,cn}
-
(i\leftrightarrow j)
\right].
$$

This term describes interference between coefficient motion and motion of the
Wannier basis.

### 7.3 External-external term inside \(P_C\)

$$
\mathcal F_{ij,mn}^{EE,C}
=
i\sum_{c\in C}
\left[
a_{i,cm}^{*}a_{j,cn}
-a_{j,cm}^{*}a_{i,cn}
\right].
$$

This term appears naturally at this intermediate stage. It will cancel against
part of the contribution from \(Q_W\).

## 8. Contribution from \(Q_W\)

The remaining piece of the occupied complement is

$$
\mathcal F_{ij,mn}^{Q_W}
=
i\langle\partial_i u_m|Q_W|\partial_j u_n\rangle
-i\langle\partial_j u_m|Q_W|\partial_i u_n\rangle.
$$

Because

$$
Q_W|u_s^W\rangle=0,
$$

the internal coefficient derivative is annihilated:

$$
Q_W
\sum_s|u_s^W\rangle\partial_iU_{sn}
=0.
$$

Hence

$$
Q_W|\partial_i u_n\rangle
=
\sum_sQ_W|\partial_i u_s^W\rangle U_{sn}.
$$

The \(Q_W\) contribution is purely external.

### 8.1 Wannier-subspace curvature

Define

$$
K_{ij,st}^{W}
=
\Omega_{ij,st}^{W}
-i[A_i^W,A_j^W]_{st}.
$$

Repeating the occupied-projector derivation for the complete \(J\)-dimensional
Wannier space gives

$$
K_{ij,st}^{W}
=
i\langle\partial_i u_s^W|Q_W|\partial_j u_t^W\rangle
-i\langle\partial_j u_s^W|Q_W|\partial_i u_t^W\rangle.
$$

Thus the curvature into states outside the Wannier space can be computed
without constructing those states explicitly. It is determined by the
Wannier-basis connection and its curl.

Rotate the ordinary Wannier curl into the Hamiltonian eigenbasis:

$$
\overline\Omega_{ij,mn}^{W}
=
\sum_{s,t}
U_{sm}^{*}\Omega_{ij,st}^{W}U_{tn}.
$$

Similarly, for any Hamiltonian-gauge indices \(p,q=1,\ldots,J\), define

$$
a_{i,pq}
=
\sum_{s,t}
U_{sp}^{*}A_{i,st}^{W}U_{tq}.
$$

Then

$$
\mathcal F_{ij,mn}^{Q_W}
=
\overline\Omega_{ij,mn}^{W}
-i\sum_{r=1}^{J}
\left[
a_{i,mr}a_{j,rn}
-a_{j,mr}a_{i,rn}
\right].
$$

Split the intermediate index \(r\) into occupied and complementary states:

$$
\mathcal F_{ij,mn}^{Q_W}
=
\overline\Omega_{ij,mn}^{W}
-i\sum_{\ell\in O}
\left[
a_{i,m\ell}a_{j,\ell n}
-a_{j,m\ell}a_{i,\ell n}
\right]
$$

$$
-
i\sum_{c\in C}
\left[
a_{i,mc}a_{j,cn}
-a_{j,mc}a_{i,cn}
\right].
$$

The connection is Hermitian:

$$
a_{i,mc}=a_{i,cm}^{*}.
$$

Therefore the final complementary sum is exactly

$$
-\mathcal F_{ij,mn}^{EE,C}.
$$

The explicit external-external term from \(P_C\) and the complementary part of
the \(Q_W\) curvature cancel.

## 9. Final gauge-covariant decomposition

After the cancellation,

$$
\mathcal F_{ij}
=
\mathcal F_{ij}^{\mathrm{int}}
+
\mathcal F_{ij}^{\mathrm{cross}}
+
\mathcal F_{ij}^{\mathrm{ext}}.
$$

The internal contribution is

$$
\mathcal F_{ij,mn}^{\mathrm{int}}
=
i\sum_{c\in C}
\left[
d_{i,cm}^{*}d_{j,cn}
-d_{j,cm}^{*}d_{i,cn}
\right].
$$

The cross contribution is

$$
\mathcal F_{ij,mn}^{\mathrm{cross}}
=
i\sum_{c\in C}
\left[
d_{i,cm}^{*}(-ia_{j,cn})
+
(-ia_{i,cm})^{*}d_{j,cn}
-
(i\leftrightarrow j)
\right].
$$

The purely external contribution is

$$
\mathcal F_{ij,mn}^{\mathrm{ext}}
=
\overline\Omega_{ij,mn}^{W}
-i\sum_{\ell\in O}
\left[
a_{i,m\ell}a_{j,\ell n}
-a_{j,m\ell}a_{i,\ell n}
\right].
$$

In block notation:

$$
\mathcal F_{ij}^{\mathrm{int}}
=
i\left[
(d_i^{co})^\dagger d_j^{co}
-
(d_j^{co})^\dagger d_i^{co}
\right],
$$

$$
\mathcal F_{ij}^{\mathrm{cross}}
=
i\left[
(d_i^{co})^\dagger(-ia_j^{co})
+
(-ia_i^{co})^\dagger d_j^{co}
-
(i\leftrightarrow j)
\right],
$$

and

$$
\mathcal F_{ij}^{\mathrm{ext}}
=
\overline\Omega_{ij}^{W,oo}
-i[a_i^{oo},a_j^{oo}].
$$

The block labels mean:

| Block | Row index | Column index | Shape |
|---|---|---|---|
| `oo` | occupied | occupied | \(M\times M\) |
| `co` | complementary | occupied | \((J-M)\times M\) |
| `oc` | occupied | complementary | \(M\times(J-M)\) |
| `cc` | complementary | complementary | \((J-M)\times(J-M)\) |

Only the total curvature is independent of how the calculation is partitioned
into internal, cross, and external pieces.

## 10. Dependence on the number of Wannier functions

The number of explicit complementary states is

$$
N_c=J-M.
$$

Both

$$
\mathcal F^{\mathrm{int}}
$$

and

$$
\mathcal F^{\mathrm{cross}}
$$

contain explicit sums over these \(J-M\) states.

The external term contains no explicit complementary-state sum after the
cancellation, but it still depends implicitly on \(J\), because changing the
Wannier space changes

$$
P_W,\quad Q_W,\quad A_i^W,\quad\Omega_{ij}^W,
$$

as well as the interpolated eigenstates and energies.

### 10.1 Minimal occupied Wannier space

If

$$
J=M,
$$

then

$$
P_C=0.
$$

There are no explicit complementary states, so

$$
\mathcal F^{\mathrm{int}}=0,
\qquad
\mathcal F^{\mathrm{cross}}=0.
$$

All curvature is carried by the motion of the occupied Wannier space:

$$
\mathcal F_{ij}
=
\overline\Omega_{ij}^{W,oo}
-i[a_i^{oo},a_j^{oo}].
$$

### 10.2 Complete point-orbital tight-binding space

At the opposite limit, suppose the finite tight-binding basis is treated as a
complete Hilbert space and its convention-I basis has

$$
A_i^W=0.
$$

Then

$$
\mathcal F^{\mathrm{cross}}=0,
\qquad
\mathcal F^{\mathrm{ext}}=0,
$$

and all curvature is obtained from the \(J-M\) explicit unoccupied states:

$$
\mathcal F_{ij}
=
\mathcal F_{ij}^{\mathrm{int}}.
$$

### 10.3 Finite first-principles Wannier space

A practical Wannier interpolation lies between these limits. Some unoccupied
directions are represented explicitly by \(P_C\); the rest remain in \(Q_W\)
and are encoded by the Wannier-basis connection and curvature.

As \(J\) increases, Hilbert-space directions may move from \(Q_W\) into \(P_C\).
The internal, cross, and external pieces then redistribute. If two Wannier
representations reproduce the same occupied projector exactly and use exact
position data, the total gauge-covariant curvature must agree:

$$
\mathcal F^{\mathrm{int}}
+
\mathcal F^{\mathrm{cross}}
+
\mathcal F^{\mathrm{ext}}
\quad\text{is independent of the partition.}
$$

In a numerical Wannier interpolation, convergence of the total depends on:

1. the frozen and outer windows;
2. the number and quality of Wannier functions;
3. localization of the Wannier functions;
4. the original Wannier90 k-mesh;
5. the accuracy of the position matrix; and
6. consistency of the Wannier gauge across any adiabatic parameter.

The individual pieces should not be expected to converge separately.

## 11. Relation to PythTB

For a Wannier90 import, PythTB uses:

| Mathematical quantity | Source or PythTB variable |
|---|---|
| \(\boldsymbol{\tau}_s\) | `seedname_centres.xyz` |
| \(H_{st}(\mathbf R)\) | `seedname_tb.dat` or `seedname_hr.dat` |
| \(r_{\mu,st}(\mathbf R)\) | `seedname_tb.dat` or legacy `seedname_r.dat` |
| \(X_{\mu,st}(\mathbf R)\) | `XI` |
| \(\mathbf b_{st}(\mathbf R)\) | `bnd` |
| \(A_{\mu,st}^{W}\) | `Aex` |
| \(\Omega_{\mu\nu,st}^{W}\) | `Om` |
| \(a_i^{co}\), \(a_i^{oo}\) | blocks of `a` |
| \(d_i^{co}\) | `d_co` |
| \(\mathcal F^{\mathrm{cross}}\) | `BX` |
| \(\mathcal F^{\mathrm{ext}}\) | `BE` |

The internal contribution is computed from the non-Abelian quantum geometric
tensor. With `non_abelian=True`, PythTB retains the full \(M\times M\) occupied
curvature:

$$
\mathcal F_{ij,mn}^{\mathrm{int}}
=
i\sum_c
\left[
d_{i,cm}^{*}d_{j,cn}
-d_{j,cm}^{*}d_{i,cn}
\right].
$$

If the model carries a position matrix and `include_external=True`, PythTB adds

$$
\mathcal F^{\mathrm{cross}}
+
\mathcal F^{\mathrm{ext}}.
$$

The implementation is in
`TBModel._berry_curvature_external` in `pythtb/tbmodel.py`.

One convention warning is important:
`W90.berry_connection_wann()` returns the plain convention-II transform

$$
A_\mu^{W,\mathrm{II}}(\mathbf k)
=
\sum_{\mathbf R}
e^{i\mathbf k\cdot\mathbf R}
\widetilde r_\mu(\mathbf R).
$$

It must not be combined directly with the convention-I
`TBModel.hamiltonian()`. The conversion

$$
A_\mu^{W,\mathrm I}
=
D^\dagger A_\mu^{W,\mathrm{II}}D-\tau_\mu
$$

must be applied. `TBModel.berry_curvature()` performs this conversion
internally.

## 12. Reduced coordinates and units

PythTB receives reduced momenta \(\boldsymbol{\kappa}\), with phase

$$
e^{2\pi i\boldsymbol{\kappa}\cdot\mathbf b_{\mathrm{red}}}.
$$

For Cartesian derivatives, PythTB converts the resulting one-forms and
two-forms using the reciprocal lattice. With

```python
cartesian=True
```

the input k-points are still reduced, but a momentum-momentum curvature is
returned in units of \(\text{Angstrom}^2\).

For a manual implementation, both the Hamiltonian and position blocks must be
divided by their Wigner-Seitz degeneracies before Fourier transformation.

## 13. Application to the axion response

For the static Chern-Simons axion angle,

$$
\theta
=
-\frac{1}{4\pi}
\int_{\mathrm{BZ}}d^3k\,
\epsilon^{ijk}
\operatorname{Tr}
\left[
\mathcal A_i\partial_j\mathcal A_k
-\frac{2i}{3}\mathcal A_i\mathcal A_j\mathcal A_k
\right].
$$

This expression requires the occupied non-Abelian connection itself and a
globally controlled gauge.

For the derivative with respect to an adiabatic parameter \(\lambda\), use the
extended coordinate space

$$
(k_x,k_y,k_z,\lambda).
$$

Then

$$
\partial_\lambda\theta
=
\frac{1}{16\pi}
\int_{\mathrm{BZ}}d^3k\,
\epsilon^{ijkl}
\operatorname{Tr}
\left[
\mathcal F_{ij}\mathcal F_{kl}
\right].
$$

The complete occupied \(M\times M\) curvature must be retained until after the
product is formed:

$$
\operatorname{Tr}(\mathcal F_{ij}\mathcal F_{kl})
\neq
\operatorname{Tr}(\mathcal F_{ij})
\operatorname{Tr}(\mathcal F_{kl}).
$$

The derivation above applies to all coordinate pairs \(i,j\), including mixed
\(k_\mu\)-\(\lambda\) planes, provided the corresponding Wannier-basis
connections and curls are known.

In the current axion-response implementation, the frozen-Wannier-gauge
approximation sets

$$
A_\lambda^W=0.
$$

The mixed external curl is then

$$
\Omega_{k_\mu\lambda}^{W}
=
-\partial_\lambda A_{k_\mu}^{W},
$$

approximated by a finite difference between consistently gauged Wannier
representations:

$$
\partial_\lambda A_{k_\mu}^{W}
\simeq
\frac{
A_{k_\mu}^{W}(\lambda+\Delta\lambda)
-A_{k_\mu}^{W}(\lambda)
}{
\Delta\lambda
}.
$$

The same number \(J\) of Wannier functions, the same orbital ordering, and a
common gauge across \(\lambda\) are essential. Otherwise the finite difference
does not represent a derivative of the same Wannier subspace.

## 14. Useful checks

A correct implementation should satisfy:

1. Hermiticity of the connection:

   $$
   A_i^{W\dagger}=A_i^W.
   $$

2. Hermiticity in occupied indices:

   $$
   \mathcal F_{ij,mn}^{*}
   =
   \mathcal F_{ij,nm}.
   $$

3. Antisymmetry in coordinate indices:

   $$
   \mathcal F_{ji}=-\mathcal F_{ij}.
   $$

4. Occupied-gauge covariance:

   $$
   \mathcal F_{ij}
   \longrightarrow
   G^\dagger\mathcal F_{ij}G.
   $$

5. Invariance of the second-Chern density:

   $$
   \operatorname{Tr}
   (\mathcal F_{ij}\mathcal F_{kl})
   \longrightarrow
   \operatorname{Tr}
   (\mathcal F_{ij}\mathcal F_{kl}).
   $$

6. Vanishing external correction in the diagonal-position approximation:

   $$
   X_i(\mathbf R)=0
   \quad\Longrightarrow\quad
   \mathcal F^{\mathrm{cross}}
   =
   \mathcal F^{\mathrm{ext}}
   =0.
   $$

7. Convergence of the total curvature and axion response with \(J\), even
   though the individual internal, cross, and external pieces may change.

## References

1. X. Wang, J. R. Yates, I. Souza, and D. Vanderbilt,
   *Ab initio calculation of the anomalous Hall conductivity by Wannier
   interpolation*, Phys. Rev. B **74**, 195118 (2006).
2. A. M. Essin, J. E. Moore, and D. Vanderbilt,
   *Magnetoelectric polarizability and axion electrodynamics in crystalline
   insulators*, Phys. Rev. Lett. **102**, 146805 (2009).
3. S. Coh, D. Vanderbilt, A. Malashevich, and I. Souza,
   *Chern-Simons orbital magnetoelectric coupling in generic insulators*,
   Phys. Rev. B **83**, 085108 (2011).
4. N. Marzari, A. A. Mostofi, J. R. Yates, I. Souza, and D. Vanderbilt,
   *Maximally localized Wannier functions: Theory and applications*,
   Rev. Mod. Phys. **84**, 1419 (2012).

