# Wannier external Berry curvature and moving centers

This note compares the external-term implementation (now in
`modules/curvature.py`) with
`references/vanderbilt_souza_2026_wannier_interpolation_notes.pdf`, especially the
external Berry-curvature terms and the mixed `(k_i, beta)` curvature.

The main conclusion is:

- The code implements the external `k-k` Berry-curvature planes consistently with
  the note, modulo the same PythTB/reduced-axis conventions.
- The code implements the full mixed `k-beta` external terms, including
  `A_beta`, all four cross terms, and `-i[a_i,a_beta]`.
- It does so in one fixed nominal convention-I embedding. Physical Wannier
  center motion remains in the residual position matrix, while a formulation
  whose nominal phase centers themselves move would require the separate terms
  derived below.

## 1. Notation

Wannier functions are labeled by orbital indices `s,t`. Lattice vectors are
`R`. The adiabatic parameter is `beta`.

The note uses Cartesian notation with phase

```math
e^{i k \cdot b_{st}(R)}.
```

The code uses reduced coordinates with

```math
e^{2\pi i k_{\rm red}\cdot b_{\rm red}}.
```

The bond vector is

```math
b_{st}(R) = R + \tau_t - \tau_s.
```

The position matrix in convention I is

```math
X_{\mu,st}(R) = \langle 0s | r_\mu - \tau_{s\mu} | Rt \rangle.
```

When `tau_s` is independent of `beta`, this is equivalent to using the raw
position matrix except for the `R=0`, `s=t` diagonal.

The external Wannier-gauge connection is

```math
A_{\mu,st}(k)
  =
  \sum_R e^{i k\cdot b_{st}(R)} X_{\mu,st}(R).
```

The code stores the three components of this object in `A_red`.

## 2. External k-k curvature in the note

From Eq. 18 of the PDF:

```math
A_{\mu,st}(k)
  =
  \sum_R e^{i k\cdot b_{st}(R)} X_{\mu,st}(R).
```

Differentiating with respect to `k_mu` gives Eq. 19:

```math
\partial_\mu A_{\nu,st}(k)
  =
  i \sum_R e^{i k\cdot b_{st}(R)}
    b_{\mu,st}(R) X_{\nu,st}(R).
```

The external curl is Eq. 21:

```math
\Omega_{\mu\nu,st}(k)
  =
  \partial_\mu A_{\nu,st}(k)
  -
  \partial_\nu A_{\mu,st}(k).
```

Therefore

```math
\Omega_{\mu\nu,st}(k)
  =
  i \sum_R e^{i k\cdot b_{st}(R)}
  \left[
    b_{\mu,st}(R) X_{\nu,st}(R)
    -
    b_{\nu,st}(R) X_{\mu,st}(R)
  \right].
```

For dual plane labels

```math
g=0,1,2 \quad \leftrightarrow \quad
(\mu,\nu)=(y,z),(z,x),(x,y),
```

define weights

```math
W_g(R)
  =
  i\left[
    b_\mu(R) X_\nu(R)
    -
    b_\nu(R) X_\mu(R)
  \right].
```

Then

```math
\Omega_g(k)
  =
  \sum_R e^{i k\cdot b(R)} W_g(R).
```

This is the `curl_red` field returned by the code.

## 3. Code mapping for the k-k fields

In `ExternalTerms.__init__`, the code starts from PythTB's real-space position
matrix:

```python
Rs = list(model._pos_r.keys())
R_int = np.array(Rs, dtype=float)
R_cart = R_int @ model.lattice.lat_vecs
deg = np.array([float(model._ham_deg.get(R, 1.0)) for R in Rs])
```

This is the list of real-space `R` vectors and Wigner-Seitz degeneracies.

It then Hermitianizes the position matrix:

```python
XI[j] = 0.5 * (model._pos_r[R] + np.conj(np.transpose(Pn, (0, 2, 1))))
```

This enforces the relation

```math
X_\mu(R) = X_\mu^\dagger(-R).
```

The code then zeroes the `R=0` diagonal:

```python
if R == (0, 0, 0):
    for ax in range(3):
        np.fill_diagonal(XI[j, ax], 0.0)
```

This is the convention-I choice

```math
X_{\mu,ss}^{I}(0) = 0,
```

meaning the centers are carried in the `tau` phase, not in the diagonal
position matrix.

Finally:

```python
XI /= deg[:, None, None, None]
```

implements the Wigner-Seitz degeneracy factor.

The code constructs the bond vectors:

```python
tau_cart = model.get_orb_vecs(cartesian=True)
bnd = (
    R_cart[:, None, None, :]
    + tau_cart[None, None, :, :]
    - tau_cart[None, :, None, :]
)
```

This is exactly

```math
b_{st}(R) = R + \tau_t - \tau_s.
```

The curl weights are:

```python
for g in range(3):
    al, be = _ALPHA_AX[g], _BETA_AX[g]
    W[:, g] = 1j * (bnd[..., al] * XI[:, be] - bnd[..., be] * XI[:, al])
```

This implements

```math
W_g(R)
  =
  i\left[
    b_\alpha(R) X_\beta(R)
    -
    b_\beta(R) X_\alpha(R)
  \right].
```

In `ExternalTerms.fields`, the code forms the Bloch sums:

```python
phR = np.exp(2j * np.pi * (kb @ self.R_int.T)).astype(self.dtype)
A = (phR @ self._XI_flat).reshape(nkb, 3, N, N)
Om = (phR @ self._W_flat).reshape(nkb, 3, N, N)
```

This gives the `e^{2 pi i k.R}` part. It then adds the convention-I intracell
phase:

```python
Dk = np.exp(2j * np.pi * (kb @ self.tau_red.T))
Dph = (np.conj(Dk)[:, :, None] * Dk[:, None, :]).astype(self.dtype)
A *= Dph[:, None]
Om *= Dph[:, None]
```

The phase factor is

```math
D_{st}(k)
  =
  e^{-2\pi i k\cdot\tau_s}
  e^{+2\pi i k\cdot\tau_t}
  =
  e^{2\pi i k\cdot(\tau_t-\tau_s)}.
```

Thus the total phase is

```math
e^{2\pi i k\cdot R}
e^{2\pi i k\cdot(\tau_t-\tau_s)}
=
e^{2\pi i k\cdot b_{st}(R)}.
```

Finally, the code converts Cartesian components to reduced-axis components:

```python
A = np.einsum("ua, kaij -> kuij", self._B_oneform, A)
Om = np.einsum("ua, kaij -> kuij", self._M_twoform, Om)
```

For one-forms:

```math
A_u^{\rm red} = B_{ua} A_a^{\rm cart}.
```

For dual two-form components:

```math
\Omega_u^{\rm red}
  =
  \left[\det(B) B^{-T}\right]_{ua} \Omega_a^{\rm cart}.
```

That is why `A` uses `B`, while `Om` uses the cofactor matrix.

## 4. Berry curvature assembly in the note

The note decomposes the gauge-covariant metric-curvature tensor as

```math
F_{\mu\nu}
  =
  \bar A_\mu^\dagger \bar A_\nu
  -
  (\bar A_\mu^E)^\dagger \bar A_\nu^E
  +
  W_{\mu\nu}^E
  -
  A_\mu^E A_\nu^E.
```

Here:

```math
\bar A_\mu = \bar A_\mu^I + \bar A_\mu^E.
```

The internal cross-gap connection is

```math
\bar A_{\mu,cv}^{I,H}
  =
  i
  \frac{
    \langle\!\langle u_c^H || \partial_\mu H || u_v^H \rangle\!\rangle
  }{
    E_v - E_c
  }.
```

The external cross-gap connection is

```math
\bar A_{\mu,cv}^{E}
  =
  (\bar U^\dagger A_\mu V)_{cv}.
```

The Berry curvature is

```math
B_{\mu\nu}
  =
  iF_{\mu\nu}
  -
  iF_{\nu\mu}.
```

Using Eq. 22 of the note, the antihermitian part of `W_{\mu\nu}` can be
replaced by the curl `Omega_{\mu\nu}`. Eq. 61 is then

```math
B_{\mu\nu}
  =
  i
  \left[
    \bar A_\mu^\dagger \bar A_\nu
    -
    (\bar A_\mu^E)^\dagger \bar A_\nu^E
  \right]
  + {\rm h.c.}
  +
  \Omega_{\mu\nu}^E
  -
  i [ A_\mu^E, A_\nu^E ].
```

Separating internal, external-internal cross, and purely external terms gives:

```math
B_{\mu\nu}
  =
  B_{\mu\nu}^{I}
  +
  B_{\mu\nu}^{X}
  +
  B_{\mu\nu}^{E}.
```

The internal Kubo term is

```math
B_{\mu\nu,mn}^{I}
  =
  i
  \sum_c
  \left[
    d_{\mu,cm}^* d_{\nu,cn}
    -
    d_{\nu,cm}^* d_{\mu,cn}
  \right],
```

with

```math
d_{\mu,cn}
  =
  \frac{
    \langle\!\langle u_c^H || \partial_\mu H || u_n^H \rangle\!\rangle
  }{
    E_n - E_c
  }.
```

The cross term is

```math
B_{\mu\nu}^{X}
  =
  i
  \left[
    d_\mu^\dagger (-i a_\nu)
    +
    (-i a_\mu)^\dagger d_\nu
    -
    (\mu \leftrightarrow \nu)
  \right].
```

The purely external occupied-space term is

```math
B_{\mu\nu}^{E}
  =
  \omega_{\mu\nu}^{oo}
  -
  i [a_\mu^{oo},a_\nu^{oo}].
```

Here

```math
a_\mu = U^\dagger A_\mu U,
\qquad
\omega_{\mu\nu}=U^\dagger \Omega_{\mu\nu} U.
```

## 5. Detailed derivation of the code formula from Eq. 57

This section derives the exact expression implemented in
`external_curvature`. It is the shortest useful path from the note's Eq. 57 to
the code.

Start from the note's Eq. 57:

```math
F_{\mu\nu}
  =
  \bar A_\mu^\dagger \bar A_\nu
  -
  (\bar A_\mu^E)^\dagger \bar A_\nu^E
  +
  W_{\mu\nu}^E
  -
  A_\mu^E A_\nu^E.
```

The full cross-gap connection is

```math
\bar A_\mu = \bar A_\mu^I + \bar A_\mu^E.
```

In the Hamiltonian gauge, Eq. 53 gives

```math
\bar A_{\mu,cv}^{I}
  =
  i
  \frac{
    \langle\!\langle u_c^H || \partial_\mu H || u_v^H \rangle\!\rangle
  }{
    E_v - E_c
  }.
```

Define the code's denominator-weighted velocity block

```math
d_{\mu,cv}
  =
  \frac{
    \langle\!\langle u_c^H || \partial_\mu H || u_v^H \rangle\!\rangle
  }{
    E_v - E_c
  }.
```

Then

```math
\bar A_\mu^I = i d_\mu.
```

Also define the Hamiltonian-gauge external connection blocks

```math
a_\mu^{co}
  =
  \bar A_\mu^E
  =
  \bar U^\dagger A_\mu V,
```

and

```math
a_\mu^{oo}
  =
  A_\mu^E
  =
  V^\dagger A_\mu V.
```

Therefore the full cross-gap connection is

```math
\bar A_\mu = i d_\mu + a_\mu^{co}.
```

Now isolate the first two terms of Eq. 57:

```math
G_{\mu\nu}
  =
  \bar A_\mu^\dagger \bar A_\nu
  -
  (\bar A_\mu^E)^\dagger \bar A_\nu^E.
```

Substitute `bar A_mu = i d_mu + a_mu`:

```math
G_{\mu\nu}
  =
  (i d_\mu + a_\mu)^\dagger
  (i d_\nu + a_\nu)
  -
  a_\mu^\dagger a_\nu.
```

Expand:

```math
G_{\mu\nu}
  =
  (-i d_\mu^\dagger + a_\mu^\dagger)
  (i d_\nu + a_\nu)
  -
  a_\mu^\dagger a_\nu.
```

Multiply out:

```math
G_{\mu\nu}
  =
  d_\mu^\dagger d_\nu
  -
  i d_\mu^\dagger a_\nu
  +
  i a_\mu^\dagger d_\nu
  +
  a_\mu^\dagger a_\nu
  -
  a_\mu^\dagger a_\nu.
```

The purely external cross-gap term cancels:

```math
G_{\mu\nu}
  =
  d_\mu^\dagger d_\nu
  -
  i d_\mu^\dagger a_\nu
  +
  i a_\mu^\dagger d_\nu.
```

Write the last two terms in the exact form used in the code:

```math
G_{\mu\nu}
  =
  d_\mu^\dagger d_\nu
  +
  d_\mu^\dagger(-i a_\nu)
  +
  (-i a_\mu)^\dagger d_\nu.
```

This is the main simplification. The full cross-gap connection is
`i d + a`, but after subtracting the pure external double count, only:

```math
d^\dagger d,
\qquad
d^\dagger(-i a),
\qquad
(-i a)^\dagger d
```

survive.

Now Berry curvature is Eq. 60:

```math
B_{\mu\nu}
  =
  iF_{\mu\nu}
  -
  iF_{\nu\mu}.
```

Applying this to the `d^\dagger d` piece gives the internal Kubo term:

```math
B_{\mu\nu}^{I}
  =
  i
  \left[
    d_\mu^\dagger d_\nu
    -
    d_\nu^\dagger d_\mu
  \right].
```

In indices:

```math
B_{\mu\nu,mn}^{I}
  =
  i
  \sum_c
  \left[
    d_{\mu,cm}^*d_{\nu,cn}
    -
    d_{\nu,cm}^*d_{\mu,cn}
  \right].
```

Applying the same antisymmetrization to the two mixed `d-a` terms gives

```math
B_{\mu\nu}^{X}
  =
  i
  \left[
    d_\mu^\dagger(-i a_\nu)
    +
    (-i a_\mu)^\dagger d_\nu
    -
    d_\nu^\dagger(-i a_\mu)
    -
    (-i a_\nu)^\dagger d_\mu
  \right].
```

This is the code's `BX`.

It remains to derive the purely occupied external term. From Eq. 57, the
occupied external part of `F` is

```math
F_{\mu\nu}^{E,oo}
  =
  W_{\mu\nu}^{E,oo}
  -
  a_\mu^{oo}a_\nu^{oo}.
```

Antisymmetrize:

```math
B_{\mu\nu}^{E}
  =
  i
  \left[
    W_{\mu\nu}^{E,oo}
    -
    a_\mu^{oo}a_\nu^{oo}
  \right]
  -
  i
  \left[
    W_{\nu\mu}^{E,oo}
    -
    a_\nu^{oo}a_\mu^{oo}
  \right].
```

Group the `W` terms:

```math
iW_{\mu\nu}^{E,oo}
-
iW_{\nu\mu}^{E,oo}
  =
  \Omega_{\mu\nu}^{E,oo},
```

using Eq. 22 of the note,

```math
\Omega_{\mu\nu}=iW_{\mu\nu}-iW_{\nu\mu}.
```

Group the connection-product terms:

```math
-i a_\mu^{oo}a_\nu^{oo}
+
i a_\nu^{oo}a_\mu^{oo}
  =
  -i[a_\mu^{oo},a_\nu^{oo}].
```

Thus:

```math
B_{\mu\nu}^{E}
  =
  \Omega_{\mu\nu}^{E,oo}
  -
  i[a_\mu^{oo},a_\nu^{oo}].
```

Putting the pieces together:

```math
B_{\mu\nu}
  =
  B_{\mu\nu}^{I}
  +
  B_{\mu\nu}^{X}
  +
  B_{\mu\nu}^{E},
```

where

```math
B_{\mu\nu}^{I}
  =
  i
  \left[
    d_\mu^\dagger d_\nu
    -
    d_\nu^\dagger d_\mu
  \right],
```

```math
B_{\mu\nu}^{X}
  =
  i
  \left[
    d_\mu^\dagger(-i a_\nu)
    +
    (-i a_\mu)^\dagger d_\nu
    -
    d_\nu^\dagger(-i a_\mu)
    -
    (-i a_\nu)^\dagger d_\mu
  \right],
```

and

```math
B_{\mu\nu}^{E}
  =
  \Omega_{\mu\nu}^{E,oo}
  -
  i[a_\mu^{oo},a_\nu^{oo}].
```

This is exactly the split used by the code:

- `berry_curvature` computes `B^I`.
- `external_curvature` computes `B^X + B^E`.

For a mixed plane `(i,beta)`, the same formulas apply if one supplies the
fourth-direction objects `d_beta`, `a_beta`, and `Omega_{i beta}`. The current
code supplies `d_beta`, but then assumes

```math
a_\beta = 0,
\qquad
\Omega_{i\beta}^{E} = -\partial_\beta A_i.
```

Under those assumptions:

```math
B_{i\beta}^{X}
  =
  i
  \left[
    (-i a_i)^\dagger d_\beta
    -
    d_\beta^\dagger(-i a_i)
  \right],
```

and

```math
B_{i\beta}^{E}
  =
  -(\partial_\beta A_i)^{oo}.
```

Those are the two lines

```python
BX = 1j * (FX_ib - FX_bi)
BE = -da_oo[:, i]
```

in `external_curvature`.

## 6. Code mapping for curvature assembly

The caller rotates velocity matrices and constructs the energy-denominator block:

```python
v_k_rot_tf = tf.matmul(
    evecs_conj_tf[None, :, :, :],
    tf.matmul(v_k_tf, evecs_t_tf[None, :, :, :]),
)

v_cond_occ_tf = tf.gather(tf.gather(v_k_rot_tf, cond_idxs, axis=-2), occ_idxs, axis=-1)
v_cond_occ_tf *= inv_delta_cond_occ
```

This is

```math
d_{\mu,cn}
  =
  \frac{
    \langle\!\langle u_c^H || \partial_\mu H || u_n^H \rangle\!\rangle
  }{
    E_n - E_c
  }.
```

The internal Kubo curvature is:

```python
q = tf.matmul(v_occ_cond_tf[:, None], v_cond_occ_tf[None, :]).numpy()
omega = 1j * (q - np.swapaxes(q, -1, -2).conj())
```

This is

```math
B_{\mu\nu}^{I}
  =
  i(q_{\mu\nu} - q_{\nu\mu}).
```

The external code receives the eigenvectors and the `d_co` block:

```python
curvature = curvature + external_curvature(
    evecs_tf.numpy(),
    v_cond_occ_tf.numpy(),
    occ_idxs,
    cond_idxs,
    A_red,
    curl_red,
    mixed=mixed_fields,
)
```

Inside `external_curvature`, the external fields are rotated:

```python
def rotate(X):
    return V_dagger[:, None] @ X @ V[:, None]

a = rotate(A_red)
curl = rotate(curl_red)
```

This implements

```math
a_\mu = U^\dagger A_\mu U,
\qquad
\omega_{\mu\nu}=U^\dagger \Omega_{\mu\nu} U.
```

Then the occupied and conduction-occupied blocks are extracted:

```python
a_co = a[:, :, cond[:, None], occ[None, :]]
a_oo = a[:, :, occ[:, None], occ[None, :]]
omz_oo = omz[:, :, occ[:, None], occ[None, :]]
```

For each k-k plane:

```python
FX_ab = adj(d_al) @ (-1j * a_be) + adj(-1j * a_al) @ d_be
FX_ba = adj(d_be) @ (-1j * a_al) + adj(-1j * a_be) @ d_al
BX = 1j * (FX_ab - FX_ba)
```

This is

```math
B_{\alpha\beta}^{X}
  =
  i
  \left[
    d_\alpha^\dagger (-i a_\beta)
    +
    (-i a_\alpha)^\dagger d_\beta
    -
    d_\beta^\dagger (-i a_\alpha)
    -
    (-i a_\beta)^\dagger d_\alpha
  \right].
```

The purely external term is:

```python
BE = omz_oo[:, g] - 1j * (ao_al @ ao_be - ao_be @ ao_al)
```

This is

```math
B_{\alpha\beta}^{E}
  =
  \omega_{\alpha\beta}^{oo}
  -
  i [a_\alpha^{oo},a_\beta^{oo}].
```

Finally:

```python
ext[al, be] = BX + BE
ext[be, al] = -(BX + BE)
```

This enforces

```math
B_{\beta\alpha} = -B_{\alpha\beta}.
```

Thus the k-k part of `modules/curvature.py` matches the note.

## 7. Fixed-center mixed k-beta formula in the note

The PDF's mixed-curvature section assumes the nominal centers `tau_s` are kept
independent of `beta`.

Define

```math
\Lambda_{st}(R) = \langle 0s | \partial_\beta | Rt \rangle.
```

Then the beta connection is Eq. 70:

```math
A_{\beta,st}(k)
  =
  \langle u_s^W | i\partial_\beta u_t^W\rangle
  =
  i \sum_R e^{i k\cdot b_{st}(R)} \Lambda_{st}(R).
```

The mixed curl is Eq. 72:

```math
\Omega_{\beta\mu,st}(k)
  =
  \partial_\beta A_{\mu,st}
  -
  \partial_\mu A_{\beta,st}.
```

Substitute the fixed-center expressions:

```math
\partial_\beta A_{\mu,st}
  =
  \sum_R e^{i k\cdot b_{st}(R)}
  \partial_\beta X_{\mu,st}(R),
```

because `b_{st}` is beta-independent, and

```math
\partial_\mu A_{\beta,st}
  =
  \partial_\mu
  \left[
    i \sum_R e^{i k\cdot b_{st}(R)} \Lambda_{st}(R)
  \right]
  =
  - \sum_R e^{i k\cdot b_{st}(R)}
  b_{\mu,st}(R)\Lambda_{st}(R).
```

Therefore:

```math
\Omega_{\beta\mu,st}(k)
  =
  \sum_R e^{i k\cdot b_{st}(R)}
  \left[
    \partial_\beta X_{\mu,st}(R)
    +
    b_{\mu,st}(R)\Lambda_{st}(R)
  \right].
```

Equivalently, with the index order used by the code:

```math
\Omega_{\mu\beta,st}(k)
  =
  -
  \sum_R e^{i k\cdot b_{st}(R)}
  \left[
    \partial_\beta X_{\mu,st}(R)
    +
    b_{\mu,st}(R)\Lambda_{st}(R)
  \right].
```

## 8. What the code computes for mixed k-beta planes

The driver uses the base Wannier centers as one fixed nominal convention-I
embedding for both endpoint models. It computes

```python
A0, curl0 = ext_base.fields(k_mesh[i:j])
A1, curl1 = ext_mode.fields(k_mesh[i:j])
dA_k = (A1 - A0) / d_beta
```

This is

```math
\partial_\beta A_\mu
  \approx
  \frac{A_\mu^{\rm mode}-A_\mu^{\rm base}}{\Delta\beta}.
```

The additional real-space input is

```math
\Lambda_{st}(R)=\langle 0s|\partial_\beta|Rt\rangle.
```

Wannier90 does not include this matrix in `*_tb.dat` or `*_r.dat`. The driver
loads it from `Y2Ir2O7_lambda_beta.npz` and prepares its Bloch interpolation:

```python
beta_external = ext_base.prepare_beta(lambda_r)
mixed_fields = beta_external.fields(k_mesh[i:j], dA_k)
```

The returned fields are

```math
A_\beta
  =
  i\sum_R e^{ik\cdot b(R)}\Lambda(R)
```
and, in the code's index order,
```math
\Omega_{i\beta}
  =
  \partial_i A_\beta - \partial_\beta A_i
  =
  -\sum_R e^{ik\cdot b(R)}
  \left[
    \partial_\beta X_i(R)+b_i(R)\Lambda(R)
  \right].
```

After rotation to the Hamiltonian gauge, `external_curvature` evaluates all
four complementary-projector cross terms,

```math
B_{i\beta}^{X}
  =
  i
  \left[
    d_i^\dagger (-i a_\beta)
    +
    (-i a_i)^\dagger d_\beta
    -
    d_\beta^\dagger (-i a_i)
    -
    (-i a_\beta)^\dagger d_i
  \right].
```

It also evaluates the full purely external term,

```math
B_{i\beta}^{E}
  =
  \Omega_{i\beta}^{oo}
  -i[a_i^{oo},a_\beta^{oo}].
```

Thus the code no longer assumes `A_beta=0`. It does assume that the endpoint
Wannier matrices and `Lambda(R)` have already been placed in one smooth orbital
gauge; a fixed set of nominal centers aligns the convention-I phases but cannot
align two otherwise unrelated Wannier gauges.

## 9. If centers move with beta

Now let

```math
\tau_s = \tau_s(\beta),
\qquad
\dot\tau_s = \partial_\beta \tau_s.
```

Then

```math
b_{st}(R,\beta)
  =
  R+\tau_t(\beta)-\tau_s(\beta),
```

and

```math
\dot b_{st}(R)
  =
  \dot\tau_t-\dot\tau_s.
```

The k-connection is still

```math
A_{\mu,st}(k,\beta)
  =
  \sum_R e^{i k\cdot b_{st}(R,\beta)}
  X^\tau_{\mu,st}(R,\beta),
```

where

```math
X^\tau_{\mu,st}(R,\beta)
  =
  \langle 0s | r_\mu-\tau_{s\mu}(\beta) | Rt\rangle.
```

The beta connection now has a phase contribution:

```math
A_{\beta,st}(k,\beta)
  =
  - k\cdot\dot\tau_s \, \delta_{st}
  +
  i\sum_R e^{i k\cdot b_{st}(R,\beta)}
  \Lambda_{st}(R,\beta).
```

The first term comes from differentiating the `tau_t` phase in the ket. The
overlap enforces `s=t`, so the result can be written with either `tau_s` or
`tau_t`.

Differentiate the k-connection with respect to beta:

```math
\partial_\beta A_{\mu,st}
  =
  \sum_R e^{i k\cdot b_{st}}
  \left[
    \partial_\beta X^\tau_{\mu,st}
    +
    i(k\cdot\dot b_{st})X^\tau_{\mu,st}
  \right].
```

Differentiate the beta connection with respect to k:

```math
\partial_\mu A_{\beta,st}
  =
  -\dot\tau_{s\mu}\delta_{st}
  -
  \sum_R e^{i k\cdot b_{st}}
  b_{\mu,st}\Lambda_{st}.
```

Therefore the full moving-center mixed curl is

```math
\Omega_{\beta\mu,st}
  =
  \partial_\beta A_{\mu,st}
  -
  \partial_\mu A_{\beta,st}
```

or

```math
\Omega_{\beta\mu,st}
  =
  \dot\tau_{s\mu}\delta_{st}
  +
  \sum_R e^{i k\cdot b_{st}}
  \left[
    \partial_\beta X^\tau_{\mu,st}
    +
    i(k\cdot\dot b_{st})X^\tau_{\mu,st}
    +
    b_{\mu,st}\Lambda_{st}
  \right].
```

The code uses the opposite orientation, so:

```math
\Omega_{\mu\beta,st}
  =
  -\dot\tau_{s\mu}\delta_{st}
  -
  \sum_R e^{i k\cdot b_{st}}
  \left[
    \partial_\beta X^\tau_{\mu,st}
    +
    i(k\cdot\dot b_{st})X^\tau_{\mu,st}
    +
    b_{\mu,st}\Lambda_{st}
  \right].
```

Using the raw position matrix

```math
X^{\rm raw}_{\mu,st}(R)
  =
  \langle 0s|r_\mu|Rt\rangle,
```

and

```math
X^\tau_{\mu,st}(R)
  =
  X^{\rm raw}_{\mu,st}(R)
  -
  \tau_{s\mu}\delta_{R0}\delta_{st},
```

we have

```math
\partial_\beta X^{\rm raw}_{\mu,st}
  =
  \partial_\beta X^\tau_{\mu,st}
  +
  \dot\tau_{s\mu}\delta_{R0}\delta_{st}.
```

Thus the moving-center formula can also be written as

```math
\Omega_{\beta\mu,st}
  =
  \sum_R e^{i k\cdot b_{st}}
  \left[
    \partial_\beta X^{\rm raw}_{\mu,st}
    +
    i(k\cdot\dot b_{st})X^\tau_{\mu,st}
    +
    b_{\mu,st}\Lambda_{st}
  \right].
```

and

```math
\Omega_{\mu\beta,st}
  =
  -
  \sum_R e^{i k\cdot b_{st}}
  \left[
    \partial_\beta X^{\rm raw}_{\mu,st}
    +
    i(k\cdot\dot b_{st})X^\tau_{\mu,st}
    +
    b_{\mu,st}\Lambda_{st}
  \right].
```

In reduced coordinates, replace

```math
i k\cdot b
```

by

```math
2\pi i k_{\rm red}\cdot b_{\rm red}.
```

Therefore the moving-center phase term becomes

```math
2\pi i
\left(
  k_{\rm red}\cdot \dot b_{\rm red}
\right)
X^\tau_{\mu}.
```

The beta-connection center term becomes

```math
A_{\beta,st}^{\rm center}
  =
  -2\pi
  \left(
    k_{\rm red}\cdot\dot\tau_{s,{\rm red}}
  \right)
  \delta_{st}.
```

## 10. Why the implementation fixes the nominal centers

Finite-differencing two endpoint connections is meaningful only if they use the
same convention-I embedding. The driver therefore passes the base centers as
`nominal_orb_vecs` to both `ExternalTerms` objects and uses those same positions
when constructing both endpoint Hamiltonians.

For the distorted model, the code does not zero the raw `R=0` position
diagonal. It subtracts the fixed nominal center:

```math
X_{i,ss}(0,\beta)
  =
  \langle 0s(\beta)|r_i|0s(\beta)\rangle
  -\tau_{s i}^{\rm nominal}.
```

Consequently, physical center motion contributes to `partial_beta X_i`. The
fixed nominal embedding removes only the avoidable beta dependence of the
Bloch phase convention. It does not freeze the Wannier functions or their
centers.

This still requires a smooth orbital gauge across beta. Reusing the same
nominal positions fixes diagonal convention-I phases, but it cannot by itself
align two independently rotated Wannier frames.

## 11. Full mixed-plane implementation

The external routine accepts:

```math
A_\beta(k),
\qquad
\Omega_{i\beta}(k).
```

After rotation:

```math
a_\beta = U^\dagger A_\beta U,
\qquad
\omega_{i\beta}=U^\dagger\Omega_{i\beta}U.
```

Then the mixed cross term should be

```math
B_{i\beta}^{X}
  =
  i
  \left[
    d_i^\dagger (-i a_\beta)
    +
    (-i a_i)^\dagger d_\beta
    -
    d_\beta^\dagger (-i a_i)
    -
    (-i a_\beta)^\dagger d_i
  \right].
```

The purely external term should be

```math
B_{i\beta}^{E}
  =
  \omega_{i\beta}^{oo}
  -
  i[a_i^{oo},a_\beta^{oo}].
```

The total external mixed curvature is

```math
B_{i\beta}^{\rm ext}
  =
  B_{i\beta}^{X}
  +
  B_{i\beta}^{E}.
```

The code evaluates these expressions directly. The former approximation is
obtained only if one deliberately supplies

```math
a_\beta = 0,
\qquad
\omega_{i\beta} = -\partial_\beta A_i,
```

which requires both:

```math
\Lambda=0,
\qquad
\dot\tau_s=0.
```

## 12. Practical options

There are two clean paths.

### Option A: fixed nominal centers

Choose a single `tau_s` set, usually the base centers or atom-centered projection
positions, and use it for both base and distorted models.

Then:

```math
\dot\tau_s = 0.
```

The mixed formula reduces to the PDF's Eq. 72:

```math
\Omega_{\beta\mu}
  =
  \sum_R e^{i k\cdot b}
  \left[
    \partial_\beta X_\mu + b_\mu\Lambda
  \right].
```

If one deliberately imposes `Lambda=0`, then:

```math
\Omega_{\mu\beta} = -\partial_\beta A_\mu.
```

This is only a special case of the formula now coded.

The important implementation detail is that `tau_s` must not be recomputed from
each endpoint's Wannier centers when forming `A_0` and `A_1`.

Keeping the nominal `tau_s` fixed does not freeze the physical Wannier centers.
Their displacement remains in the diagonal `R=0` part of the residual position
matrix `X`.

### Option B: moving centers

Let `tau_s(beta)` move with the Wannier centers.

Then the mixed curvature must include:

```math
\dot\tau_s,
\qquad
\dot b_{st},
\qquad
A_\beta,
\qquad
\Lambda.
```

Even if one assumes `Lambda=0`, the center term remains:

```math
A_{\beta,st}^{\rm center}
  =
  -k\cdot\dot\tau_s\,\delta_{st}
```

and

```math
\Omega_{i\beta}^{\rm center}
  =
  -\partial_i A_\beta^{\rm center}
  =
  -\dot\tau_{s i}\delta_{st}
```

in Cartesian note notation. In reduced coordinates, insert the corresponding
`2 pi` factors as described above.

The implementation chooses Option A. Option B would require the explicit center
phase terms above and is not mixed into the fixed-embedding path.

## 13. Validation status

The existing validator checks the `k-k` external terms against PythTB:

```python
ref_tot = model_base.berry_curvature(
    k_test, occ_idxs=occ, non_abelian=True, include_external=True
)
ref_int = model_base.berry_curvature(
    k_test, occ_idxs=occ, non_abelian=True, include_external=False
)
ref_ext = ref_tot - ref_int
```

and compares this to `external_curvature`.

That validates the `k-k` block. The validator also supplies a nonzero,
off-diagonal `Lambda(0)` and compares every mixed block against a direct
assembly of the four cross terms, the ordinary mixed curl, and
`-i[a_i,a_beta]`. It additionally checks antisymmetry in direction indices and
Hermiticity in occupied-band indices.
