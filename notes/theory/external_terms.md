# External (position-matrix) terms in the axion response

This note documents, in math, how the **external** Berry-curvature terms are
built and added to the **internal** (Kubo) terms in the axion-response pipeline.
The production code follows the same boundaries as the derivation:

| derivation object | code |
|---|---|
| Wannier inputs: the $H(R)/A(R)$ pair as `_hr.dat` / `_r.dat` | [`wannier/`](../../wannier/README.md) (stage 1); PythTB reads the files |
| $A$, $\operatorname{curl}\Omega$, the $d$ matrix, and full $\Omega$ | [`modules/curvature.py`](../../modules/curvature.py) |
| $c_2$ and the BZ integral $\mathrm{d}\theta/\mathrm{d}\beta$ | [`modules/axion.py`](../../modules/axion.py) |
| Y2Ir2O7 calculation choices | [`calculations/run_axion.py`](../../calculations/run_axion.py) |

The position-matrix algebra lives in
[`modules/curvature.py`](../../modules/curvature.py). The construction reproduces pythtb's
`TBModel._berry_curvature_external` ($B^X+B^E$, validated against postw90's
$J_0+J_1$) and adds the mixed $(k,\beta)$ planes for the adiabatic phonon
coordinate.

---

## 0. What we are computing

The magnetoelectric (axion) response is the **second Chern number** of the 4D
parameter space $(k_x,k_y,k_z,\beta)$, where $\beta$ is the frozen-phonon
amplitude (in Å; the code's `d_beta`). Writing the non-Abelian Berry curvature
of the occupied subspace as $\Omega_{\mu\nu}$ (a $4\times4$ antisymmetric array
of $n_{\text{occ}}\times n_{\text{occ}}$ matrices, indices
$\mu,\nu\in\{k_x,k_y,k_z,\beta\}=\{0,1,2,3\}$), the response density is the
Pontryagin (second-Chern) density, in the code's normalization

$$
c_2(k,\beta) = \frac{1}{16\pi}\,\epsilon^{\mu\nu\rho\sigma}\,
\operatorname{Tr}\!\big[\Omega_{\mu\nu}\,\Omega_{\rho\sigma}\big],
$$

and the response is its BZ integral at fixed $\beta$,

$$
\mathrm{d}\theta[\beta]\ =\ \int_{\mathrm{BZ}} \mathrm{d}^3k\; c_2(k,\beta).
$$

In code this contraction is owned by
[`axion.c2_density` and `axion.dtheta`](../../modules/axion.py),

```python
c2_density[:, beta] = einsum("ijkl,ij...mn,kl...nm->...", epsilon, Omega, Omega) / (16*pi)
dtheta[beta]        = sum_k c2_density * d3k
```

The $\epsilon^{\mu\nu\rho\sigma}$ runs over all four axes and double-counts each
plane pair; the $1/16\pi$ prefactor is the normalization chosen to match the
$\theta$ (radians) convention of `axion.py`. Because
$\epsilon^{\mu\nu\rho\sigma}$ needs all four distinct axes, every
surviving term carries exactly one $\beta$ leg — the response is the mixed
$k$–$\beta$ curvature, as it must be.

**The key point of this note:** the curvature splits as

$$
\boxed{\;\Omega_{\mu\nu} = \Omega^{\text{int}}_{\mu\nu} + \Omega^{\text{ext}}_{\mu\nu}\;}
$$

and [`curvature.omega`](../../modules/curvature.py) sums the pieces
(`omega_internal + omega_cross + omega_external`).

---

## 1. Internal (Kubo) curvature — the baseline

In the Hamiltonian (band) gauge, with occupied states $m,n$ and empty states
$c$, velocity matrix elements $v_{\mu} = \langle u | \partial_\mu H | u\rangle$
and energy denominators $D_{mc}=E_m-E_c$, the internal non-Abelian curvature is
the standard interband Kubo sum:

$$
\Omega^{\text{int}}_{\mu\nu,\,mn}
= i\sum_c \frac{v_{\mu,mc}\,v_{\nu,cn} - v_{\nu,mc}\,v_{\mu,cn}}{D_{mc}\,D_{nc}}.
$$

In code ([`curvature.connection` and `omega_internal`](../../modules/curvature.py)) this is formed from the
conduction–occupied "velocity over gap" block

$$
d_{\mu,cn} \equiv \frac{v_{\mu,cn}}{E_n - E_c},
\qquad
V_{\mu,mc} \equiv \frac{v_{\mu,mc}}{E_m - E_c},
$$

as $q_{\mu\nu} = V_\mu\,d_\nu$ (matmul over $c$) and
$\Omega^{\text{int}}_{\mu\nu} = i\,(q_{\mu\nu}-q_{\mu\nu}^\dagger)$. The block
$d_{\mu,cn}$ (`Connection.d`) is **reused** by the external terms, so they cost no
extra diagonalization.

For $\mu$ or $\nu=\beta$, the velocity is the finite-difference of the two
Wannierized Hamiltonians
([`axion.beta_terms`](../../modules/axion.py)):

$$
v_\beta = \partial_\beta H \approx \frac{H_{\text{mode}} - H_{\text{base}}}{d_\beta}.
$$

The internal term alone is what postw90 would call the Wannier-interpolated
curvature with the position operator taken as $i\partial_k$ only — it neglects
the off-diagonal (interband) matrix elements of the true position operator that
the finite spread of the Wannier functions carries. Those are the external
terms.

---

## 2. External position-matrix fields $A(k)$ and $\bar\Omega(k)$

The external contribution comes from the **off-diagonal Wannier position
matrix**

$$
X_a(R)_{nm} = \langle 0n |\, r_a \,| Rm\rangle,\qquad a\in\{x,y,z\},
$$

The production path rebuilds this matrix from the `.chk/.mmn` link overlaps,
with the right-orbital phase applied before the finite-difference neighbor sum,
and then applies the same orbital-dependent real-space mapping as $H(R)$.
`prefix_r.dat` is retained for centre/provenance checks and for an explicit
legacy comparison; it is not a finite-mesh substitute for the `.mmn` route.

Following convention I (orbital centers live in the $\tau$ phases, not in the
on-site position), the code Hermitianizes the resulting matrix and removes the
$R=0$ diagonal:

$$
X^{I}_a(R) = \frac{1}{\deg(R)}\cdot\tfrac12\Big(X_a(R) + X_a(-R)^\dagger\Big),
\qquad X^{I}_a(0)_{nn}=0,
$$

with $\deg(R)$ the Wigner–Seitz degeneracy. Two $k$-space fields are built from
it by Bloch summation ([`curvature.fields`](../../modules/curvature.py)):

### (a) External Berry connection (a one-form, $N\times N$ per axis)

$$
A_a(k)_{nm} = \Big[\sum_R e^{2\pi i\,k\cdot R}\, X^{I}_a(R)_{nm}\Big]\;
e^{\,2\pi i\,k\cdot(\tau_m-\tau_n)},
$$

where the last factor is the convention-I intra-cell phase
([`curvature.fields`](../../modules/curvature.py)). It is then converted to **reduced**
one-form components with the reciprocal-lattice matrix $B$
with the stored one-form transform:
$A_u^{\text{red}} = B_{ua}\,A_a$.

### (b) External curl $\bar\Omega(k)$ (a two-form, dual components $g=(yz,zx,xy)$)

Using the Cartesian bond vectors $b_{nm}(R) = R + \tau_m - \tau_n$
([`curvature.position_terms`](../../modules/curvature.py)), define the $k$-independent curl weights

$$
W_g(R) = i\big(b_\alpha\, X^{I}_\beta - b_\beta\, X^{I}_\alpha\big),
\qquad (\alpha,\beta) = \text{axes of plane } g,
$$

so that

$$
\bar\Omega_g(k)_{nm} = \Big[\sum_R e^{2\pi i\,k\cdot R}\, W_g(R)_{nm}\Big]\;
e^{\,2\pi i\,k\cdot(\tau_m-\tau_n)}.
$$

The dual components of this antisymmetric two-form convert to reduced axes with
the **cofactor** matrix $M = \det(B)\,B^{-\top}$, not with $B$ itself.

Both fields are returned by `curvature.fields(model, pos, k)` in the Wannier gauge;
[`curvature.connection`](../../modules/curvature.py) supplies the eigenvectors $U$ and
rotates them to the Hamiltonian gauge, $a = U^\dagger A\,U$ and
$\bar\omega = U^\dagger \bar\Omega\,U$.

---

## 3. External curvature in the $k$–$k$ planes: $B^X + B^E$

With the rotated connection $a_\mu$ (blocks $a^{co}$ = conduction–occupied,
$a^{oo}$ = occupied–occupied) and the rotated curl $\bar\omega$, the external
curvature for a $k$–$k$ plane $(\mu,\nu)$ is assembled by
[`omega_cross` and `omega_external`](../../modules/curvature.py):

$$
\Omega^{\text{ext}}_{\mu\nu} = B^X_{\mu\nu} + B^E_{\mu\nu}.
$$

**Cross term $B^X$** (conduction–occupied mixing of the Kubo block $d$ with the
external connection $a$):

$$
B^X_{\mu\nu} = i\Big[
\big(d_\mu^{co}\big)^{\!\dagger}(-i\,a_\nu^{co})
+\big(-i\,a_\mu^{co}\big)^{\!\dagger} d_\nu^{co}
\;-\;(\mu\leftrightarrow\nu)\Big].
$$

**Curl + commutator term $B^E$** (purely occupied–occupied):

$$
B^E_{\mu\nu} = \bar\omega_{\mu\nu}^{oo} - i\big[a_\mu^{oo},\,a_\nu^{oo}\big].
$$

The result is antisymmetrized, $\Omega^{\text{ext}}_{\nu\mu} = -\Omega^{\text{ext}}_{\mu\nu}$.
Adding $B^X+B^E$ to $\Omega^{\text{int}}$ reproduces pythtb's
`berry_curvature(include_external=True, non_abelian=True)` and postw90's
$J_0+J_1$ (checked by
[`check_vs_pythtb.ipynb`](../../validation/external_curvature/check_vs_pythtb.ipynb)).

---

## 4. The mixed $k$–$\beta$ planes

The novel part is the adiabatic $\beta$ direction, sampled by **two
Wannierizations** (base and mode). Under the **frozen-Wannier-gauge
approximation** the external connection along $\beta$ vanishes,

$$
a_\beta = 0,
$$

because the trial orbitals are held fixed (the rigidly-shifted base MLWFs) — the
$\beta$ dependence enters only through $H$ and through the $k$-connection, not
through a new position matrix along $\beta$. The mixed curl is then supplied by
the $\beta$-derivative of the $k$-connection
([`axion.beta_terms`](../../modules/axion.py)):

$$
\bar\Omega_{i,\beta} = -\,\partial_\beta A_i,
\qquad
\partial_\beta A_i \approx \frac{A_i^{\text{mode}} - A_i^{\text{base}}}{d_\beta}
\;\equiv\; \mathrm{d}A_i,
$$

computed once per $k$-batch and shared between the two $\beta$ endpoints by
[`axion.dtheta`](../../modules/axion.py). With $a_\beta=0$, the cross and curl
terms collapse to, inside [`omega_cross` and `omega_external`](../../modules/curvature.py):

$$
B^X_{i\beta} = i\Big[\big(-i\,a_i^{co}\big)^{\!\dagger} d_\beta^{co}
- \big(d_\beta^{co}\big)^{\!\dagger}(-i\,a_i^{co})\Big],
\qquad
B^E_{i\beta} = -\,\big(\mathrm{d}A_i\big)^{oo}
$$

(the commutator $[a_i,a_\beta]$ drops since $a_\beta=0$; $\mathrm{d}A_i$ is
rotated to the Hamiltonian gauge before taking the occupied block). Again
$\Omega^{\text{ext}}_{i\beta} = B^X_{i\beta}+B^E_{i\beta}$ and
$\Omega^{\text{ext}}_{\beta i} = -\Omega^{\text{ext}}_{i\beta}$.

---

## 5. Assembly, the two endpoints, and the gauge check

Putting it together, on each $k$-batch and for each $\beta\in\{0,1\}$:

1. Build the 4D velocity set $v_\mu = (v_{k_x},v_{k_y},v_{k_z},v_\beta)$ with
   $v_\beta = (H_{\text{mode}}-H_{\text{base}})/d_\beta$.
2. Compute $\Omega^{\text{int}}_{\mu\nu}$ (Kubo) from $v_\mu$ and the gaps.
3. Add $\Omega^{\text{ext}}_{\mu\nu}$: the $k$–$k$ planes from §3 using the
   fields $(A,\bar\Omega)$ of the model **at this $\beta$** (base fields for
   $\beta{=}0$, mode fields for $\beta{=}1$), and the $k$–$\beta$ planes from §4
   using the shared $\mathrm{d}A$.
4. Contract into the second-Chern density $c_2$ and integrate over $\mathrm{d}^3k$.

This yields **two** estimates of the response, `dtheta[0]` and `dtheta[1]`,
anchored at the base and mode endpoints respectively:

$$
\mathrm{d}\theta[\,\beta\,] = \int_{\mathrm{BZ}}\mathrm{d}^3k\; c_2\big(k,\beta\big).
$$

They are derivatives evaluated at two different $\beta$ values, so they need
not be identical. Their difference contains physical variation of the response
with $\beta$, finite-$d_\beta$ error, and any residual gauge discontinuity. A
shared rigid-shift gauge removes the avoidable gauge-rotation part, but endpoint
agreement by itself is not a pass/fail correctness gate.

---

## 6. Conventions and caveats

- **Convention I** (`cartesian=False`, reduced $k$-axes): orbital centers are
  carried by the $\tau$ phases $e^{2\pi i k\cdot(\tau_m-\tau_n)}$, and the $R=0$
  diagonal of $X^I$ is zeroed. One-forms convert with $B$; two-form dual
  components with the cofactor $M=\det(B)B^{-\top}$.
- **Frozen Wannier gauge** ($a_\beta=0$) is an *approximation*: it assumes the
  trial-orbital gauge does not vary with $\beta$. It is consistent with using
  the rigidly-shifted base MLWFs and projection-only Wannierization, where the
  gauge is anchored to fixed trial functions.
- **Validation**: the $k$–$k$ external terms are checked against pythtb's
  `berry_curvature(include_external=True, non_abelian=True)` (itself validated
  against postw90 $J_0+J_1$) in `check_vs_pythtb.py`. The mixed
  $k$–$\beta$ terms are the new physics and are validated indirectly by the
  endpoint agreement of §5.
- **Requires** compatible `.chk`, `.mmn`, `.nnkp`, `_hr.dat`, and `_r.dat`
  inputs. The `.mmn` supplies the production $A(R)$; `_r.dat` supplies a
  centre/provenance check and the PythTB position container that is replaced
  during loading.
