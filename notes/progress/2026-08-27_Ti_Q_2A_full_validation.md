# Ti_Q_2A Wannier-curvature validation

Date: 2026-08-27

## Decision

The implementation checks pass.  The trial-03/trial-04 curvature difference
is not caused by a different occupied DFT projector, a corrupt `.mmn`, the
inverse Fourier transform, the orbital-aware R assignment, analytic
derivatives, the internal/cross/external formulas, or a broken crystal
symmetry.

It is a finite-source-mesh interpolation ambiguity.  The two Wannierizations
give the same occupied subspace and the same gauge-invariant Wilson loops on
the 8x8x8 first-principles mesh, but they define different derivatives of that
subspace between and even *at* those sample points.  Curvature is derivative
information; equality of a projector at isolated nodes does not fix it.

This means:

* The compact, provenance-checked `A(R)` cache can replace the large `.mmn`
  **as stored input to our present interpolation scheme**.
* Ordinary `_r.dat` must not replace it.
* `postw90` with `transl_inv_full=true` and `use_ws_distance=true` is
  implemented consistently, but it is a second finite-difference
  discretization, not an exact reference.
* Production axion accuracy is not yet certified from the 8x8x8 data.  A
  denser *ab-initio* overlap mesh and the real beta-connection gate remain.

## Checks and results

### 1. Raw first-principles occupied geometry

For each checkpoint, let `C_a(k)` be the 100-by-38 occupied frame in the common
DFT-band space and `P_a(k)=C_a(k) C_a(k)^dagger`.

| Check | Result | Status |
|---|---:|---|
| max Frobenius `P_03-P_04`, all 512 k points | 3.463e-14 | PASS |
| max projector entry difference | 7.661e-15 | PASS |
| minimum principal-overlap singular value | 0.9999999999999963 | PASS |
| gauge-aligned raw occupied-link relative error | 4.692e-15 max | PASS |
| gauge-aligned polar-link relative error | 4.579e-15 max | PASS |
| max occupied-link singular-value difference | 7.883e-15 | PASS |
| max trial difference in polar Wilson flux | 1.762e-14 rad | PASS |

All uncompressed and gzipped MMN copies have the same decompressed SHA-256:

`017bb5f49f946cfd39a764f7f2700f971ad30206895c5848b1686e84a85e7646`

Thus the raw overlap input and the discrete occupied quantum geometry are
identical.  The approximately 1e-13 reverse-link residual is common to both
trials and is the expected consequence of the 12-digit MMN text format.

Artifacts:

* `validation/inputs/check_raw_mmn_wilson.py`
* `calculations/SrTiO3/Ti_Q_2A/raw_mmn_wilson_comparison/`

### 2. MMN to A(R), H(R), and orbital-aware R mapping

All 14 required interpolation-mechanics checks pass for both trials.

| Check | trial 03 | trial 04 | Status |
|---|---:|---:|---|
| inverse FFT `A(k)->A(R)->A(k)`, relative | 2.081e-16 | 1.902e-16 | PASS |
| raw-to-pair-R equality at source nodes, relative | 2.081e-16 | 1.902e-16 | PASS |
| regenerated versus cached A(R) | exactly zero | exactly zero | PASS |
| H(R) versus checkpoint H(k), max | 5.271e-5 eV | 5.332e-5 eV | PASS |
| off-grid H/A/curl Hermiticity, relative | <=9.52e-16 | <=1.18e-15 | PASS |
| reciprocal covariance, relative | <=1.44e-15 | <=1.82e-15 | PASS |
| analytic dH versus central difference, relative | 2.19e-10 | 3.83e-10 | PASS |
| analytic curl versus central difference, relative | 6.56e-10 | 8.98e-10 | PASS |

The derivative errors decrease as `h^2`, as required for the central
difference.  The H/checkpoint residual is below 7.5e-5 eV and is consistent
with the six-decimal `_hr.dat` text precision.

Two expected diagnostics are important:

1. The raw WYSV Eq. (44) finite-difference `A(k)` is not exactly Hermitian:
   its relative anti-Hermitian residual is 2.46% in trial 03 and 6.96% in
   trial 04.  Hermitian projection changes it by 1.23% and 3.48%.  Wannier90
   explicitly documents this in `get_oper.F90` and takes the same Hermitian
   part.
2. Pair-dependent R reassignment preserves the 8x8x8 nodes but changes the
   off-grid continuation of A by 27.6% and 37.0%.  This is the alias freedom
   that a source-mesh convergence test must resolve.  It is not an arithmetic
   failure.

Artifacts:

* `validation/inputs/check_interpolation_inputs.py`
* `calculations/SrTiO3/Ti_Q_2A/interpolation_validation/`

### 3. Curvature formulas and independent postw90 reference

The analytic covariance suite now covers:

* curvature moving exactly between internal and external representations;
* a moving frame with internal, cross, and external terms all nonzero;
* a non-Abelian occupied commutator;
* the corresponding mixed `(k_i,beta)` cases, including nonzero internal,
  cross, and external terms.

All 24 repository tests pass.  The analytic total-curvature agreements are at
approximately 1e-12 or better.

A fresh local build of Wannier90 4.0.1 was also used as an independent
reference with

```text
transl_inv_full = true
use_ws_distance = true
write_aa_r = true
```

Our evaluator, after reading postw90's full A(R), agrees with postw90 on the
same generic line as follows:

| | all-component relative L2 | max absolute |
|---|---:|---:|
| trial 03 | 5.272e-5 | 7.275e-5 A^2 |
| trial 04 | 6.959e-8 | 4.242e-7 A^2 |

The larger trial-03 residual is consistent with our use of the printed
six-decimal H(R), while postw90 reconstructs H from its checkpoint/eigenvalue
data.  Trial 04 has no internal term because all 38 Wannier bands are occupied,
so it is essentially insensitive to that H rounding.

Postw90 itself finds a 1.518% all-component trial difference on the generic
line.  Our evaluator supplied with the same postw90 A(R) finds 1.520%.  The
trial dependence therefore exists independently of our curvature evaluator.

Artifacts:

* `validation/unit/test_berry_curvature_covariance.py`
* `validation/external_curvature/check_vs_postw90_line.py`
* `calculations/SrTiO3/Ti_Q_2A/postw90_reference_comparison/`

### 4. Direct source-node curvature and symmetry

At all 512 original mesh nodes, the two continuous interpolants give:

| Component | relative L2 trial difference | max absolute difference |
|---|---:|---:|
| yz | 12.708% | 0.4878 A^2 |
| zx | 12.708% | 0.4878 A^2 |
| xy | 38.370% | 0.0180 A^2 |
| all | 12.717% | 0.4878 A^2 |

The small xy component has a large relative percentage because its scale is
small.  The occupied eigenvalues agree within 7.29e-5 eV.

Each trial separately obeys every surviving symmetry at generic off-grid
points:

| Symmetry | trial 03 relative residual | trial 04 relative residual |
|---|---:|---:|
| time reversal | 2.67e-10 | 2.35e-10 |
| C4z | 8.44e-9 | 1.03e-9 |
| Mx | 8.04e-9 | 1.10e-9 |

The discrepancy is therefore symmetry-preserving.

Artifacts:

* `validation/symmetry/check_trial_curvature_geometry.py`
* `calculations/SrTiO3/Ti_Q_2A/trial_curvature_geometry/`

### 5. Integrated continuum curvature versus the raw Wilson flux

For every original plaquette, 4x4 Gauss-Legendre quadrature was used to
evaluate

`I_mn = integral_P Tr Omega_mn dk_m dk_n`.

The finite-link Wilson relation is `arg det W_mn = -I_mn` plus finite-link
corrections.  The order-4 results are:

| Metric | xy | yz | zx |
|---|---:|---:|---:|
| trial 03 continuum versus raw flux, relative L2 | 6.584% | 1.735% | 1.735% |
| trial 04 continuum versus raw flux, relative L2 | 40.804% | 9.289% | 9.289% |
| trial 03 versus trial 04 continuum integrals | 39.329% | 10.715% | 10.715% |

The order-3 to order-4 relative change in the dominant yz/zx integrals is
2.99e-5 for trial 03 and 1.05e-5 for trial 04.  Hence the approximately 10.7%
trial difference is converged with respect to the numerical quadrature.

Every raw and interpolated two-torus slice has zero total flux to round-off
(raw <=1.45e-14 rad, interpolants <=1.06e-15), so the Chern/symmetry closure
check also passes.

Artifacts:

* `validation/external_curvature/check_vs_wilson_flux.py`
* `calculations/SrTiO3/Ti_Q_2A/plaquette_flux_interpolant_comparison/`

### 6. Mixed curvature and Born effective charge

Fresh mixed-response runs and the stored Sr endpoint give:

| Check | Result | Status |
|---|---:|---|
| zero displacement, max electronic response | exactly 0 | PASS |
| tau=centres versus tau=0 total, max difference | 7.13e-15 | PASS |
| redistribution among internal/cross/external under tau change | 0.3281 | expected |
| Ti `Z*_zz`, nk=16 | 7.30099 vs DFPT 7.28829 (+0.174%) | PASS |
| Sr `Z*_zz`, nk=16 | 2.55120 vs DFPT 2.55488 (-0.144%) | PASS |
| Ti nk=12 to 16 change | 1.67e-4 | PASS |
| Sr nk=12 to 16 change | 1.26e-5 | PASS |
| Ti/Sr forbidden transverse components | <=2.75e-10 | PASS |
| occupied-only curvature versus centre formula | <=4.30e-6 | PASS |
| occupied-only MMN versus postw90 Z* | <=2.02e-8 | PASS |
| calibrated deep-band null | 0.381 < 0.5 | PASS |
| raw QE acoustic-sum residual | 0.00164 | PASS |
| QE ASR-applied residual | 1.0e-5 | PASS |

For the 48-WF calculation, postw90-minus-MMN changes `Z*_Ti,zz` by 0.03196;
for Sr it changes it by 0.000349.  This is a discretization-sensitivity
diagnostic, not a failed identity.  In the occupied-only model the difference
cancels to about 1e-8, proving that it lives in the off-diagonal/empty-sector
position information.

Artifact:

* `validation/observable/audit_mixed_born_response.py`
* `calculations/SrTiO3/mixed_born_validation/`

## Why equal projectors do not imply equal Omega here

The raw test proves `P_03(k_l)=P_04(k_l)` only at the discrete points `k_l`.
Curvature depends on derivatives,

`Omega_ij(k) = -i Tr P(k) [partial_i P(k), partial_j P(k)]`,

or equivalently on the internal, cross, and external derivative terms in the
Wannier representation.  Infinitely many smooth functions pass through the
same set of discrete projectors while having different derivatives there.

Trial 03 interpolates a 48-dimensional space and obtains part of the answer
from virtual transitions to ten in-model empty states.  Trial 04 interpolates
only the 38 occupied states; its internal and cross terms vanish and all
curvature is carried by the external position matrix.  Gauge covariance says
the totals would agree if the two representations supplied exact, mutually
consistent H, A, and derivatives.  A finite 8x8x8 overlap stencil does not
guarantee that consistency.  The raw Wilson loop constrains an integrated
lattice holonomy, while Eq. (44) followed by Fourier interpolation chooses one
particular local continuum connection.

## Remaining production gates

These cannot be completed from the files currently in the repository.

1. **Ab-initio overlap-mesh convergence.**  Repeat both trial sets on at least
   10x10x10 and preferably 12x12x12 NSCF/MMN meshes.  Repeat the raw-link,
   source-node, line, and plaquette tests.  A reasonable release criterion is
   below 1% trial dependence in the dominant curvature components and stability
   of the final axion observable at the chosen scientific tolerance.
2. **Displacement truncation.**  Generate a +0.005 A Ti endpoint in the same
   fixed-projection gauge as +0.010 A.  Compare both finite differences and use
   `2 response(0.005)-response(0.010)` for the O(delta-u^2) extrapolate.
3. **Real beta connection for axion response.**  No `Lambda(R)` is present.
   The Born charge cannot validate it because the Brillouin-zone integral of
   `Tr partial_k A_beta` vanishes and commutator traces vanish.  A non-Abelian
   axion integrand does not enjoy that cancellation.  Obtain the projection
   connection from Sternheimer/wavefunction-response data or construct and
   validate an explicit endpoint parallel transport before calling the axion
   result exact.
4. **Full interpolated acoustic sum rule (recommended).**  Add an O-displaced
   endpoint.  The available Sr and Ti endpoints pass their individual DFPT
   benchmarks, but an all-sublattice Wannier ASR cannot be formed without O.

Until gate 1 passes, use trial spread as an interpolation error bar.  Until
gate 3 passes, frozen `A_beta=0` results are controlled approximations rather
than the exact projection-connection response.

## Reproduction commands

From the repository root, with the `axion` Python environment:

```bash
python validation/inputs/check_raw_mmn_wilson.py
python validation/inputs/check_interpolation_inputs.py
python validation/external_curvature/check_vs_wilson_flux.py
python validation/symmetry/check_trial_curvature_geometry.py
python validation/observable/audit_mixed_born_response.py
python -m unittest discover -s tests -v
```

`check_vs_postw90_line.py` additionally takes the two directories containing
fresh `SrTiO3-curv.dat` and `SrTiO3_r_full.dat` exports.
