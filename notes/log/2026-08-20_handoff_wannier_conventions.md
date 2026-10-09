# Handoff: Wannier interpolation conventions, SrTiO3 Z* and the axion path

Status as of 2026-08-20. Written for an agent picking this up cold.

## 1. The one idea that explains almost everything

A finite `N1 x N2 x N3` k-mesh determines the real-space matrices `H(R)` and
`A(R)` **only up to the aliasing class** `[R] = {R + sum_i n_i N_i a_i}`. Every
choice of representative reproduces the ab-initio data *exactly* on the source
mesh and they differ only *between* those points.

So the choice is **not a gauge** - it is a choice of **interpolating function**.
Consequences:

* Using one choice **consistently** for `H(R)` and `A(R)` is necessary for them
  to define one interpolant, but it is not sufficient for symmetry covariance.
  The representative-selection rule must itself commute with the space-group
  action. Under `g`, `tau_s -> tau_s' + L_s`, hence
  `R -> gR + L_t - L_s`; centre-independent `min |R|` need not map a chosen
  representative to the chosen representative of its image. The physical-bond
  criterion `|R + tau_n - tau_m|` is point-group invariant, so its minimiser
  does map to the minimiser.
* **Mixing** two choices is fatal. The Berry-curvature formula assumes `A(k)` is
  the connection *of the eigenstates of* `H(k)`. If `H` is interpolated one way
  and `A` another, that relation holds only on the source mesh and fails
  everywhere else. The result is not the curvature of any Hamiltonian, so it
  respects no symmetry.
* It is **invisible at the ab-initio k-points**, which is why it survived so
  long. It only appears on the dense interpolation mesh used for BZ integrals.

Wannier90's `_hr.dat` / `_r.dat` minimise `|R|` (centre-independent).
WannierBerri (and wannier90 with `use_ws_distance = .true.`) minimise the
physical separation `|R + tau_n - tau_m|`, which differs per orbital pair.

**These `.win` files have `use_ws_distance = false`**, so the shipped
`wsvec.dat` encodes the identity (all `ndeg = 1`, shift `0,0,0`). Do not waste
time applying it - it is a no-op.

## 2. Environment

```
python : /Users/treycole/miniforge3/envs/axion/bin/python   (has pythtb 2.0.0
         editable from ~/Repos/pythtb, wannierberri 1.8.0 editable from
         ~/Repos/wannier-berri, scipy, spglib 2.7.0)
repos  : ~/Repos/axion-phonon, ~/Repos/pythtb, ~/Repos/wannier-berri,
         ~/Repos/wannier90
```
`/opt/homebrew/bin/python3` does **not** have wannierberri. Always use the axion env.

## 3. Bugs found and fixed (all verified)

### 3a. PythTB embedding was a silent no-op  -- FIXED
`TBModel.lattice` returns `copy.copy(self._lattice)` (`pythtb/tbmodel.py:233`),
so `model.lattice.orb_vecs = tau` mutates a throwaway. `ExternalTerms` received
`tau` while `H` and `v` kept the Wannier centres -> a mismatch **linear in tau**.
That is the true origin of the "defect 1 tau artefact" recorded in earlier notes;
there is no Fourier-truncation tau artefact.

Fix: build with `W90.model(orb_vecs=...)` plus a post-build assertion. After the
fix everything is **exactly tau-invariant** (tau = 0, 1x, 2x, 3x, 4x centres and
a random non-collinear tau all give identical answers). The int/cross/ext *split*
is not invariant - it shifts by +-0.3281 between tau = 0 and centres while the
total is fixed - so always state tau when quoting a decomposition.

Applied in `validation/observable/born_effective_charge.py:~835` and
`calculations/Y2Ir2O7/phonon/axion_response_vs_nk.py:~330` and
`check_vs_pythtb.py:~111`. In the axion script the embedding is now set
**unconditionally**, not only under `include_external`, because
`v_beta = (h_mode - h_base)/d_beta` needs the shared tau even with external
terms off.

**`calculations/MnBi2Te4/phonon/axion_response_nk_sweep.py:262-273` still has the
bug** - it builds both models with default `orb_vecs` and never attempts a shared
embedding. It also uses **no external terms at all**. The paper's headline
`dQ theta ~ 10 rad/A` for MBT comes from that script and needs rerunning. *(2026-09-24: the script is
deleted; the paper's numbers are in flux and will be recomputed on the new pipeline.)*

### 3b. H(R) and A(R) were in different WS conventions  -- FIXED
We took `A(R)` from WannierBerri but left `H(R)` as pythtb's `_hr.dat`. Located
by swapping one ingredient at a time on the displaced trial_03 line:

| what | max abs Tr Omega_xy | corr w/ WB |
|---|---|---|
| our fields + our d_cv | 4.33e-03 | 0.929 |
| WB formula fed with OUR fields | 4.33e-03 (identical) | assembly is FINE |
| our A(k), Omega(k) vs WB's | rel 0.00% (3e-6 / 3e-9) | fields are FINE |
| pythtb H(R) vs WB H(R) | rel 0.0105%, max 3.3e-4 eV | <- the difference |
| ALL-WB inputs + WB formula | 6.90094e-03 | 1.000000, max abs d 1.8e-10 |

Fix, **entirely in-house**: `wannier/wannier_io.py` now has
`ws_reassign()` and `apply_ws_reassignment(w90, mp_grid)`. Call it on the `W90`
object **before** `.model()`; it rewrites `ham_r` with `deg = 1`, which keeps a
later `pos_r_from_mmn` self-consistent. Driver flags `--ws-reassign` /
`--no-ws-reassign` (default on), outputs tagged `_nows` when off.
Prune zero images or the R-set balloons 729 -> 15625 and every Bloch sum slows ~20x.

Verified against WannierBerri: rel **0.0004%**, max abs 7.06e-7 (that residual is
`_hr.dat`'s 6-decimal printing).

## 4. Results after the fixes

**Every symmetry check collapses to numerical zero for all three sets.**

| set | H(R) conv | PT zeros | mirror kx=1/2 | forced zeros |
|---|---|---|---|---|
| a | min abs R | 6.13e-02 | 4.63e-02 | 1.07e-02 |
| **a** | **pair-dep** | **9.94e-09** | **5.24e-08** | **5.43e-08** |
| b | either | 1.58e-09 | 7.99e-10 | 7.69e-10 |
| c | min abs R | 3.23e-03 | 1.85e-03 | 8.50e-04 |
| **c** | **pair-dep** | **1.03e-08** | **1.13e-07** | **5.41e-08** |

So the symmetry violations were **our convention bug, not the Wannier sets**.
Set (c)'s physical curvature also grows 4.33e-3 -> 6.81e-3, matching wberri's
6.90e-3; the old value was 37% too small.

**Z\* barely moves** (measured, not assumed): trial_02 7.37013 -> 7.37013
identical, trial_03 7.30088 -> 7.30099, trial_01 7.42721 -> 7.43250. DFPT is
7.28829. All Born-charge conclusions stand: (a) +1.98%, (b) +1.12%, (c) +0.17%.

Why Z\* is insensitive: BZ integration over a symmetric mesh annihilates exactly
the error component that transforms incorrectly under the point group - which is
what the symmetry residual measures. **Residual quality and integrated accuracy
are different error components.** A clean residual does not certify Z\*, and vice
versa.

**David's question (invariance w.r.t. Wannier content) still fails, and now
cleanly.** Sets (a) and (b) span the *same* occupied manifold (both n_occ = 30,
same 10 excluded bands), so Tr Omega_xy must be identical. After the fix:
max abs a = 2.396e-02 vs max abs b = 3.365e-03, **ratio 7.1x, corr -0.95**. Both
are symmetry-clean, so this is genuine interpolation quality; set (a)
(unsupported Sr d - `Sr-sp_r.upf` has no d channel - spread 5.69 A^2) is the
outlier.

The generic-line figure was regenerated on 2026-08-20 at
`calculations/SrTiO3/berry_curvature_generic_line.png`. All five marked forced
zeros are now numerical noise: max residuals are 5.43e-08 for (a), 7.69e-10
for (b), and 5.41e-08 A^2 for (c). The physical maxima are 2.395779e-02,
3.364779e-03, and 6.812123e-03 A^2 respectively.

A later strict comparison kept both pipelines internally consistent and used
identical 729-vector `H/A` R sets: (i) centre-aware `.mmn` plus pair-dependent
`H`, and (ii) direct `_r.dat + _hr.dat` with no reassignment. Direct files still
fail the informative symmetry tests: on the undisplaced generic PT-null path,
sets (a,c) have maxima 2.208e-01 and 4.057e-02 A^2, versus 5.53e-09 and 8.57e-09
for the centre-aware path. Set (b) is zero in both for the structural reason
below and is not a validation. Full Born-charge convergence, mirror-plane and
generic-path data are in
`calculations/SrTiO3/convention_analysis/analysis_report.md`.

## 5. postw90 disagreement (open, do not cite)

A postw90 `kpath_task = curv` run exists at
`.../185116bc889/postw90/185136bc889/` on the identical line. Column 4 is
**minus** Omega_z in Ang^2. Its inputs are correct (num_wann 48, right `.chk`
timestamp 17Aug2026 16:31:57, right `.mmn`, `fermi_energy = 11.0` inside a
2.95 eV gap with exactly 38 states below at every k). On the undisplaced cell
where P+T force Tr Omega_xy = 0:

```
postw90                1.08e-01     <- violates
ours, _r.dat           4.06e-02
ours, .mmn             9.04e-04
WannierBerri own calc  1.08e-08     <- satisfies
```

postw90 is **uncorrelated** with ours (+0.055 / -0.065) and matches no single
piece of our decomposition. Cause unidentified. Two independent implementations
satisfy the constraint and postw90 does not.

## 6. In-house A(R): COMPLETE (2026-08-20)

Goal: drop the WannierBerri run-time dependency. `H(R)` is **done** (section 3b).
`A(R)` genuinely needs the `.mmn` - it cannot be recovered from `_r.dat`, and the
WS reassignment does **not** fix it (tested: 6.72% -> 6.72%). A convention
hypothesis was also tested and rejected: the best `_r.dat` vs wberri agreement is
same-convention-on-both (1.12% in A(k)), so it is a genuine construction
difference, not a phase convention.

`wannier/wannier_io.py` now implements the
whole route in process. WannierBerri is no longer imported or launched at run
time.

| piece | status |
|---|---|
| `read_chk` | **written, validated EXACT** - v_matrix, centres, kpt_red, real_lattice all max abs d = 0.000e+00 vs WannierBerri |
| `bshell_weights` | **written, validated EXACT** - 12.36209079, matches wberri, completeness 0.00e+00 |
| `read_mmn` | **done** - keeps the on-disk `(ik, ib)` order, uses the Wannier90/WannierBerri matrix transpose, and does not sort by neighbour index (which changes order with k) |
| A(R) assembly | **done, validated to machine precision** - 729/729 R-vectors; max abs d `1.11e-16` A and relative d `1.44e-16` vs `System_w90(...).get_R_mat('AA')` on base trial 03 |

Spec for the assembly, matching wberri's defaults
(`transl_inv_JM=False`, `transl_inv_MV=False`) - read
`chk.get_AABB_q_ib` and `system_w90.sum_matrix_b`:

```
M^W(k,b) = V^dag(k) M^ab(k,b) V(k+b)          # V = v_matrix from .chk
A_a(k)_mn = i * sum_b w_b * b_a * M^W_mn(k,b) # NO delta_mn (it cancels: sum_b w_b b = 0)
                                              # multiply each link by exp(i b.tau_n)
                                              # plain formula on the diagonal too
                                              #   (MV log only if transl_inv_MV)
A(R) = FFT to real space, then the SAME pair-dependent WS reassignment as H(R)
```

The formatted `recip_lattice` block from `.nnkp` is load-bearing for bit-level
agreement: WannierBerri builds its b-vectors from that explicitly rounded block,
not from the higher-precision `.chk` lattice and not by inverting the separately
rounded `.nnkp` real lattice. The latter choices leave 1e-9-scale differences in
`A(R)`.

Cold parse/assembly is 34.3 s for the 1.1 GB trial-03 `.mmn`; a provenance-checked
cache hit is 0.25 s. The real PythTB bridge gives 729 R-vectors, exact pair
Hermiticity, b-shell completeness error 0, and preserves the existing intentional
replacement of the R=0 diagonal by the `.chk` Wannier centres. Regression tests
are in `validation/unit/test_wannier_io.py`.

## 7. Traps that cost time - do not repeat

* **Do not infer by elimination from a weak test.** "It's Fourier truncation" was
  concluded three times and was wrong three times. The block-Frobenius-norm
  covariance test is *necessary but not sufficient* - blind to a wrong rotation
  inside a block. Use the **singular-value spectrum** per centre-block (also
  D-independent, much stronger): it exposed H(R) at 1.4e-05 where the norm test
  said 2.3e-10.
* **Point-group audits need the pair-dependent lattice shift.** Under g the WF at
  tau_s maps to `g.tau_s = tau_{s'} + L_s`, so the R index shifts by `L_t - L_s`.
  Assuming `R -> gR` gives garbage (~2.0 violations).
* `min_hopping_norm = 0` improves E(k)/v(k) covariance 6x but changes **nothing**
  downstream (curvature, Z\*). Use it, expect nothing.
* Set (b) passes every symmetry test **structurally** (occupied-only -> no
  conduction block -> int and cross vanish identically, and the full-space trace
  depends only on the manifold). It cannot fail. Never cite it as validation.
* `--self-test` is a smoke test only: `model_mode = model_base` makes every
  finite difference *identically* zero, so it returns bitwise 0 and tests almost
  nothing. The Wannier-centre cross-check (4.3e-06 for set b) is the real one.
* `System_w90(..., wannier_centers_from_chk=False)` is **broken in wberri 1.8** -
  it calls `chk.get_AA_q(mmn, ...)` but the signature is
  `get_AA_q(bkvec, mmn, kptirr, weights_k)`.
* wberri's `symmetrize` is hard-disabled in 1.8 (`TODO FIXME`), so its results are
  *not* symmetrised - the 1e-8 it achieves is earned.
* Scripts that fork a multiprocessing pool need `if __name__ == "__main__":`.
* Piping a long run through `| tail` buffers everything; log to a file instead.

## 8. Open items, in priority order

1. Fix and rerun **MnBi2Te4** (section 3a) - the paper's headline number.
2. Y2Ir2O7 dtheta/dQ: embedding fixed and rerun to 12^3. Endpoint pair
   dtheta[0]/dtheta[1] at nk=12: -0.001433 / -0.001940. **This is NOT a failed
   gate** - they are dtheta/dbeta at two different beta, so they agree only if
   theta(beta) is linear. The half-difference is more nk-stable (7.7%) than the
   mean (15.5%), so it looks like genuine curvature d2theta/dbeta2 ~ -0.0505.
   Settle it with a second amplitude (only `Q_0.01A` exists).
3. `Lambda(R) = 0` (frozen gauge) is still untested and is the last unverified
   physical assumption. Invisible on occupied-only sets (it enters as
   `Tr d_i A_beta`, a total derivative that integrates to zero), so it must be
   tested on set (c).
4. Latent: `ExternalTerms` converts A with `_one_form = B` but the curl with
   `_two_form = det(B) inv(B).T`. For cubic both equal `(2pi/a)^2` so SrTiO3 and
   Y2Ir2O7 are unaffected; **MnBi2Te4 is rhombohedral** and they differ. Check
   before trusting any non-cubic result.
