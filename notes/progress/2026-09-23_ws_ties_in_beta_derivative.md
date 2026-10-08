# Wigner–Seitz ties in the β derivative, and `shared_ws_centres` (2026-09-23)

Follows [2026-09-20_jaemo_ndegen_yio_mode1.md](2026-09-20_jaemo_ndegen_yio_mode1.md) ("base and mode do not
share a Wigner-Seitz rule"). Equation numbers are those of
[berry_curvature_derivation.md](../berry_curvature_derivation.md). Rendered version with the worked example:
https://claude.ai/artifact/9FEMzkP6fMjTtBwDXX3hup

## 1. What X subtracts

`X_{μ,st}(R) = <0s| r_μ − τ_{sμ} |Rt>` (8). Differentiating (6) naturally gives the ket's position `R + τ_t`;
swapping it for the bra's `τ_s` changes the element by a constant times `<0s|Rt> = δ_{R0} δ_{st}`, where
`b_ss(0) = 0`. So `X` is the plain position matrix except `X_ss(0) = r̄_s − τ_s`. With both models on the
base `τ` (`curvature.py:233`), `∂_β X_ss(0) = ∂_β r̄_s` is the centre-motion channel. `R = 0, s = t` has the
unique zero-length bond, so it is never part of a tie.

## 2. Where the image choice enters, and A3′

A mesh fixes only class sums `O[R] = Σ_L O(R + L)`; the interpolation hands each class to images with weights
`w_st(R)` (ours: `min |R + τ_t − τ_s|`, ties split equally). Plain Bloch sums (`H(k)`, `A(k)`) are
image-independent **on the ab-initio mesh only**; `b`-weighted objects (`∂_k H`, the curl (9)) depend on the
images everywhere.

§8.2 differentiates (8) at fixed `b_st(R)`. That needs `τ` fixed (A3) **and** the image weights fixed:

> **A3′.** `τ_s` and the image weights `w_st(R)` of every aliasing class are the same for both structures of
> the finite difference.

The in-house fold used each run's own centres (`wannier_io.py` → `apply_orbital_dependent_R_mapping`, and
`checkpoint["wannier_centres"]` in `_connection_to_A_R`), while the Bloch phase uses the shared base `τ`. A tie
split (½, ½) in the base and (1, 0) in the mode adds `(h / 2Δβ)(e^{ik·b1} − e^{ik·b2})` to `∂_β H(k)` (and the
same to `∂_β A(k)`): zero on the ab-initio mesh, `O(1/Δβ)` off it. In the frozen gauge the β leg of (39) uses
only plain Bloch sums (`d_β` from `∂_β H`, `Ω_{βl} = ∂_β A_l`), so the tie term vanishes exactly at ab-initio
k-points. `dtheta` samples `k = j/nk` from Γ and the YIO `mp_grid` is 8, so **nk = 4 and 8 are tie-free; 6,
10 and 12 are not**. Existing `.mmn` sweep: −0.0056/−0.0057 (nk 8) vs −0.0062/−0.0067 (nk 10, 12), which is
suggestive but not conclusive.

## 3. `load_wannier_pair(shared_ws_centres=True)`

Folds `H(R)` and `A(R)` of **both** structures on the base centres (one array for all four folds; it is also
the embedding `τ`). `A(R)` is re-folded from the cache by `remap_A_R_to_centres`; this is exact because in the
`.mmn` route the images of a class are copies (the link phase is R-independent). Then `ws_pair_consistency`
must report zero changed entries (asserted). Unit test: `validation/unit/test_shared_ws_centres.py`.
Driver: `AXION_POSITION_SOURCE=mmn AXION_SHARED_WS=1 python run_axion.py`.

### When it is valid

The image choice is a choice of interpolant, not physics. The mesh-determined content (class sums) is identical
under any rule, so a fixed rule loses no physical signal; it only removes the `1/Δβ` change-of-interpolant term.
All rules converge together as the ab-initio mesh grows (tied bonds sit at the Wigner–Seitz boundary, length
~ N a / 2, where hoppings decay exponentially). It needs:

1. **One lattice.** Same cell ⇒ same aliasing classes. Checked (`atol 1e-8`).
2. **One-to-one Wannier functions.** WF `s` of the mode is the continuation of WF `s` of the base (same
   projections, order, rigid-shift gauge). Required by any finite difference; only the count is checked.
3. **Small centre motion.** By the triangle inequality the base rule picks a mode bond at most
   `|Δr̄_s| + |Δr̄_t|` longer than the mode's own shortest image. YIO mode 1: `max |Δr̄| ≈ 0.015 Å` against
   tie-bond lengths of order half the 8-cell period. Negligible. The shift is logged.
4. **A `.mmn`-family `A(R)`.** Refused for `rfull` and `rdat_ndegen_applied`: under `transl_inv_full`, images
   differ by `(−1)^{n_i}`, so class sums cannot be redistributed; that rule must be fixed inside Wannier90
   (`plot.F90:198`).

Why the **base** centres: the base is the high-symmetry reference, its ties are exact, and equal splitting keeps
its interpolant covariant under m-3m′. The mode's group (−4′3m′) is a subgroup, so the same rule is covariant
for the mode too. The base-end estimate uses exactly the base's natural interpolant.

### When it is not (or what it does not fix)

* Different cells for base and mode (strain, or a zone-boundary mode computed in a cell other than the base's):
  the classes differ and there is no shared rule. Refused.
* Wannier functions that change identity between the runs (a centre jumping, disentanglement picking a different
  state): the base rule may then put a mode hopping on a genuinely long bond. The finite difference is invalid
  in that case anyway; the logged centre shift and `Λ ≠ 0` are the diagnostics.
* It does **not** settle the `.mmn` vs `rfull` interpolant-family gap (factor ~2 in dθ). It only makes each
  family consistent across β.
* The mode-end estimate uses a rule that is not the mode's own shortest-bond rule. That is an `O(Δβ)` effect,
  the same order as the endpoint spread itself.

## Result (2026-09-23): YIO mode 1, `.mmn`, frozen gauge, unpatched QE

`AXION_POSITION_SOURCE=mmn AXION_SHARED_WS=1 python run_axion.py`: 27 min, peak RSS 8.6 GB. Output
`.../181565bc889/nk_sweep_mmn_shared_ws/`. `ws_pair_consistency`: 0 / 6,436,812 entries changed (was 92,706,
5.5%). dθ in rad/Å (base end / mode end); "own" = `nk_sweep_ext_ws_mmn` (Aug 24):

| nk | on ab-initio mesh? | own centres | shared base centres | base-end change |
|---|---|---|---|---|
| 4 | yes | +0.003768 / +0.000572 | +0.003786 / +0.002231 | +0.5% |
| 6 | no | −0.006597 / −0.007130 | −0.005831 / −0.006074 | −12% |
| 8 | yes | −0.005554 / −0.005669 | −0.005552 / −0.005640 | 0.0% |
| 10 | no | −0.006237 / −0.006756 | −0.005293 / −0.005445 | −15% |
| 12 | no | −0.006224 / −0.006746 | −0.005281 / −0.005429 | −15% |

* The base end is unchanged at nk = 4 and 8 and moves by 12–15% at nk = 6, 10, 12, as §2 predicts.
* The mode end can move on-mesh too (nk = 4). Correction to §2: the β-derivatives are image-independent on-mesh,
  but the mode-end estimate also uses the mode's own `∂_k H` and `Ω_kk`, which are `b`-weighted and so depend on
  the mode's images even at mesh points. The base end has no such dependence (its rule is the same in both runs).
* Converged: nk 10 → 12 changes by 0.2%. The endpoint spread at nk = 12 drops from 8.4% to 2.8%, and the
  nk = 8 → 10 jump (−0.0056 → −0.0062) is gone (−0.0056 → −0.0053).
* **The tie inconsistency was a 15–18% error in the `.mmn` dθ** (−0.0062/−0.0067 → −0.0053/−0.0054 at nk = 12).
* Jae-Mo (`rdat_ndegen_applied`, own-centre rule, cannot be re-folded) is −0.00377/−0.00341. It is now
  0.71×/0.63× the `.mmn` value, up from 0.61×/0.51×. The remaining gap mixes the interpolant-family difference
  with Jae-Mo's own, still uncorrected, tie inconsistency.

Still on the unpatched QE DFT (see [2026-09-21_qe_new_ns_nc_bug.md](2026-09-21_qe_new_ns_nc_bug.md)).

## Jae-Mo's formula, re-implemented from the `.chk` (2026-09-23)

`position_matrix_transl_inv_full` (`wannier/position_matrix.py`) transcribes Wannier90's `hamiltonian_get_rmn`
(`transl_inv_full` + `write_ndegen_applied`) from the checkpoint's `m_matrix`: per b, k-space phase
`exp(i b·(r_i+r_j)/2)` with the run's own centres → class sums → fold → `i w_b b exp(−i b·R/2)` at the final
image → sum b; `R = 0` diagonal = centres. The fold can use shared centres, which a file cannot. Source
`chk_transl_inv_full` in `load_wannier_pair` / `run_axion.py`. No `.mmn` needed.

**1:1 against Jae-Mo's `_r.dat`** (`validation/inputs/compare_chk_transl_inv_full.py`, results in
`validation/results/chk_transl_inv_full/`): SrTiO3 trial04, YIO base and YIO mode all **100.000% of components
equal after rounding to six decimals** (33,980,672 rows each for YIO), same R set, max difference 7.07e-7
(= half a print step in both parts), Hermiticity 7e-15.

**1:1 in dθ, and a print-precision finding.** Own-centre dθ at nk = 6:

| input | base end | mode end |
|---|---|---|
| Jae-Mo's files (Sep 20, and rerun through today's pipeline) | −0.0045095 | −0.0045006 |
| ours, full precision | −0.0046468 | −0.0046070 |
| ours, A(R) rounded to 6 decimals | −0.0045106 | −0.0045015 |

Rounding reproduces the files to 2e-4. **`_r.dat`'s `F12.6` print costs ~3% of dθ**: rounding errors (≤ 5e-7 Å)
are independent in base and mode and are divided by Δβ = 0.01. Full precision is the better number. The stock
`_hr.dat` is also printed at 6 decimals (eV) and feeds every route; building H(R) from `.chk` + `.eig` at full
precision is the analogous follow-up (not done).

**Also fixed:** H(R) was folded on PythTB's centres (from text files) while A(R) was folded on the checkpoint's;
they differ by 5e-9 Å, which flips 20 near-tie H(R) entries at the 1e-5 Å tolerance edge (max 2.4e-4 eV) — a
small CLAUDE.md-#1 violation. `load_wannier_pair` now folds H(R) on the checkpoint centres (from the A(R) cache
for `.mmn`). Effect on dθ at nk = 6: 2e-7, negligible, but the rule is now one array per structure. The `.mmn`
shared-centre sweep above used one (PythTB) array for all four folds, so it was already internally consistent.

## Shared rule with Jae-Mo's formula, and a correction to §2 (2026-09-23)

`AXION_POSITION_SOURCE=chk_transl_inv_full AXION_SHARED_WS=1`, nk 4–10, full precision. dθ (base / mode):

| nk | `.mmn` own | `.mmn` shared | tif own | tif shared |
|---|---|---|---|---|
| 4 | +0.00377 / +0.00057 | +0.00379 / +0.00223 | +0.0384 / +0.0411 | +0.0143 / +0.0110 |
| 6 | −0.00660 / −0.00713 | −0.00583 / −0.00607 | −0.00465 / −0.00461 | −0.00372 / −0.00450 |
| 8 | −0.00555 / −0.00567 | −0.00555 / −0.00564 | −0.00369 / −0.00168 | −0.00330 / −0.00413 |
| 10 | −0.00624 / −0.00676 | −0.00529 / −0.00545 | −0.00392 / −0.00355 | −0.00281 / −0.00365 |

**Correction to §2.** "The tie term vanishes on the ab-initio mesh" holds only when the images of a class are
copies (`.mmn`, H(R)). Under `transl_inv_full` each image carries `exp(−i b·R/2)` per b, and tied images differ by
`(−1)^{n_i}`: a (½, ½) split can cancel at mesh k while a (1, 0) split does not. So for the tif family the
own-centre tie error is O(1)/Δβ in `∂_β A(k)` **at every k, on-mesh included**. That is why the tif base end moves
at nk = 4 and 8. The shared rule matters more for tif, not less; own-centre tif (and Jae-Mo's files) carry it
everywhere.

**Open.** tif-shared is not converged at nk = 10 (base end −0.00372 → −0.00330 → −0.00281) and its endpoints differ
by 20–30% (`.mmn`-shared: 3%). Endpoint spread is a linearity test, not validity, but it is a warning; needs
nk ≥ 12. The family gap persists: tif/`.mmn` = 0.53/0.67 at nk = 10, shared. Settling it needs a denser
ab-initio mesh.

## Both formulas from one `.chk`: the gap is the formula (2026-09-23)

`position_matrix_mmn_formula` applies the production `.mmn` formula to the checkpoint's `m_matrix` (which is
`V^† M_raw V`, what `connection_from_mmn` forms from the `.mmn`). Source `chk_mmn`. Against the cached
`.mmn` A(R): same 1097 R, relative L2 9e-7, max 4e-7 Å (`validation/inputs/compare_mmn_formula_chk_vs_cache.py`,
base only; origin of the 1e-6 not pinned down). dθ, shared centres: identical to the cache run to the printed
digits at nk 6–10 (nk 4 differs by 6e-6 absolute; plausibly the `.chk`-vs-PythTB centre switch at edge ties).

Matched comparison (same `.chk`, H(R), shared rule, code; only the A(R) formula differs):

| nk | `.mmn` formula | tif formula | tif / `.mmn` |
|---|---|---|---|
| 6 | −0.00583 / −0.00607 | −0.00372 / −0.00450 | 0.64 / 0.74 |
| 8 | −0.00555 / −0.00564 | −0.00330 / −0.00413 | 0.59 / 0.73 |
| 10 | −0.00529 / −0.00545 | −0.00281 / −0.00365 | 0.53 / 0.67 |

The family gap is the formula alone: k-space phase (`r_j` vs `(r_i+r_j)/2`), the real-space `exp(−i b·R/2)`,
and Hermiticity (explicit pairing vs by construction). `.mmn` is nearly nk-converged (−4.7%, then −0.2% to nk 12);
tif is still drifting (−11%, −15%). Next: tif to nk 12–14; toggle the three ingredients one at a time.

## All six at nk = 4 and 8 (2026-09-23)

dθ (base / mode). Per-run and his-files rows run through today's pipeline, logs in the session scratchpad.

| | nk = 4 | nk = 8 |
|---|---|---|
| shared: `.mmn` formula (`.chk`) | +0.003780 / +0.002226 | −0.005552 / −0.005640 |
| shared: tif formula (ours) | +0.014267 / +0.010972 | −0.003299 / −0.004134 |
| per-run: `.mmn` formula (`.chk`) | +0.003766 / +0.000573 | −0.005554 / −0.005670 |
| per-run: tif, full precision | +0.038367 / +0.041148 | −0.003687 / −0.001678 |
| per-run: tif, A(R) rounded to 6 dp | +0.038752 / +0.041575 | −0.003517 / −0.001520 |
| per-run: Jae-Mo's files | +0.038744 / +0.041573 | −0.003518 / −0.001519 |

* 1:1: rounded-ours vs his files agree to 2e-4 / 4e-5 (nk 4) and 3e-4 / 5e-4 (nk 8); with nk 6 (2e-4) the
  re-implementation reproduces his dθ at every nk tested. His print rounding costs 1% (nk 4), 5% / 10% (nk 8).
* `.mmn` formula: per-run vs shared agree on-mesh to 0.03% / 0.5% (nk 8) and 0.4% base (nk 4); the mode end
  moves at nk 4 through the mode's own `b`-weighted k-derivatives. The `.chk` rebuild matches the Aug 24 cache run
  to ~1e-4.
* tif formula: per-run vs shared differ by O(1) even on-mesh (nk 4 base +0.0384 → +0.0143; nk 8 mode −0.00168 →
  −0.00413), the image-dependent `(−1)^{n_i}` phases. Per-run tif (and his files) are tie-dominated; only the
  shared rows compare the formulas fairly (nk 8: 0.59 / 0.73 of `.mmn`).
