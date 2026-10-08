# #17 — The segni fix (#16) stalls the YIO base SCF at ~3e-5 Ry

> **Status: OPEN** (cause confirmed 2026-10-07, candidate fix in testing).
> - **Cause:** the #16 patch itself. From the same stalled state with the same settings, the unpatched binary
>   converges in 153 iterations and the patched binary does not converge in 338.
> - **Candidate fix:** `paw_mag_regularize.patch`, smooth |m| in the local-frame PAW GGA. Cluster
>   `~/q-e-segnireg` (branch `fix/paw-mag-regularize`), δ set at run time by `QE_PAW_MAG_DELTA`. Tests R1-R3
>   (δ = 1e-5, 1e-4, 1e-3 e/bohr³) running.
> - **Production meanwhile:** base 207546 (stalled at 3.1e-5 Ry) is used as the best available base; its
>   eigenvalue asymmetry is 0.10 meV worst case (bug #16 file, production table).

**Found:** 2026-10-06/07. **Affects:** YIO base (m-3m′) SCF with the segni-fixed pw.x (`e270b05b`). The three
modes (−4′3m′) converged with the same binary (mode_1 and mode_3 plain β 0.03, mode_2 after switching to Davidson).

## Symptom

Every base SCF with the patched binary plateaus at an estimated accuracy of 2e-5 to 2e-4 Ry and never reaches
`conv_thr = 1e-7`:

| run | mixing | diagonalization | result |
|---|---|---|---|
| 207546 | plain β 0.03 | CG | ~3e-5 Ry after 500 iterations (used for production) |
| 209024 | local-TF β 0.05 | CG | 2.5-6.7e-5 after 88 |
| 209042 | local-TF β 0.05 | Davidson | plateau |
| 209049 | plain β 0.1 | Davidson | 5e-5 to 2e-4 after 70 |
| T1c 209057 | plain β 0.03 | CG | 2e-5 to 1e-4 after 338 (killed) |

The slow mode is the small O(48f) moment drifting by ~7%. Before the patch the same input converged (192344 in
252 iterations, 193145 in 411). Davidson fixed mode_2's stall but not base's, so CG band noise isn't the cause
for base.

## The decisive test (`base/u_3.0/stall_test_2026-10-07/`)

T1 and T1c restart from the identical stalled state (207546's SCF state at 14:32, `startingpot = 'file'`) with
identical inputs (plain β 0.03, `mixing_ndim` 20, CG, 48 cores). Only the binary differs.

| run | pw.x | result |
|---|---|---|
| **T1** (209056) | unpatched `c41c409` | **converged in 153 iterations**, −4624.80687790 Ry, \|m\|abs 3.32 μB |
| **T1c** (209057) | patched `e270b05b` | 338 iterations, still 2e-5 to 1e-4 Ry, \|m\|abs 3.37 μB |

So the patch makes the base SCF map non-convergent; the starting state and the mixing are not the problem.

## Suspected mechanism

With the patch and `lsign = .false.`, the PAW one-centre GGA always uses the local frame,
n↑,↓ = ½(n ± |m|). `|m|` has a cone-shaped kink at m = 0, and inside the PAW sphere this kinked function is
projected onto a finite Y_lm set (`PAW_rad2lm`, lmax = `i%l`) before it is differentiated (`PAW_gradient`), so the
truncation rings. The field goes back along m/|m| (`compute_pot_nonc`), a direction set by noise wherever |m| is
small. Near m = 0 the map m → (n↑, n↓, B) is therefore not differentiable, and a linear mixer can't settle the
slow mode. Base has the highest symmetry, so more of its oxygen magnetization sits near zero than in the modes.
Not yet proven: this is the hypothesis the regularization tests.

The plane-wave part (`gradcorr.f90`) uses the same local frame without the Y_lm truncation, and it converges.

## Candidate fix

`calculations/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/paw_mag_regularize.patch`, on top of
`paw_segni_lsign.patch`. Only in the local frame (`lsign = .false.`):

- `compute_rho_spin_lm`: |m| → |m̃| = |m|² / √(|m|² + δ²). Smooth (analytic in m), ≈ |m|²/δ near zero,
  ≈ |m| − δ²/(2|m|) far from it.
- `compute_pot_nonc`: B = vs · ∂|m̃|/∂m = vs · m · (|m|² + 2δ²) / (|m|² + δ²)^{3/2}, the exact derivative, so the
  energy and potential stay consistent. With δ = 0 it reduces to the old vs · m/|m|.
- δ from the environment variable `QE_PAW_MAG_DELTA` (e/bohr³, default 0 = unchanged behaviour), printed once in
  the output. `compute_drho_spin_lm` (linear response only) is unchanged.

Depends only on |m|, so it stays exactly rotation-invariant and keeps #16's symmetry. Tests R1-R3 rerun T1c with
δ = 1e-5, 1e-4, 1e-3 e/bohr³.

## Remaining

1. R1-R3: does any δ converge, and how far does the total energy move against δ?
2. If one converges: check the eigenvalue symmetry of its NSCF (should stay at the #16 level), then decide
   whether to redo base production with it.
3. If none converges: the kink hypothesis is wrong. Next: dump the PAW-sphere |m| of the slow O(48f) mode.

## Update 2026-10-08

**The |m| regularization does not fix it.** R1-R3 (δ = 1e-5, 1e-4, 1e-3 e/bohr³) each ran 500 iterations at about
3e-5 Ry, with |m|abs 3.35. The kink hypothesis is not supported.

**LDA converges the base.** `lda_test_2026-10-08/`, all from scratch with the segni-fixed binary:

| run | result |
|---|---|
| YP (PBE control) | still wandering 2e-7 to 1.5e-5 Ry at 72+ iterations, \|m\|abs 3.32 → 3.40 |
| YL1 (`input_dft = 'pz'`, `rel-pbe` PAW) | converged in 21 iterations |
| YL2 (`rel-pz` PAW) | converged in 22 iterations |

Under LDA: Ir moments 0.283 μB, four equal, exactly along ⟨111⟩; gap 0.68 eV on the SCF mesh.

This matches MnBi2Te4's noncollinear-GGA problem (bug #18, `../progress/2026-10-08_mbt_noncollinear_gga_lda.md`
§12). **Open:** whether YIO moves to LDA+U (check the full-mesh gap at U = 3 first).

