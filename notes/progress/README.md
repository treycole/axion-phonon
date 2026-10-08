# Bug and progress log

One dated record of every bug found in this project: what it was, why it was a bug, how it was fixed, and
where it stands. Newest is last. Problems are numbered **#1–#17** in the order found; a fix, refinement or
open item cites the number of the problem it belongs to. Progress milestones are unnumbered. Each number links
to that bug's detailed file in [`../bugs/`](../bugs/). Each entry points to the note (in this folder, `YYYY-MM-DD_<topic>.md`) or memory record holding the evidence. Two dates
(2026-07-27, and the Born-charge item on 2026-09-21) come from when the finding was recorded in memory, not from
a dated note.

**Status at 2026-10-06:** the cause of the YIO symmetry floor is found and fixed (a QE bug, 2026-10-05), and
all four YIO structures are being recomputed with the patched pw.x. Every YIO number computed before then,
including dθ, is superseded.

## 2026-07-27

**[#1](../bugs/01_wannier_diffuse_projections_gap_collapse.md) — Bug: YIO Wannier model closed the gap (diffuse empty projections), and n_occ was wrong.**
- **What:** dθ came out at ±thousands and jumped with nk. The interpolated minimum gap was 0.006–0.017 eV,
  against ~0.72 eV for a healthy model.
- **Why it was a bug:** the projections included empty, diffuse shells (Y 4d, Ir 6s/6p). These can't be
  one-shot Wannierized: their spreads were 9–24 Å² against ~1 Å² for real WFs. They injected spurious
  long-range hoppings into `_hr.dat` that nearly closed the gap, and the Kubo curvature ∝ 1/gap² blew up.
  Separately, n_occ = 104 was carried over from a model without Y 4p / O 2s, which puts the chemical potential
  inside the O 2p manifold.
- **Fix:** "Option A" projections, `Y: p`, `Ir: d`, `O: s;p` (176 WF, all spreads ≲ 1.1 Å²), and
  **n_occ = 156**.
- **Source:** memory `wannier-gap-collapse-Os`.

## 2026-08-19

**[#2](../bugs/02_pythtb_embedding_noop.md) — Bug: the PythTB embedding was a silent no-op in the axion scripts.**
- **What:** `model.lattice.orb_vecs = tau` did nothing, because `TBModel.lattice` returns a copy.
- **Why it was a bug:** H and v kept each run's own Wannier centres while the external terms used `tau`. That
  gives a mismatch linear in τ. In the dθ finite difference, base and mode didn't share an embedding, which
  violates assumption A3 of the derivation. One symptom: the two endpoint estimates of dθ had opposite signs
  (e.g. 173199: −0.32 / +0.30).
- **Fix:** build with `W90.model(orb_vecs=tau)` plus a post-build assertion. Results are then exactly
  τ-invariant.
- **Status:** every dθ computed before this is invalid. The MBT driver behind the paper's ~10 rad/Å still had
  the bug; it was deleted 2026-09-24, and MBT will be recomputed on the new pipeline.
- **Source:** `2026-08-20_handoff_wannier_conventions.md` §3a; memory `dtheta-dq-gauge-pitfall`.

## 2026-08-20

**[#3](../bugs/03_hr_ar_ws_convention_mismatch.md) — Bug: H(R) and A(R) used different Wigner–Seitz conventions.**
- **What:** A(R) came from WannierBerri's representative rule and H(R) from `_hr.dat`'s.
- **Why it was a bug:** a finite k-mesh fixes H(R) and A(R) only up to the aliasing class, and the
  representative rule is a choice of interpolant. Mixing two rules gives a curvature that belongs to no
  Hamiltonian. It is invisible at the ab-initio k-points and shows up only on the dense mesh. It produced
  spurious symmetry violations: CLAUDE.md pitfall #1.
- **Fix:** `ws_reassign()` / `apply_ws_reassignment()` re-fold H(R) with the same rule before `.model()`.
  Agreement with WannierBerri is rel 0.0004%.
- **Source:** `2026-08-20_handoff_wannier_conventions.md` §1, §3b.

## 2026-08-27

**[#4](../bugs/04_stock_rdat_position_matrix.md) — Bug (a Wannier90 limitation): stock `_r.dat` cannot carry the position matrix.**
- **Why:**
  - The `6F12.6` print format zeroes 52% of components.
  - No link phase is applied, and the sum over neighbours happens before the file is written, so the
    information is unrecoverable.
  - `use_ws_distance` never reaches the writer.
- **Fix:** build A(R) from `.mmn`, with the `_AA_cache/*.npz` replacing 40.5 GB of `.mmn` by 2.6 GB.
  Production later moved to `.mmn`-free routes (2026-09-13 onward).
- **Source:** `2026-08-27_rfull_vs_mmn_position_matrix.md` §0; memory `rdat-position-matrix-unusable`.

## 2026-09-10

**[#5](../bugs/05_pr702_rdat_ignores_use_ws_distance.md) — Bug (Wannier90): PR 702's `_r.dat` silently ignores `use_ws_distance`.**
- **What:** PR #702 (`transl_inv_full`) fixes the link phase. Its `_r.dat` equals postw90's A(R) to 4.4e-16,
  but at `use_ws_distance = F` whatever the `.win` says. The file is byte-identical with the flag T or F.
- **Why it was a bug:** that's the wrong Wigner–Seitz representative for symmetry work. It is still 318,000×
  worse than `.mmn` on a symmetry-forced zero. The print format contributes only ~1% of that.
- **Fix:** see 2026-09-13 (Jae-Mo's branch).
- **Source:** `2026-09-10_pr702_transl_inv_full_rdat.md`.

## 2026-09-13

**Fix for #4 and #5: Jae-Mo Lihm's `write_ndegen_applied` branch.** With `transl_inv_full = T, use_ws_distance = T,
write_ndegen_applied = T`, the representative problem is fixed. The output matches postw90's `write_aa_r` at
`use_ws_distance = T` bit-exactly (after rounding) on two test systems.
- **Source:** `2026-09-13_jaemo_transl_inv_full_ws_distance.md`.

## 2026-09-16

**Progress:** the Jae-Mo route (`source="rdat_ndegen_applied"`) is closed on the production SrTiO3 system.
It is symmetry-exact and agrees with an independently built postw90 `rfull` export to ~1e-8.
- **Source:** `2026-09-16_jaemo_rdat_production_closure.md`.

## 2026-09-20

**[#6](../bugs/06_base_mode_ws_tie_inconsistency.md) — Bug found: base and mode don't share a Wigner–Seitz rule.**
- **What:** on the first YIO base/mode pair through the Jae-Mo route, 5.5% of |ΔH| came from ties split
  differently in base and mode. This affected both position-matrix routes.
- **Status:** open on this date; diagnosed and fixed on 2026-09-23 (below, #6).
- **Source:** `2026-09-20_jaemo_ndegen_yio_mode1.md`.

## 2026-09-21

**[#7](../bugs/07_qe_new_ns_nc_symmetrizer.md) — Bug (QE 7.5/7.6): `new_ns_nc` symmetrizes the noncollinear Hubbard occupations with the wrong atom.**
- **What:** it builds atom na's ns from its *image* `irt(isym, na)` instead of the *pre-image*
  `irt(invs(isym), na)`.
- **Why it was a bug:** that's harmless for self-inverse operations. It's wrong for every operation of order
  3, 4 or 6 that permutes Hubbard atoms (28 of 48 in YIO). `symmetry_with_labels = .true.` with four Ir
  species is what lets such operations appear. The result was Ir moments differing by ~9% and tilted off
  their axes, and a Hubbard spin channel that was too weak.
- **Fix:** a two-line patch (`calculations/qe_patches/`), verified in a real pw.x run: the four Ir moments
  become identical.
- **Status:** all YIO DFT reran on QE 7.6 + this patch (2026-09-22 to 09-29). The YIO symmetry floor
  remained; its cause was #16.
- **Source:** `2026-09-21_qe_new_ns_nc_bug.md`; memory `qe-new-ns-nc-bug`.

**[#8](../bugs/08_born_charge_ionic_term.md) — Bug (bookkeeping): Born effective charge ionic term.**
- **What:** the ionic term is not simply `z_valence`.
- **Why it was a bug:** bands below `dis_win_min` that sit on the displaced atom move rigidly with it and
  contribute −1 each. SrTiO3/Ti: Z^ion,eff = 12 − 8 = +4, n_occ = 30, not 40. Getting it wrong moved Z* from
  7.1 to 15.1.
- **Fix:** Z^ion,eff = z_valence − (excluded electrons on the displaced atom). trial_02 stands at 7.370;
  trial_03 moved to 6.665.
- **Source:** memory `born-charge-ionic-accounting` (recorded on this date).

## 2026-09-22

**[#9](../bugs/09_mbt_t2_rounded_cell_and_labels.md) — Bug (input): MnBi2Te4 T2⁻ cells lose C3, and `symmetry_with_labels` triggers `new_ns_nc`.**
- **Why it was a bug:**
  - The T2⁻ inputs used a `CELL_PARAMETERS bohr` cell rounded to 8–9 digits. The ~5e-5 bohr error is above
    QE's 1e-6 tolerance, so C3 was lost (6 or 12 operations instead of the full set).
  - With the flag on, T2⁻ gains two order-6 operations swapping Mn1 ↔ Mn2, exactly the `new_ns_nc` trigger.
- **Fix:** precise ångström cell with positions = base + Q·mode; flag on only with the patched QE.
- **Source:** memory `mbt-symmetry-with-labels`.

## 2026-09-23

**#6, diagnosed and fixed: WS image weights must also be shared across the β finite difference (assumption A3′).**
- **What:** the in-house fold used each run's own Wannier centres, while the Bloch phase used the shared base
  τ.
- **Why it was a bug:** a tie split (½, ½) in base and (1, 0) in mode adds an `O(1/Δβ)` term to ∂_β H off the
  ab-initio mesh. The derivation's β derivative assumes fixed τ *and* fixed image weights.
- **Fix:** `load_wannier_pair(shared_ws_centres=True)` folds both structures on base's centres.
- **Source:** `2026-09-23_ws_ties_in_beta_derivative.md`.

**[#10](../bugs/10_hr_ar_folded_on_different_centres.md) — Bug (small): H(R) and A(R) folded on different centre arrays.** H(R) used PythTB's centres from text files
and A(R) the checkpoint's. They differ by 5e-9 Å, enough to flip 20 near-tie H(R) entries (max 2.4e-4 eV).
- **Fix:** one array per structure (the checkpoint's).

**[#11](../bugs/11_rdat_print_precision.md) — Finding:** `_r.dat`'s 6-decimal print costs ~3% of dθ. Full-precision A(R) from the `.chk` is preferred.

## 2026-09-27

**[#12](../bugs/12_qe_paw_vrad_grad_uninitialized.md) — Bug (QE 7.4–7.6): PAW `v_rad` / `g_rad` uninitialized in `compute_pot_nonc`.**
- **What:** the noncollinear-GGA PAW one-centre XC routine allocates these arrays and never zeroes them.
  `g_rad` is *accumulated into* for fully relativistic PAW.
- **Why it was a bug:** the result depends on leftover memory. It caused MBT's iteration-1 blow-up and the
  −2 μB net moment: Bi and Te start at m = 0, where the garbage enters fully. The zeroing lines were dropped in
  QE commit 4ce366b33 (2024-03).
- **Fix:** fixed upstream (QEF #863, `87ac195cf`, 2026-07-10). The cluster `~/q-e` build (pw.x md5 `c41c409`)
  has it.
- **Status:** it does not cause the YIO symmetry breaking (that was #16; runs A/B, 2026-10-05). It does shift the YIO SCF
  state (C vs A: 7.6 → 4.8 meV).
- **Source:** memory `qe-paw-nc-uninit-vrad`.

## 2026-09-29

**Progress:** first YIO base/mode pair on QE 7.6 + the `new_ns_nc` fix. H(R)-only symmetry residuals fell
2–4×, but the total got worse, and the floor appeared to sit in A(R). dθ (mode 1) ≈ −0.0175 rad/Å, about 5×
the unpatched value. That number is now superseded (2026-10-05). Two diagnostics of the A(R) floor came back
negative (centre symmetrization; symmetrized BZ integral).
- **Source:** `2026-09-29_yio_qe76_new_ns_nc_pair.md`.

## 2026-09-30

**[#13](../bugs/13_ws_tie_tolerance_too_tight.md) — Bug: the Wigner–Seitz tie tolerance (1e-5 Å) was too tight.**
- **What:** near-ties whose two candidate images differ by < 1e-3 Å, only because of noise in the centres,
  weren't treated as ties.
- **Why it was a bug:** one image got all the weight instead of a symmetric split. That broke base's
  inversion: 14% of such near-ties were inversion-inconsistent.
- **Fix:** `_DEFAULT_WS_DISTANCE_TOLERANCE` raised to 1e-3 Å. Base inversion residual 12–14% → 0.3%. Refined
  on 2026-10-05 (#13, refined).
- **Source:** `2026-10-05_yio_ws_tolerance_full_group.md` §1.

## 2026-10-04

**[#14](../bugs/14_mbt_stale_binary.md) — Bug (operational): MBT u_3.0 was running a stale, buggy pw.x.**
- **What:** u_3.0 still had the old binary (md5 `1541022`) after the 09-28 update, which had only touched u_4.0.
- **Why it was a bug:** the run's total magnetization drifted to (0, −0.45, −8.33) μB, while u_4.0 on the fixed
  binary stayed at 0. About 1.5 h went into auditing QE's symmetrizer, which was correct.
- **Fix:** copied the fixed pw.x in and resubmitted.
- **Lesson:** md5sum pw.x first, and monitor the total magnetization on every check.
- **Source:** memories `qe-paw-nc-uninit-vrad`, `monitor-magnetization-and-check-binaries`.

## 2026-10-05

**[#15](../bugs/15_mbt_u3_base_rounded_cell.md) — Bug (input): MBT base u_3.0 also used a rounded bohr cell (same defect as #9).**
- **What:** QE found 4 symmetry operations instead of 12. That input also had ecutrho 240 against u_4.0's 500.
- **Fix:** u_3.0 rebuilt as an exact copy of u_4.0, with only the U value changed (job 202971).
- **Source:** memory `mbt-symmetry-with-labels`.

**#13, refined: the tolerance's physical window.** Base's Wannier centres are symmetric to 5.7e-4 Å, so
noise-broken ties have gaps ≤ 2.05e-3 Å. Genuinely distinct images start at 2.12e-2 Å, with nothing in
between. Any tolerance in [2.05e-3, 2.12e-2) gives the identical fold. The 1e-3 default still misses 34k
forced ties, which carried what was left of the inversion residual.
- **Recommendation:** 5e-3 (not yet applied).
- **Result:** at 5e-3 inversion is clean (0.18%), but the rotations stay at ~1.5% total / 26% external, so
  what's left isn't tie-breaking.
- **Source:** `2026-10-05_yio_ws_tolerance_full_group.md`.

**[#16](../bugs/16_qe_paw_segni_stale_ux.md) — Bug (QE, all versions checked): the PAW noncollinear GGA uses a stale `ux` without the `lsign` guard. This is
the root cause of the YIO symmetry floor.**
- **What:**
  - `compute_ux` sets `ux` to the first magnetic atom's starting moment (Ir1's, (−1,−1,−1)). When the moments
    turn out not to be parallel it sets `lsign = .false.`, but it never clears `ux`.
  - The smooth-grid GGA (`compute_rho.f90`) correctly uses `segni = 1` when `lsign` is false.
  - The PAW one-centre GGA (`paw_onecenter.f90::compute_rho_spin_lm`, line 2340, and
    `compute_drho_spin_lm`, line 2554) always uses `segni_rad = SIGN(1, m·ux)`.
- **Why it was a bug:** inside each PAW sphere, n↑ and n↓ are swapped wherever m points away from Ir1's axis.
  They jump across the surface m ⊥ ux, and the GGA differentiates those jumps, giving a spurious potential.
  Because `ux` doesn't rotate with the crystal, symmetry-equivalent atoms get different potentials. Only −3m′
  about [111] survives.
- **Effects in YIO:**
  - DFT eigenvalues break m-3m′ by 4.8–7.6 meV.
  - Ir1's on-site levels are 12 meV off the other three Ir.
  - The curvature residual of ~1.5% total (≈ 26% external) tracks it, with corr ≈ 0.9 for all parts.
- **What it corrects:** A(R) was never a separate source, and the 09-29 "floor moved into A(R)" reading was an
  artifact of relative normalization.
- **Why it was missed:** the SCF symmetrizes ρ, m, ns and becsum, which hides it. Only an NSCF with `nosym`
  exposes it.
- **How it was isolated** (cluster NSCF tests, `base/u_3.0/test_nscf_sym_2026-10-05`):
  - the PAW `v_rad` fix doesn't remove it;
  - it isn't heap- or process-count-dependent;
  - it isn't Hubbard or ortho-atomic;
  - LDA removes it (0.0 meV), because LDA has no gradients;
  - moving `ux` moves it;
  - the stored ρ, m, ns and becsum are all exactly symmetric.
- **MnBi2Te4 is not affected:** collinear Mn ±z gives `lsign = .true.` with `ux = z` (it prints "Fixed
  quantization axis for GGA"), and every MBT operation maps z to ±z.
- **Source:** `2026-10-05_yio_symmetry_floor_dft_origin.md`; the plain-language explanation with GGA examples is in
  `2026-10-06_segni_fix_validated_and_reruns.md` §2.

## 2026-10-06

**Fix for #16 validated: the PAW `segni_rad` `lsign` guard** (`2026-10-06_segni_fix_validated_and_reruns.md`).
- **Patch:** `calculations/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/paw_segni_lsign.patch`.
  It sets `segni_rad = SIGN(1, m·ux)` only `IF (lsign)`, and 1 otherwise, as `compute_rho.f90` does.
- **Build:** cluster `~/q-e-segnifix`, branch `fix/paw-segni-lsign`, commit `4b480ac68`; pw.x md5
  **`e270b05b`**.
- **Test:** run C2 (same fixed SCF save and k-points as C) takes the worst eigenvalue asymmetry from
  **4.83 meV to 46 μeV** over all 47 operations, with the Ir1 pattern gone.

**Progress:** production reruns submitted with the patched pw.x. Only the binary changed; each structure folder
has a `NOTE`.

| structure | job | started |
|---|---|---|
| base/u_3.0 | 207546 | 11:15 |
| mode_1 | 207547 | 11:34 |
| mode_2 | 207548 | 11:15 |
| mode_3 | 207549 | 11:15 |

To free a 64-core node, MBT u_4.0 job 205075 was killed for stagnation: 82 iterations, SCF accuracy flat at
1.4–2.8e-3. Its output is saved in `u_4.0/`. MBT u_3.0 (202971) is still converging slowly.

**Progress (evening):** with the patched pw.x, mode_1 and mode_3 SCFs converged. The base and mode_2 SCFs
stagnated (base hit 500 iterations at ~3e-5 Ry). mode_2 was killed and restarted with `local-TF` mixing, β = 0.1
(job 209009), and its stale runs (~884 GB) were deleted. The base NSCF continues for the DFT symmetry check only.

**Open:**
1. (#16) Wannier stage on the new DFT: base, then each mode's rigid-shift `.amn`, wannier90, and the
   `restart = plot` job.
2. (#16) `.eig` and curvature symmetry checks on the new runs, then dθ.
3. (#13) Apply the 5e-3 WS tolerance.
4. (#16) Delete the old buggy YIO runs (~3 TB) once the new ones validate; needs the user's OK.
5. (#12, #16) Update CLAUDE.md pitfall 7 and its "correct versions" list. **Done 2026-10-07.**
6. (#16) Decide whether to report the `segni` bug to QEF.

## 2026-10-07

**Progress: the 10×10×10 ab-initio mesh doesn't reduce the YIO symmetry floor** (on the pre-segni DFT;
`2026-10-07_validation_methods_and_nk10_mesh.md` §1).
- **Run:** base 192344, new 10×10×10 NSCF + Wannierization (202714; its final `wannier90.x` on 64 ranks was killed
  by `h_vmem` and rerun on one rank as 202972, 67 min, 29.6 GB).
- **Result:** full m-3m′ table, median total `chk` route 1.57% → **3.60%**, Jae-Mo files 8.10% → 11.9%. Only
  the cross part improves. The band residual stays at ~5 meV on every mesh and route, which is #16's DFT
  asymmetry.
- **Side result:** the `.chk` route is 2-5× more symmetric than the Jae-Mo `_hr`/`_r` files at both meshes
  (inversion 0.2% vs 13-15%).
- **Open:** redo on the segni-fixed base before deciding on a mode-1 10×10×10 run.

**Fix: the gzipped-checkpoint reader didn't handle Fortran subrecords** (same note, §2). The 10×10×10 `m_matrix`
record (3.96 GB) exceeds 2 GiB, so ifort splits it. `_SequentialFortranReader.read_record` now joins subrecords,
and `m_matrix` is a zero-copy `complex128` view. Checks: bit-identical on 8×8×8; the 10×10×10 centres match
`_centres.xyz` to 5e-9 Å.

**Documentation:** how every symmetry check is done (the `new_ns_nc` validation chain, the eigenvalue check and
why every gk lies on a Γ-centred grid, the on-site levels, the C2 validation, the curvature metric), and why the
SCF density is symmetric while the NSCF eigenvalues are not (same note, §3-4). CLAUDE.md pitfall 7 and the
"correct versions" list now include #16.

**[#17](../bugs/17_segni_fix_base_scf_stall.md) — Bug: the segni fix (#16) stalls the YIO base SCF.**
- **What:** every base SCF with the patched pw.x plateaus at 2e-5 to 2e-4 Ry (CG or Davidson, plain or local-TF
  mixing); the modes converge.
- **Test:** from the same stalled state and inputs, the unpatched binary converged in 153 iterations, the patched
  one not in 338. So the patch causes it.
- **Suspected why:** the local-frame |m| has a kink at m = 0, which the PAW sphere's finite Y_lm expansion and
  its m/|m| field direction turn into a non-differentiable SCF map where the magnetization is small.
- **Candidate fix (testing):** smooth |m| → |m|²/√(|m|²+δ²) with the matching field, δ from `QE_PAW_MAG_DELTA`
  (`paw_mag_regularize.patch`, `~/q-e-segnireg`); runs R1-R3 at δ = 1e-5, 1e-4, 1e-3 e/bohr³.
- **Meanwhile:** production uses base 207546 (stalled at 3.1e-5 Ry, worst eigenvalue asymmetry 0.10 meV).

**Decision: production A(R) is Jae-Mo's formula from the `.chk` (`transl_inv_full_from_chk`, shared WS).** It is now the
default of `run_axion.py` and `base_full_group_run.py`. From 09-30 to 10-06 both had silently used `right_centre_from_chk`
(the `.mmn` formula from the `.chk`). Redone with Jae-Mo's formula on the pre-segni pair (#13): the 1e-5 → 1e-3 tolerance
change moves dθ by ≤ 2% (mode 1 −0.01753, mode 2 +0.04517, mode 3 +0.1689 rad/Å at 1e-3, base end); the formula choice moves
it by up to 6%. Old-base full-group symmetry with Jae-Mo's formula: median total 2.02% (1e-3), 1.94% (5e-3).

## 2026-10-08

**[#18](../bugs/18_qe_noncollinear_gga_spurious_state.md) — Bug: noncollinear PBE has a spurious lower-energy MBT state.**
- **What:** every MBT PBE+SOC base SCF on the new QE stalled. From a random start (the default) they drift into a
  layer-asymmetric state with |m|abs ≈ 13 μB (vs 10.2) and 3.5 μB of transverse interstitial magnetization. It
  sits 55 mRy *below* the symmetric state, with the whole gain in the grid XC (−285 mRy), and never converges.
- **Why it is a bug:**
  - QE drops MBT's anti-translation {1′|½,½,½} ("supercell"), and the random start breaks the layer equivalence.
  - The same physics in collinear PBE without SOC (C1) converges symmetric; noncollinear (C2) finds the low state.
  - LDA (E1/E2) has no such state.
  - Mechanism: QE's fixed-axis GGA spin density ½(n ± sign(m_z)|m|) **jumps** wherever m_z crosses zero at finite
    |m| (84k grid links; +48.5% total variation over m_z). This is not the |m| kink.
- **Fix:** MBT moves to **LDA+U** (`rel-pz` PAW; Te generated with `ld1.x` from the official recipe). Base 209653
  converged in 12 iterations, symmetric to 3e-5 μB. `startingwfc = 'atomic'` also works in PBE (test D, 18
  iterations), but can't protect the T2⁻ modes.
- **Also found:**
  - The production u_4.0 cell is the exact CIF geometry; 152188's cell is the rounded one.
  - 152188 converged because its leftover restart flag set ethr = 1e-5 and its random start happened to be
    symmetric.
  - The QE 7.2 → 8.0dev code audit found no regression; noncollinear +U now runs as kind 0 Dudarev, equivalent at
    J = 0.
- **YIO:** the #17 |m|-kink regularization (R1-R3, δ to 1e-3) did **not** fix the base stall.
  - **LDA converges the YIO base** in 21-22 iterations (210526/210527): Ir 0.283 μB, exactly along ⟨111⟩;
    0.68 eV SCF-mesh gap.
  - **PBE (210525) wanders**, with |m|abs creeping up from 3.32 to 3.40. Same family as #18.
- **Why the atomic start converges:** the SCF preserves the symmetric subspace, and the atomic start lies in it. In
  PBE the symmetric state is a saddle, so it works only for the symmetric base, not for the T2⁻ modes.
- C2 and B2 were killed after their results were recorded.
- **Source:** `2026-10-08_mbt_noncollinear_gga_lda.md`.

