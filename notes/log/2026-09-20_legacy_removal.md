# Removing the earlier pipeline (2026-09-20)

`modules/berry_curvature.py`, `modules/axion_angle.py` and the old Y2Ir2O7 driver
`axion_response_vs_nk.py` are deleted. `modules/` is now `curvature.py` and `axion.py`
(plus two unrelated utilities). Everything that imported them was ported and, where the
data still exists, rerun on the new code and compared with what the old code had saved.

**None of the removed files were in git.** The only copy is
`tmp/legacy_archive_2026-09-20/` (git-ignored; see its `README.txt`).

## What changed

| | before | after |
|---|---|---|
| `modules/` | `berry_curvature.py`, `axion_angle.py`, `curvature.py`, `axion.py` | `curvature.py`, `axion.py` |
| libs | `WannierConnection` objects threaded through `evaluate_trial`, `trace_omega_xy`, `build_model`, `mixed_curvature_traces` | a `PositionTerms` (`pos`) is threaded instead; `fields`/`connection`/`omega` do the rest |
| notebooks | 9 called the old classes directly; 4 more threaded the connection object through the libs | all ported (variables `connection`/`fields` renamed `pos`/`f` to avoid clashing with the new function names) |
| unit tests | 56 | 49: minus 5 covariance and 3 pipeline tests (ported or duplicated in `test_curvature.py`) and 3 legacy-equivalence tests; plus 4 that need no legacy code |
| `axion.py` | `_prepare_lambda`, `_parametric_connection` private | `prepare_lambda`, `parametric_connection` public (the Born-charge check needs Lambda(R)) |
| `check_pr702_rdat_symmetry` | monkeypatched TensorFlow with numpy | shim removed; `curvature.py` is numpy only |

The four new unit tests replace the legacy comparisons with independent ones: the Bloch sums
against a brute-force loop over R (from the raw stored blocks), `A_beta` Hermitian, `dA_beta`
equal to a central difference of `A_beta`, and the frozen-gauge / Lambda `beta_terms`.

**Kept on purpose:** `calculations/MnBi2Te4/phonon/axion_response_nk_sweep.py`. It produced the
paper's MnBi2Te4 number (`paper/main.tex`, `d_Q theta ~ 10.0 rad/A`), which
[`2026-08-20_handoff_wannier_conventions.md`](2026-08-20_handoff_wannier_conventions.md) (3a) flags as needing a rerun. It
now runs on a verbatim copy of the old module in `calculations/MnBi2Te4/phonon/frozen_legacy/`.
Delete both when the MnBi2Te4 structures are regenerated. It has only diagonal Wannier centres
(these runs kept no `_r.dat`), so it could not be ported to the `_r.dat` route.

*2026-09-24:* both deleted, with their `nk_sweep_shared_embedding*` outputs. The paper's numbers are in
flux (everything will be recomputed), so reproducing the old 10 rad/A was no longer worth keeping them.
See `implementation_map.md` for the old-path analysis code that stays to be ported.

## Old versus new

**Direct comparison, legacy code still present.** Eigenvalues of every Hermitian plane matrix
`Omega_{xy,yz,zx}` (gauge invariant, scale 1 to 8 A^2), at random k, for the production band group
and for the full space, on the real models the notebooks build (SrTiO3 sets a, b, c, sr_b, sr_c;
conventions `mmn_pair`, `rdat_direct`, `rdat_pair`, `rfull_pair` where the inputs exist):
**worst relative difference 1.2e-9.** (Traces are useless here: for an occupied-only set they are
zero to 1e-16 and the comparison is noise against noise, CLAUDE.md pitfall 5.)

**Notebook reruns against the results the old code saved.** Only quantities above 1e-5 are
compared; those below are noise-floor residuals (symmetry-forced zeros at 1e-9 to 1e-16) whose
relative difference means nothing.

| notebook | compared with | worst relative difference |
|---|---|---|
| `check_vs_pythtb` (YIO, 176 WF, 156 occ) | PythTB itself | k-k: **5.7e-17** at scale 0.61; mixed k-beta block, hand assembly and Hermiticity (5e-15): all assertions pass |
| `born_effective_charge`, `rdat` source, trial_03 (48 WF, 38 occ) | saved `born_charge_Ti_Q_01A_mhn0_nows.npz` | **4e-16**, every nk from 4 to 16, int, cross, ext and total separately. `Z*_zz` = 7.57196 (+3.9% vs DFPT 7.28829) |
| `compare_position_matrix_sources` | `results/position_matrix_sources` | 4.9e-13 (2587 entries) |
| `check_convention_symmetry` | `results/convention_symmetry` | 2.6e-7 (23 curve arrays) |
| `check_pr702_rdat_symmetry` | `results/pr702_rdat_symmetry` | 6.1e-7; its own hard-coded `BASELINE` assertion passes |
| `check_rdat_ndegen_applied_symmetry` | `results/rdat_ndegen_applied_symmetry` | 1.6e-6 (23 leaves) |
| `plot_curvature_generic_line` | nothing saved before | reran; new output in `results/curvature_generic_line` |
| `check_magnetic_symmetry` (YIO, 176 WF) | `results/magnetic_symmetry` | median 1.9e-6, worst **1e-4** (mode-1 `cross`); explained below |

`check_interpolation_inputs` cannot run whole (below), so its ported functions were run on an
available model instead: A(k) against its central difference 5.5e-9, analytic `dH` 1.8e-9, off-grid
covariance of `H`, `A`, `curl` 4e-16, 9e-16, 5e-15, the same regime as the saved trial03/trial04 values.

The old saved tables also contain rows the new run cannot reproduce: the four PR #702
`_r.dat` variants in `pr702_rdat_symmetry` (files gone). They are untouched by the port.

## What could not be rerun (inputs, not the port)

* `.chk` of Ti_Q_2A trial03 is missing: `compare_trial_curvature`, `check_trial_curvature_geometry`,
  `check_vs_wilson_flux`, `check_interpolation_inputs`. The two libs they share with the notebooks that
  did run (`evaluate_trial`, `symmetry_checks`) were exercised by `compare_position_matrix_sources` and
  `check_rdat_ndegen_applied_symmetry`.
* the SrTiO3 `.mmn` files are deleted: `born_effective_charge` with its default `mmn` source (the
  `rdat` source above was used instead).
* `check_vs_postw90_line` and `check_vs_postw90_full_mesh` need external postw90 output (`TRIAL*_POSTW90`
  and `KPOINTS_NPZ` are `None`). Their edit is the same three lines as everywhere else.
* `calculations/diagnostics/2026-09-20_ndegen_pair/refold_driver.py` was repointed at
  `run_axion.py` (it now patches `pythtb.W90`). Not rerun: it is a parked 25 minute experiment.
* `meetings/2026-08-31.ipynb` still imports the removed module. It is a dated meeting record, so it was
  left alone; its saved outputs are intact.

## Findings along the way

1. **A latent trap the port removes.** The old `WannierConnection(model)` with no `embedding=`, on a
   model built with `orb_vecs = 0`, assumed the stored centres equalled `orb_vecs`, zeroed the R = 0
   diagonal and gave an O(1) wrong curvature (relative 1.2 against the tau = centres value). Passing
   `embedding=orb_vecs.copy()` was correct. This is the `rdat_pair_tau0` convention of
   `check_convention_symmetry_lib`. **No saved result used it** (checked: no `tau0` key in any file under
   `validation/results`). `position_terms` always subtracts `model.orb_vecs`, and in the rerun tau = 0 and
   tau = centres agree to 3.5e-13 on all 14 curves with signal, the exact invariance the physics requires.
2. **Checklist item 13 in `berry_curvature_decomposition.md` still said the two dtheta endpoints "must agree"**
   as a validity diagnostic. Corrected: the endpoint spread is a linearity measure, there is no right
   endpoint, only a sign flip is diagnostic (as already stated in the `axion.py` docstring).
3. **Most of `validation/`, `wannier/` and the new modules are untracked.** Do not rely on `git checkout`
   there. Consider committing.
4. The reruns that failed in 2 s failed while loading inputs or checking config, before any curvature
   was computed, so they say nothing about the port either way.

## The 1e-4 on YIO mode 1, and why it is the old code's precision

`check_magnetic_symmetry` agrees with the old saved numbers to a median of 1.9e-6 but to only 1e-4 in a few
mode-1 `cross` residuals (which are themselves differences of O(0.5) terms). The old notebook built its
connection at `complex64` (the `WannierConnection` default); the new one runs in `complex128`. The first
guess, that this is single precision, did **not** survive a naive test: rerunning the new code at all-`complex64`
matched the old numbers *worse* (2.2e-4), because it also degrades the eigensolve, which the old code ran in double.

The settling experiment ran the old code itself (the frozen copy) on the mode-1 z axis, 17 points:

| comparison | internal | cross | external |
|---|---|---|---|
| old @ `complex128` vs new @ `complex128` | 4.6e-12 | 2.3e-12 | 7.5e-14 |
| old @ `complex64` vs old @ `complex128` | 0 | 3.9e-5 | 4.0e-6 |
| new @ `complex128` vs saved old curve | 4.6e-12 | 3.9e-5 | 4.7e-6 |
| old @ `complex64` rerun vs its own saved curve | 2e-12 | 1.6e-5 | 4.5e-6 |

At equal precision the two implementations agree to 1e-12. The old-vs-new gap is exactly the old code's own
`complex64` error, which does not even reproduce its own saved curve better than 1.6e-5. The new numbers are
the accurate ones. Nothing in the YIO symmetry conclusions moves (residuals 1e-2 to 2, four digits agree).

## Where things are now

`notes/code.md` (function to equation table), `validation/README.md` (rewritten pipeline
diagram and unit-test table), `notes/theory/external_terms.md` and `berry_curvature_decomposition.md` (code
references moved to `curvature.py`/`axion.py`; the math is unchanged).
