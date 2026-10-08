# #13 — Wigner–Seitz tie tolerance too tight (1e-5 Å)

> **Status: PARTIALLY FIXED.**
> - **2026-09-30:** `_DEFAULT_WS_DISTANCE_TOLERANCE` raised from 1e-5 to **1e-3 Å**
>   (`wannier/position_matrix.py:501`).
> - **2026-10-05:** the physically allowed range was found to be **[2.05e-3, 2.12e-2) Å**. 1e-3 is just below
>   it. **5e-3 recommended, not yet applied**; dθ was measured at 1e-3 on 09-30 (pre-segni DFT): mode 1 −0.01856/−0.01864, mode 2 +0.04387/+0.04558, mode 3 +0.1697/+0.1787 rad/Å at nk 12, `right_centre_from_chk`, shared WS. The tolerance was confirmed 2026-10-07 by rebuilding the base gap (matches 1e-3 to 1e-7 eV, not 1e-5). Not yet measured at 5e-3. **With Jae-Mo's formula (`transl_inv_full_from_chk`, shared WS; 2026-10-07), 1e-5 → 1e-3:** mode 1 −0.01752/−0.01761 → −0.01753/−0.01760, mode 2 +0.04593/+0.04758 → +0.04517/+0.04668, mode 3 +0.1684/+0.1775 → +0.1689/+0.1779 rad/Å. The tolerance moves dθ by ≤ 2%; the formula choice (tif vs right_centre) moves it up to 6%. Outputs `nk_sweep_tif_shared_ws_tol1e-3/` and (mode 2) `nk_sweep_tif_shared_ws_tol1e-5/`; logs `calculations/diagnostics/2026-10-05_ws_tolerance_window/dtheta_old_tif_*.log`.

**Found:** 2026-09-30, refined 2026-10-05. **Affects:** the in-house folds of H(R) and A(R)
(`_map_R_by_orbital_positions`), used by every `.mmn`/`.chk` route.

## Symptom
- **Base inversion residual 12-14%:** at 1e-5, YIO base's curvature broke inversion by 12-14% (total).
- **Concentrated in near-ties:** 14% of (orbital pair, R) entries whose two nearest images were within 1e-3 Å of
  each other were inversion-inconsistent.
- **Wider gaps were consistent:** they showed 0% inconsistency.

## Root cause
- **How the fold works:** it gives the closest image all the weight, splitting equally only among images whose
  distances are within `tolerance` of the minimum.
- **Why ties break:** the Wannier centres are symmetric only to ~5.7e-4 Å (numerical noise in the
  Wannierization). Bonds that symmetry makes exactly tied therefore differ by up to ~2e-3 Å.
- **The consequence:** at 1e-5 those ties weren't recognized, so one image got 100% where symmetry requires a
  split.

## Why it was a bug
The fold must be covariant under the crystal's symmetry. Breaking a symmetric tie arbitrarily makes H(R) and A(R)
non-covariant: a choice of interpolant that the symmetry forbids. It looks like a physical symmetry violation.

## Fix
- **2026-09-30:** default raised to 1e-3. Base inversion residual 12-14% → ~0.3% (total, max over k).
- **2026-10-05:**
  - **Measured the window:** all ~34M (pair, R) entries × 27 images give a noise cluster of excess distances up
    to 2.046e-3 Å, then **nothing** until 2.125e-2 Å.
  - **The upper side is safe:** entries that aren't ties even for atom-pinned centres have gaps ≥ 0.32 Å.
  - **One fold for the whole window:** any tolerance in [2.05e-3, 2.12e-2) gives the identical fold.
  - **What 1e-3 still misses:** 34k forced-tie entries, which carried what was left of the inversion residual.
  - **At 5e-3:** inversion is clean (0.18% total, 0.9% external).

## Verification

| tolerance | full m-3m′ group median total | median external | inversion total |
|---|---|---|---|
| 1e-5 | 6.55% | 89.6% | 10.5% |
| 1e-3 | 1.57% | 25.9% | 0.24% |
| 5e-3 | 1.55% | 25.6% | 0.18% |

What remained at 5e-3 on the rotations (~1.5%) was **not** tie-breaking. It was #16, in the DFT.

## Remaining
- **Change the default to 5e-3**, the middle of the window. It moves dθ, which should be measured on the new DFT.
- **The window is base's:** shared-centre pairs fold the mode on base's centres, so it should carry over.

**Sources:** [`../progress/2026-10-05_yio_ws_tolerance_full_group.md`](../progress/2026-10-05_yio_ws_tolerance_full_group.md);
scripts in `calculations/diagnostics/2026-09-29_centre_symmetrization/` and `.../2026-10-05_ws_tolerance_window/`.
