# WS tie tolerance in YIO base: it matters a lot, but it isn't what's left of the symmetry floor (2026-10-05)

Follows [`2026-09-29_yio_qe76_new_ns_nc_pair.md`](2026-09-29_yio_qe76_new_ns_nc_pair.md) §1, which
left the A(R) symmetry floor open. Writes up three diagnostic runs from Sep 30 to Oct 1 that were
never put in a note (logs only in a job tmp dir), and adds today's full-group runs at the old
tolerance and at the top of the physically allowed range.

**Bottom line.** The WS tie tolerance in `_map_R_by_orbital_positions` (used for both the H(R) and A(R)
folds) controls much of the floor. Raising it from 1e-5 to 1e-3 Å cuts the residual about 4x for
**every** operation of base's m-3m' group, not just inversion: median total 6.6% → 1.6%, median
external 90% → 26%. But the physically allowed range is **[2.05e-3, 2.12e-2) Å**, and any value in it
gives the identical fold (§4). At 5e-3, every tie that symmetry forces is caught. Inversion then falls to
0.18% total / 0.9% external, while the other 46 operations stay at **1.55% total / 26% external
(median)**, essentially the 1e-3 values. **The rest of the floor isn't a tie-breaking problem.**

## 0. Setup

- Base: `data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889` (the
  QE 7.6 + `new_ns_nc` pair from the 09-29 note), `n_occ = 156`, `mp_grid = 8 8 8`.
- Position source `right_centre_from_chk` (`position_matrix_mmn_formula` on the `.chk`'s stored M
  matrices; no `.mmn` file, which is deleted). H(R) and A(R) are both folded on the `.chk` centres by
  `_map_R_by_orbital_positions` at the same tolerance. The module default
  `_DEFAULT_WS_DISTANCE_TOLERANCE` (`wannier/position_matrix.py:501`) is **now 1e-3 Å**, was 1e-5. The
  file is untracked, so git doesn't record when the default changed; its mtime is 2026-09-30 19:40.
  The Fortran-side `ws_distance_tol` in every YIO `.win` is still Wannier90's default (1e-5); the
  `jaemo` route (fork-written `_r.dat`) never goes through the Python fold, so its numbers in
  `validation/results/magnetic_symmetry_jaemo_qe76/` can't change with this setting.
- Metric is the same as `check_magnetic_symmetry.ipynb`: 7 generic k-points (rng 20260825), relative
  residual = max |Ω(gk) − predicted| / max |Ω|, of the trace-vector of each curvature part.
- Magnetic point group from `derive_point_groups.py` (from the DFT atom positions and the converged
  AIAO Ir moments, not recalled): base m-3m' (48 ops), modes 1-3 -4'3m' (24 ops). `point_groups.npz`
  stores only the rotations. The Fd-3m operations carry quarter-lattice fractional translations at
  the `.win` origin; these drop out of Ω(k), but §4's centre check has to restore them.

Files:
- Sep 30 – Oct 1 scripts: `calculations/diagnostics/2026-09-29_centre_symmetrization/`. Their logs,
  `point_groups.npz` and `derive_point_groups.py` are copied into its `logs/` (originally under a
  `~/.claude/jobs/` tmp dir).
- Today's scripts, logs and `base_full_group_tolerances.npz`:
  `calculations/diagnostics/2026-10-05_ws_tolerance_window/`.

## 1. Why the tolerance was widened (Sep 30)

`check_ws_tie_consistency.py`: the closest-image choice at 1e-5 Å was inversion-*inconsistent*
14% of the time when the two nearest images are within 1e-3 Å of each other. These are near-ties the
tight tolerance doesn't treat as ties, so one image gets all the weight instead of a symmetric split.
Wider candidate gaps were 0% inconsistent.

`check_tie_tolerance_widened.py` (`logs/tie_tolerance_widened.log`), base inversion residual, max
over on-mesh points, total curvature:

| route | 1e-5 | 1e-3 | 1e-2 |
|---|---|---|---|
| `right_centre_from_chk` | 0.120 | 0.0037 | 0.0021 |
| `transl_inv_full_from_chk` | 0.141 | 0.0034 | 0.0023 |

The external part goes from 0.75-1.23 to 0.02-0.03. This is the result that read as "the
symmetry problem is fixed". It was only inversion.

## 2. The full group at 1e-3 (Oct 1, reproduced 2026-10-05)

`full_symmetry_tables.py` (`logs/full_symmetry_tables.log`). A base-only rerun today
(`base_full_group_widened.py`, fresh build from the `.chk`, pythtb `external` @ `9db1288`, newer
than the Oct 1 run) reproduces the Oct 1 base table **digit for digit, all 48 rows**. The S4z control
(not a symmetry of base) gives ~199% (`logs/base_symmetry_widened.log`), so the test still tells a real
symmetry from noise. The modes (-4'3m', Oct 1) look like base: total median 1.46% (mode 1), 1.30%
(mode 2), 2.33% (mode 3); external 4.8-34%.

## 3. Full group at 1e-5, 1e-3 and 5e-3 (2026-10-05)

`base_full_group_tolerances.py` passes the tolerance explicitly to both folds. Comparison table:
`compare_tol.py`. Median over the 47 non-identity operations (max in brackets):

| part | 1e-5 | 1e-3 | 5e-3 |
|---|---|---|---|
| internal (H(R) only) | 2.25% (3.92) | 1.06% (1.41) | 1.06% (1.31) |
| cross | 36.5% (51.6) | 16.0% (23.4) | 16.3% (21.3) |
| external | 89.6% (153) | 25.9% (36.8) | 25.6% (33.6) |
| total | 6.55% (12.3) | 1.57% (2.42) | 1.55% (1.99) |

By operation class (median over the class), total and external:

| class | priming | n | total 1e-5 / 1e-3 / 5e-3 | external 1e-5 / 1e-3 / 5e-3 |
|---|---|---|---|---|
| inversion | unprimed | 1 | 10.5 / 0.24 / **0.18** | 151 / 3.3 / **0.9** |
| C2 ⟨110⟩ | primed | 6 | 5.73 / 1.03 / 0.90 | 91.5 / 19.0 / 17.2 |
| I·C2 ⟨110⟩ | primed | 6 | 7.85 / 1.06 / 0.89 | 91.2 / 19.2 / 17.3 |
| C3 ⟨111⟩ | unprimed | 8 | 6.05 / 1.54 / 1.51 | 86.3 / 26.0 / 25.9 |
| I·C3 (S6) | unprimed | 8 | 6.63 / 1.58 / 1.52 | 92.2 / 26.3 / 25.9 |
| C4 | primed | 6 | 6.42 / 1.58 / 1.55 | 85.0 / 28.6 / 26.2 |
| I·C4 (S4) | primed | 6 | 6.57 / 1.58 / 1.56 | 84.7 / 28.4 / 26.4 |
| C2 ⟨100⟩ | unprimed | 3 | 4.92 / 1.91 / 1.76 | 80.8 / 30.4 / 30.5 |
| I·C2 ⟨100⟩ | unprimed | 3 | 8.51 / 1.83 / 1.73 | 131 / 32.3 / 31.2 |

- **1e-5 → 1e-3** helps every class about 4x (external about 3.5x). It's not just inversion.
- **1e-3 → 5e-3** fixes the rest of inversion and changes everything else by only a few percent (median
  per-operation ratio 0.91 total, 0.96 external; 7 of 47 operations get slightly worse).
- At 1e-5, g and I·g differ (e.g. C2⟨100⟩ 4.9% vs I·C2⟨100⟩ 8.5%) because inversion itself is broken.
  From 1e-3 on they agree, so what's left is set by the proper rotation. ⟨110⟩ two-folds are the
  cleanest at ~0.9% total, ⟨100⟩ two-folds the worst at ~1.75%.

## 4. The physical tolerance window

The tolerance merges images whose distances `|R + L + τ_r − τ_l|` (`L` an 8×8×8 supercell vector)
differ by ≤ tol. This is only legitimate for images that symmetry makes exactly equidistant and that
differ only through noise in the centres.

- **How symmetric the centres are** (`centre_defect.py`): set-wise defect `max_w min_{w',L} |g c_w + t_g
  − c_w' − L|` over all 48 operations, with each operation's fractional translation `t_g` recovered from
  the atoms (atom basis then maps to 3e-15 Å). **Max 5.7e-4 Å** (all O functions; Y and Ir centres sit
  exactly on their atoms, while O centres are up to 0.069 Å off-site, symmetrically). A tie forced by
  symmetry can therefore be broken by at most ~4 × 5.7e-4 ≈ 2.3e-3 Å.
- **Distribution of excess distances** (`gap_histogram.py`, `all_image_window.py`, all 34M
  (left, right, R) entries × 27 images): a noise cluster from ~1e-7 up to **2.046e-3 Å**, then
  **nothing** until **2.125e-2 Å**, at any image rank. The cluster's top sits just inside the 2.3e-3 bound.
- **Upper side** (`tie_window.py`): entries that aren't even ties for atom-pinned centres have gaps
  ≥ 0.32 Å. The 2e-2 to 0.1 Å population consists of ties for atom-pinned centres that are genuinely broken
  by the symmetric O off-site displacements, so they must not be merged.

So the allowed range is **[2.05e-3, 2.12e-2) Å**, and every tolerance in it gives the same fold for
base. 1e-3 is just below it: it misses 34 k forced-tie entries in (1e-3, 2.05e-3]. Those turned out to
carry the remaining inversion residual. 1e-5 misses ~1.7 M. Sep 30's 1e-2 test was inside the range.

## 5. What this leaves

- **What remains is not WS tie-breaking** (traced to the DFT NSCF Hamiltonian: [`2026-10-05_yio_symmetry_floor_dft_origin.md`](2026-10-05_yio_symmetry_floor_dft_origin.md)). At 5e-3 every forced tie is split symmetrically, and no
  further tolerance changes the fold. The ~1.5% total / ~26% external rotation residual must be in the
  matrix elements H(R), A(R) themselves, or in how they transform (orbital mixing under rotations).
  Inversion only maps each orbital to ± itself, and it's now clean. That fits the 09-29 note §3's
  orbital-pairing issue and the `zaxis`/`xaxis`-free projections, but neither has been tested.
- **The production default should probably be 5e-3**, the middle of the window instead of just below
  it. Not changed: it would also move dθ, which hasn't been measured at either 1e-3 or 5e-3 (next
  point). The window is base's. Shared-centre pairs fold the mode on base's centres, so it should carry
  over, but that hasn't been checked for the modes' own centres.
- **dθ at the new tolerance.** The 09-29 dθ (−0.0175 rad/Å, `transl_inv_full_from_chk`) was computed
  before the default changed. A Sep 30 `right_centre_from_chk` + shared-WS rerun (job log
  `mode1_rerun.log`, output `nk_sweep_chkmmn_shared_ws/`, saved 20:11) gives **−0.01856 / −0.01864**
  at nk = 12 (modes 2, 3: +0.0439/+0.0456 and +0.1697/+0.1787). It's ~6% from the 09-29 value, but the
  log doesn't record the tolerance and the run may have started before the 19:40 change, so it isn't
  a clean before/after.
  **Settled 2026-10-07:** rebuilding base H(k) at both tolerances, the 09-30 runs' stored base gaps match
  1e-3 to 1e-7 eV (nk 4 and 6 minimum-gap points) and miss 1e-5 by 7e-5 to 3e-4 eV, so all three 09-30 runs
  used 1e-3. They still change the A(R) source at the same time (tif → right_centre), so the comparison
  isn't a pure tolerance effect.
