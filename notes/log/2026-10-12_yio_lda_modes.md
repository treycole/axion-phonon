# 2026-10-12 — Y₂Ir₂O₇ LDA modes 1-3: SCF/NSCF checks and the Wannier chain

The user chose LDA+U for YIO (2026-10-12): the base curvature symmetry is exact there with the bug #22 fix, and the
PBE mode-2 SCF never converged ([log 10-11 §5](2026-10-11_production_binary_ws_cleanup.md)). This note covers the
three Γ₂⁻ mode runs in LDA and their Wannierization.

## 1. Setup

`Y2Ir2O7/phonon/mode_{1,2,3}/fromYL2_2026-10-12/` on the cluster (`phonon/setup_modes_fromYL2.sh`):
- **Inputs:** the YL2 base `scf.in`/`nscf.in` with only ATOMIC_POSITIONS replaced (Q = 0.01 Å; mode 1 and 3 positions
  from `fromYP_2026-10-08`, mode 2 from `Q_0p01_from_Q0p005_2026-10-11`). The SCF adds `tprnfor` and
  `startingpot = 'file'` (YL2 density, Hubbard ns, PAW becsum).
- **Run:** pw.x `db4b99b4`; the NSCF (512 k, `-nk 8`) runs only if the SCF converges.
- **Jobs:** 212755 / 212756 / 212757.

## 2. Results

| mode (job) | atoms moved | iterations | E − E(YL2) | −½ F·Δr | ratio | SCF correction | ops |
|---|---|---|---|---|---|---|---|
| 1 (212755), Y | 1-4 | 13 | +1.6133 mRy | +1.6142 mRy | 0.999 | 2.6e-4 | 24 |
| 2 (212756), Ir | 5-8 | 15 | +1.8093 mRy | +1.8096 mRy | 1.000 | 1.2e-4 | 24 |
| 3 (212757), O 48f | 11-22 | 11 | +0.5528 mRy | +0.5486 mRy | 1.008 | 6.7e-4 | 24 |

- **Base:** E(YL2) = −4499.85691362 Ry. F_base·Δr = 0, since the mode is not totally symmetric.
- **Magnetization:** total 0 and |m|abs 3.21 for all three, the same as the base.
- **Ir moments:** |m_i| = 0.2326 μB per Cartesian component (base 0.2330/0.2333).
- **Moment pairs:** the printed Ir1,2 vs Ir3,4 values differ by ≤ 3e-5 in the modes and by 3e-4 in the base, which has
  48 operations. The difference comes from the moment printout, not from the symmetry.
- **Symmetry:** QE finds the 24 operations of −4′3m′ in every mode.
- **NSCF:** 512 k, JOB DONE; HOMO/LUMO 10.443/11.128 eV.

### Window edges

The base 212615 edges are in eigenvalue gaps that are open at every k, for every mode:

| | dis_froz_max 14.396 | dis_win_max 24.522 |
|---|---|---|
| base 210590 | gap 14.374..14.418 | 24.502..24.542 |
| mode 1 | 14.373..14.416 | 24.502..24.542 |
| mode 2 | 14.371..14.405 (9 meV) | 24.502..24.542 |
| mode 3 | 14.373..14.413 | 24.501..24.542 |

Below each edge, every mode has the same band counts as the base at every k: 176-178 bands below 14.396 and 260-268
below 24.522. `dis_win_min` = −14 sits in a 1.8 eV gap. So the mode `.win` keeps the base windows unchanged
(`win_gaps.py` in the session scratchpad; the parse also counted the 300 occupation numbers, which are subtracted
here).

## 3. Wannier chain (`phonon/setup_w90_modes_YL2.sh`)

**Prep jobs 212777 / 212779 / 212781**, run in each mode's run folder:
- rename the 13 MB start tarball to `start_from_YL2.tar.gz`;
- gunzip the NSCF save in place (the untarred `tmp/` the Julia code reads);
- tar the save as the run's one `tmp.tar.gz`.

**W90 jobs 212778 / 212780 / 212782**, each held on its prep job:
1. untar that `tmp.tar.gz` on the node;
2. `wannier90 -pp`;
3. `generate_amn.jl --no-displacement` from the base 212615 `.chk`/`.nnkp` (now gunzipped) and the base save
   210590;
4. pw2wannier90 `8f5368ba -nk 8` (mmn and eig only);
5. wannier90 `26c8b506 -np 1`.

**Inputs and files:**
- `.win` = the base `.win` with each mode's atoms_cart. Checked: only the displaced atoms differ.
- There is no tarball copy in the w90 folders. The PBE chain left one in every w90 folder.

## 4. Next

- **Mode symmetry check:** `archive/diagnostics/2026-10-12_yio_lda_modes/mode_full_group.py`. It checks the 24
  operations of −4′3m′ at 7 generic k with H from `.chk`, no hopping cut, tol 5e-3, folded on the base centres.
- **Then dθ:** with `run_axion.py`, after the bug #22 production decision.

## 5. The w90 jobs failed: no Julia on the compute nodes (16:00)

**Failure:** all three w90 jobs (212778 / 212780 / 212782) stopped at `julia: command not found`.

**Cause:** the compute nodes mount a different `/mnt/home` from the login node. The nodes see `fs:/homel` (938 GB);
the login node sees `192.168.101.200:/homel` (1.5 TB). `~/.juliaup` exists only on the login node.

**Earlier chains:**
- The PBE mode chain failed the same way (210722, 210723), so no PBE mode Wannierization ever existed.
- The September pairs (e.g. 192347) generated the `.amn` on the login node (`generate_rigid_shift_amn.sh`) before the
  job ran.

**Fix:** `phonon/amn_login_YL2.sh`, run under nohup on the login node, does one mode after another:
1. generates the `.amn` with `/mnt/home/bc889/.juliaup/bin/julia`;
2. resubmits the w90 job without the Julia line (`send_job_w90.sh`).

The prep jobs had worked: each run folder has its 50.6 GB `tmp.tar.gz` and a gunzipped `tmp/`.

## 6. Production switch: H(R) from the `.chk`, no hopping cut (user, 17:00)

The bug #22 fix is now the production default:
- **`wannier/wannier_io.py`:** new `hamiltonian_classes_from_checkpoint`, plus
  `load_wannier_pair(hamiltonian_source=None → "chk")`. The exception is `rdat_ndegen_applied`, which keeps
  `"hr_dat"`. The model options set `min_hopping_norm=None`.
- **`run_axion.py`:** settings `H_SOURCE = "chk"` and `MIN_HOPPING_NORM = None`. The old settings write to a folder
  with an `_hrdat`/`_cut` suffix, and both settings go into the results metadata.
- **Docstring examples** in `modules/axion.py` and `modules/curvature.py` use `min_hopping_norm=None`.

**Checks on base 212615 (`load_wannier_pair`, 12 mesh k, all frozen states):**
- Class sums differ from `_hr.dat` by at most 3.4e-6 eV, which is print rounding.
- With the `.chk` H, the Wannier eigenvalues match the DFT frozen states to **1.5e-13 eV**; the `_hr.dat` H gives
  3.5e-5 eV.
- The base/mode WS images are identical (0 entries changed).
- `validation/unit`: 55 passed.

**Environment:** `~/Repos/pythtb` is now on `dev` (a4e6608), which has no `W90.model(orb_vecs=...)`. Run production
code with the `external` worktree first on the path (`PYTHONPATH=<worktree>`); the session's worktree is
`scratchpad/pythtb-external` at 9db1288.

## 7. PBE modes 1 and 3 Wannierized for the PBE vs LDA comparison (user, 18:15)

**What the earlier PBE symmetry checks used:** the base Wannierizations (210600, then 212616 with edges in gaps). A
base needs only pw2wannier90's own projections, with no Julia. The PBE **mode** Wannierizations never happened: their
rigid-shift `.amn` needs Julia, which failed (§5).

**Mode 2 is left out:** it has no converged PBE SCF.
- 210597 crashed.
- 210623 (local-TF restart) ended at 6e-5 Ry, but its NSCF ran anyway, so that NSCF is not usable.
- 212707 stalled ([log 10-11 §5](2026-10-11_production_binary_ws_cleanup.md)).

**Window edges:** the base 212616 edges (dis_froz_max 14.858, dis_win_max 25.155) are in global gaps of the mode 1
(210596) and mode 3 (210598) NSCFs. The closest band is 16 meV from an edge, and the band counts match the base at
every k.

**Chain (`phonon/setup_w90_modes_YP_frozfix.sh`):**
- **Tarballs:** each failed chain's 48 GB NSCF-save tarball was moved into its run folder as the run's one
  `tmp.tar.gz`. The 13 MB start tarball is kept as `start_from_YP.tar.gz`.
- **Prep 212983 / 212984:** gunzip each run's `tmp/` save.
- **`.amn`:** `phonon/amn_login_YP.sh` makes it on the login node after the LDA ones finish, then submits
  `send_job_w90.sh`.
- **Folders:** `mode_{1,3}/fromYP_2026-10-08/w90_trial02_no-displacement_frozfix/`.
- **Base folder:** the base `.chk` was gunzipped, and the `.nnkp`/`.eig` copied into 212616bc889. Julia and the
  `.chk` H need them there.
- **Not removed:** each run folder still has a 48 GB `tmp_bands` copy, which the cleanup script deletes.

## 8. MBT base Wannier bands (`analysis/bands.ipynb`)

The notebook now runs with `BANDS_OUT = None`. It draws:
- the Wannier bands on the `bands.in` path;
- the base NSCF's `.eig` (all 300 bands) as dots at the 8 path points that are on the 4×4×4 mesh;
- the `.win` windows as lines.

There is no DFT band-structure run for the LDA base, only the 9-k Γ–L check 210652.

**MBT LDA base 210651** (E_F = 7.406 eV):
- The Wannier bands span about −13.5 to +4.5 eV relative to E_F, i.e. ≤ ~12 eV absolute.
- The frozen window ends at +2.6 (10 eV absolute); the outer window at +12.6 (20 eV).
- A `dis_win_max` of 13 eV absolute (+5.6) is above the top Wannier band everywhere on the path.

## 9. MBT outer window ~15 eV (user, 18:30)

**User's decisions:**
- Compare PBE and LDA on YIO modes 1 and 3 only.
- Put the MBT outer window top near 15 eV, so the NSCF can use fewer bands.

**What the base `.eig` gives** (210622/210651, 4×4×4, 300 bands; 56 bands below −7 eV):

| top | bands below, per k | in [−7, top], max |
|---|---|---|
| 13 | 184-192 | 136 |
| 15 | 196-222 | 166 |
| 20 | 266-288 | 232 |

**Window edges:**
- **15.047:** the centre of the global gap 15.020..15.073 (53 meV).
- **`dis_froz_max = 10`:** sits 5 meV from a band (gap 9.995..10.023), so it moves to 10.009.
- **`dis_win_min = −7`:** in a 1.1 eV gap.

**Is the base Wannierization affected?** Yes: the base subspace changes. The run is projection-only and the frozen
window holds 104-112 states per k (160-168 below 10 eV minus the 56 semicore states) for 116 WFs. The remaining WFs
come from the outer-window states, and the outer window then contains 140-166 states instead of 210-232. So the base
needs a new wannier90 run, which is cheap. The modes must use the same windows, with every edge in a gap of each
mode NSCF.

**Test job 212985** (`w90_trial_sp_d_recentred/window_o15_f10_2026-10-12`): the recentred `.amn` plus the 210622
`.mmn`/`.eig`, wannier90 only. The 2026-10-09 window scan (11.5-13 eV) used the pre-recentring `.amn` (#20), so its
results don't carry over.

**NSCF band count:** about 240 for the denser mesh, to check against the denser mesh's own band counts.

**Fix history:** before/after for each fix is now one table in `notes/bugs/README.md`.

## 10. MBT: 15 eV window test, then 8×8×8 NSCFs (user, 19:15)

**Window test 212985** (wannier90 only, same `.amn`/`.mmn`/`.eig` as 210651):

| | outer top 15.047, frozen top 10.009 | outer top 20, frozen top 10 (210651) |
|---|---|---|
| Ω total | 345.5 Å² | 382.2 Å² |
| largest spread | 5.67 Å² | 6.24 Å² |
| path gap, indirect / direct | 105 / 189 meV | 104 /185 meV |
| 12³ grid gap, indirect / direct | 119 / 198 meV | 118 / 198 meV |

- **DFT gap on the 4×4×4 mesh:** 198 meV.
- **Bands** (`analysis/bands.ipynb`, both runs drawn): identical below the frozen top, with small
  differences in the top empty Wannier bands.
- **Conclusion:** the narrower window loses nothing and localizes the functions better. The remaining
  indirect-gap error is the 4×4×4 off-mesh error.

**User's decision: 8×8×8 mesh.**

**NSCFs** (`MnBi2Te4/setup_nscf_k8_lda.sh`), in `nscf_k8_2026-10-12/` next to each SCF run:
- base 209653 → 212986;
- T₂⁻ modes 1-4 (SCFs 212708-11) → 212987-212990.

**Inputs:** each structure's own `nscf.in`, with only `nbnd` 300 → 240 and K_POINTS → the 512-k list in Wannier90's
order. Each starts from its own SCF state (density, xml, occupations, PAW), with pw.x `db4b99b4`, `-nk 8` and 48
slots.

**Expected cost:** about 1.5 days each (the base 4×4×4 NSCF took 5h49m on 64 processes with 300 bands) and
~250-350 GB per structure with the tarball. The cluster is already over target (~5 TB vs 2 TB) until the cleanup runs.

**Correction (23:35, measured):** that estimate was wrong. Each pool takes ~6,300 s per k point (Davidson, 240
bands, 6 processes per pool), and each pool has 64 k points, so a run needs ~112 h. The jobs were submitted with
h_rt = 96 h, so all five would be killed at ~86% of k points unless the limit is raised (queue maximum 248 h).
**Resubmitted 2026-10-13 00:50 (user OK):** the five jobs were killed after 5.5 h (3 of 64 k points per pool) and
resubmitted unchanged except for `h_rt=248:00:00`, the queue maximum: base → 212997, modes 1-4 → 212998-213001. On
the kill, the epilog copied each job's partial scratch (~20 GB, 98 GB total) back into its run folder. Those folders
were moved to `/mnt/wkdk/bc889/axion/_killed_2026-10-13/` so the new jobs don't stage them, then deleted (96 GB, user
OK).

**Next, after each NSCF:**
1. Check that band 240 stays above 15.047 eV at every k, and that 15.047 / 10.009 / −7 are in gaps of every NSCF.
2. Tar the save as the run's one `tmp.tar.gz` and gunzip `tmp/`.
3. Base: pw2wannier90 (amn, mmn, eig) → wannier90, with the recentred `.win` plus num_bands 240, mp_grid 8 8 8 and the
   new windows.
4. Modes: pw2wannier90 (mmn, eig) and the `--no-displacement` `.amn` from the base `.chk` on the login node →
   wannier90.
5. Add the pinned-position option (#21) to `run_axion.py` → mode symmetry checks → dθ.

## 11. Pinned positions in production; methods draft and meeting note (user, 21:45)

**Pinning (#21) is a production option.**
- **Code:** `load_wannier_pair(..., pin_to_atoms=True)`, set by `PIN_TO_ATOMS` in `run_axion.py`. Off for YIO, on
  for MBT; the output folder gets `_pinned`, and the metadata records `pinned_to_atoms`.
- **Positions:** each base Wannier centre is replaced by its nearest atom image (atoms from the base `.win`,
  `atoms_cart` or `atoms_frac`), in the H and A folds of both structures and in τ. The A(R) formula itself keeps
  the run's own centres.
- **Requirements:** `shared_ws_centres=True` and a `.chk` A(R) source. A centre more than 1 Å from every atom is
  refused.
- **H(R) class representative:** now the shortest member of each class (it was the first member). The fold searches
  ±1 supercell around the representative, and in MBT's skewed rhombohedral cell that covers every candidate image
  only near the origin. The MBT diagnostic already did this.
- **Unit tests:** 59 pass (4 new in `validation/unit/test_shared_ws_centres.py`).
- **Real-data check** (`archive/diagnostics/2026-10-11_ws_positions/pinned_loader_check.py`, `.log`):

  | | result |
  |---|---|
  | MBT 210651, production loader, centres | bands break −3m′ by 64.042 meV (12 ops × 7 generic k) |
  | MBT 210651, production loader, pinned | **0.000 meV**, the same as `mbt_eig_symmetry_fold.py` |
  | YIO 212615, first vs shortest representative | 184/512 classes change member; folded H identical (0.0 eV) |

- **To check on 8×8×8:** the longest MBT hopping is 55.5 Å on 4×4×4 and roughly doubles on 8×8×8, close to the
  loader's `max_distance = 125` Å.

**Write-ups.**
- **Methods:** `paper/methods_draft.md` is the draft of the paper's computational sections. It is not in the `.tex`
  yet, and has a "what changes relative to main.tex" table.
- **Meeting note:** `meetings/2026-10-12.md` explains the symmetry tests (operations, the three tests, what the
  residual percentage means), the YIO and MBT bug stories, and the open questions.

**YIO mode Wannier jobs:**
- the `.amn` for all three modes is done (login node, `--no-displacement`);
- w90 212982 (mode 1), 212994 (mode 2) and 212995 (mode 3) are running;
- the PBE `.amn` script (`amn_login_YP.sh`) waits for the LDA log's "all done", which appeared at 21:31.

## 12. Mode Wannier models and their symmetry (2026-10-13)

**Wannier90:** w90 212982 / 212994 / 212995 finished. Every `.amn` was made with `--no-displacement` (checked in each
`generate_amn.log`). Total spread 190.68 / 190.38 / 190.33 Å² (base 193.76). Downloaded to
`data/Y2Ir2O7/phonon/GM2-/mode{1,2,3}/Q_0.01A/output/YL2_21275{5,6,7}/w90_trial02_no-displacement/<job>bc889/`.

**DFT eigenvalues** (`.eig`, 512 k, lowest 260 bands) under the 24 operations of −4′3m′: ≤ 0.001 meV for every mode,
the print precision, the same as the base.

**Curvature table** (`archive/diagnostics/2026-10-12_yio_lda_modes/mode_full_group.py`: H from `.chk`, no cut,
tol 5e-3, `transl_inv_full_from_chk`, fold and τ on the base 212615 centres), median (max) over 23 operations, %:

| mode | internal | cross | external | total |
|---|---|---|---|---|
| 1 (212982) | 0.008 (0.008) | 0.222 (0.227) | 0.011 (0.014) | 0.000 (0.001) |
| 1, own centres | 0.009 (0.009) | 0.220 (0.233) | 0.058 (0.072) | 0.000 (0.001) |
| 2 (212994) | 0.008 (0.009) | 0.207 (0.214) | 0.004 (0.007) | 0.000 (0.001) |
| 3 (212995) | 0.008 (0.008) | 0.209 (0.213) | 0.013 (0.018) | 0.000 (0.000) |
| base 212615 (48 ops, same settings) | 0.008 (0.008) | 0.220 (0.224) | 0.012 (0.020) | 0.000 (0.000) |

- **The modes reproduce the base.** The total is exact. The cross part is 0.21-0.23% on the 16 operations that
  don't keep the z axis and ≤ 0.009% on the 7 that do.
- **Why the parts are not exact:** the split into internal / cross / external depends on τ, and τ is not exactly
  covariant. With spinor projections along z, the operations that don't keep z mix same-site Wannier functions whose
  centres differ by up to ~9e-4 Å; the measured centre defect is 5e-4 Å on exactly those operations (both base and
  mode centres; ≤ 1e-7 Å on the others). The total does not depend on τ.
- **Not the shared-centre fold:** folding the mode on its own centres leaves the cross part at 0.22%.
- **The base's "0.000%"** (status page, bug #22) is the total. `symmetry_comparison_table.md` has always shown the
  cross part at 0.22%.

## 13. dθ/dQ, LDA modes (2026-10-13)

`archive/code/run_axion.py` (now pointing at base 212615 and mode 1 212982): `transl_inv_full_from_chk`, SHARED_WS,
H from `.chk`, no cut, unpinned, frozen gauge (A_β = 0), Q = 0.01 Å, n_occ 156. Output in each mode folder's
`nk_sweep_tif_shared_ws/`; modes 2 and 3 ran from copies of the driver with only MODE changed.

**Mode 1 (212982), rad/Å** (base-end, mode-end estimate; minimum gap 0.956 eV):

| nk | base end | mode end |
|---|---|---|
| 4 | −0.01953 | −0.01945 |
| 6 | −0.01828 | −0.01841 |
| 8 | −0.01761 | −0.01778 |
| 10 | −0.01799 | −0.01811 |
| 12 | −0.01798 | −0.01810 |

Converged to 0.05% between nk 10 and 12. The two estimates differ by 0.7%; their mean is **−0.0180 rad/Å**.

**Modes 2 (212994) and 3 (212995), rad/Å** (base end / mode end; minimum gap 0.954 / 0.953 eV):

| nk | mode 2 | mode 3 |
|---|---|---|
| 4 | +0.02972 / +0.03127 | +0.06069 / +0.06256 |
| 6 | +0.05357 / +0.05419 | +0.06553 / +0.06641 |
| 8 | +0.05288 / +0.05351 | +0.06563 / +0.06656 |
| 10 | +0.05383 / +0.05450 | +0.06502 / +0.06596 |
| 12 | +0.05384 / +0.05450 | +0.06502 / +0.06596 |

**Summary (nk 12, mean of the two estimates):**

| mode | atoms | dθ/dQ (rad/Å) | estimates differ by |
|---|---|---|---|
| 1 | Y | −0.0180 | 0.7% |
| 2 | Ir | +0.0542 | 1.2% |
| 3 | O 48f | +0.0655 | 1.4% |

- All three are converged to ≤ 0.05% between nk 10 and 12.
- Q is the amplitude convention of the SMODES displacement files (open question: per-atom displacement or norm of
  the full displacement vector; meeting 10-12 §8).

## 14. `_r.dat` no longer needed for the `.chk` sources (2026-10-13, user OK to delete)

The laptop disk reached 97-99% full. With `transl_inv_full_from_chk` / `right_centre_from_chk` the A(R) comes from
the `.chk`, but the loader still demanded `_r.dat`: PythTB's A(R) was a placeholder replaced before use, and its
R=0 diagonal served as the check that the `.chk` belongs to the same run.

- **Change:** `load_wannier_pair` no longer requires `_r.dat` for the two `.chk` sources (it requires
  `_centres.xyz` instead). `_position_matrix_for_model(..., reference_centres=)` does the same-run check against
  `_centres.xyz` (8 decimals; it matches the `.chk` to 5e-9 Å on all three modes) when the model has no `_r.dat`.
  The `.mmn` and `rdat` sources are unchanged.
- **Check:** mode 1 at nk 4 without `_r.dat` gives −0.01953139 / −0.01945369, identical to every printed digit;
  unit tests 59 passed.
- **Deleted:** `_r.dat` (3.3 GB each) of LDA modes 1-3. The PBE mode downloads (212996, 213002) skip it. The base
  folders keep theirs.


## 15. PBE vs LDA, modes 1 and 3 (2026-10-13)

The PBE modes 1 (212996) and 3 (213002) were Wannierized against the PBE base 212616 (`_frozfix` windows) and run with
the production settings of §13 (`transl_inv_full_from_chk`, SHARED_WS, H from `.chk`, no cut, unpinned, frozen gauge,
Q = 0.01 Å, n_occ 156). Scripts: `scratchpad` copies of `run_axion.py`; logs `run_axion_2026-10-13.log` in each mode
folder and `archive/diagnostics/2026-10-12_yio_lda_modes/pbe_mode{1,3}_full_group_*.log`.

**Symmetry (−4′3m′, 24 ops, 7 k), median (max), excluding the identity:**

| run | internal | cross | external | total |
|---|---|---|---|---|
| PBE base 212616 (m-3m′) | | 0.40% | 0.25% | 0.01% (0.02%) |
| PBE mode 1 212996 | 0.012% | 0.395% | 0.209% | 0.011% (0.016%) |
| PBE mode 3 213002 | 0.012% | 0.358% | 0.143% | 0.011% (0.015%) |
| LDA mode 1 212982 | 0.008% | 0.222% | 0.011% | 0.000% (0.001%) |

The PBE modes match the PBE base, so they are as good as PBE gets. PBE's 0.01% total floor is the residual
symmetry breaking of the PBE DFT (#16 fixed gives 46 μeV there vs 1.2 μeV in LDA).

**dθ/dQ, rad/Å** (base end / mode end):

| nk | PBE mode 1 | PBE mode 3 |
|---|---|---|
| 4 | −0.01760 / −0.01754 | +0.05465 / +0.05607 |
| 6 | −0.01775 / −0.01784 | +0.05874 / +0.05942 |
| 8 | −0.01731 / −0.01744 | +0.05851 / +0.05924 |
| 10 | −0.01752 / −0.01762 | +0.05810 / +0.05884 |
| 12 | −0.01752 / −0.01761 | +0.05810 / +0.05884 |

Converged to 0.04% (mode 1) and 0.01% (mode 3) between nk 10 and 12; the estimates differ by 0.5% / 1.3%.

**Comparison (nk 12, mean of the two estimates):**

| mode | LDA (§13) | PBE now | PBE/LDA | PBE before the fixes (09-30/10-05, tol 1e-3) |
|---|---|---|---|---|
| 1 (Y) | −0.01804 | −0.01756 | 0.973 | −0.01753 |
| 3 (O 48f) | +0.06549 | +0.05847 | 0.893 | +0.1689 |

- **Functional:** LDA vs PBE moves mode 1 by 2.7% and mode 3 by 11%, same sign. The answer is not sensitive to
  the functional at the level that matters for the paper.
- **The 2.6× drop in mode 3 was the fixes, not the switch to LDA:** PBE alone went from +0.1689 to +0.0585 (2.9×).
  Mode 1 barely moved under the fixes (−0.01753 → −0.01756). The pre-fix runs had the segni bug (#16), the
  rounded `_hr.dat` + hopping cut (#22) and the 1e-3 WS tolerance (#13); which of these carried mode 3 is not
  separated (not needed: all are superseded).
- PBE mode 2 has no converged SCF (§7), so there is no PBE number for it.

## 16. Γ₂⁻ dynamical matrix and physical modes, LDA (2026-10-13, first pass)

`archive/diagnostics/2026-10-12_yio_lda_modes/dynmat/gm2_dynmat.py` (inputs copied from the cluster SCFs).

- **Q convention, as applied:** Δr = Q·v with v the raw SMODES vector (entries ±1). Y and Ir move 0.0173 Å each
  (0.01 per Cartesian component), the 48f O 0.01 Å each; |Δr| = 0.0346 Å for every mode. It is not the normalized
  vector the paper and methods draft say.
- **Force constants** Φ_ij = −v_i·F(Q_j)/Q (eV/Å², Y, Ir, O): diagonal 439 / 492 / 149, off-diagonal Y-Ir 30,
  Y-O 44, Ir-O −63. Φ − Φᵀ is 0.3% of Φ (one-sided differences are O(Q²): no cubic Γ₂⁻ invariant).
- **Modes** (Φc = ω²Mc, cᵀMc = 1 amu; dθ/dq = Σ c_i dθ/dQ_i with the §13 values; signs of c arbitrary):

| ω (cm⁻¹) | Y / Ir / O share of \|Δr\|² | dθ/dq (rad Å⁻¹ amu⁻½) | dθ per Å of \|Δr\| |
|---|---|---|---|
| 229 | 0.04 / 0.67 / 0.29 | −0.00205 | −0.0240 |
| 328 | 0.80 / 0.02 / 0.18 | +0.00122 | +0.0108 |
| 471 | 0.01 / 0.00 / 0.99 | −0.00427 | −0.0176 |

- **Caveat:** the base is the experimental structure (a = 10.19 Å, x = 0.212932), not relaxed in LDA+U. Its forces
  have no Γ₂⁻ component by symmetry, but Φ is evaluated there. The base SCF did not print forces.
