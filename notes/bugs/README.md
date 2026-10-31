# Bugs

One file per bug found in this project, numbered in the order found. Each file opens with its **status**, then
describes the symptom, the root cause, why it was a bug, the fix, how the fix was verified, and what remains. The
evidence is in the dated notes in [`../log/`](../log/). The old date-by-date bug log (through 2026-10-08, with the
milestones between bugs) is kept in [`../archive/bug_log_to_2026-10-08.md`](../archive/bug_log_to_2026-10-08.md).

| # | Status | Found | Bug |
|---|---|---|---|
| [1](01_wannier_diffuse_projections_gap_collapse.md) | **FIXED** | 2026-07-27 | YIO Wannier model closed the gap (diffuse empty projections); n_occ 104 → 156 |
| [2](02_pythtb_embedding_noop.md) | **FIXED** | 2026-08-19 | PythTB embedding was a silent no-op; base and mode didn't share τ |
| [3](03_hr_ar_ws_convention_mismatch.md) | **FIXED** | 2026-08-20 | H(R) and A(R) used different Wigner–Seitz conventions |
| [4](04_stock_rdat_position_matrix.md) | **FIXED** (worked around) | 2026-08-27 | Stock `_r.dat` cannot carry the position matrix |
| [5](05_pr702_rdat_ignores_use_ws_distance.md) | **FIXED** (Jae-Mo's branch) | 2026-09-10 | PR 702's `_r.dat` silently ignores `use_ws_distance` |
| [6](06_base_mode_ws_tie_inconsistency.md) | **FIXED** (in-house folds); open for the Fortran `_r.dat` route | 2026-09-20 | Base and mode didn't share a Wigner–Seitz rule (A3′) |
| [7](07_qe_new_ns_nc_symmetrizer.md) | **FIXED** | 2026-09-21 | QE `new_ns_nc` symmetrizes Hubbard ns with the wrong atom |
| [8](08_born_charge_ionic_term.md) | **FIXED** | 2026-09-21 | Born effective charge: wrong ionic term |
| [9](09_mbt_t2_rounded_cell_and_labels.md) | **FIXED** (exact cell in the 10-11 T₂⁻ reruns; labels harmless with #7 patched) | 2026-09-22 | MBT T2⁻ inputs: rounded bohr cell loses C3; labels trigger #7 |
| [10](10_hr_ar_folded_on_different_centres.md) | **FIXED** | 2026-09-23 | H(R) and A(R) folded on slightly different centre arrays |
| [11](11_rdat_print_precision.md) | **FIXED** (mitigated) | 2026-09-23 | `_r.dat` 6-decimal print costs ~3% of dθ |
| [12](12_qe_paw_vrad_grad_uninitialized.md) | **FIXED** | 2026-09-27 | QE PAW `v_rad`/`g_rad` uninitialized (noncollinear GGA) |
| [13](13_ws_tie_tolerance_too_tight.md) | **FIXED** (YIO; default 5e-3 since 2026-10-11) | 2026-09-30 | WS tie tolerance too tight; 5e-3 Å sits in YIO's empty gap (1.7e-3 to 2.1e-2) |
| [14](14_mbt_stale_binary.md) | **FIXED** | 2026-10-04 | MBT u_3.0 ran a stale, buggy pw.x |
| [15](15_mbt_u3_base_rounded_cell.md) | **FIXED** | 2026-10-05 | MBT base u_3.0 rounded bohr cell (same defect as #9) |
| [16](16_qe_paw_segni_stale_ux.md) | **FIXED** (production reruns in progress) | 2026-10-05 | QE PAW GGA uses stale `ux` without `lsign` guard: root cause of the YIO symmetry floor |
| [17](17_segni_fix_base_scf_stall.md) | **OPEN** (candidate fix in testing) | 2026-10-07 | The segni fix (#16) stalls the YIO base SCF at ~3e-5 Ry |
| [18](18_qe_noncollinear_gga_spurious_state.md) | **WORKED AROUND** (MBT → LDA); open for YIO | 2026-10-08 | Noncollinear PBE has a spurious lower-energy, layer-asymmetric MBT state; random starts fall into it and stall |
| [19](19_qe_sym_ftrans_trev_energy.md) | **FIXED** (production pw.x `db4b99b4` since 2026-10-11; MBT T₂⁻ reruns 212708-11) | 2026-10-09 | QE 7.5CN/7.6 `irreducible_BZ` unfolds k with time reversal and drops −k stars for the MBT T₂⁻ group (13 k instead of 16): E 0.3 mRy too low, \|m\|abs 10.27. One-line patch |
| [20](20_w90_centre_wrap_long_cell.md) | **FIXED** (recentred sites) | 2026-10-09 | wannier90 centres wrap for sites at crystal x > 2/3 in the long MBT cell, so the WS images are wrong; gap −346 meV |
| [21](21_mbt_wf_centres_not_c3_covariant.md) | **FIXED** 2026-10-12 (`PIN_TO_ATOMS` in `run_axion.py`) | 2026-10-11 | MBT p_x/p_y WF centres sit off-axis, so a WS fold and τ on WF centres break C3: interpolated bands 62 meV asymmetric; pinned: exact |
| [22](22_hopping_cut_and_hr_rounding.md) | **FIXED** 2026-10-12 (production: H from `.chk`, no hopping cut) | 2026-10-11 | PythTB `min_hopping_norm=1e-5` and `_hr.dat`'s 1e-6 eV print break the spin-mixing operations: the whole YIO 0.5-0.9% internal floor; without them 0.000% |
| [23](23_generate_amn_no_paw_overlap.md) | **OPEN** (fix in progress: `wannier/paw_overlap.py`) | 2026-10-28 | `generate_amn.jl` overlaps PAW pseudo-wavefunctions without the augmentation term: the self-test `.amn` is O·V, not V (MBT \|A−V\|/\|V\| 0.22, Mn d norm 0.66); explains the self-test Ω_I shift; Λ = 0 no longer exact |

**Status meanings:**
- **FIXED:** the fix is in the code, builds and inputs used now.
- **(worked around) / (mitigated):** the defect still exists in the tool, but production doesn't depend on it.
- **PARTIALLY FIXED / OPEN:** see the file's "Remaining" section.

## What each fix changed (before → after)

Assembled 2026-10-12 from the bug files and the YIO base symmetry table
(`archive/diagnostics/2026-10-05_ws_tolerance_window/symmetry_comparison_table.md`). The curvature numbers are
the YIO base's relative residual under m-3m′ at 7 generic k, as the median over the 47 non-identity operations
(internal / external / total), with A(R) from `transl_inv_full_from_chk`. Each fix was measured on the runs that
existed at the time, so a row's "before" is the previous row's "after" only when the run is the same.

**Y₂Ir₂O₇**

| date | fix | measured on | before → after |
|---|---|---|---|
| 09-21 | [#7](07_qe_new_ns_nc_symmetrizer.md) QE `new_ns_nc` | DFT, old QE 7.5 runs → 192344 | Ir \|m\| spread 2.5% → 0.00%, tilts 0.7° → 0°; min gap 0.72 → 1.04 eV; H(R) symmetry residual 2-4× lower (the 4-14% curvature floor stayed: #16) |
| 09-30 | [#13](13_ws_tie_tolerance_too_tight.md) WS tie tolerance 1e-5 → 1e-3 Å | old PBE 192344 | base inversion 12-14% → ~0.3%; full group total 6.55% → 1.57% (external 89.6% → 25.9%) |
| 10-06 | [#16](16_qe_paw_segni_stale_ux.md) QE PAW `segni` | DFT eigenvalues | m-3m′ breaking 4.83 meV → 46 μeV (PBE); LDA YL2 1.18 μeV |
| 10-10 | window edges in eigenvalue gaps ([log 10-08 §12.8](../log/2026-10-08_mbt_noncollinear_gga_lda.md)) | PBE 210600 → 212616, tol 1e-3 | total 4.04% → 0.67%, external 85.9% → 16.4% (210600's dis_froz_max sat 0.16 meV from a band) |
| 10-11 | [#13](13_ws_tie_tolerance_too_tight.md) tolerance 1e-3 → 5e-3 Å | LDA 212615 (PBE 212616) | external 17.8% → 0.13% (16.4% → 0.22%); total 0.68% → 0.52% (0.67% → 0.68%) |
| 10-12 | [#22](22_hopping_cut_and_hr_rounding.md) H(R) from `.chk`, no hopping cut | LDA 212615 (PBE 212616), tol 5e-3 | internal 0.54% → 0.01% (0.71% → 0.01%); total 0.52% → **0.00%** (0.68% → 0.01%) |

Earlier pipeline fixes (#2 embedding, #3/#10 H and A on one WS rule, #4/#5/#11 `_r.dat`) changed dθ itself, not a
symmetry residual. See their files.

**MnBi₂Te₄**

| date | fix | before → after |
|---|---|---|
| 09-27 | [#12](12_qe_paw_vrad_grad_uninitialized.md) PAW `v_rad`/`g_rad` | uninitialized noncollinear GGA potential → zeroed (upstream 87ac195cf); moved the MBT SCF state |
| 10-04/05 | [#14](14_mbt_stale_binary.md), [#15](15_mbt_u3_base_rounded_cell.md)/[#9](09_mbt_t2_rounded_cell_and_labels.md) stale binary, rounded cell | 4 symmetry operations → 12 (C3 restored by the exact hex-CIF cell) |
| 10-08 | [#18](18_qe_noncollinear_gga_spurious_state.md) noncollinear PBE → LDA+U | \|m\|abs drifting 10.3 → 12.6-13.5 μB, layers asymmetric → base 209653 converged in 12 iterations, NSCF symmetric to 0.043 μeV |
| 10-09 | [#20](20_w90_centre_wrap_long_cell.md) projection sites recentred | interpolated path gap −346 meV (centres 27 Å off) → +104 meV indirect, 185 meV smallest direct (DFT 198 meV on the mesh; the rest is the 4×4×4 off-mesh error) |
| 10-11 | [#19](19_qe_sym_ftrans_trev_energy.md) QE `irreducible_BZ` | T₂⁻ modes: 13 k → 16 k; E 0.2-0.3 mRy below the base → above it (E vs −½F·Δr 1.01-1.24); \|m\|abs 10.27 → 9.99 |
| 10-11/12 | [#21](21_mbt_wf_centres_not_c3_covariant.md) atom-pinned fold (production option 10-12) | interpolated bands break C3 by 62 meV → 0.000 meV |
