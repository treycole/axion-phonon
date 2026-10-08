# Bugs

One file per bug found in this project, numbered as in the dated log
[`../progress/README.md`](../progress/README.md). Each file opens with its **status**, then describes the symptom,
the root cause, why it was a bug, the fix, how the fix was verified, and what remains.

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
| [9](09_mbt_t2_rounded_cell_and_labels.md) | **OPEN** (fix known) | 2026-09-22 | MBT T2⁻ inputs: rounded bohr cell loses C3; labels trigger #7 |
| [10](10_hr_ar_folded_on_different_centres.md) | **FIXED** | 2026-09-23 | H(R) and A(R) folded on slightly different centre arrays |
| [11](11_rdat_print_precision.md) | **FIXED** (mitigated) | 2026-09-23 | `_r.dat` 6-decimal print costs ~3% of dθ |
| [12](12_qe_paw_vrad_grad_uninitialized.md) | **FIXED** | 2026-09-27 | QE PAW `v_rad`/`g_rad` uninitialized (noncollinear GGA) |
| [13](13_ws_tie_tolerance_too_tight.md) | **PARTIALLY FIXED** | 2026-09-30 | WS tie tolerance too tight; 5e-3 recommended, not applied |
| [14](14_mbt_stale_binary.md) | **FIXED** | 2026-10-04 | MBT u_3.0 ran a stale, buggy pw.x |
| [15](15_mbt_u3_base_rounded_cell.md) | **FIXED** | 2026-10-05 | MBT base u_3.0 rounded bohr cell (same defect as #9) |
| [16](16_qe_paw_segni_stale_ux.md) | **FIXED** (production reruns in progress) | 2026-10-05 | QE PAW GGA uses stale `ux` without `lsign` guard: root cause of the YIO symmetry floor |
| [17](17_segni_fix_base_scf_stall.md) | **OPEN** (candidate fix in testing) | 2026-10-07 | The segni fix (#16) stalls the YIO base SCF at ~3e-5 Ry |
| [18](18_qe_noncollinear_gga_spurious_state.md) | **WORKED AROUND** (MBT → LDA); open for YIO | 2026-10-08 | Noncollinear PBE has a spurious lower-energy, layer-asymmetric MBT state; random starts fall into it and stall |

**Status meanings:**
- **FIXED:** the fix is in the code, builds and inputs used now.
- **(worked around) / (mitigated):** the defect still exists in the tool, but production doesn't depend on it.
- **PARTIALLY FIXED / OPEN:** see the file's "Remaining" section.
