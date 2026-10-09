# Investigation log

Dated notes, one per question settled. Each is correct for its date and is **not edited afterwards**: when a later
note corrects an earlier one, the earlier one gets a one-line `> Superseded by …` banner and nothing else. For what
is true *now*, read [`../README.md`](../README.md); for each bug's final state, [`../bugs/README.md`](../bugs/README.md).

New note: name it `YYYY-MM-DD_<topic>.md`, open with a "Bottom line", and add one row here.

| date | note | outcome | bugs |
|---|---|---|---|
| 08-20 | [handoff_wannier_conventions](2026-08-20_handoff_wannier_conventions.md) | An R-representative is a choice of interpolant (§1); PythTB embedding no-op and H(R)/A(R) WS mismatch found and fixed; traps list (§7) | #2, #3 |
| 08-27 | [rfull_vs_mmn_position_matrix](2026-08-27_rfull_vs_mmn_position_matrix.md) | Stock `_r.dat` cannot carry A(R); `.mmn` can be replaced by an `_AA_cache` | #4 |
| 08-27 | [Ti_Q_2A_full_validation](2026-08-27_Ti_Q_2A_full_validation.md) | SrTiO₃ trial-03/04 curvature difference is finite-mesh interpolation ambiguity, not a bug | |
| 08-28 | [Ti_Q_2A_nk12_A_R_comparison](2026-08-28_Ti_Q_2A_nk12_A_R_comparison.md) | In-house `.mmn` A(R) vs postw90 `transl_inv_full` at 12³: conventions compared | |
| 09-10 | [pr702_transl_inv_full_rdat](2026-09-10_pr702_transl_inv_full_rdat.md) | PR #702 `_r.dat` matches postw90 exactly, but at the wrong WS representative | #5 |
| 09-13 | [jaemo_transl_inv_full_ws_distance](2026-09-13_jaemo_transl_inv_full_ws_distance.md) | Jae-Mo's `write_ndegen_applied` fixes it on small test systems | #4, #5 |
| 09-16 | [jaemo_rdat_production_closure](2026-09-16_jaemo_rdat_production_closure.md) | …and on production SrTiO₃: symmetry-exact, matches postw90 `rfull` to 1e-8. **Closed** | #5 |
| 09-20 | [jaemo_ndegen_yio_mode1](2026-09-20_jaemo_ndegen_yio_mode1.md) | First YIO base/mode pair through the Jae-Mo route; base/mode WS-tie inconsistency found | #6 |
| 09-20 | [legacy_removal](2026-09-20_legacy_removal.md) | Old `berry_curvature`/`axion_angle` pipeline deleted; validation ported; what can't be rerun | |
| 09-21 | [yio_symmetry](2026-09-21_yio_symmetry.md) | YIO curvature obeys the magnetic group only to a 4-14% floor. *Floor explained by 10-05* | |
| 09-21 | [qe_new_ns_nc_bug](2026-09-21_qe_new_ns_nc_bug.md) | QE `new_ns_nc` symmetrizer uses the wrong atom; two-line patch | #7 |
| 09-23 | [ws_ties_in_beta_derivative](2026-09-23_ws_ties_in_beta_derivative.md) | WS image weights must be shared across the β finite difference (A3′); `shared_ws_centres` | #6, #10, #11 |
| 09-29 | [yio_qe76_new_ns_nc_pair](2026-09-29_yio_qe76_new_ns_nc_pair.md) | First pair with the `new_ns_nc` fix: H(R) floor down 2-4×, dθ ≈ −0.0175 rad/Å. *Superseded by 10-05* | #7 |
| 10-05 | [yio_ws_tolerance_full_group](2026-10-05_yio_ws_tolerance_full_group.md) | WS tie tolerance 1e-5 → 1e-3 Å cuts residuals ~4×; what's left isn't tie-breaking | #13 |
| 10-05 | [yio_symmetry_floor_dft_origin](2026-10-05_yio_symmetry_floor_dft_origin.md) | **Root cause of the YIO floor: QE PAW `segni` bug** (stale `ux`, no `lsign` guard) | #16 |
| 10-06 | [segni_fix_validated_and_reruns](2026-10-06_segni_fix_validated_and_reruns.md) | Patch validated (4.83 meV → 46 μeV); how the bug works; production reruns submitted | #16 |
| 10-07 | [validation_methods_and_nk10_mesh](2026-10-07_validation_methods_and_nk10_mesh.md) | How every symmetry check is done; a 10³ mesh makes the floor worse, not better | #16, #17 |
| 10-08 | [mbt_noncollinear_gga_lda](2026-10-08_mbt_noncollinear_gga_lda.md) | Noncollinear PBE has a spurious MBT state; **MBT moves to LDA+U**; YIO PBE base converges from a fresh start | #17, #18 |

All dates are 2026.
