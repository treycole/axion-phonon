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
| 10-11 | [production_binary_ws_cleanup](2026-10-11_production_binary_ws_cleanup.md) | pw.x `db4b99b4` production (k lists unchanged except MBT T₂⁻); WS tolerance 5e-3; MBT per-WF centres break C3 (pinned fixes); YIO symmetry table; cleanup plan | #13, #19, #21 |
| 10-12 | [yio_lda_modes](2026-10-12_yio_lda_modes.md) | YIO LDA modes 1-3 converged (E vs F within 1%, 24 ops); base windows sit in gaps of every mode NSCF; mode Wannier chain | #22 |
| 10-14 | [chern_simons_theta](2026-10-14_chern_simons_theta.md) | MBT θ = π (parities); the old SCDM 3-form was singular at Γ. Exact-connection CS 3-form reproduces YIO dθ/dQ (c₂) to ≤ 0.25% for all three modes at nk 10; frozen gauge ≲ 0.3% | — |
| 10-14 | [mbt_t2m_dynmat](2026-10-14_mbt_t2m_dynmat.md) | MBT T₂⁻ 4×4 dynamical matrix (LDA+U, F − F(0)): 78 / 124 / 140 cm⁻¹; lowest mode = rigid antiphase septuple-layer (interlayer) motion, unresolved (−20 ± 13 cm⁻¹) | — |
| 10-15 | [mbt_rigid_layer_mode](2026-10-15_mbt_rigid_layer_mode.md) | Rigid-layer SCFs 213463/213464 (Q 0.03/0.05 Å): soft T₂⁻ mode stable at **18.0 cm⁻¹**, E ∝ Q² to 0.1%, E vs F 0.5%; Q = 0.01 runs had the wrong sign along c. Basis reruns 213555-213558 (Q 0.03 Å, conv_thr 1e-10): E/F 1.000-1.001 (Mn 1.24 was SCF noise), **17.8 / 82.7 / 128.7 / 143.3 cm⁻¹** from forces alone, rigid check 2% in stiffness | — |
| 10-16 | [mbt_interlayer_stiffness](2026-10-16_mbt_interlayer_stiffness.md) | Corrected two-gap counting: K_gap = k_Q/8, about 2.9–3.0 × 10¹⁹ N/m³; phonon frequencies unchanged. Roughly half the experimental chain-model stiffness, close to earlier DFT quoted by Bartram et al. | — |
| 10-17 | [yio_amplitude_series](2026-10-17_yio_amplitude_series.md) | YIO Γ₂⁻ E(Q), F(Q), θ(Q) series: 3 basis modes × 0.01/0.02/0.03/0.05 Å at conv_thr 1e-10 + base (213998-214010) | — |
| 10-24 | [cluster_cleanup_plan](2026-10-24_cluster_cleanup_plan.md) | Cluster inventory (10.9 TB): what to delete (superseded, ~4.3 TB), what to slim to outputs (~3.1 TB), production decisions (~1.5 TB) | — |

| 11-03 | [status_through_11-02](2026-11-03_status_through_11-02.md) | The notes/README status, production settings and "Next" history as of 2026-11-02, moved verbatim when the README was rewritten for the wrap-up | — |

All dates are 2026.
