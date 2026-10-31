# axion-phonon — status

**Goal.** dθ/dλ, the linear response of the axion angle θ to a phonon amplitude λ, from a Brillouin-zone integral
of Tr[ΩΩ] on Wannier-interpolated DFT Hamiltonians, for two magnetic insulators:
- **MnBi₂Te₄** (θ = π), zone-boundary T₂⁻ modes;
- **Y₂Ir₂O₇** (θ = 0, AIAO), zone-centre Γ₂⁻ (A₂ᵤ) mode.

The precise statement is in [`../paper/main.tex`](../paper/main.tex). The motivation is David's 2024 proposal,
[`../README.md`](../README.md).

## Where things stand (2026-10-15)

**Restructure (2026-10-29/30, [`restructure/design.md`](restructure/design.md)):** production code is now
`../axion_phonon/` and the cluster steps `../pipeline/` ([`code.md`](code.md)); the old `modules/`, `wannier/`,
`run_axion.py` and the dated diagnostics are in `../archive/` (paths in these notes were rewritten). Steps 1-3 done
(design, package, pipeline dry run), the code side of step 6, and step 4 (MnBi₂Te₄ tree, cluster and `data/`, 10-30:
`MnBi2Te4/README.md`; 1,575 GB of superseded runs and 497 GB of duplicate tarballs in its `archive/` for the user to
delete). Step 5 (Y₂Ir₂O₇) next.

Latest meeting notes (progress 10-21 → 10-27: MBT dθ/dQ for all four patterns, k convergence, O(Q²) estimates, window study, symmetry with pinned positions, c₂ maps per pattern, PBE vs LDA bands): [`meetings/2026-10-27.md`](meetings/2026-10-27.md). Before that (final YIO numbers, θ(Q), MBT phonons and tube, bands, c₂ maps): [`meetings/2026-10-21.md`](meetings/2026-10-21.md). Previous, with the symmetry tests and every bug explained: [`meetings/2026-10-12.md`](meetings/2026-10-12.md).
Methods draft for the paper (not yet in the `.tex`): [`../paper/methods_draft.md`](../paper/methods_draft.md).
Full rewrite of the paper, `main.tex` kept unchanged (2026-10-18): [`../paper/main2.tex`](../paper/main2.tex) with [`../paper/ref2.bib`](../paper/ref2.bib) (every reference checked against Crossref / the publisher; MBT dθ and figures marked TBD).
Paper strategy (MBT lead, claims to verify, venue): [`../paper/pitch.md`](../paper/pitch.md); detection mechanics and YIO estimates: [`../paper/detection_draft.md`](../paper/detection_draft.md).
Beginner-level explainer of the symmetry selection rules (HTML, published at https://claude.ai/artifact/4v4SbB4tx8eCo4GtAKtuaZ): [`../paper/selection_rules_explainer.html`](../paper/selection_rules_explainer.html).

**YIO (LDA) dθ is done and valid (2026-10-13). MBT production (2026-10-30, window f10.962_w12.810_k8, PAW-overlap trial functions, `data/MnBi2Te4/production/`): patterns 1-4 at Q = 0.01 Å −1.95890 / 3.81959 / 6.38695 / −4.54479 rad/Å (base end); pattern 3 Q → 0 (Q² fit through 0.01 and 0.03 Å) 6.390 (base end) / 6.379 (structure end).** The earlier MBT values (window f10.278_w15.169, Julia `.amn`) were −1.985 / 3.833 / 6.427 / −4.543. Every earlier value was computed on DFT
with a since-fixed QE bug. The pipeline itself (Wannier → Ω → dθ) is validated.
- **YIO:** bugs [#7](bugs/07_qe_new_ns_nc_symmetrizer.md) and [#16](bugs/16_qe_paw_segni_stale_ux.md).
- **MBT:** the paper's ~10 rad/Å came from a driver with bug [#2](bugs/02_pythtb_embedding_noop.md) and is invalid.
  MBT is being redone in LDA+U; the T₂⁻ modes also had [#19](bugs/19_qe_sym_ftrans_trev_energy.md).

**Interlayer-stiffness audit (2026-10-16):** the per-gap conversion is k_Q/8, giving about 2.9–3.0 × 10¹⁹ N/m³;
the phonon frequencies below are unchanged. Experimental comparison and its chain-model limitation:
[log 10-16](log/2026-10-16_mbt_interlayer_stiffness.md).

| | MnBi₂Te₄ | Y₂Ir₂O₇ |
|---|---|---|
| **Functional** | LDA+U, U = 4 eV ([#18](bugs/18_qe_noncollinear_gga_spurious_state.md)) | **LDA+U**, U = 3 eV, `rel-pz` PAW (user, 2026-10-12). PBE (YP) superseded |
| **Base DFT** | 209653; symmetric to 0.04 μeV | YL2 210527 (22 iterations) + NSCF 210590 |
| **Wannier** | base recentred (210651), 4×4×4: off-mesh error ~0.1 eV near Γ; **per-WF centres break C3 (interpolated bands 62 meV); atom-pinned positions make them exact** ([#21](bugs/21_mbt_wf_centres_not_c3_covariant.md); `PIN_TO_ATOMS = True` in `run_axion.py` since 10-12). The base's traced curvature is PT-forced zero (1e-6), so it can't test symmetry. Denser mesh: decision pending | base 212615 (window edges in gaps). WS tol 5e-3 (default) + H from `.chk`, no hopping cut (#22, production since 2026-10-12): curvature symmetry **exact** (0.000% total, all 47 operations) |
| **Modes** | **T₂⁻ 1-4 on `db4b99b4` converged** (212708-212711, Q = 0.01 Å, 16 k, \|m\|abs 9.99; these carry dθ). **T₂⁻ frequencies 17.8 / 82.7 / 128.7 / 143.3 cm⁻¹**, all stable, from force-constant SCFs 213555-213558 (Q = 0.03 Å, conv_thr 1e-10): E vs F 1.000-1.001 (the Mn 1.24 at Q = 0.01 was SCF noise), Φ symmetric to 0.07%; 17.8 is the rigid interlayer mode, checked against rigid-layer SCFs 213463/213464 (18.2, 2% in stiffness) ([log 10-15](log/2026-10-15_mbt_rigid_layer_mode.md); notebook `analysis/phonon_frequencies.ipynb`); 8×8×8 NSCFs (nbnd 240) for base + modes: 212997 (base), 212998-213001 (modes 1-4), h_rt 248 h (log 10-12 §10); 10-17 10:25: all at k 61/64 per pool, due ~10-17 16:00-18:00; 1-2 unconverged eigenvalues at k 23 and/or 37 per pool in each (check they lie above the window) | **LDA modes 1-3 converged** (212755-212757, Q = 0.01, `db4b99b4`): E vs −½F·Δr 0.999 / 1.000 / 1.008, 24 ops, \|m\|abs 3.21; NSCF done, base windows in gaps of each. Mode Wannier 212982/212994/212995 (`--no-displacement`): DFT eigenvalues symmetric to 1 μeV, curvature total **0.000%** (max 0.001%) over 24 ops, parts as in the base ([log 10-12 §12](log/2026-10-12_yio_lda_modes.md)). dθ/dQ (nk 12, mean of both estimates; conv_thr 1e-7, frozen top 14.396 eV, atomic base: **superseded 10-19**): mode 1 −0.0180, mode 2 +0.0542, mode 3 +0.0655 rad/Å (§13). **Current (10-19): −0.017145 / 0.054737 / 0.063766 rad/Å** (Q → 0, conv_thr 1e-10, series windows 12.422 / 22.644, self-test base, nk 10; log 10-17). PBE check: mode 1 −0.0176, mode 3 +0.0585 (PBE/LDA 0.97 / 0.89; symmetry total 0.011%, as the PBE base) (§15) |
| **Last dθ** | **production 10-30 (window f10.962_w12.810_k8, Q = 0.01 Å): −1.95890 / 3.81959 / 6.38695 / −4.54479 rad/Å** (base end; `data/MnBi2Te4/README.md`). Before: −1.985 / 3.833 / 6.427 / −4.543 (10-22, window f10.278_w15.169) | LDA modes 1-3 (2026-10-13), above |

Full detail: [`log/2026-10-11_production_binary_ws_cleanup.md`](log/2026-10-11_production_binary_ws_cleanup.md);
earlier [`log/2026-10-08_mbt_noncollinear_gga_lda.md`](log/2026-10-08_mbt_noncollinear_gga_lda.md) (MBT §11, YIO §12).

## Production settings

- **QE:** pw.x md5 **`db4b99b4`** (cluster `~/q-e-kfix`, branch `production/axion-phonon`): new_ns_nc, segni and
  irreducible_BZ fixes. **Wannier90:** 3.1.0 for the `.chk`, projection only. Details and why:
  [`howto/qe_wannier_builds.md`](howto/qe_wannier_builds.md).
- **A(R):** `transl_inv_full_from_chk` with `SHARED_WS = True`. State the A(R) source with every number.
- **H(R):** from the `.chk` + `.eig`, with no PythTB hopping cut (bug [#22](bugs/22_hopping_cut_and_hr_rounding.md), 2026-10-12).
  Keep the `.eig` next to every `.chk`.
- **Fixed embedding τ:** shared across base and mode. **WS tie tolerance 5e-3 Å** (default since 2026-10-11,
  [#13](bugs/13_ws_tie_tolerance_too_tight.md)); both materials use atom-pinned positions (`PIN_TO_ATOMS = True`; MBT since 10-12, [#21](bugs/21_mbt_wf_centres_not_c3_covariant.md); YIO since 10-26, which changes its dθ/dQ by ≤ 2×10⁻⁶ relative).
- **Per material:** [`materials/MnBi2Te4.md`](materials/MnBi2Te4.md), [`materials/Y2Ir2O7.md`](materials/Y2Ir2O7.md).

## Open decisions

1. ~~MBT Wannier mesh~~ decided 2026-10-12: 8×8×8, outer window top 15.047 eV, nbnd 240 (log 10-12 §9-10).
2. **Cleanup:** every cluster run folder is now marked (`axion/INDEX.md`, `_STATUS`; ~5.2 TB SUPERSEDED, log 10-17). `scripts/cluster_cleanup_2026-10-11.py` (keep list updated 10-17; dry run first) frees ~5-6 TB; the automatic run was blocked
   by the permission check, so it waits for the user to run it. Log 10-11 §6.
3. **Upstream:** drafts in [`upstream/`](upstream/README.md); all QE branches are on the GitLab fork (2026-10-11);
   the merge requests are not opened yet.

## Next

1. **YIO:** done (LDA dθ for modes 1-3, PBE comparison for 1 and 3; [log 10-12 §12-15](log/2026-10-12_yio_lda_modes.md)).
   Independently confirmed by the Chern–Simons 3-form difference to ≤ 0.25% at nk 10 ([log 10-14](log/2026-10-14_chern_simons_theta.md)).
   Open: the Q convention; optional interpolated-band symmetry for the modes. **Symmetry table on the 1e-10 runs (10-26,
   user):** self-test base 214249 + ampwin10 modes 214147/214151/214155: total 0.000% (max 0.001%), internal 0.012%,
   cross 0.47-0.52% (1e-7 runs: 0.22%), DFT eigenvalues ≤ 1.35 μeV; paper table updated
   (`../archive/diagnostics/2026-10-26_yio_symmetry_conv1e-10/`). **YIO switched to atom-pinned positions (10-26, user):** images and τ from the 214249 atoms. Symmetry (`symmetry_table_pinned.log`): total 0.000% (max 0.001%), internal 0.000%, cross ≤ 0.018%, external ≤ 0.012% (base centres: cross 0.47-0.53%); paper table, §II reference-position paragraph and settings table updated. dθ/dQ: c₂ series rerun pinned (`../archive/diagnostics/2026-10-26_yio_pinned/run_pinned_series.py`, output `nk_sweep_tif_shared_ws_selftestbase_pinned/`): every end within 2×10⁻⁶ of the base-centre value, Q → 0 limits unchanged to 6 digits (−0.017145 / 0.054737 / 0.063766). 3-form θ(Q) pinned on the cluster (`.../Y2Ir2O7/phonon/cs_theta_series_pinned_2026-10-26/`, jobs 217907-217932; all 24 done; `assemble_pinned.py`: Q → 0 intercepts identical to the base-centre ones in every printed digit (−0.017144 / 0.054737 / 0.063791, nk 12). Fig. 5 now reads the pinned c₂ and θ(Q) (re-rendered 10-26, same numbers). WS ties with pinned positions are exact (≤ 1e-12 Å; next distinct distance 0.39 Å YIO, 0.68 Å MBT; `tie_window_pinned.log`), so the 5e-3 tolerance no longer matters; paper §II updated, and its settings table now lists conv_thr 1e-10 for both materials. Earlier pinned check on the 1e-7 runs: `../archive/diagnostics/2026-10-22_yio_pinned_symmetry/`. **Self-test base (10-19):** base and mode
   ends agree at Q → 0 for all three modes (gap c₀ ≤ 0.08%): −0.017145 / 0.054737 / 0.063766 rad/Å (nk 10, series
   windows; log 10-17 "YIO self-test base: all modes"); these move ~0.5-1% from the atomic-base values.
   **θ(Q) from the Chern–Simons 3-form (10-20, job 214446, nk 10 = nk 12):** Q → 0 intercepts −0.017144 / 0.054737 /
   0.063791 rad/Å, = c2 to ≤ 0.04%; θ linear in Q to 0.6 / 3.1 / 1.6% at 0.05 Å (log 10-17 "θ(Q) from the
   Chern–Simons 3-form"). Paper amplitude figure: θ_CS(Q) with the line Q·dθ/dQ|₀.
   **Amplitude series running** (E(Q), F(Q), θ(Q); 3 modes × 0.01-0.05 Å, conv_thr 1e-10, jobs 213998-214005, 214011-214013, 214015-214016):
   SCF + NSCF done; E(Q), F(Q) harmonic to ≤ 1.5% at 0.05 Å (`analysis/amplitude_series.ipynb`); series windows 12.422 / 22.644 eV (bands cross the production edges at larger Q), base 214064 → amn/w90 214065-214088 → c₂ sweep ([log 10-17](log/2026-10-17_yio_amplitude_series.md)).
2. **MBT:** 8×8×8 NSCFs done (10-17) → chain: base prep 214042 → base Wannier 214043 → mode amn 214056/58/60/62 → mode Wannier 214057/59/61/63; 8×8×8 windows −7 / 10.278 / 15.169 eV (the 4×4×4 edges sit on bands; log 10-17);
   mode symmetry checks (pinned) → dθ with `PIN_TO_ATOMS = True` (first check the longest hopping against `max_distance` 125 Å).
   10-18: base 214043 and mode 3 (214061) done, modes 1/2/4 due 10-19 morning; Q = 0 self-test base 214245. **Mode 3 dθ
   (provisional, mode SCF 1e-7): 6.37 / 6.47 rad/Å (base / mode end), converged to ~0.1%** by in-plane refinement
   around the small-gap k_z line through Γ (uniform meshes fall ~1/nk²: 48 → 2.7 rad/Å at nk 4 → 16; cluster job 214347,
   `../archive/diagnostics/2026-10-19_mbt_refined_bz/`, log 10-17 "MBT c2 integrand: a small-gap tube").
   Final pairs: the same tube scheme (half-width 3, m = 4, 8 k_z planes, 7,928 k per pair) on the cluster.
   **10-19 rebuild at conv_thr 1e-10 (user):** base + modes 1-4 at Q = 0.01 and 0.03 Å, SCF+NSCF on dko64m4 (214290-214298),
   Wannier chain with Q = 0 self-test base (214299-214317); folders `*_conv1e-10_k8_2026-10-19/`,
   script `../archive/setup_templates/MnBi2Te4/modes/setup_mbt_conv10_rebuild.sh`. Base SCF converged (|m|abs 9.99,
   13 k); each 512-k NSCF takes ~29 h on dko64m4 (0.44 h per k per pool; 10-19 18:06: base 37/64, modes 1, 2 34/64; wave 2 = modes 3, 4 and mode 3 Q 0.03 ~10-20 07:00 → ~10-21 12:00). The old 1e-7 k8 chains become
   superseded when it finishes (~1.3 TB to delete, with OK). 10-20: base, modes 1, 2 NSCFs done (|m|abs 9.99, 512 k);
   base Wannier 214300 and self-test 214301 done (Ω_I 325.5 / 344.5 Å²); mode 1, 2 Wannier 214303/214305 in pw2wannier90
   (10-21 00:20, due ~01:30). Wave 2 NSCFs (mode 3, 4 Q 0.01, mode 3 Q 0.03: 214293/214294/214297) at k 33-37 of 64 per
   pool on 10-21 00:15, ~0.48 h per k → done ~10-21 13:00-16:00, then .amn (~2 h) + Wannier (~8 h). **dθ chain scheduled (10-20, user):** per pair stage → c2 array (10 tasks) → assemble,
   held on the self-test 214301 and that mode's w90 job: mode 3 Q 0.01 214625-27, mode 1 214628-30, mode 2 214631-33,
   mode 4 214634-36, mode 3 Q 0.03 214637-39; tube scheme as above; result in
   `MnBi2Te4/dtheta_tube_conv1e-10_k8_2026-10-19/<pair>/dtheta_tube.log`; script
   `../archive/setup_templates/MnBi2Te4/modes/setup_mbt_dtheta_tube.sh`.
   **Mode 1 done (10-21 11:53, 214630): dθ/dQ −1.9881 / −1.9831 rad/Å** (base / mode end; m 4, 8 planes, box
   half-width 3; Q = 0.01 Å, not yet Q → 0). Box half-width 2 → 3: 0.6%; in-plane m 2 → 4: 8% (mode 3's 1e-7 run had
   3.2% there and 0.04% from m 4 → 8; m 8 not run for mode 1). Chunks, klist and log copied to
   `data/MnBi2Te4/archive/julia_amn/p1/dtheta/`.
   **Mode 2 done (10-21 18:06, 214633): dθ/dQ +3.8411 / +3.8929 rad/Å** (same scheme; ends 1.3% apart; min direct gap
   192.3 / 189.7 meV). Copied to `.../mode2/Q_0.01A/lda/dtheta_tube_conv1e-10_214633/`.
   **Convergence not shown for these modes:** at m 4, 4 → 8 planes moves mode 1 by 5.5% and mode 2 by 8.4%;
   m 2 → 4 at 8 planes moves them 8.7% / 5.8%. Mode 3's 1e-7 run had the same-sized steps (4 → 8 planes 5.5%, m 2 → 4
   3.2%) and then converged (8 → 32 planes 0.015%, m 4 → 8 0.04%; log 10-17 "MBT c2 integrand"), so 8 planes and m 4
   are probably enough, but only 16 planes and m 8 for these modes would show it.
   **Mode 3 Q 0.03 stopped (10-21 14:47):** the `.amn` job 214314's window check failed (`WINDOW_FAIL`): dis_win_max
   15.169 eV is 3.1 meV from band 223 at k 293 (same band count as the base at every k); w90 214315 exited, and the
   dθ chain 214637-39 ran anyway (SGE releases holds on exit 100) and wrote an all-NaN `dtheta_tube.log` (deleted).
   **Rerun (10-21 22:01, user OK'd a 3 meV margin):** amn 214976 → w90 214977 → stage 214978 → c2 214979 → assemble
   214980 (stage.sh reads the mode from `214977bc889`; failed logs in `mode3_Q0p03/logs/failed_2026-10-21/`).
   **Convergence run (10-21 22:10, user):** every pair extended to in-plane level m 8 and 16 k_z planes (53,488 k; the
   lists nest, so only 45,560 new k per pair, merged with the production chunks), code
   `/mnt/wkdk/bc889/axion/code/2026-10-21` = 2026-10-19 + pythtb 5776792 (11.5× faster per k), first checked by
   validation 214981 (100 production points of mode 1). c2 / assemble: mode 1 214982/3, mode 2 214984/5, mode 3
   214986/7, mode 4 214988/9, mode 3 Q 0.03 214990/1 (each held on its production assemble job); result
   `<pair>/convergence_M8_N16/dtheta_convergence.log`; scripts `../archive/diagnostics/2026-10-21_mbt_tube_convergence/`.
   **Result, modes 1 and 2 (10-21 22:40; validation: 100 production points to 5.8e-7):** m 8, 16 planes −1.9858 /
   −1.9808 (mode 1), 3.8376 / 3.8893 rad/Å (mode 2); vs production (m 4, 8 planes) −0.12% / −0.09%: 8 → 16 planes
   0.02%, m 4 → 8 0.10% / 0.07%. Box half-width 2 → 3 still 0.6% / 0.8% at m 8 (the cells outside the box are coarse
   points; ring 3 moved 0.6% / 0.9% from m 1 to m 4), so the box was the largest remaining k error: see the box test below. Copied to `data/.../mode{1,2}/Q_0.01A/lda/dtheta_tube_conv1e-10_*/convergence_M8_N16/`. Assemble jobs failed
   on OpenBLAS memory (1 slot, threads unset; fixed in every assemble.sh, rerun by hand).
   **Mode 3 Q 0.01 done (10-21 ~22:30, 214627): 6.4344 / 6.6027 rad/Å** (production scheme; ends 2.6% apart; 1e-7
   provisional was 6.37 / 6.47); m 2 → 4 3.4%, 4 → 8 planes 5.5%; convergence run 214986 running. Copied to
   `data/.../mode3/Q_0.01A/lda/dtheta_tube_conv1e-10_214627/`.
   **Box test (10-21 23:00, user):** box half-width 5 (m 4, 8 planes; 16,568 k, 8,640 new per pair; list
   `box_H5/`): c2 / assemble mode 1 214996/7, mode 2 214998/9, mode 3 215000/1, mode 4 215002/3 (hold 214636), mode 3
   Q 0.03 215004/5 (hold 214980); tables for half-widths 5-2 in `<pair>/box_H5/dtheta_convergence.log`.
   **Result (10-21 23:20):** half-width 3 → 4 → 5 (m 4, 8 planes, base end): mode 1 −1.98812 → −1.98765 → −1.98701
   (0.06%), mode 2 3.84110 → 3.83800 → 3.83671 (0.11%), mode 3 6.43441 → 6.43143 → 6.43119 (0.05%). With m 8 / 16
   planes, the production scheme is within 0.2% of the converged k integral for every mode so far (mode 3's m 8,
   16-plane value 6.43067, −0.06%). Copied to `data/.../mode{1,2,3}/.../box_H5/`.
   **All pairs done (10-22 ~01:00).** The production assemble of mode 3 Q 0.03 (214980, no `OPENBLAS_NUM_THREADS` in its
   assemble.sh) and the m 8 / 16-plane assembles of mode 4 and mode 3 Q 0.03 failed on OpenBLAS memory; all k present,
   rerun on the login node 10-22. **Converged values** = box half-width 5 (m 4, 8 planes) + [(m 8, 16 planes) − (m 4,
   8 planes)] at half-width 3 (the two refinements touch different k, so the corrections add), base / mode end, rad/Å:
   mode 1 −1.9847 / −1.9797, mode 2 3.8332 / 3.8849, mode 3 6.4274 / 6.5956, mode 4 −4.5426 / −4.5306 (Q = 0.01 Å),
   mode 3 Q 0.03 6.4058 / 8.0205. Each is 0.11-0.21% below the half-width 3, m 4, 8-plane value. Mode 4: min direct gap
   192 / 191 meV, ends 0.26% apart. **Mode 3 Q → 0:** the base end hardly moves (6.427 → 6.406), the end difference
   grows 9.6× from Q 0.01 to 0.03 (∝ Q²), so it is not a k or window problem; Q² fits give 6.430 (base end) and
   6.418 (mode end). **Why (10-22):** T₂⁻ is parity-odd, so inversion maps Q → −Q and θ → −θ: θ − π is odd in Q and
   dθ/dQ and every gap are even, changing as Q². Mode 3 shrinks the tube gap most (3.2 / 23.7 meV at Q 0.01 / 0.03 on
   k with base gap < 220 meV; modes 1 / 2 / 4 move the minimum 0.3 / 2.6 / 0.8 meV at 0.01), and 95% of the end difference
   comes from k with base gap < 260 meV (82% of dθ), in every k_z plane. Copied to `data/.../mode4/Q_0.01A/lda/dtheta_tube_conv1e-10_214636/` and
   `data/.../mode3/Q_0.03A/lda/dtheta_tube_conv1e-10_214980/` (each with `convergence_M8_N16/` and `box_H5/`), and
   mode 3 Q 0.01's `convergence_M8_N16/`.
   **Paper and couplings (10-22, user):** main2 MBT TBDs filled except the displaced-structure symmetry table, the MBT
   window/mesh test and the Fig. 5 MBT panels (10-23: symmetry table and window test filled, items 5-6; mesh test and
   Fig. 5 remain). Values = undistorted-structure estimate at Q = 0.01 Å, k-converged
   (pattern 3: within 0.05% of its Q → 0 limit; the mean is 1.3% above). The methods paragraph now says both estimates
   are O(Q²) for an axion-active pattern (symmetry) and that the Q → 0 limit is reported (was: "differ at O(Q); report
   the mean"). Normal modes (`../archive/diagnostics/2026-10-22_mbt_mode_couplings/`): ∂θ/∂q = −0.123 / +0.011 /
   +0.394 / −0.065 rad Å⁻¹ amu⁻½ (17.8 / 82.7 / 128.7 / 143.3 cm⁻¹), per Å −1.45 / +0.12 / +3.66 / −0.61; the 82.7 mode
   is a cancellation (+0.013 with the mean). Axion charges Mn −0.99, Bi +0.96, Te1 +1.61, Te2 −1.14 rad/Å.
   **c₂ figures (10-22):** paper draft `../analysis/figures/c2_distribution.ipynb` now on the 1e-10 mode 3 tube
   (214627), layout unchanged; every symmetry mode: `../analysis/c2_distribution_modes.ipynb` →
   `data/MnBi2Te4/studies/c2_distribution_modes.png`, `data/Y2Ir2O7/phonon/GM2-/c2_distribution_modes.png`.
   **PBE vs LDA bands (10-22, user):** PBE+U DFT bands on the LDA paths: MBT 215098 (from bug #18 test D 209052, the
   symmetric PBE state), YIO 215099 (from YP 210525); `../archive/setup_templates/MnBi2Te4/base/setup_dft_bands_pbe.sh`,
   folders `dft_bands_pbe_2026-10-22/` next to those runs; plot `../analysis/bands_pbe_vs_lda.ipynb`
   → `data/pbe_vs_lda_bands.png`; outputs in `data/<material>/base/pbe/dft_bands_pbe_<job>/`. Both done (pw.x
   db4b99b4; MBT 6.3 h, YIO 23 min). Gaps on the path, LDA+U → PBE+U: YIO indirect 685 → 774 meV, min direct 955 →
   1054 (rigid ~+0.1 eV upshift of the Ir t2g conduction bands, same dispersion); MBT indirect 163 → 165 meV, direct
   at Z 198 → 167 meV (PBE gap nearly direct at Z; LDA's VBM is the off-Γ valence maximum on Γ–L). Valence bands of
   both overlay within ~0.1 eV. MBT direct gap along Γ–Z (the inverted k_z line, where dθ/dQ sits): LDA 218 / 226
   (midway) / 198 meV, PBE 221 / 199 / 167 meV, so ~8% smaller on average and equal at Γ. dθ/dQ per k_z plane (8
   planes of the 1e-10 tubes, undistorted end): 10-15% each, every mode; Z is not special (14% for mode 3).
   10-21 21:30: mode 3 / mode 4 Q 0.01 w90 214307 (since 14:15, n126) / 214309 (since 18:05, n108) running; their dθ
   chains 214625-27 / 214634-36 held on them.
   **Keep the c2 distribution (user, 10-21):** each pair's `chunks/chunk_*.npz.part.npz` (idx into `klist.npz`, c2 per k
   at both ends, direct gap) plus the tube's `klist.npz` are the k-resolved integrand for paper Fig. 7c (the plane
   across the k_z line through Γ). Small; do not delete in cluster cleanups; copy to `data/` when the pairs finish.
3. **Band figures (10-20, user): DFT vs Wannier-interpolated bands of both bases.** DFT bands on each base's `.win`
   kpoint_path from the 1e-10 SCF state (`../archive/setup_templates/MnBi2Te4/base/setup_dft_bands.sh`, both
   materials): MBT 214640 (104 k, from 214290; `.../u_4.0/conv1e-10_k8_2026-10-19/dft_bands_2026-10-20/`), YIO 214641
   (208 k, from 213998; `.../conv1e-10_2026-10-17/dft_bands_2026-10-20/`). Wannier bands: YIO self-test 214249 has
   wannier90's `_band.dat`; MBT's comes with the rebuild (214300 / self-test 214301). Overlay in
   `../analysis/bands.ipynb` (10-20: production `.chk` model by default, `PIN_TO_ATOMS` for MBT; checks the
   cells and k lists agree, x in Å⁻¹ as Wannier90). **YIO done:** DFT 214641 vs self-test 214249 in its frozen
   window: max 2.3 meV, median 0.22 meV over 208 k (old 14.396 window 214125: 235 meV). The 16 Ir e_g-like bands above
   the frozen window are off DFT by +0.47 eV on average (projection only, Y 4d not in the set; log 10-17 "Band figures").
   **MBT done (10-21):** DFT 214640 vs self-test 214301 (pinned, `.chk` H): within E_F ± 1 eV max 13 meV, median ~1 meV;
   gap 191.8 meV vs DFT 194.8 meV (min on the path, next to T); 0.05 meV at the 12 path points on the mesh. Errors
   above E_F + 2.75 eV reach 249 meV (bands within ~0.1 eV under the frozen top 10.278 eV, off-mesh; log 10-17
   "Band figures: MBT"). **Open (10-21):** near Γ off the Γ–Z line (inside the dθ tube) the bands bordering the gap are
   off by up to 10 meV between mesh points and the direct gap by −2%: 8×8×8 interpolation, effect on dθ not measured;
   tests: MBT window test, dθ vs Wannier mesh (log 10-17 "Does the conduction-band misfit matter?"). Window test done
   10-23 (item 5): ≤ 1% at frozen tops ≥ VBM+3 eV, up to 17% below; the mesh test is still open.
   **Paper figures (10-21, user):** `../analysis/figures/fig<n>.ipynb`, one per figure of main2 → `../paper/figures/fig<n>.pdf`; fig3 (bands, copy of `bands.ipynb`) and fig5 (YIO ΔE(Q), 3-form θ(Q) vs the four-curvature slope, intercepts to ≤ 0.04%) done, the rest stubs listing their data. Earlier: last cell of `bands.ipynb` → Fig. 3 +
   Sec. IV E; YIO windows there now the series windows. Per-panel energy ranges down to the lowest occupied Wannier
   bands (YIO O 2s −22.6 eV, MBT Te 5s −13.2 eV rel. VBM).
   **Paper YIO numbers updated (10-21):** tables V (charges −0.0025 / +0.0079 / +0.0053 rad/Å), VI (couplings −0.00204 /
   +0.00118 / −0.00415), appendix D, abstract, summary and the response paragraph now use the production values
   (series windows, 1e-10, self-test base, Q → 0; `../archive/diagnostics/2026-10-21_yio_mode_couplings/`). Band
   physics paragraphs added for both materials (YIO edges: `2026-10-21_yio_band_edges/`; MBT inversion at Γ:
   `2026-10-21_mbt_band_inversion/`) and the MBT integrand tube (mode 3 provisional numbers).
4. **YIO window study (10-19, user):** 5 new (froz, outer) sets + production, Q = 0.01 Å, self-test bases (214318-214343,
   `.../conv1e-10_2026-10-17/nscf_from_213998/w90_trial02_base/window_study_2026-10-19/`); dθ locally by
   `../archive/diagnostics/2026-10-19_yio_window_study/run_window_study.py`. **Done (10-19 13:36):** frozen top
   10.8-13.1 eV and outer top 20.0-24.5 eV agree with the series windows to ≤ 0.75% (all modes); frozen top 14.396
   (production) is the outlier: Ω_I 117 → 159 Å², mode 1 +5.7%, mode 2 −1.1%, mode 3 +0.3%. Use the series windows
   (log 10-17 "YIO window study: result").
5. **MBT window study (10-22, user):** mode 3 Q 0.01, 4 frozen tops at outer 15.169 + 3 outer tops at frozen 10.278,
   each a new base w90 (pw2wannier reused from 214300), self-test base (julia --no-displacement) and mode w90 (from
   214307), then the production tube chain (7928 k). Setup
   `../archive/setup_templates/MnBi2Te4/modes/setup_mbt_window_study.sh`; cluster
   `/mnt/wkdk/bc889/axion/MnBi2Te4/window_study_2026-10-22/<set>/{base,selftest,mode_3,tube}` (jobs 217302-217343;
   window checks pass in base and mode NSCF); logs, `.win` and window checks copied to
   `data/MnBi2Te4/archive/phonon/T2-/mode3/Q_0.01A/lda/window_study_2026-10-22/`. **Done (10-23):** dθ/dQ, m = 4, 8 planes, undistorted / displaced end
   (production f10.278/w15.169 = VBM+3.01/+7.90 eV: 6.43441 / 6.60273 rad/Å, self-test Ω_I 344.5 Å²):

   | frozen / outer top (eV, VBM+) | dθ/dQ (rad/Å) | vs production | Ω_I self-test (Å²) |
   |---|---|---|---|
   | 7.809 / 15.169 (0.54 / 7.90) | 6.25613 / 6.42691 | −2.77 / −2.66% | 255.5 |
   | 8.344 / 15.169 (1.08 / 7.90) | 7.53708 / 7.72007 | +17.14 / +16.92% | 304.1 |
   | 9.564 / 15.169 (2.30 / 7.90) | 5.85244 / 6.01050 | −9.04 / −8.97% | 378.4 |
   | 10.790 / 15.169 (3.52 / 7.90) | 6.48061 / 6.63847 | +0.72 / +0.54% | 319.7 |
   | 10.278 / 12.814 (3.01 / 5.55) | 6.49655 / 6.66733 | +0.97 / +0.98% | 299.6 |
   | 10.278 / 13.696 (3.01 / 6.43) | 6.42837 / 6.59762 | −0.09 / −0.08% | 311.6 |
   | 10.278 / 15.412 (3.01 / 8.14) | 6.41475 / 6.58206 | −0.31 / −0.31% | 349.4 |

   Plateau (≤ 1%) only at frozen tops ≥ VBM+3.0 eV; lower frozen tops move dθ by −9 to +17%, not monotonic, not
   tracking Ω_I. All windows hold every occupied band in the frozen window, so the models agree on the mesh and the
   spread is off-mesh interpolation in the tube (cf. the 10 meV band errors near Γ, item 3). Paper §IV E updated.
   **Mesh test 10³ (user, 10-26) running:** `../archive/setup_templates/MnBi2Te4/modes/setup_mbt_k10_mesh.sh`, folders `<1e-10 run>/../nscf_k10_2026-10-26/{nscf,w90}`. NSCFs 218391 (base, n109) and 218395 (mode 3, n108), 1000 k, ~2.5 days; then base w90 + self-test 218396 and mode 3 w90 218397 (held). The first submission (217903-217906) was lost when n108/n109 went down (~10-26 afternoon, `failed 19`); its input folders are in `nscf_k10_2026-10-26_lost_node_reboot/`. After the outage **n113 had no `/mnt/wkdk`** (218392 failed at launch), so these jobs exclude `!n138&!n113` (n113 mounted again 10-27). Progress 10-28 01:40: k 70/125 per pool at ~1610 s per k, so the NSCFs end ~10-29 02:30 EDT. The lost run left 82 GB scratch folders on both nodes (n109 `/scratch/217903bc889`, n108 `/scratch/217904bc889`), for the user to delete. **Window variants at 10³ (user, 10-28):** `setup_mbt_k10_windows.sh`, `window_study_k10_2026-10-28/`: frozen tops 8.344 (+17% at 8³) and 9.564 (−9%), outer 15.169; per set base + self-test (218671 / 218684, hold 218396), mode 3 (218672 / 218685, hold 218397), tube c2 → sum (218683 / 218696); plus the production-window tube on 218396/218397 (sum 218707). Window edges are checked against the 10³ NSCFs inside the jobs (WINDOW_FAIL → new edges).
   **10-29: NSCFs done (NSCF_OK), but the production windows fail the 10³ window check.** dis_froz_max 10.278 eV is
   0.0 meV from an eigenvalue (k 38, band 163), and dis_win_max 15.169 is 0.5 meV from one (k 18, band 207). WINDOW_FAIL
   from base w90 218396 cascaded: 218397, the 4 variant w90 jobs and 30 c2 tasks failed. The c2 tasks and sums were
   deleted 10-29; the cluster admin flagged the failures.
   - **Edge search (10-29, `archive/setup_templates/MnBi2Te4/modes/window_search.py`, valid on the 8³ and 10³ NSCFs
     of base + mode 3 at once):** no single 10³-only edge works at 8³, and nothing above 13.92 eV works for the outer top.
     | edge | candidates (eV) | notes |
     |---|---|---|
     | frozen top | 10.698 (gap 20.9 meV wide; frozen states 110-116); 10.082 (0.4 meV wide) | 10.698 is VBM + 3.43, in the 8³ plateau |
     | outer top | 13.915; 13.344 | the outer top changed dθ by ≤ 1% at 8³ |
     | variant 8.344 | 8.395 (14 meV wide) | |
     | variant 9.564 | none within 0.5 eV | |
     **Chosen (user, 10-29):** 10.698 / 13.915 (reference) and 8.395 / 13.915, on both meshes.
   - **Mesh test with windows valid on both meshes** (`setup_mbt_windows_k8_k10.sh`, `windows_k8_k10_2026-10-29/`):
     corrected `.amn` throughout, base end = atomic-trial base. Submitted one stage at a time.
     - Stage 1 (10-29 11:45):
       - k8 base w90: 219090 / 219091;
       - k10 pw2wannier90: base 219092, mode 219093;
       - k10 PAW overlaps: 219094.
       - all OK; base Ω_I 279.195 (10.698 / 13.915) and 284.719 Å² (8.395 / 13.915).
     - Stage 2 (10-29 22:30), all OK: k8 mode 3 219252 / 219250 (Ω_I 279.217 / 284.772 Å², +0.008% / +0.019% of
       the base); k10 base 219253 / 219251 (Ω_I 290.958 / 304.494 Å²).
     - Stage 3 (10-29 23:40): k10 mode 3 219266 / 219254; k8 tubes c2 219267-76 + sum 219277 (10.698) and
       219255-64 + sum 219265 (8.395); k8 c₂ maps on k_z = 0, π/c (10.698) 219278-87.
       Results (10-30 01:00), all OK: k10 mode 3 Ω_I 290.986 / 304.540 Å² (+0.010% / +0.015% of the base);
       **8³ tube dθ/dQ (base / mode end, m 4, 8 planes): 10.698 / 13.915 → 6.29562 / 6.45394 rad/Å (−2.1% vs the
       production window's 6.43227); 8.395 / 13.915 → 7.07453 / 7.24244 (+10.0%).**
     - Stage 4 (10-30 01:05): k10 tubes, 10.698: c2 219299-308 + sum 219309; 8.395: c2 219288-97 + sum 219298.
       **All 20 c2 jobs ran out of memory** building the models (the A(R) fold over 27 000 R x 116² at 10³ needs
       more than 8 slots x 3.6 GB; it fit at 8³); the sums found 0 chunks and failed with exit 1. Resubmitted
       10-30 01:35 with 24 slots (86 GB), nothing else changed: 8.395 c2 219310-19 + sum 219320; 10.698 c2
       219321-30 + sum 219331. `pipeline/6_c2.sh` notes SLOTS = 24 for 10³.
     - **Result (10-30 ~02:00): 10³ tube dθ/dQ, 10.698 / 13.915 → 6.23783 / 6.38785 rad/Å (−0.9% / −1.0% vs 8³);
       8.395 / 13.915 → 6.92349 / 7.09004 (−2.1%).** k_z-plane c₂ maps: 10 chunks (219278-87). Everything moved into
       the new tree (`MnBi2Te4/README.md`): `dtheta/T2-_p3_Q0.01_k{8,10}_f{10.698,8.395}_w13.915/`. Decision 6
       (10.698 / 13.915 as production if 8³ and 10³ agree): open, for the user.
     - **Why the window matters (10-30, `../analysis/studies/2026-10-30_window_parts/`, jobs 219357-66; log in
       `data/MnBi2Te4/p3/scf/nscf_k8/studies/window_bands_and_parts/assemble.log`):** in all 10 windows the cross + external
       terms change dθ/dQ by ≤ 0.8%: the window dependence is in the internal (Kubo) curvature. Frozen tops below
       VBM + 3 eV misplace the conduction bands 0-3 eV above the VBM by 120-380 meV (median; up to 575 meV), and those
       are the windows whose dθ/dQ jumps (5.78-7.56). At ≥ VBM + 3.0 eV the bands agree to ~10 meV (median) and dθ/dQ
       is 6.30-6.49 rad/Å; the remaining 3% tracks the outer top (projection-chosen bands > VBM + 3.4 eV).
     - **Frozen-top ceiling:** 116 WFs allow a frozen top of at most 10.962 eV (VBM + 3.69; valid on all six 8³
       structures, `../analysis/studies/2026-10-30_frozen_top/window_search.log`). Valid outer tops on all six: 12.810,
       13.046, 13.697, 14.020, 15.171 (none in 14.02-15.17 or above 15.171; production 15.169 is < 5 meV from a band in
       one of the other patterns). Scan at the ceiling, 10-30: base wannier90 219737-41 (`base/wannier_k8_f10.962_w*`),
       then pattern 3 and window_parts.
   **PAW overlap fix, bug [#23](bugs/23_generate_amn_no_paw_overlap.md) (user, 10-28):** `wannier/paw_overlap.py`
   (S validated against QE to 1e-14). First production pair: MBT pattern 3, Q 0.01, 8³.
   - Folder: `MnBi2Te4/paw_overlap_2026-10-28/mode3_Q0p01_k8/`.
   - Jobs: overlaps 218716-23 → assemble 218724 (Λ of the production pair, corrected `.amn`) → wannier90 mode 218725
     and base check 218726 → tube 218727-36 → sum 218737.
   - **Result:** the corrected pair gives 6.43227 / 6.60007 rad/Å (−0.03%, −0.04%).
     - Λ of the corrected pair: 1e-10 Å⁻¹.
     - Λ of the old pair: 0.03 Å⁻¹.
     - The corrected self-test reproduces the base (Ω_I 325.534103932 vs ...937).
   - **Corrected window study at 8³ (user, 10-28):** `paw_overlap_2026-10-28/window_study_k8/`, assemble 218804 →
     7 sets (wannier90 + tube, sums 218816 ... 218888). The base end of each is that window's atomic-trial base.
   **Open:** dθ vs Wannier mesh (tube dθ on the 10³ models once they finish), and the
   interpolated bands near Γ for each window set (cheap, from the sets' `.chk`). Not started.
6. **MBT displaced-structure symmetry (10-22, user):** `../archive/diagnostics/2026-10-22_mbt_mode_symmetry/`
   (job 217301, 3.4 h; models as the production tube chains, pinned). **Done:** all four modes, −3′m′ (12 ops): DFT
   .eig ≤ 0.06 μeV (162 bands below 10.278; ≤ 21 μeV over all 240), interpolated bands ≤ 0.1 μeV; T_ab = Re Tr[F_a F_b]
   (Tr F is PT-forced zero here too: {−1′|τ½} survives) at 7 random k: internal 0.000%, cross ≤ 0.020%, external
   ≤ 0.14%, total ≤ 0.0065%, all in the C₃ operations, same as the undistorted structure (10-18); at 5 tube k: total
   0.0042%. Primed elements repeat their unprimed partners' test (T is PT-even): 5 independent ops. Controls: unprimed
   inversion (broken at O(Q)) 1.0-1.9% random / 0.5-1.3% tube; C4z 12% / 2%. Paper: Table "symmetry-mbt" and the PΘ
   sentence corrected (it said only the undistorted structure's Tr F vanishes).

---

## Where to find things

| folder | what | edit policy |
|---|---|---|
| this file | current status | rewritten at each meeting |
| [`theory/`](theory/) | [`derivation.md`](theory/derivation.md) (**canonical**); [`decomposition.md`](theory/decomposition.md) (symbol map, the two silent failures, the 13-point validation checklist); [`external_terms.md`](theory/external_terms.md); [`wannier_external_omega_derivation.md`](theory/wannier_external_omega_derivation.md) (checks the code against David & Ivo's notes); [`rigid_shift.md`](theory/rigid_shift.md); [`axion_selection_rules.md`](theory/axion_selection_rules.md) (which phonons and magnons move θ; computed) | living |
| [`code.md`](code.md) | the three code stages, function ↔ equation table | living |
| [`materials/`](materials/) | per-material setup: structure, irreps, Wannier set, current runs | living |
| [`howto/`](howto/) | cluster use, QE/Wannier90 builds, rigid-shift commands, magnetic symmetry by hand (Bilbao) | living |
| [`bugs/`](bugs/README.md) | one file per bug (#1-#21), status table in its README | updated when status changes |
| [`upstream/`](upstream/README.md) | drafts of the bug reports and merge requests for QE, Wannier90, PythTB | drafts |
| [`log/`](log/README.md) | dated investigation notes, one-line index in its README | **append-only** |
| [`archive/`](archive/) | superseded derivation drafts, the old bug log, early YIO run notes | frozen |

Outside `notes/`:
- [`../references/`](../references/README.md): the PDFs (literature and David's handwritten notes).
- [`meetings/`](meetings/): one file per group meeting.
- [`../validation/README.md`](../validation/README.md): the check suite.
- [`../archive/calculations_README.md`](../archive/calculations_README.md): how the drivers and analysis notebooks are run.
- `../data/` (not in git): the runs. Y2Ir2O7 is reorganized (2026-10-13): `base/{lda,pbe}/`,
  `phonon/GM2-/mode*/Q_0.01A/{lda,pbe}/`, superseded runs in `archive/` with a README; the run table is in
  [`materials/Y2Ir2O7.md`](materials/Y2Ir2O7.md). Templates and cluster workflow scripts:
  `../archive/setup_templates/Y2Ir2O7/`. MnBi2Te4 the same way (2026-10-13): `base/lda/`,
  `phonon/T2-/mode*/Q_0.01A/lda/`, `archive/` (every PBE run); run table in
  [`materials/MnBi2Te4.md`](materials/MnBi2Te4.md), templates in `../archive/setup_templates/MnBi2Te4/`.
  SrTiO3 is not reorganized yet.

**Conventions.**
- Each fact has one home. Status lives here, a bug's story in `bugs/`, the evidence in `log/`; elsewhere, link to
  them.
- A log note is never rewritten: correct it with a new note, and add a `> Superseded by …` banner to the old one.
