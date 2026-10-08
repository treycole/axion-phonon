# axion-phonon

## Goal

Compute the linear response ∂θ/∂λ of the Chern-Simons axion angle θ to a
phonon amplitude λ, as a Brillouin-zone integral of a gauge-invariant
4-curvature `Tr[ΩΩ]`, evaluated on Wannier-interpolated tight-binding
Hamiltonians from DFT. Applied to two magnetic insulators with contrasting
axion topology:

- **MnBi₂Te₄** — θ=π, protected by inversion + a composite antiunitary
  symmetry; axion-active phonons are the zone-boundary `T₂⁻` mode.
- **Y₂Ir₂O₇** — θ=0 (correlation-driven, all-in-all-out order), a
  non-topological counterpoint; axion-active phonon is the zone-center
  `Γ₂⁻` (`A₂ᵤ`) mode, IR- and Raman-silent.

`paper/main.tex` has the precise current statement. `README.md` at repo root
is David Vanderbilt's original 2024 project proposal (motivation), not a
software readme.

## Cluster practices

- When a file is stale in the cluster, and has large size, delete it. If we don't need it, we don't need to keep it. The cluster should be at less the 2 TB of data.
- We need one good DFT run per structure (undistorted and each distortion from the phonon modes for MBT and YIO). In this ouput folder should be a tmp.tar.gz for Quantum ESPRESSO restarts and pw2wannier.x. There should be an untarred tmp/ folder for the Julia code that generates the .amn of the distorted structure. There should be one .mmn per DFT run, since it doesn't depend on trial wave functions. Extra copies of these files should not be kept around.
- Be sure that we are using the correct versions of Quantum ESPRESSO and Wannier90. The cluster has a few different versions of QE and Wannier90, and we need to be sure that we are using the correct ones. The correct versions are:
  - Quantum ESPRESSO 7.6.0 (with the patch for the new_ns_nc bug and with the PAW initialization bug fixed); for Y2Ir2O7 also the PAW segni/`lsign` fix (bug #16): cluster `~/q-e-segnifix`, pw.x md5 `e270b05b`. MnBi2Te4 runs use the same binary.
  - Wannier90 3.1.0 (with the patch for the rdat_ndegen_applied bug)
- **Functionals (since 2026-10-08):**
  - **MnBi2Te4: LDA+U**, U = 4 eV on Mn 3d, `rel-pz` PAW. Base: `MnBi2Te4/base/soc_lda/u_4.0/` (job 209653).
  - **Pseudopotentials:** Mn `spn` 0.3.1 and Bi `dn` 1.0.0 are the official pslibrary files. Te `rel-pz-n` 1.0.0 is
    not on the QE site, so it was generated with `ld1.x` from the official `rel-pbe` recipe; see
    `MnBi2Te4/base/soc_lda/pseudo_gen/`.
  - **No PBE for MBT with SOC** (pitfall 8).
  - **Y2Ir2O7 is still PBE+U** (U = 3). LDA converges its base where PBE stalls
    (`Y2Ir2O7/phonon/base/u_3.0/lda_test_2026-10-08/`); whether to switch is open.
- **Monitor every SCF:** check total *and* absolute magnetization at each check. `md5sum pw.x` before debugging
  anything at source level.


## Where to look

Read in this order before touching code on an unfamiliar part of this repo:

1. **[`notes/INDEX.md`](notes/INDEX.md)** — map of every physics note: which
   derivation is canonical, which are superseded drafts, and the chronological
   investigation log that explains *why* current defaults are what they are.
2. **[`notes/implementation_map.md`](notes/implementation_map.md)** — the code: three
   stages (`wannier/` writes `_hr.dat`/`_r.dat` → `modules/curvature.py` builds
   Omega → `modules/axion.py` does the β finite difference and c2), a
   function ↔ equation table, and what is legacy.
3. **[`validation/README.md`](validation/README.md)** — the check suite.
   Every check is a Jupyter notebook (`.ipynb`), run cell-by-cell with the
   `axion` conda kernel; `RUN_COMPUTATION = False` (default) reloads the last
   saved run from `results/<check>/` instead of recomputing.
4. **[`calculations/README.md`](calculations/README.md)** — how the code in
   `calculations/` is laid out and run. The user's convention: every script and
   notebook starts with a settings block of plain variables (full paths into
   `data/`, n_occ, amplitude, options with their choices listed), and there are no
   command-line flags or automatic run detection. Derive a path only when it is
   always predictable (`<run>/jaemo/`, `<scf>/<prefix>.scf.out`, the output
   `nk_sweep_*` folder). Some run folders carry a plain-text `NOTE`.

`modules/` (production code), `calculations/` (`run_axion.py` for dθ/dQ of a
base/mode pair set at the top of the file, `analysis/` one notebook per task, `setup/`, `diagnostics/`,
`qe_patches/`), `data/` (git-ignored: the QE/Wannier90 runs and everything computed
from them, e.g. `nk_sweep_*` dθ results next to the mode run, `_AA_cache/`),
`validation/` (the check suite), `notes/` (theory, investigation log, and per-material
notes in `notes/MnBi2Te4/`, `notes/Y2Ir2O7/`), `paper/` (the writeup; its numbers are in
flux until the reruns).

## Environment

Use the `axion` conda env (`conda activate axion`) for everything that
touches Wannier data — it has pythtb on the **`external`** branch (needed for
`W90.pos_r`/`_pos_r`; `main` and other topic branches raise `no '_pos_r'`),
plus scipy, matplotlib, nbformat/nbclient/nbconvert. `validation/unit`
(pytest) is the only part of this repo that doesn't need it.

## The eight things that cause the most silent wrong answers

1. **An R-representative choice is a choice of interpolant, not a gauge.** A
   finite k-mesh fixes `H(R)`/`A(R)` only up to the aliasing class
   `[R] = {R + Σ n_i N_i a_i}`. `H(R)` and `A(R)` must use the *same*
   representative-selection rule or the result isn't the curvature of any
   Hamiltonian — and this is invisible at the ab-initio k-points, only
   showing up on the dense interpolation mesh. See
   `notes/progress/2026-08-20_handoff_wannier_conventions.md` §1.
2. **Stock `_r.dat` cannot carry the position matrix, and PR #702 alone
   doesn't fix that** — but Jae-Mo Lihm's fork now does, under specific
   conditions. Three tiers, don't conflate them:
   - **Stock wannier90.x**: unusable. `6F12.6` rounding destroys >50% of
     components, and the per-link phase is baked in before the file is
     written, so it's unrecoverable — not just imprecise. `use_ws_distance`
     is silently ignored either way.
   - **PR #702** (`transl_inv_full=T`): fixes the link-phase problem, matches
     postw90's `_r_full.dat` bit-exact — but `use_ws_distance` *still* never
     reaches the writer, so it's frozen at the wrong WS representative and is
     still **318,000x off** on the symmetry-forced-zero test. Bit-exact at
     the file level is not sufficient evidence it's usable.
   - **Jae-Mo's `write_ndegen_applied` branch**, with all three flags
     together (`use_ws_distance=T`, `transl_inv_full=T`,
     `write_ndegen_applied=T`): genuinely fixes it. Validated symmetry-exact
     on the production SrTiO3 trial (matching `.mmn`'s own ~1e-9 residuals)
     and cross-checked against an independently-built postw90 `rfull` export
     to ~1e-8 — two separate Fortran implementations agreeing. This is
     `source="rdat_ndegen_applied"` in `wannier/position_matrix.py`, and it
     needs no `.mmn` at all.
   - **Caveat, not a bug**: `rdat_ndegen_applied` still differs from `.mmn`
     by ~6-38% pointwise at current mesh density — that's the same known
     `rfull`-family-vs-`mmn`-family interpolant gap as postw90's own
     `write_aa_r` (§1's aliasing-class point), confirmed by the two
     constructions landing on the identical answer. Not evidence either is
     wrong. **Production (2026-10-07) is
     `transl_inv_full_from_chk`**: Jae-Mo's formula rebuilt from the `.chk`
     overlap matrices, with `SHARED_WS = True` (his own files can't be folded on
     the base's WS images). It is the default in `calculations/run_axion.py`
     and the symmetry scripts. `right_centre_from_chk` is the `.mmn` formula
     from the same `.chk` (the default 09-30 to 10-06, by accident); state
     the A(R) source with every number.
   See `notes/progress/2026-08-27_rfull_vs_mmn_position_matrix.md` §0 (why stock fails),
   `notes/progress/2026-09-10_pr702_transl_inv_full_rdat.md` (why PR #702 alone still
   fails), `notes/progress/2026-09-13_jaemo_transl_inv_full_ws_distance.md` (the fix,
   validated on small test systems), and
   `notes/progress/2026-09-16_jaemo_rdat_production_closure.md` (production-scale
   closure).
3. **The embedding τ must be identical across both structures/endpoints of a
   finite difference** (assumption A3 of the derivation) — build both PythTB
   models with `W90.model(orb_vecs=tau)` using that same array
   (`position_terms` reads it from the model; `axion.dtheta` checks it). Otherwise the
   cross/external terms pick up a spurious artifact linear in τ (a past bug:
   `TBModel.lattice` returns a copy, so setting `model.lattice.orb_vecs`
   directly is a silent no-op). See `notes/berry_curvature_decomposition.md`
   §2 and `notes/progress/2026-08-20_handoff_wannier_conventions.md` §3a.
4. **`A_beta = 0` (frozen Wannier gauge) is a valid explicit choice, never a
   silent default to trust unexamined.** Pass `Lambda_R=...` to
   `axion.dtheta` for the full mixed k-β treatment; omitting it means
   every term containing `A_beta` is skipped, not multiplied by zero. Know
   which one a given run is doing. See `notes/berry_curvature_decomposition.md`
   §2.
5. **Occupied-only Wannier sets pass every curvature symmetry test
   structurally** (no conduction block ⇒ internal and cross vanish
   identically) — this cannot fail and is never evidence the pipeline is
   correct. See `notes/progress/2026-08-20_handoff_wannier_conventions.md` §7.
6. **Point-group audits need the pair-dependent lattice shift.** Under a
   space-group operation `g`, a Wannier center at `τ_s` maps to
   `τ_{s'} + L_s`, so the real-space index shifts as `R → gR + L_t − L_s` —
   assuming `R → gR` gives garbage. See `notes/progress/2026-08-20_handoff_wannier_conventions.md`
   §7.

7. **Every Y2Ir2O7 DFT run through 2026-09-28 carries a QE 7.5 bug**, `new_ns_nc` (noncollinear
   DFT+U) symmetrising the Hubbard occupations with the rotation applied in the wrong orientation,
   wrong for every operation of order 3, 4 or 6 that permutes Hubbard atoms (`symmetry_with_labels
   = .true.` with four Ir species is what makes such operations appear). Result: the Ir moments
   differ by ~9% and tilt off their 3-fold axes, a Hubbard spin channel that is too weak, and a
   curvature symmetry floor that was 4-14% and flat across H(R)/A(R). Fix: `calculations/qe_patches/`
   (two lines), or `nosym = .true.`, or drop `symmetry_with_labels`. First real base/mode pair with
   the fix: `notes/progress/2026-09-29_yio_qe76_new_ns_nc_pair.md` — H(R) residuals fell 2-4x and dtheta moved
   ~5x, but the floor didn't close. **Its cause is a second QE bug, #16** (found 2026-10-05): the PAW
   one-centre noncollinear GGA uses a stale `ux` (Ir1's moment) without the `lsign` guard, so the
   NSCF Hamiltonian breaks m-3m′ by ~5 meV (patched `~/q-e-segnifix`, pw.x md5 `e270b05b`: 46 μeV);
   the Wannier/WS/A(R) layers were never the floor. Every YIO run before the segni fix is superseded.
   See `notes/bugs/07_qe_new_ns_nc_symmetrizer.md` and `notes/bugs/16_qe_paw_segni_stale_ux.md`. A separate QE bug,
   uninitialized `v_rad`/`g_rad` in noncollinear-SOC-PAW-GGA runs (found 2026-09-27), predates this
   pair too, but a direct comparison run (192344 vs 193145) showed it does **not** move YIO
   meaningfully (iteration-1 energy to 2e-5 Ry, Ir moment +1.5%, unlike MBT where it was large) — not
   a reason to hold off quoting YIO. See `notes/progress/2026-09-21_qe_new_ns_nc_bug.md` and
   `notes/progress/2026-10-06_segni_fix_validated_and_reruns.md`.
   - **The segni-fixed binary stalls the YIO base SCF at ~3e-5 Ry (bug #17).** The modes converge; the base does
     not. Smoothing the |m| kink in the PAW GGA (runs R1-R3, δ up to 1e-3) did **not** fix it.
   - **Production base meanwhile:** 207546, stalled at 3.1e-5, worst eigenvalue asymmetry 0.10 meV.
   - **Same family as pitfall 8.** In `lda_test_2026-10-08/`, LDA converges the YIO base in 21-22 iterations (Ir
     0.283 μB exactly along ⟨111⟩, 0.68 eV SCF-mesh gap) while PBE wanders, |m|abs creeping from 3.32 to 3.40.
   - **Open:** whether YIO moves to LDA+U.

8. **Noncollinear PBE has a spurious lower-energy MnBi2Te4 state (bug #18); MBT uses LDA+U.** With SOC the run
   must be noncollinear.
   - **The state:** QE's noncollinear GGA admits a layer-asymmetric state 55 mRy *below* the physical one, with
     3.5 μB of transverse interstitial magnetization. It never converges (stuck at 1e-5 to 1e-3 Ry).
   - **Symptoms:** |m|abs climbs from about 10.2 to 13 μB, and the Mn1-layer and Mn2-layer moments stop mirroring.
   - **Why nothing prevents it:** the symmetry that forbids it is the anti-translation {1′|½,½,½} (time reversal ×
     half the doubled cell, which maps the Mn1 layer onto the Mn2 layer). QE can't enforce it: it says "This is a
     supercell" and drops it.
   - **Why a random start finds it:** the default random part of the starting wavefunctions breaks the layer
     equivalence by about 4e-3 in iteration 1.
   - **Mechanism:** the fixed-axis GGA spin density ½(n ± sign(m_z)|m|) jumps wherever m_z changes sign at finite
     |m|. It's a jump, not the |m| kink.
   - **Proofs:** collinear PBE without SOC converges symmetric while noncollinear PBE doesn't; LDA converges from
     any start in 12 iterations.
   - **Rules for every MBT run:**
     - LDA+U (`rel-pz`) and `startingwfc = 'atomic'`. Also `diago_thr_init = 1.0d-5`, which QE 7.2's reference
       152188 got by accident from a leftover restart flag.
     - The exact hex-CIF cell (a = 4.33360, c = 81.85200 Å, atoms at u·c, as in the u_4.0 input). The
       rhombohedral-CIF cell of 151625/152188 is rounded by up to 1e-6 Å.
     - Check layer mirroring and ∫|m⊥| on new runs, with
       `calculations/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py`.
   - **Why it matters for the phonons:** the T2⁻ modes break the anti-translation physically, so in PBE no
     symmetric start can protect them. LDA has no such state at all.
   - **Not the cause:** the QE 7.2 → 8.0dev code audit found no regression. Noncollinear +U now runs as `kind = 0`
     Dudarev rather than Liechtenstein, which is equivalent at J = 0.
   - See `notes/progress/2026-10-08_mbt_noncollinear_gga_lda.md` and
     `notes/bugs/18_qe_noncollinear_gga_spurious_state.md`.

If you hit a curvature discrepancy that isn't one of these, check
`notes/berry_curvature_decomposition.md`'s 13-point validation checklist and
which `validation/` notebook covers each point before assuming new code is
wrong.
