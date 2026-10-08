# Notes index

A map of `notes/`, not a summary of physics. If you're picking this project up
cold, read this file, then `berry_curvature_derivation.md`, then
`implementation_map.md`, then `../validation/README.md` — in that order —
before touching code.

The project goal itself isn't documented here: `../README.md` is David
Vanderbilt's original 2024 project proposal (why this is worth computing), and
`paper/main.tex`'s abstract is the current, precise statement of what's being
computed for which two materials (MnBi₂Te₄, θ=π; Y₂Ir₂O₇, θ=0).

---

**[`progress/`](progress/)** — the dated investigation log and bug record. Every note there is named
`YYYY-MM-DD_<topic>.md`, so the folder lists in chronological order. **[`progress/README.md`](progress/README.md)**
is the numbered bug log (#1–#16, one header per date): what each bug was, why it was a bug, the fix, the status,
and which dated note holds the evidence. Start there for "where does the project stand".

**[`bugs/`](bugs/)** — one file per bug (`NN_<slug>.md`, numbered as in the log), each opening with its status
(FIXED / PARTIALLY FIXED / OPEN), then the symptom, root cause, why it was a bug, the fix, verification and what
remains. [`bugs/README.md`](bugs/README.md) is the status table.

---

## 1. Theory — the derivation, and its two retired drafts

**[`berry_curvature_derivation.md`](berry_curvature_derivation.md)** is
*canonical*. It derives the gauge-covariant non-Abelian Berry curvature of an
occupied subspace from a finite Wannier representation, split into
internal/cross/external pieces, through §8 (the phonon direction β) and §9
(the axion response itself). Every other note and the implementation
(`modules/curvature.py`, `modules/axion.py`) cites this one by equation number.

It went through two earlier full drafts, both now in
[`archive/`](archive/) with a superseded banner — keep them only if you need
to see how the derivation evolved, never cite them:

1. `archive/wannier90_pythtb_nonabelian_curvature.md` (2026-07-29) — first draft.
2. `archive/external_derivation.md` (2026-08-03) — notation cleanup of draft 1.
3. `berry_curvature_derivation.md` (2026-08-04) — **current**, one day later.

Two companions extend the canonical derivation rather than re-deriving it:

- **[`external_terms.md`](external_terms.md)** — the math of just the
  external terms, with an explicit table mapping each derivation object to
  the module that implements it (`wannier/`, `curvature.py`, `axion.py`,
  the Y₂Ir₂O₇ driver `run_axion.py`).
- **[`berry_curvature_decomposition.md`](berry_curvature_decomposition.md)**
  — the full symbol map (derivation ↔ [VS] ↔ code), **the two things that
  fail silently** (Λ optionality, fixed-embedding τ), and the 13-point
  validation checklist with which `validation/` script covers which point.
  **Read this before debugging a curvature discrepancy** — it will usually
  tell you which of the two silent-failure modes you're looking at.

One more compares the code against outside literature rather than deriving
from scratch:

- **[`wannier_external_omega_derivation.md`](wannier_external_omega_derivation.md)**
  — checks `modules/curvature.py` against
  `Wannier_interpolation_of_geometric_quantities-2.pdf` (an external
  reference note, not ours), for the external k-k planes and the mixed
  k-β planes specifically. Its §12 ("Practical options") is where the
  fixed-vs-moving-center choice is actually argued, not just asserted.

---

## 2. Code architecture

**[`implementation_map.md`](implementation_map.md)** — one page: the three
stages (`wannier/` → `modules/curvature.py` → `modules/axion.py`), a
function ↔ equation table, how to use it yourself, and what is legacy. Start
here for "what calls what."

---

## 3. Investigation log (`progress/`) — chronological, not living reference

Each of these settles one specific open question, is dated, and explicitly
opens or closes the question raised by its neighbor. Treat them as a lab
notebook: correct for the date it was written, not necessarily reflecting the
current default (check `../modules/*.py` or `../validation/_paths.py` for
what's actually live today). Read in this order if you're reconstructing how
a current default was decided:

1. **[`progress/2026-08-20_handoff_wannier_conventions.md`](progress/2026-08-20_handoff_wannier_conventions.md)**
   (2026-08-20) — the aliasing-class/interpolant framing (§1, "the one idea
   that explains almost everything"), two real bugs found and fixed (silent
   PythTB embedding no-op; H(R)/A(R) convention mismatch), and a "traps that
   cost time" section (§7) worth rereading whenever something looks like
   Fourier-truncation error.
2. **[`progress/2026-08-27_rfull_vs_mmn_position_matrix.md`](progress/2026-08-27_rfull_vs_mmn_position_matrix.md)**
   (2026-08-27) — closes HANDOFF's open questions. Establishes `_r.dat`
   cannot carry the position matrix at all (§0), and that `.mmn` can be
   dropped from disk in favor of a `_AA_cache/*.npz` (§0b).
3. **[`progress/2026-08-27_Ti_Q_2A_full_validation.md`](progress/2026-08-27_Ti_Q_2A_full_validation.md)**
   — same day, the validation verdict: the trial-03/trial-04 curvature
   difference is a finite-mesh interpolation ambiguity, not a bug.
4. **[`progress/2026-08-28_Ti_Q_2A_nk12_A_R_comparison.md`](progress/2026-08-28_Ti_Q_2A_nk12_A_R_comparison.md)**
   — whether the in-house and postw90 position-matrix constructions use
   compatible conventions, at 12×12×12.
5. **[`progress/2026-09-10_pr702_transl_inv_full_rdat.md`](progress/2026-09-10_pr702_transl_inv_full_rdat.md)**
   — Wannier90 PR #702's `_r.dat` matches postw90's `A(R)` exactly, but at
   the wrong Wigner-Seitz representative for symmetry work.
6. **[`progress/2026-09-13_jaemo_transl_inv_full_ws_distance.md`](progress/2026-09-13_jaemo_transl_inv_full_ws_distance.md)**
   — extends #5; Jae-Mo Lihm's `write_ndegen_applied` branch fixes exactly
   the representative problem PR #702 left open, validated on two small
   Wannier90 test-suite systems. Left the production-scale rerun open.
7. **[`progress/2026-09-16_jaemo_rdat_production_closure.md`](progress/2026-09-16_jaemo_rdat_production_closure.md)**
   — closes #6 on the real production SrTiO3 system: `rdat_ndegen_applied`
   (`wannier/position_matrix.py`) is symmetry-exact and cross-validated
   against an independently built postw90 `rfull` export to ~1e-8. **Closed.**
8. **[`progress/2026-09-20_jaemo_ndegen_yio_mode1.md`](progress/2026-09-20_jaemo_ndegen_yio_mode1.md)**
   — first real base/mode *pair* through `rdat_ndegen_applied` (YIO mode 1),
   from the existing `.chk`s with `restart = plot`, no DFT. The route works and
   `H(k)` matches production, but dθ is ~0.5–0.6× the `.mmn` value, and the new
   `ws_pair_consistency` check finds a 5.5% base/mode Wigner-Seitz tie
   inconsistency in `H(R)` that the `.mmn` route shares. Also records why `A(R)`
   images are not copies under `transl_inv_full` (per-b `exp(-i b.R/2)`). **Open**:
   what the tie inconsistency does to dθ (next steps listed in the note). Also records the regression of the rewritten
   `curvature.py`/`axion.py` against this run.

`validation/README.md` is the living counterpart to this chain — it's what
currently runs and passes, not what was learned along the way.

9. **[`progress/2026-09-20_legacy_removal.md`](progress/2026-09-20_legacy_removal.md)**
   — `berry_curvature.py`, `axion_angle.py` and the old YIO driver deleted; the
   validation libs, notebooks and unit tests ported to `curvature.py`/`axion.py`
   and rerun where the data still exists. Old-versus-new agreement is 1e-9 on plane
   eigenvalues of real models and 1e-16 on the mixed-plane Born charge; the primary
   PythTB gate on YIO passes at 6e-17. Records which checks cannot be rerun (inputs
   gone), the frozen MnBi2Te4 driver behind the paper's number (deleted 2026-09-24), and a latent
   `embedding=` trap for tau = 0 models that the port removes. **Nothing removed
   was in git**; the archive is `tmp/legacy_archive_2026-09-20/`.

10. **[`progress/2026-09-21_yio_symmetry.md`](progress/2026-09-21_yio_symmetry.md)**
    — does the YIO curvature obey the magnetic point group at base and displaced? To a floor of
    4 to 14% in both structures and both position-matrix routes (`.mmn` and Jae-Mo's), against
    about 190% for a non-symmetry. The floor is common to the two routes (shared H(R)) but **not**
    common to the two structures, and the 0.01 A inversion breaking sits at or below it. Corrects
    the notebook's mode-1 z-axis expectation (all components vanish).

11. **[`progress/2026-09-21_qe_new_ns_nc_bug.md`](progress/2026-09-21_qe_new_ns_nc_bug.md)**
    — root cause of the YIO symmetry floor: a bug in QE 7.5's noncollinear DFT+U symmetriser (`new_ns_nc`
    pairs the rotation with the wrong atom for operations of order >= 3). Evidence ported from QE's own
    routines, a two-line patch, and a real `pw.x` run (identical raw data, patched vs unpatched) in which
    the four Ir moments become identical. **Every YIO DFT run so far is affected** and should be rerun
    with the patch, `nosym`, or without `symmetry_with_labels`.

12. **[`progress/2026-09-23_ws_ties_in_beta_derivative.md`](progress/2026-09-23_ws_ties_in_beta_derivative.md)**
    — the derivation's β derivative needs the Wigner-Seitz image weights fixed as well as τ (A3′). The
    in-house fold used each run's own centres, which gives an `O(1/Δβ)` term off the ab-initio mesh
    (zero at nk = 4, 8 for an 8³ grid). Adds `load_wannier_pair(shared_ws_centres=True)` (`.mmn` only)
    and says when it is valid.

13. **[`progress/2026-09-29_yio_qe76_new_ns_nc_pair.md`](progress/2026-09-29_yio_qe76_new_ns_nc_pair.md)**
    — first real base/mode pair on QE 7.6 + the `new_ns_nc` fix (the later PAW `v_rad`/`g_rad`
    fix is confirmed **not** to matter for YIO, unlike MBT). H(R)-only symmetry residuals fall
    2-4x, but the total gets *worse* (dominated by A(R)); **corrects #10's "4-14% floor, common
    to both routes"** — the floor moved into A(R) and is no longer flat. `dtheta` (mode 1,
    `transl_inv_full_from_chk` + shared centres) converges to ≈ -0.0175 rad/Å, about 5x the unpatched
    pair's value. Two negative diagnostics on the A(R) floor's origin: a centre-symmetrization
    test (contaminated by same-site orbital-degeneracy mismatching for 15% of WFs) and a
    symmetrized-BZ-integral test (mathematically guaranteed to be a null result, not
    informative). `num_iter > 0` is **ruled out** as a next test — incompatible with the
    external-curvature formula and not even reliably symmetry-preserving. **Cause still open.**

14. **[`progress/2026-10-05_yio_ws_tolerance_full_group.md`](progress/2026-10-05_yio_ws_tolerance_full_group.md)**
    — the in-house WS tie tolerance (`_DEFAULT_WS_DISTANCE_TOLERANCE`, both H(R) and A(R) folds).
    1e-5 → 1e-3 Å cuts base's residual ~4x for **every** m-3m' operation (median total 6.6% → 1.6%,
    external 90% → 26%). The physical range is [2.05e-3, 2.12e-2) Å (centres symmetric to 5.7e-4 Å;
    empty band in the image-distance gaps), and any value in it gives the same fold. At 5e-3 inversion
    is clean (0.18% total) but the rotations stay at ~1.5% total / 26% external, so **what's left of
    the floor isn't tie-breaking**. Corrects the reading "the tolerance fixed the symmetry floor". dθ
    at the new tolerance not cleanly measured.

15. **[`progress/2026-10-05_yio_symmetry_floor_dft_origin.md`](progress/2026-10-05_yio_symmetry_floor_dft_origin.md)**
    — **root cause of the YIO symmetry floor: a QE bug.** `paw_onecenter.f90::compute_rho_spin_lm` uses
    `segni_rad = SIGN(1, m·ux)` without the `lsign` guard that `compute_rho.f90` has, and `compute_ux`
    leaves `ux` at the first magnetic atom's (Ir1's) moment. So the PAW GGA swaps up/down by the sign of
    m·(Ir1 axis), and only −3m′ about [111] survives. NSCF tests on the cluster: the PAW `v_rad` fix,
    Hubbard and ortho-atomic are not the cause; LDA → 0.0 meV; moving `ux` moves the breaking. Stored ρ,
    m, ns and becsum are exactly symmetric. The curvature residual tracks the DFT asymmetry (corr ≈ 0.9),
    so A(R) isn't a separate source (corrects #13). **Patch validated 2026-10-06**: `~/q-e-segnifix` pw.x
    (md5 e270b05b) brings the DFT asymmetry from 4.83 meV to 46 μeV over all 47 ops. All YIO DFT must be redone with it.

16. **[`progress/2026-10-06_segni_fix_validated_and_reruns.md`](progress/2026-10-06_segni_fix_validated_and_reruns.md)**
    — the segni fix validated (C2: 4.83 meV → 46 μeV), a plain-language explanation of the bug with worked GGA
    examples (local frame vs fixed axis, the 1D antiferromagnet and 2D rotating-moment cases), why MnBi2Te4 is
    unaffected, and the production reruns of base and modes 1-3 (jobs 207546-207549).

17. **[`progress/2026-10-07_validation_methods_and_nk10_mesh.md`](progress/2026-10-07_validation_methods_and_nk10_mesh.md)**
    — a 10×10×10 ab-initio mesh makes the (pre-segni) base symmetry **worse**, not better (median total
    1.57% → 3.60%, `.chk` route); the ~5 meV band residual is mesh-independent (#16). The `.chk` route is 2-5× more
    symmetric than the Jae-Mo files. Also: a Fortran-subrecord fix to the gzipped `.chk` reader, why the SCF
    density is symmetric while the NSCF eigenvalues aren't, and how every symmetry check is done (gk on a
    Γ-centred grid, on-site levels, ns/ρ covariance, run C2, the curvature metric).

18. **[`progress/2026-10-08_mbt_noncollinear_gga_lda.md`](progress/2026-10-08_mbt_noncollinear_gga_lda.md)**
    — **MnBi2Te4 moves to LDA+U** (bug #18). Noncollinear PBE (forced by SOC) has a spurious layer-asymmetric MBT
    state, 55 mRy *below* the physical one, with 3.5 μB of transverse interstitial m. QE drops the anti-translation
    that forbids it, and random starts fall in and stall. The fixed-axis GGA spin density ½(n ± sign(m_z)|m|) jumps
    where m_z crosses zero at finite |m|, and the XC gain (−285 mRy) is all grid XC. Collinear PBE (no SOC) and LDA
    converge symmetric in 12-13 iterations. Also: the settings audit and why QE 7.2's 152188 converged; the QE
    7.2 → 8.0dev code audit (no regression); the exact CIF cell; the full run table; the LDA trade-offs and the DMC
    U shift (LDA+U ≈ 4.4 vs PBE+U ≈ 3.5-4 eV); the generated Te `rel-pz` PAW; the YIO connection (the #17
    |m|-kink regularization did not fix the stall; YIO LDA tests running).

---

## 4. Workflow notes

- **[`rigid_shift.md`](rigid_shift.md)** — Jae-Mo's Julia tool for rigidly
  shifting MLWFs from the undistorted structure into trial orbitals for a
  distorted one, so the two Wannierizations share a gauge and the finite
  difference across them is physical rather than noise. Cluster commands
  included.
- **[`cluster_tips.md`](cluster_tips.md)** — QE restart mode (tar the `tmp/`
  save directory, the cluster doesn't copy it), which nodes have 48 cores,
  the `rsync` pattern for pulling results back.
- **[`MnBi2Te4/`](MnBi2Te4/)** and **[`Y2Ir2O7/`](Y2Ir2O7/)** — per-material working notes (structure,
  SMODES, k-paths, Wannierization trials, progress log), moved here from `calculations/` on 2026-09-24.
  Older than the investigation log above; file names they mention (e.g. `SMODES.out.txt`) may have moved
  into `data/<material>/`.

---

## 5. Reference PDFs (external literature, uncited elsewhere in `notes/`)

`2312.10769v1.pdf`, `Beta.pdf`, `beta-derivatives.pdf`, `jh1y-x9nz.pdf` sit in
`notes/` alongside `Wannier_interpolation_of_geometric_quantities-2.pdf`
(the one actually cited, by §1). Not indexed further here — open them
directly if you need them.

---

## 6. Related, outside `notes/`

| Doc | What it holds |
| --- | --- |
| [`implementation_map.md`](implementation_map.md) | Code architecture (see §2 above) |
| [`../calculations/README.md`](../calculations/README.md) | Layout of `calculations/` and `data/`; the settings-block convention; running `run_axion.py` |
| [`../validation/README.md`](../validation/README.md) | The check suite: what each notebook proves, how to run one, `RUN_COMPUTATION` toggle |
| [`../CLAUDE.md`](../CLAUDE.md) | Top-level orientation for an agent picking this up cold |
| [`../paper/main.tex`](../paper/main.tex) | The current precise statement of the physics goal and results |
