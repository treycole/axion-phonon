# axion-phonon — status

**Goal.** dθ/dλ, the linear response of the axion angle θ to a phonon amplitude λ, from a Brillouin-zone integral
of Tr[ΩΩ] on Wannier-interpolated DFT Hamiltonians, for two magnetic insulators:
- **MnBi₂Te₄** (θ = π), zone-boundary T₂⁻ modes;
- **Y₂Ir₂O₇** (θ = 0, AIAO), zone-centre Γ₂⁻ (A₂ᵤ) mode.

The precise statement is in [`../paper/main.tex`](../paper/main.tex). The motivation is David's 2024 proposal,
[`../README.md`](../README.md).

## Where things stand (2026-10-08)

**No dθ number is currently valid.** Every earlier value was computed on DFT with a since-fixed QE bug. The
pipeline itself (Wannier → Ω → dθ) is validated; the DFT is being redone.
- **YIO:** bugs [#7](bugs/07_qe_new_ns_nc_symmetrizer.md) and [#16](bugs/16_qe_paw_segni_stale_ux.md).
- **MBT:** the paper's ~10 rad/Å came from a driver with bug [#2](bugs/02_pythtb_embedding_noop.md) and is invalid.
  MBT is being redone in LDA+U.

| | MnBi₂Te₄ | Y₂Ir₂O₇ |
|---|---|---|
| **Functional** | LDA+U, U = 4 eV (switched 10-08, [#18](bugs/18_qe_noncollinear_gga_spurious_state.md)) | PBE+U, U = 3 eV; LDA under test |
| **Base DFT** | 209653 converged in 12 iterations; NSCF eigenvalues symmetric to 0.04 μeV | YP (210525) converged from a fresh start; NSCF 210589 running |
| **Wannier** | base: 92-WF trial submitted (210619) | base + 3 modes: held chain on YP (210599-210606) |
| **Modes** | T₂⁻ on LDA: not started | Γ₂⁻ modes 1-3 rerun from YP's state (210596-210598) |
| **Last dθ** | none valid | pre-fix −0.0175 / +0.045 / +0.169 rad/Å (modes 1-3), **superseded** |

Full detail: YIO in [`log/2026-10-08_mbt_noncollinear_gga_lda.md`](log/2026-10-08_mbt_noncollinear_gga_lda.md) §12;
MBT in §11 of the same note.

## Production settings

- **QE:** segni-fixed pw.x, md5 `e270b05b`. **Wannier90:** 3.1.0 for the `.chk`, projection only. Details and why:
  [`howto/qe_wannier_builds.md`](howto/qe_wannier_builds.md).
- **A(R):** `transl_inv_full_from_chk` with `SHARED_WS = True`. State the A(R) source with every number.
- **Fixed embedding τ:** shared across base and mode. WS tie tolerance: see [#13](bugs/13_ws_tie_tolerance_too_tight.md).
- **Per material:** [`materials/MnBi2Te4.md`](materials/MnBi2Te4.md), [`materials/Y2Ir2O7.md`](materials/Y2Ir2O7.md).

## Open decisions

1. **YIO: PBE (base YP) or LDA+U?** The NSCF gap/symmetry comparison is running (210587 PBE, 210588 LDA).
2. **YIO WS tolerance:** whether to apply 5e-3 Å ([#13](bugs/13_ws_tie_tolerance_too_tight.md)).
3. **Cleanup:**
   - Delete the old buggy YIO runs (~3 TB) once the new ones validate.
   - Delete the superseded MBT PBE runs.
4. **Report upstream?** The `segni` bug (#16) and the noncollinear GGA state (#18) to QE.

## Next

1. **YIO:** check each mode's branch (ΔE, forces, Ir |m|) → Wannier → symmetry checks → dθ.
2. **MBT:** finish the base Wannierization; build the T₂⁻ structures on LDA (exact cell); then the same chain.

---

## Where to find things

| folder | what | edit policy |
|---|---|---|
| this file | current status | rewritten at each meeting |
| [`theory/`](theory/) | [`derivation.md`](theory/derivation.md) (**canonical**); [`decomposition.md`](theory/decomposition.md) (symbol map, the two silent failures, the 13-point validation checklist); [`external_terms.md`](theory/external_terms.md); [`wannier_external_omega_derivation.md`](theory/wannier_external_omega_derivation.md) (checks the code against David & Ivo's notes); [`rigid_shift.md`](theory/rigid_shift.md) | living |
| [`code.md`](code.md) | the three code stages, function ↔ equation table | living |
| [`materials/`](materials/) | per-material setup: structure, irreps, Wannier set, current runs | living |
| [`howto/`](howto/) | cluster use, QE/Wannier90 builds, rigid-shift commands | living |
| [`bugs/`](bugs/README.md) | one file per bug (#1-#18), status table in its README | updated when status changes |
| [`log/`](log/README.md) | dated investigation notes, one-line index in its README | **append-only** |
| [`archive/`](archive/) | superseded derivation drafts, the old bug log, early YIO run notes | frozen |

Outside `notes/`:
- [`../references/`](../references/README.md): the PDFs (literature and David's handwritten notes).
- [`../meetings/`](../meetings/): one file per group meeting.
- [`../validation/README.md`](../validation/README.md): the check suite.
- [`../calculations/README.md`](../calculations/README.md): how the drivers and analysis notebooks are run.

**Conventions.**
- Each fact has one home. Status lives here, a bug's story in `bugs/`, the evidence in `log/`; elsewhere, link to
  them.
- A log note is never rewritten: correct it with a new note, and add a `> Superseded by …` banner to the old one.
