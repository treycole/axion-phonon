# axion-phonon — status

**Goal.** dθ/dλ, the linear response of the axion angle θ to a phonon amplitude λ, from a Brillouin-zone integral
of Tr[ΩΩ] on Wannier-interpolated DFT Hamiltonians, for two magnetic insulators:
- **MnBi₂Te₄** (θ = π), zone-boundary T₂⁻ modes;
- **Y₂Ir₂O₇** (θ = 0, AIAO), zone-centre Γ₂⁻ (A₂ᵤ) modes.

The paper draft is [`../paper/draft/main_draft.tex`](../paper/draft/main_draft.tex) (how to compile: [`../paper/README.md`](../paper/README.md)).
**To reproduce any figure, table or number of the paper, start at [`../REPRODUCE.md`](../REPRODUCE.md).** The motivation is
David's 2024 proposal, [`../README.md`](../README.md). The history of everything below, up to 2026-11-02, is in
[`log/2026-11-03_status_through_11-02.md`](log/2026-11-03_status_through_11-02.md).

## Where things stand (2026-11-03)

**Both materials are done at production settings.** The settings, trial sets, the response of every pattern, the
axion charges, the normal modes and every convergence test (with its status) are in **[`../RESULTS.md`](../RESULTS.md)**;
each number there is in the paper with a `%` comment naming its data.

**Validation of the mixed (k, u) curvature:** SrTiO₃ Ti Born charge 7.404 (all-valence) against DFPT 7.289
(`data/SrTiO3/README.md`, `analysis/validation/sto_born_charge.ipynb`).

**Paper figures** (`../analysis/figures/README.md`): Figs. 1, 3, 5, 6, 7 and the c₂ distribution are made from
production data; Figs. 2 (selection rules), 4 (workflow) and 8 (detection) are schematic placeholders. Remaining
`TBD` in the draft: the Zenodo DOI.

**Running (11-03/04):**
- The Y₂Ir₂O₇ U scan (`../data/Y2Ir2O7/studies/u_scan/README.md`): θ′ of pattern 3 against the gap at U = 2.75, 2.50,
  2.25, 2.00 eV. U = 2.75: c₂ (221057-60); U = 2.50: pw2wannier90 and the pattern-3 NSCF.
- Model convergence (decided 11-03, after the MnBi₂Te₄ window spread): MnBi₂Te₄ production window on the 10³ Wannier
  mesh; larger trial sets for both materials (Y₂Ir₂O₇ `Yd`, 216 WFs; MnBi₂Te₄ a probe set first); projectability of
  the production sets. Status rows in `../RESULTS.md`.

**Cluster:** `/mnt/wkdk/bc889/axion/` holds 1.8 TB (target < 2 TB) after deleting the large files of the MnBi₂Te₄
`archive/` on 11-03 (list in `MnBi2Te4/archive/deleted_large_files_2026-11-03.txt`). One SCF/NSCF save per structure
in `dft/`.

## Production settings

- **QE:** pw.x md5 **`db4b99b4`** (`/mnt/wkdk/bc889/axion/bin/qe_db4b99b4/`; branch `production/axion-phonon`, commit
  `fa9b19c8f`): the new_ns_nc, PAW v_rad, segni and irreducible_BZ fixes. **Wannier90:** 4.0.3 (develop
  v4.0.1-160-gc4a924a5, wannier90.x md5 `26c8b506`), projection only, `mpirun -np 1`. Details:
  [`howto/qe_wannier_builds.md`](howto/qe_wannier_builds.md).
- **A(R):** `transl_inv_full_from_chk` with `SHARED_WS = True`. **H(R):** from the `.chk` + `.eig`, no hopping cut
  (bug [#22](bugs/22_hopping_cut_and_hr_rounding.md)). Keep the `.eig` next to every `.chk`.
- **Fixed embedding τ:** shared across base and displaced structure, pinned to the atoms (`wannier.pair_models`).
- **Trial functions of a displaced structure:** the base's Wannier functions through the PAW overlaps of step 4.
- **SCF conv_thr 1e-10** for every structure ([conv_thr for dθ](log/2026-10-17_yio_amplitude_series.md)).
- **Per material:** [`../analysis/mbt/MnBi2Te4.md`](../analysis/mbt/MnBi2Te4.md), [`../analysis/yio/Y2Ir2O7.md`](../analysis/yio/Y2Ir2O7.md) (setup), and the data READMEs [`../data/MnBi2Te4/README.md`](../data/MnBi2Te4/README.md), [`../data/Y2Ir2O7/README.md`](../data/Y2Ir2O7/README.md) (every run).

## Open

1. **Paper:** Figs. 2, 4, 8; the Zenodo deposit (DOI); the U-scan result for the "proximity to band inversion" paragraph.
2. **Upstream:** drafts in [`upstream/`](upstream/README.md) (QE patches for bugs #7, #16, #19); merge requests not opened.

## Where to find things

| folder | what | edit policy |
|---|---|---|
| this file | current status | rewritten at each meeting |
| [`theory/`](theory/) | [`derivation.md`](theory/derivation.md) (**canonical**); [`paw_overlap.md`](theory/paw_overlap.md) (the cross-structure PAW overlap of pipeline step 4); [`decomposition.md`](theory/decomposition.md) (symbol map, the two silent failures, the 13-point validation checklist); [`external_terms.md`](theory/external_terms.md); [`wannier_external_omega_derivation.md`](theory/wannier_external_omega_derivation.md) (checks the code against David & Ivo's notes); [`rigid_shift.md`](theory/rigid_shift.md); [`axion_selection_rules.md`](theory/axion_selection_rules.md) (which phonons and magnons move θ; computed) | living |
| [`code.md`](code.md) | the three code stages, function ↔ equation table | living |
| `../analysis/{mbt,yio}/*.md` | per-material setup: structure, irreps, Wannier set (moved from `notes/materials/`) | living |
| [`howto/`](howto/) | cluster use, QE/Wannier90 builds, rigid-shift commands, magnetic symmetry by hand (Bilbao) | living |
| [`bugs/`](bugs/README.md) | one file per bug (#1-#23), status table in its README | updated when status changes |
| [`upstream/`](upstream/README.md) | drafts of the bug reports and merge requests for QE, Wannier90, PythTB | drafts |
| [`log/`](log/README.md) | dated investigation notes, one-line index in its README | **append-only** |
| [`archive/`](archive/) | superseded derivation drafts, the old bug log, early YIO run notes | frozen |

Outside `notes/`:
- [`meetings/`](meetings/): one file per group meeting.
- [`../tests/README.md`](../tests/README.md): the test suite. The older notebook suite (old code only): [`../archive/validation/`](../archive/validation/README.md).
- [`../archive/calculations_README.md`](../archive/calculations_README.md): how the drivers and analysis notebooks are run.
- `../data/` (not in git): the runs, one README per material, the same tree as the cluster (`dft/`, `windows/`,
  `phonons/`, `studies/`): MnBi2Te4, Y2Ir2O7, and SrTiO3 (validation only).
- [`../REPRODUCE.md`](../REPRODUCE.md): paper item → code → data, and how to rerun everything.

**Conventions.**
- Each fact has one home. Status lives here, a bug's story in `bugs/`, the evidence in `log/`; elsewhere, link to
  them.
- A log note is never rewritten: correct it with a new note, and add a `> Superseded by …` banner to the old one.
