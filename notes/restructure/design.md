# Restructuring design (agreed 2026-10-29)

Agreed by the user on 2026-10-29, with one change: file readers live in `axion_phonon/io/`, apart from the physics
(section 3). Section 7 records the decisions; section 8 is the order of work. The plan:
- **Code:** one small package for the production calculation, in the order the calculation runs.
- **Runs:** one folder tree per material, identical layout for every structure. A displaced structure differs from
  the base only in its atomic positions.
- **Old work:** archived, never deleted.

## 1. Rules

1. **One production run per structure.** Other attempts go to `archive/`; you delete them by hand.
2. **Every structure folder has the same layout** (base and every displaced pattern, both materials). A script that
   works for one structure works for all.
3. **Every result lives in exactly one place** on the cluster, mirrored at the same relative path in local `data/`.
4. **The production code is only the production path.** Alternatives that were tried and rejected (other
   position-matrix sources, the Julia `.amn`, the self-test pairing) move to `archive/` with a note saying why.
5. **Settings stay as plain variables at the top of each script** (your convention): full paths, no CLI flags, no
   registry, no auto-detection. Physical parameters live in the run's own input files (`.scf.in`, `.win`), the files
   QE and Wannier90 read.

## 2. The calculation, in order

This is the dependency chain. Every number in the paper comes out of the last box.

```text
 per structure (base, and each pattern at amplitude Q)
 ─────────────────────────────────────────────────────
 1 scf          pw.x                       → charge density, Hubbard ns
 2 nscf_k<N>    pw.x on the Wannier mesh   → Bloch states (tmp/ save), eigenvalues
 3 pw2wannier   pw2wannier90.x             → .mmn, .eig (base also: .amn onto atomic orbitals)
 4 mode_amn     axion_phonon.paw           → .amn onto the base Wannier functions  (displaced structures only;
                                             reads both saves + the base .chk; replaces generate_amn.jl)
 5 wannier      wannier90.x (projection    → .chk  (V = U_dis U: the Wannier functions)
                only)

 per pair (base + one displaced structure)
 ─────────────────────────────────────────
 6 dtheta       axion_phonon.wannier       → H(R), A(R) of both models (shared fold, atom-pinned τ)
                axion_phonon.curvature     → Ω
                axion_phonon.axion         → c2 on the tube k-list → ∂θ/∂Q
```

Step 4 is the only one that exists for displaced structures alone. The base's step 5 must finish before any
pattern's step 4 starts.

## 3. Repository layout

```text
axion-phonon/
  axion_phonon/              the production code (target ~1,200 lines in all)
    io/                      file reading and writing only; no physics
      qe.py                  QE: save XML (cell, atoms), wfcN.dat, UPF, pw.x outputs (eigenvalues, k list)
      wannier90.py           Wannier90: .chk, .eig, .win, .amn (read and write)
    paw.py                   PAW metric S, overlaps between two structures, the mode .amn          (step 4)
    wannier.py               H(R), A(R) (transl_inv_full), WS fold, pinning → the two PythTB models (step 6)
    curvature.py             Ω and its parts                                                       (step 6)
    axion.py                 β terms, c2, BZ integral                                              (step 6)
    tube.py                  tube k-list and its quadrature                                        (step 6)
  pipeline/                  one script per step, the same for every structure (section 5)
  analysis/                  notebooks: figures and paper numbers (from calculations/analysis, calculations/figures)
  validation/                unit tests (no data) and checks (notebooks)
  notes/  paper/  references/            unchanged
  archive/                   everything else, with a README saying what each piece was and why it was retired
  data/                      git-ignored; mirrors the cluster tree (section 4)
```

**The `io/` rule.** The physics files (`paw.py`, `wannier.py`, `curvature.py`, `axion.py`, `tube.py`) never open a
file. Their functions take arrays (the overlap matrices M_mn(k,b), V(k), eigenvalues, β projectors, wavefunction
coefficients) and return arrays, so each one reads as the equations it implements. The pipeline scripts and notebooks
call `io/` to get the arrays, then call the physics. Nobody needs to read `io/` to follow the calculation; it is
checked once, by the unit tests (round trips, and agreement with what QE and Wannier90 print).

Where today's code goes:

| today | lines | goes to |
|---|---|---|
| `modules/curvature.py`, `modules/axion.py` | 849 | `axion_phonon/` unchanged |
| `wannier/paw_overlap.py` | 487 | readers (`read_wfc`, `read_structure`, `read_upf`) → `io/qe.py`; the production `exact` overlap and the mode `.amn` product → `paw.py` |
| `wannier/wannier_io.py` | 1,058 | the `.chk`/`.eig`/`.win` readers → `io/wannier90.py`; `hamiltonian_classes_from_checkpoint`, `pin_centres_to_atoms`, and the `transl_inv_full_from_chk` branch of `load_wannier_pair` → `wannier.py`; the `.mmn`/`_r.dat`/postw90 branches → `archive/` |
| `wannier/position_matrix.py` | 1,945 | `position_matrix_transl_inv_full`, the orbital-dependent R fold and the A(R) helpers it uses (~300 lines) → `wannier.py`; the five unused sources → `archive/` |
| `modules/qe_read_file.py`, `modules/utils.py` | 560 | readers that production or analysis uses → `io/qe.py`; the rest → `archive/` |
| `calculations/diagnostics/2026-10-19_mbt_refined_bz/` (`make_klist_tube.py`, `c2_klist.py`, `assemble_tube.py`) | | `tube.py` + `pipeline/` |
| `calculations/run_axion.py` | 130 | `pipeline/` (the dθ step) |
| `calculations/setup/*/template/*.sh` (~15 one-off setup scripts) | | replaced by `pipeline/`; originals → `archive/` |
| `calculations/diagnostics/` (31 dated folders) | | `archive/diagnostics/` as is (each is evidence for a log entry) |
| `calculations/analysis/`, `calculations/figures/` | | `analysis/` |
| `wannier/check_outputs.py`, `run_plot_stage.sh` | | `archive/` unless still used |

The split between `wannier/` and `modules/` goes away. It was meant as "read Wannier90's files" vs "physics", but
`wannier/` grew to build the models too, so the boundary no longer means anything.

**Acceptance test for the code step:** recompute c2 at about 100 tube k-points of the corrected pattern-3 pair with
the new package. The values must equal the cluster chunks (the ∂θ/∂Q = 6.43227 / 6.60007 rad/Å run) to ~1e-12, and
the unit tests must pass. Nothing is archived until this passes.

## 4. Run layout (cluster, mirrored in `data/`)

```text
/mnt/wkdk/bc889/axion/<material>/
  README.md                     what each structure and pair is; the paper numbers and the folder each comes from
  base/  T2-_p1_Q0.01/  ...     structure folders, all with the same inside, nested by dependency (user, 10-30):
    scf/                          SCF: inputs, outputs, OK
      bands/                        pw.x bands from this density
      nscf_k8/                      NSCF on 8^3, from this SCF; + tmp/ (untarred save)
        pw2wannier/                   .mmn, .eig (+ the base's atomic .amn)
        overlap/                      PAW overlaps with the base NSCF (displaced structures)
        wannier_<window>/             one wannier90 run per window (no k tag: the parent NSCF has it)
          dtheta/                       dtheta/dQ with the base's wannier_<window> (c2 chunks, dtheta.log)
        studies/                      studies over several Wannier runs of this NSCF
      nscf_k10/ ...
  phonons/                      force-constant SCFs (Q = 0.03 Å), rigid layer, amplitude series
  studies/                      studies across structures
  archive/                      everything superseded, with a list; you delete by hand
```

- **Structure name:** `<irrep>_p<pattern>_Q<amplitude in Å>`. "Pattern" means the symmetry-adapted displacement we
  apply (SMODES basis, today's `mode_3`). It is not a phonon eigenvector, so the name stops clashing with the normal
  modes (17.8, 82.7, ... cm⁻¹).
- **Window label:** `f<frozen top>_w<outer top>`, e.g. `f10.698_w13.915`.
- **Markers and dates:** each step folder holds its job folder(s) `<job>bc889/`, a `NOTE` (what, when, inputs) and
  `OK` or `FAIL`. Dates go in `NOTE`, not in folder names.
- **Y₂Ir₂O₇:** the same, with `GM2-_p1_Q0.01` ... and its amplitude series (`Q0.02`, `Q0.03`, `Q0.05`) as further
  structure folders.

Local `data/<material>/` has the same tree. It holds only what was copied back (`.chk`, `_hr.dat`, `.win`, c2 chunks,
logs), never the saves.

## 5. Pipeline scripts

One script per step in `pipeline/`. Each has a settings block whose main variable is the structure folder (or the
pair, for steps 4 and 6). The script reads the input files there, submits one job, and returns. The rules from the
10-29 cluster failures are built in:

- **Inputs:** the job copies its own inputs into node scratch and checks them; it does not rely on the prolog.
- **Markers:** it writes `OK` or `FAIL` and exits 1 on failure (never 100, so no error state).
- **Order:** a step refuses to submit unless its inputs' `OK` markers exist. Nothing is queued behind a step that
  has not succeeded. Running the chain means running step 1, then step 2 when step 1 is `OK`, and so on.
- **Binaries:** the md5s of the binaries and pseudopotentials are checked, as now.

You can run everything yourself with: set the structure folder, run the step, wait for `OK`.

## 6. What maps where (current production → new tree)

Proposed. The final list is made by an inventory script and checked by you before anything moves.

**MnBi₂Te₄** (all LDA+U, U = 4, conv_thr 1e-10; today under `MnBi2Te4/base/soc_lda/u_4.0/conv1e-10_k8_2026-10-19/`
and `MnBi2Te4/phonon_lda/T2-/mode_<p>/Q_<Q>_conv1e-10_k8_2026-10-19/`):

| new | scf + nscf_k8 | wannier_k8 (production windows 10.278 / 15.169) | other |
|---|---|---|---|
| `base/` | 214290 | 214300 (atomic trials); self-test 214301 → archive | nscf_k10: 218391; dft_bands_2026-10-20 |
| `T2-_p1_Q0.01/` | 214291 | 214303 (Julia `.amn`: redo with step 4) | |
| `T2-_p2_Q0.01/` | 214292 | 214305 (redo) | |
| `T2-_p3_Q0.01/` | 214293 | corrected 218725 (`paw_overlap_2026-10-28/`); 214307 → archive | nscf_k10: 218395 |
| `T2-_p3_Q0.03/` | 214297 | 214315 (redo) | |
| `T2-_p4_Q0.01/` | 214294 | 214309 (redo) | |
| `dtheta/` | | tubes in `dtheta_tube_conv1e-10_k8_2026-10-19/` (+ `box_H5`, `convergence_M8_N16`), corrected pair `paw_overlap_2026-10-28/mode3_Q0p01_k8/tube` | |
| `phonons/` | force-constant SCFs 213555-213558, rigid-layer 213463/213464 | | |

**Y₂Ir₂O₇** (LDA+U, U = 3; today under `Y2Ir2O7/phonon/base/u_3.0/lda_test_2026-10-08/YL2_lda_relpz/` and
`Y2Ir2O7/phonon/mode_<p>/`):

| new | scf + nscf | wannier (series windows 12.422 / 22.644) |
|---|---|---|
| `base/` | SCF 213998, NSCF 214124 | 214162 (atomic trials); self-test 214249 → archive after the PAW fix |
| `GM2-_p<1,2,3>_Q<0.01..0.05>/` | amplitude series 213999-214016 | 214147-214158 (Julia `.amn`: redo with step 4) |

**Studies** (window studies, the 10³ mesh test, PAW checks) become clearly labelled subfolders of `dtheta/` or
`archive/`, never new top-level folders. Today's top-level study folders (`window_study_*`, `paw_overlap_*`,
`windows_k8_k10_*`, `w90_k*_pipeline`, `dtheta_refined_*`) are sorted into those two places.

**Size:** the cluster holds 5.7 TB against the 2 TB target. Production keeps one untarred save per structure and
mesh. The bulk of the rest (old attempts, PBE comparisons, duplicate `tmp.tar.gz` next to an untarred `tmp/`) goes to
`archive/` for you to delete.

## 7. Decisions (agreed 2026-10-29)

1. **Structure names:** `<irrep>_p<pattern>_Q<amplitude>`, e.g. `T2-_p3_Q0.01`.
2. **Package name:** `axion_phonon/`, with readers in `axion_phonon/io/`.
3. **Saves:** keep only the untarred `tmp/`; drop `tmp.tar.gz`. CLAUDE.md is updated in step 6 of section 8.
4. **The self-test:** not a pipeline step. With the PAW fix it equals the base; it stays as a check in `validation/`.
5. **The remaining Julia `.amn` runs** (MBT p1, p2, p4, p3 at Q 0.03; all YIO) are redone with step 4 as they move
   into the new tree.
6. **Windows:** 10.698 / 13.915 becomes MBT production if the 8³/10³ mesh test agrees; otherwise it stays a study.

## 8. Order of work, with stop points

1. **This page:** done 2026-10-29.
2. **Code:** built 2026-10-29, `axion_phonon/` (1,830 lines with docstrings). Acceptance (`validation/restructure/`):
   - step 4, `accept_paw.py`: new `paw.py` vs `wannier/paw_overlap.py` on the local QE test saves, difference 0.0;
   - step 6, `accept_wannier.py` (+ `.log`), corrected pattern-3 pair (214300 / 218725), 100 tube k, complex128:
     new vs old A(k), curl identical, H 6e-16, dH/dk 4e-15, c2 4e-14 (H's class order differs); same 941 R vectors;
   - new code in production complex64 vs the cluster chunk (218727): c2 summed over the 100 k 2359.9810 vs 2359.9849
     (1.7e-6), single-precision level, as the old code vs the cluster;
   - unit tests `validation/unit/test_axion_phonon_{core,paw}.py`: 12 pass.
   Checked by you 2026-10-29.
3. **Pipeline:** built 2026-10-29, `pipeline/` (README there), copied to the cluster as `code/2026-10-29/` (with
   `MD5SUMS`). Steps 1-7 as section 2, with three refinements: pw2wannier90 is its own step and folder
   (`pw2wannier_k<N>/`, one `.mmn` per structure and mesh); the PAW overlaps are their own step (`overlap_k<N>/`, no
   window dependence), and step 5 makes the mode `.amn` from them; the BZ sum is step 7 (login node, seconds).
   - Dry run (`validation/restructure/dryrun_pipeline.sh` + `.log`), MBT base and T2- pattern 3, 8³, 10.278 / 15.169,
     in `/tmp` on the login node with production job folders as finished steps: the pattern-3 `nscf.in`, `.win`
     (atoms 0.0 Å from production) and `.nnkp` the scripts derive are identical to production's; the window checks
     pass; the tube k list is identical; step 7 on the production chunks (218727-218736) prints 6.43227 / 6.60007.
   - Job drivers run locally as the jobs run them, PythTB at the cluster's commit
     (`validation/restructure/accept_pipeline.py` + `.log`): `dtheta.py` on the production chunks 6.43227 / 6.60007;
     `c2.py` on 100 tube k of the corrected pair against the cluster chunk, sums agree to 1.7e-6 (single precision,
     as in step 2); `overlap.py` S-orthonormality 1.2e-14; `mode_amn.py` from the production overlap chunks and base
     `.chk` gives the 218725 `.amn` (28.5 M numbers; 178 differ by one unit in the 12th printed decimal, BLAS rounding).
   - Found and fixed on the way: `wannier90.x -pp` on the login node needs the Intel modules; the epilog gzips the
     overlap chunks (167 MB), so step 5 unzips them.
   Waiting for your check before step 4.
4. **MnBi₂Te₄ tree:** inventory done 2026-10-30 (`scripts/mbt_tree_inventory_2026-10-30.py`, read-only on the login
   node; list `mbt_inventory_2026-10-30.md` here; generated `scripts/mbt_moves_2026-10-30.sh`, not run). 3,000 GB now;
   929 GB production after the moves; 497 GB duplicate `tmp.tar.gz`; 1,575 GB superseded runs left in `archive/`
   (1,063 GB the 1e-7 `kfix` NSCFs, 322 GB old base 209653, 112 GB the first LDA pattern runs, 58 GB the Julia-.amn
   window study, 20 GB tests). Approved by the user with three choices (Julia-.amn tubes kept in `dtheta/` with a
   `_julia-amn` suffix; pattern 3's 8 overlap chunks gathered into one job folder; the 10³ NSCFs kept as tarballs).
   - **Done 2026-10-30** after the mesh test: `scripts/mbt_moves_2026-10-30.sh` (100 `mv`, exit 0; log
     `MnBi2Te4/archive/mbt_moves_2026-10-30.log`). Every production `OK` points at an existing job folder (the mesh
     test's empty markers were filled with their job ids). `MnBi2Te4/README.md` on the cluster and in `data/`.
   - **Local `data/MnBi2Te4/`** mirrored the same day (`scripts/local_mbt_tree_2026-10-30.py`, 69 moves): runs at the
     cluster's relative paths, the production pair's `.chk` and the c2 / overlap chunks moved there from the acceptance
     folder, local test data in `checks/`, the rest in `archive/` at its old path. 78 path references rewritten (plus
     hand edits in `amplitude_series.ipynb`, `fig1.ipynb`, two acceptance scripts); all ten notebooks that read MBT
     data rerun without error. Text backup `data/backup/text_before_local_mbt_2026-10-30.tar.gz`.
   - An inventory script lists every folder: production → its new place, everything else → `archive/`.
   - You check the list.
   - The moves are `mv` on the same filesystem, so they are instant and reversible.
   - Then the material README.
5. **Y₂Ir₂O₇ tree:** the same.
6. **Analysis:** point the notebooks and the paper-number scripts at the new paths, re-render one figure as a check,
   then update `notes/README.md`, `notes/code.md` and CLAUDE.md.
   - Code side done 2026-10-30 (ahead of 4-5, which wait for the mesh test; user): `modules/`, `wannier/`,
     `run_axion.py`, the 31 diagnostics folders, the setup templates and legacy drivers moved to `archive/` (README;
     backup `data/backup/code_before_restructure_2026-10-30.tar.gz`); `calculations/analysis` → `analysis/`,
     `calculations/figures` → `analysis/figures/`, setup notebooks → `setup/`, `qe_patches/` to the top; 250 path
     references in 105 notes and scripts rewritten. New: `axion_phonon/chern_simons.py` (the 3-form check),
     `wannier.run_model` (one run's model, H only), the QE band/PDOS readers in `io/qe.py`. Every notebook in
     `analysis/` imports only `axion_phonon`; `theta.ipynb` reruns to the same parities and η₄ = 2; validation
     notebooks run against `archive/code` via `_paths.py`. Unit tests: 33 pass. Paths into `data/` are unchanged.
   - `chern_simons.py` (`validation/restructure/accept_chern_simons.py` + `.log`, YIO 212615 / mode 1 212982): θ at
     nk 4 bit for bit equal to the 10-14 code at the same k batch (1e-12 with the new default batch, summation
     order); ∂θ/∂Q at nk 6 on the log's base-centre positions −0.018277 rad/Å (log −0.018277). On atom-pinned
     positions (production) it is −0.018233 at nk 6, 0.24% away: a different interpolant, inside the nk-6 mesh error
     (c₂ and the 3-form differ by 0.4% there; 0.09% at nk 10).

The running jobs (mesh test stage 1, window-study reruns) finish where they are. Their folders are sorted in step 4.
Nothing new is launched until step 3 is done.
