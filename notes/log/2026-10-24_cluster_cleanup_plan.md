# Cluster cleanup plan (2026-10-24)

Inventory of `/mnt/wkdk/bc889` on 2026-10-24, made with `du`, `find -size +100M` and a `du` of every QE scratch folder.
Nothing was running and **nothing has been deleted.** The user deletes by hand.

## Where the space is

- **Total: 10.9 TB, against the 2 TB target in CLAUDE.md.**
  - `axion/Y2Ir2O7`: 6.7 TB.
  - `axion/MnBi2Te4`: 3.5 TB.
  - `axion/SrTiO3`: 0.43 TB.
  - `envs/`, `julia/`, `bin/`, `code/` and `pseudo/`: under 5 GB together. All of these stay.
- **By content:**
  - 7.0 TB is QE scratch folders (`tmp/`, `tmp_keep/`, `tmp_bands/`, wave functions).
  - About 1.3 TB is `tmp.tar.gz` archives.
  - About 0.2 TB is `.mmn`/`.amn` files.
  - Everything else (outputs, `.chk`, `_hr.dat`, logs) comes to under 0.5 TB.
- **`axion/INDEX.md` and the `_STATUS` files date from 10-17.** They predate the conv_thr 1e-10 rebuilds, which are
  production now, so some folders marked PRODUCTION there are superseded.

Sizes below are in GiB, as `du -h` prints them.

## Tier 1: delete the whole folder

These are wrong or superseded, and no paper number needs their files. Total about 4.3 TB.

**MnBi2Te4**
- **`MnBi2Te4/phonon` (935):** PBE T1+ and T2- runs, wrong because of bug #18.
- **`MnBi2Te4/phonon_lda/T2-/mode_{1,2,3,4}/Q_0p01` (4 × 28):** run on the unpatched binary (bug #19). The bug's evidence
  is kept in the `Q_0p002_*` test folders.
- **`MnBi2Te4/base/wo_soc` (0.3):** not used.
- **`MnBi2Te4/base/soc_lda/u_4.0/209653bc889/w90_trial_sp_d` (14):** the centre-wrap run (bug #20). Its `tmp.tar.gz` is
  a hard link of the one in `..._recentred/`, so deleting it frees about 1 GB until that copy goes too (tier 2).
- **Keep:** `MnBi2Te4/base/soc` (2). It holds the bug #18 diagnosis runs, and test D is the PBE state behind the 10-22
  PBE bands.

**Y2Ir2O7** (pre-segni-fix PBE runs and tests; every number from them was replaced by LDA)
- **`Y2Ir2O7/wosoc` (131): moved to tier 2 (user, 2026-10-26).** Kept for history: the inputs, `scf`/`nscf`/`bands`
  outputs, the 61 PDOS files and the projwfc output (about 360 MB). Its `tmp/`, `tmp_bands/`, two `tmp.tar.gz` and
  executables go.
- **`mode_1/` (803):** `192345bc889` 518, `207547bc889` 161, `178365bc889` 62, `207547_nscf_save_2026-10-07` 62.
- **`mode_2/` (443):** `fromYP_2026-10-08` 154 (its SCF never converged), `nscf_nk_test_2026-10-06` 153,
  `Q_0p005_fromYP_2026-10-09` 105, `209040bc889` 30, `Q_0p01_from_Q0p005_2026-10-11` 1, and the three small folders
  `forces_209040*` and `wannier_templates_from_192346`.
- **`mode_3/` (1367):** `178377bc889` 510, `192347` 360, `174937bc889` 313, `nscf_david_nk8_2026-10-07` 119,
  `207549bc889` 65.
- **`base/u_3.0/` (553):** `192344bc889` 417, `207546bc889` 128, `stall_test_2026-10-07` 5, `209049bc889`,
  `test_newqe_48`, `test_nscf_sym_2026-10-05`, and `lda_test_2026-10-08/YL1_lda_pbepp` 1.

## Tier 2: delete the wave functions, keep the outputs

These runs are superseded, but the paper cites them or they are the record behind a result. Delete their `tmp*/`
folders, `*.tar.gz`, `.mmn`/`.amn` and `UNK*`. Keep `*.out`, `*.wout`, `.chk`, `.eig`, `.win`, `_centres.xyz` and
`NOTE`. Total about 3.1 TB.

| folder | frees | why the outputs stay |
|---|---|---|
| `MnBi2Te4/phonon_lda/T2-/mode_{1,2,3,4}/Q_0p01_kfix_2026-10-11` | 4 × 265 | conv_thr 1e-7 (superseded by the 10-19 rebuild); the 1e-7 tube pilot behind the paper's "first calculation" numbers |
| `MnBi2Te4/base/soc_lda/u_4.0/209653bc889` (incl. `nscf_k8_2026-10-12`, `tmp/`, `tmp_bands/`) | 331 | the original LDA base SCF; the 4×4×4 θ = π parity Wannier run; the 10-18 base symmetry test (undistorted row of the MBT symmetry table) |
| `Y2Ir2O7/phonon/mode_{1,2,3}/fromYL2_2026-10-12` | 3 × 161 | Wannier runs 212982/212994/212995: the paper's YIO symmetry table and the pinned check |
| `Y2Ir2O7/phonon/base/u_3.0/lda_test_2026-10-08/YL2_lda_relpz/nscf_from_210527` | 237 | 212615: the base row of the YIO symmetry table |
| `.../lda_test_2026-10-08/YP_pbe` and `Y2Ir2O7/phonon/mode_{1,3}/fromYP_2026-10-08` | 244 + 2 × 207 | the PBE comparison quoted in the paper (3% and 11%); keep `210525bc889` (the source of the 10-22 PBE bands) and `dft_bands_pbe_2026-10-22` |
| `SrTiO3` | 429 | early Born-charge validation, not in the paper; the validation notebooks reload from the repo's `results/` |

Commands for one folder `F` (dry run first, then delete):

```bash
find F -type d \( -name 'tmp*' -o -name '*.save' \) -prune -print
find F -type f \( -name '*.tar.gz' -o -name '*.mmn' -o -name '*.mmn.gz' -o -name '*.amn' -o -name '*.amn.gz' -o -name 'UNK*' \) -print
# after checking the lists:
find F -type d \( -name 'tmp*' -o -name '*.save' \) -prune -exec rm -rf {} +
find F -type f \( -name '*.tar.gz' -o -name '*.mmn' -o -name '*.mmn.gz' -o -name '*.amn' -o -name '*.amn.gz' -o -name 'UNK*' \) -delete
```

Every `tmp*` folder on the cluster is QE scratch (118 `tmp`, 25 `tmp_bands`, 12 `tmp_keep`), so the pattern catches
nothing else.

## Tier 3: production runs. Decide first

After tiers 1 and 2 about 3.4 TB remains, still above 2 TB.

**MBT 1e-10 runs:** base `conv1e-10_k8_2026-10-19` (214290), modes 1-4 at Q 0.01 and mode 3 at Q 0.03.
- Each run holds an 83 GB `tmp.tar.gz` and an 86 GB untarred `tmp/` with the same save.
- CLAUDE.md asks for both. The `tmp/` alone serves pw2wannier, the `.amn` generation and a new NSCF from its charge
  density.
- Dropping the six archives frees about 500 GB. This needs the user's OK, because it changes the CLAUDE.md rule.
- MBT is still active (mesh test, window studies), so the `tmp/` folders stay.

**YIO amplitude series:** `mode_{1,2,3}/Q_0p0{1,2,3,5}_conv1e-10_2026-10-17`.
- Each holds a 108 GB `tmp_keep/`, the NSCF save, which is needed only to project a new `.amn`.
- The production numbers, the Chern-Simons θ(Q) series and a switch to pinned positions need only the `.chk` and
  `.eig`, and those stay.
- Deleting the wave functions of Q 0.02, 0.03 and 0.05 (9 runs) frees about 970 GB.
- Keep the charge density and xml, so that an NSCF can be redone:

```bash
find .../tmp_keep -name 'wfc*' -print   # check first; only the wave-function files go
```

- Keeping Q 0.01 and the 1e-10 base `conv1e-10_2026-10-17` (240) leaves one complete YIO set.

**Staged copies of `.chk` files (about 70 GB):**
- MBT `dtheta_tube_conv1e-10_k8_2026-10-19/<pair>/{base,mode}/` and `window_study_2026-10-22/*/tube/{base,mode}/`:
  uncompressed copies of `.chk`/`_hr.dat` whose gzipped originals sit in the Wannier job folders. Keep `chunks/`,
  `klist.npz` and the logs (Fig. 7c).
- MBT `window_study_2026-10-22/*/base/ref/`.
- YIO `window_study_2026-10-19/*/base/ref/`, `ref_prodwin/`, `ref_ampwin/` and `selftest_ampwin_2026-10-18/base_ref/`.

After all three tiers about 1.9 TB remains, or about 1.5 TB if every YIO `tmp_keep` goes.

## What is still needed, and for what

| need | files |
|---|---|
| MBT mesh test (new NSCF on 10³/12³) | the 1e-10 base and mode 3 `tmp/` (charge density) |
| MBT window or trial studies | the 1e-10 `tmp/` of the base and modes, plus the base `.mmn` and `.eig` |
| YIO pinned switch, symmetry tables, θ(Q) | `.chk` + `.eig` of the 1e-10 self-test base 214249, the ampwin10 mode runs, 212615 and 212982/94/95 |
| MBT frequencies | `forces_2026-10-09`, `forces_conv1e-10_2026-10-17`, `T2-/mode_*/Q_0p03_2026-10-15`, `rigid_layer` (all small) |
| Paper Fig. 7c | MBT tube `chunks/` and `klist.npz` |

After the cleanup, rerun `scripts/cluster_mark_status_2026-10-17.sh` with the updated list, so that `INDEX.md` and the
`_STATUS` files match.

## Progress (2026-10-26)

- **Total 8.9 TB** after the user's first deletes (YIO `mode_3` superseded runs, `wosoc` scratch).
- **Copied to `data/` before deleting** (text outputs only; `data/Y2Ir2O7/archive/README.md`): mode 3 192347 SCF/NSCF
  (#16) and 207549 (#17); mode 1 207547 (#17), the trial_01 projection-only model of 178365 (#1), 192345 bands; the PBE
  mode force SCFs 210611 / 210615 (log 10-08 §12.5) into `phonon/GM2-/mode{1,3}/Q_0.01A/pbe/forces_<job>/`.
- **`fromYL2_2026-10-12` and `fromYP_2026-10-08` can go whole** (every mode): their outputs, `.chk`, `.eig`, `_hr.dat`
  are in `data/`, and the YIO symmetry table is rebuilt on the 1e-10 runs
  (`archive/diagnostics/2026-10-26_yio_symmetry_conv1e-10/`). `setup_amplitudes_YL2.sh` now reads its Q = 0.01
  positions from `Q_0p01_conv1e-10_2026-10-17` (identical).
- **Hidden in the earlier listing:** each `fromYL2`/`fromYP` job folder also holds a 50.6 GB `tmp.tar.gz` (the `tmp*`
  prune pattern matched files too). The 13 MB `tmp.tar.gz` in the `Q_*` and `from*` submit folders are the start
  densities.
- **YIO symmetry table rebuilt on 1e-10 (10-26):** paper table updated (cross 0.47-0.52%, total unchanged), so
  `fromYL2` has no remaining user.
- **Mode 2 copied (10-26):** every superseded PBE folder's text outputs (the #17 stall chain, 209040, fromYP, the
  Q = 0.005 tests, 212707) into `data/Y2Ir2O7/archive/mode2/`; only the input templates `wannier_templates_from_192346`
  were not copied.
