# Validation of the Wannier-interpolated axion response

Everything that checks `modules/` is here. Each directory is one **layer of the
argument**, ordered so that a check only depends on layers above it. If a check
fails, the layers above it tell you whether the cause is upstream.

```
unit/                synthetic models, no first-principles data, seconds
inputs/              H(R) and A(R) themselves, before any curvature
external_curvature/  Omega_ext against pythtb, postw90, WannierBerri, and raw Wilson loops
symmetry/            does the curvature transform correctly
observable/          does it reproduce a measured number (Z*)
results/             what every check last wrote
```

---

## The pipeline being validated

```
W90 run directory  (.chk .mmn .nnkp _hr.dat _r.dat _wsvec.dat _centres.xyz)
  |
  |  wannier/  ................................ stage 1.  Current route: run_plot_stage.sh makes Wannier90
  |    (Jae-Mo's fork) write _hr.dat / _r.dat; check_outputs.py checks them; PythTB reads them.
  |    The in-house .mmn route below is what most notebooks in this suite still use:
  |
  |  wannier/wannier_io.py .................... file readers + the pair loader
  |    read_wannier_checkpoint, read_mmn_overlaps, load_wannier_pair
  |
  |  wannier/position_matrix.py ................ A(R), one entry point
  |    build_model_with_position_matrix(w90, ..., source=...)
  |      source = mmn | rfull | rdat_pair | rdat_direct
  |      -> A(R) = <0i|r|Rj>, on the same R-set as H(R)   (+ _AA_cache/*.npz)
  v  model carrying H(R) and A(R); WannierPair for a base/mode pair
  |
  |  modules/curvature.py ..................... all k-space algebra (stage 2)
  |    position_terms(model)        -> X(R), curl weights, the fixed tau
  |    fields(model, pos, k)        -> H(k), dH/dk, A(k), curl A(k)
  |    connection(f, n_occ, beta)   -> d, a, A_vv, Omega_vv   (band gauge)
  |    omega_internal / omega_cross / omega_external -> the three pieces of Omega
  |                                    (cross + external = B^X + B^E, the position-matrix terms)
  |    omega(c)                     -> Omega = internal + cross + external
  v  Omega_mu,nu   (D x D antisymmetric, n_occ x n_occ; D = 3 for k only, 4 with beta)
  |
  |  modules/axion.py ......................... the beta direction and c2 (stage 3)
  |    beta_terms(...), dtheta(...)  ->  c2 = (1/16pi) eps Tr[Omega Omega]  ->  integral d3k
  v  dtheta/dbeta
```

The **external** terms are the piece this suite exists for. The internal (Kubo)
curvature only needs `H(R)`; the external terms need the off-diagonal Wannier
position matrix `A(R)`, and every way of getting that matrix wrong produces a
plausible-looking number. The math is in
[`notes/theory/external_terms.md`](../notes/theory/external_terms.md); the symbol map and the
list of things that fail silently are in
[`notes/theory/decomposition.md`](../notes/theory/decomposition.md).

---

## Start here

Everything except `unit/` is now a Jupyter notebook (`.ipynb`), not a script you
run with flags — see [Running them](#running-them) below for the one-time kernel
setup. If you are picking this up cold and want to know whether the external
terms are right, in this order:

```bash
python -m pytest validation/unit -q
```

then open `external_curvature/check_vs_pythtb.ipynb` and run it top to bottom
with `RUN_COMPUTATION = True`, then `observable/born_effective_charge.ipynb`.

The pytest suite is seconds and needs no data. `check_vs_pythtb` is the primary
correctness gate and takes a few minutes. `born_effective_charge` is the only
check against a physically measured number; read its first markdown cell before
committing to a run.

---

## What each check proves

### `unit/` — no first-principles data, runs in seconds

| Script | Proves |
| --- | --- |
| `test_curvature.py` | `modules/curvature.py` and `modules/axion.py` on analytic models and fixed-seed fakes: each piece of `Omega` and the total against exactly solvable moving frames (gauge *covariant*, not invariant, including the beta versions), Hermiticity and antisymmetry, the Kubo equation, eq 38 = eq 39, the Bloch sums against a brute-force sum over R, `A_beta` Hermitian and equal to its own k derivative, and `dtheta = 0` for a k-independent model. |
| `test_wannier_io.py` | `.mmn` parsing, finite-difference weights, the `A(R)` cache round-trip, and `postw90` position-matrix installation. |
| `test_pipeline_end_to_end.py` | The file boundaries: what `load_wannier_pair` refuses, and the layout `save_dtheta_results` writes. |
| `test_position_matrix_compare.py` | The comparison helper in `inputs/validate_transl_inv_full_lib.py` reports the right element and magnitude. |

### `inputs/` — `H(R)` and `A(R)`, before any curvature

These isolate the layer *upstream* of `Omega`, so a failure here means the
curvature code is being fed something wrong rather than computing it wrong.

| Notebook | Proves |
| --- | --- |
| `check_interpolation_inputs.ipynb` | `H(R)` reproduces the checkpoint Hamiltonian on the source mesh; the `A(k)` inverse FFT is exact; mesh-node values are invariant under the orbital-aware R remapping; the on-disk `A(R)` cache matches; analytic `dH` and `curl(A)` match central differences. |
| `check_raw_mmn_wilson.ipynb` | Two Wannierizations span the *same* occupied subspace on the source mesh, using raw `.mmn` link matrices and Wilson loops only — no `A(k)`, no `A(R)`, no interpolant. Establishes what trial-to-trial differences downstream *cannot* be blamed on. |
| `compare_position_sources.ipynb` | How the in-house `.mmn` route and postw90's `_r_full.dat` differ: whether it is a choice of aliasing representative (leaves `A(k)` untouched on the mesh) or a different discretisation (changes it everywhere). Only the latter can move a converged observable. |
| `validate_transl_inv_full.ipynb` (lib: `validate_transl_inv_full_lib.py`) | Every complex component of `<0i\|r\|Rj>` agrees between a Wannier90 export and the `.mmn` construction, after putting both in one convention. `compare_matrices` is also imported by `unit/test_position_matrix_compare.py`, so its logic lives in the `_lib.py` module, not the notebook. |
| `compare_pr702_rdat.ipynb` (lib: `compare_pr702_rdat_lib.py`) | Wannier90 PR #702's `_r.dat` is postw90's `A(R)`, rounded to `6F12.6` — bit-exact after rounding, not just close. `read_matrix`/`fortran_round` are also imported by `compare_jaemo_ndegen_rdat.ipynb`. |
| `compare_jaemo_ndegen_rdat.ipynb` | Keyed (not positional) comparison for Jae-Mo's `write_ndegen_applied` `_r.dat`, whose folded R list can't be matched row-for-row the way `compare_pr702_rdat` does. |

### `external_curvature/` — the gates on `Omega_ext`

| Notebook | Proves |
| --- | --- |
| `check_vs_pythtb.ipynb` | **The primary gate.** On the real Y2Ir2O7 production pair: our batched `Omega_ext` equals pythtb's `berry_curvature(include_external=True) - (include_external=False)`, which is itself validated against postw90 `J0+J1`. Also assembles the mixed k–beta block by hand with a nonzero `A_beta` and checks all four cross terms, the curl, and `-i[A_i, A_beta]`. Nothing is persisted to `results/` — it is pure assertions, so its notebook has no cached-load path. |
| `check_curvature_vs_wannierberri.py` (a script, quick) | The rewritten `modules/curvature.py` against WannierBerri's non-Abelian block `formula.covariant.Omega.nn` via `wb.evaluate_k(..., formula=..., iband=...)`, both fed from the same Wannier90 run (SrTiO3 trial04: `_tb.dat` for WannierBerri, `_hr.dat`/`_r.dat` for us). Eigenvalues of every plane agree to <= 3e-8 relative, for the total, the internal piece and cross + external separately, with the same sign; `H(k)` agrees to 4e-14 eV. trial04 is occupied-only, so band *subsets* (lowest 30, 20) are used to give a non-empty complement and exercise the internal and cross terms. k planes only: WannierBerri has no beta direction. |
| `check_vs_postw90_line.ipynb` | Our curvature against postw90's `kpath_task = curv` at postw90's own k-points, on identical `H(R)`/`A(R)` inputs. |
| `check_vs_postw90_full_mesh.ipynb` | `.mmn` and postw90 position matrices compared on complete source meshes, peak-normalized so zeros of the reference do not produce singular percentages. |
| `check_vs_wilson_flux.ipynb` | The continuum `integral Tr Omega` over each source plaquette, by Gauss–Legendre quadrature, against the raw `.mmn` Wilson-loop phase. Separates finite-link discretisation error from interpolation error. |

Note what is **not** covered: `check_vs_pythtb.ipynb` tests the k–beta assembly's
*algebra*, but the mixed planes have no published implementation to compare
numbers against. Their accuracy rests on `observable/` and on the analytic
moving-frame tests (`unit/test_curvature.py`). The two endpoint estimates of
`dtheta` are a linearity test, not an accuracy measure. See §10.4 of
[`notes/theory/derivation.md`](../notes/theory/derivation.md).

### `symmetry/` — does the curvature transform correctly

| Notebook | Proves |
| --- | --- |
| `check_magnetic_symmetry.ipynb` | Y2Ir2O7 base and mode 1: the Cartesian dual trace curvature obeys `omega(Rk) = det(R) R omega(k)` for unitary operations and `omega(-Rk) = -det(R) R omega(k)` for antiunitary ones, plus the invariant-line consequences (base z axis: all components vanish; mode 1 z axis: **all** components vanish, C2z gives `omega_x = omega_y = 0` and the antiunitary S4z T removes `omega_z`; [111]: `omega_x = omega_y = omega_z`; [1,-1,0]: `omega_x = omega_y`). Two calibrations: S4z **without** T is a control that must be violated (relative residual near 2, against 0.04 to 0.13 for the real symmetries), and `mode1_inversion` shows the 0.01 A displacement breaks inversion below the floor. `POSITION_SOURCE` = `mmn`, `rdat` or `jaemo` (Jae-Mo's `<run>/jaemo/`, read as written). The residual floor is the shared H(R)'s, not the curvature code or the position-matrix route: [`notes/log/2026-09-21_yio_symmetry.md`](../notes/log/2026-09-21_yio_symmetry.md). |
| `check_convention_symmetry.ipynb` (lib: `check_convention_symmetry_lib.py`) | Four internally consistent interpolation conventions (`mmn_pair`, `rdat_direct`, `rdat_pair`, `rfull_pair`) evaluated on paths where symmetry forces `Tr Omega_xy = 0`. Asserts `H(R)` and `A(R)` share an R-set before trusting any of them. Also imported by `check_pr702_rdat_symmetry.ipynb`. |
| `check_trial_curvature_geometry.ipynb` (lib: `check_trial_curvature_geometry_lib.py`) | Per-interpolant: how far apart two trials' curvatures are at the source nodes, and whether each separately obeys the surviving crystal symmetries *between* nodes. Also imported by `check_rdat_ndegen_applied_symmetry.ipynb`. |
| `compare_trial_curvature.ipynb` (lib: `compare_trial_curvature_lib.py`) | The occupied traced curvature for two Ti_Q_2A trials on a generic line. Imported as a library by `check_vs_wilson_flux.ipynb`, `check_trial_curvature_geometry.ipynb`, `check_rdat_ndegen_applied_symmetry.ipynb`, and `compare_position_matrix_sources.ipynb`. |
| `check_pr702_rdat_symmetry.ipynb` | Whether PR #702's `_r.dat` changes the *symmetry* verdict from `check_convention_symmetry`, not just the file-level residuals `compare_pr702_rdat` measures. |
| `check_rdat_ndegen_applied_symmetry.ipynb` | Same symmetry gate for Jae-Mo's `write_ndegen_applied` `_r.dat` route (no `.mmn` read at all). |
| `compare_position_matrix_sources.ipynb` | Unlike `compare_trial_curvature` (two different Wannierizations), this holds the Wannierization fixed and swaps only the `A(R)` source, isolating what a position-matrix choice alone can move. |

### `observable/` — against a measured number

| Notebook | Proves |
| --- | --- |
| `born_effective_charge.ipynb` (lib: `born_effective_charge_lib.py`) | `Z*(Ti)` in cubic SrTiO3: QE DFPT versus the BZ integral of the mixed `(k,u)` curvature. **The only quantitative check on the mixed planes.** Splits `B = B_int + B_cross + B_ext` and integrates each separately; only the sum is physical. Note `Z^ion,eff` is *not* `z_valence` — semicore bands outside the Wannier window ride with the displaced atom. Also imported by `audit_mixed_born_response.ipynb`. |
| `audit_mixed_born_response.ipynb` | Reads the NPZ artifacts above and separates executable pass/fail checks from accuracy gates that would need new first-principles data. Does not recompute the BZ integrals. |
| `summarize_convention_results.ipynb` | Consolidates the Born-charge runs across trials and conventions into one table and plot. |
| `plot_curvature_generic_line.ipynb` | `Tr[Omega_xy]` for the three SrTiO3 trials on the generic line `k(t) = (t, 0.30, 0.15)`. |

---

## Running them

`unit/` is still plain pytest, needs no data, and runs in seconds:

```bash
python -m pytest validation/unit -q
```

Every other check is a Jupyter notebook (`.ipynb`), meant to be opened and run
cell by cell, not invoked with flags. One-time setup:

```bash
conda activate axion
python -m ipykernel install --user --name axion   # registers the "axion" kernel
jupyter lab   # or: jupyter notebook, or open the .ipynb in VS Code / Cursor
```

Open a notebook, pick the **axion** kernel, and run top to bottom. Each one
follows the same shape:

1. **Imports** — a `sys.path` bootstrap that walks up to find `validation/`
   (no `__file__` inside a notebook), then the module's imports.
2. **Config** — every former CLI flag is a plain variable here instead, with
   its old default. Edit these directly rather than passing `--flag value`.
3. **`RUN_COMPUTATION = False`** — the check's headline setting. `False` (the
   default) skips straight to reloading the last saved run from
   `results/<check name>/` and plotting it — the fast path for "what did we
   last measure?". `True` rebuilds everything from the first-principles data
   named in `_paths.py`, and re-saves. A few checks (`check_vs_pythtb` is the
   main one — it's pure assertions) never persisted anything, so their
   "loaded" path is a markdown note rather than a fake reload — read the
   first cell if a notebook's load path looks like a no-op.
4. **The computation, unrolled** — what used to be a single `main()` is now
   one cell per stage (build model → connection → curvature → compare →
   save), each guarded by `if RUN_COMPUTATION:`. Run them and the intermediate
   objects — `model`, `connection`, the raw curvature tensors, not just a
   pass/fail summary — stay in the kernel as normal variables you can
   `.shape`, plot, or index into from a new cell below.
5. **Plot** — rendered inline (`plt.show()`), same figure the old script
   used to save as a PNG.

**Shared logic lives in `<name>_lib.py`, next to the notebook that owns it.**
A handful of checks are also imported by other checks (`compare_trial_curvature`
by four others, `compare_pr702_rdat` by `compare_jaemo_ndegen_rdat`,
`check_convention_symmetry` by `check_pr702_rdat_symmetry`,
`check_trial_curvature_geometry` by `check_rdat_ndegen_applied_symmetry`,
`born_effective_charge` by `audit_mixed_born_response`) or, in one case, by a
`unit/` pytest test (`validate_transl_inv_full`'s `compare_matrices`, needed by
`unit/test_position_matrix_compare.py`). For each of those, the plain functions
and constants live in a `*_lib.py` module with no argparse and no `__main__`;
the notebook imports from it just like every other importer does. If you need
one of these functions from a script instead of a notebook, `import` the
`_lib` module directly — nothing about it is notebook-specific.

**Data locations are written out once.** Every run folder a check reads is a
full path under `data/`: either a constant in [`_paths.py`](_paths.py)
(`STO_BASE_RUN`, `YIO_MODE1`, ...) or a `DATA / "..."` path in the check's own
config cell. To repoint a check at a new run, edit that path. Missing inputs raise
a readable error via `_paths.require`.

**pythtb branch matters.** Everything that touches `A(R)` needs a pythtb build
with the Wannier position-matrix API (`W90.pos_r`, `TBModel._pos_r`, read from
`prefix_r.dat`). That lives on pythtb's **`external`** branch; `main` and the
topic branches do not have it, and a model built from them fails with a
`no '_pos_r'` message naming this. `validation/unit` does not need it. Check
with `git -C <pythtb> branch --show-current`. The `axion` conda environment is
already pinned to the right branch.

**Results** go to `results/<check name>/`, not next to the data. That
directory is what `RUN_COMPUTATION = False` reloads.

The `.mmn` files these notebooks stream are multi-gigabyte and are not in
version control; several checks read the much smaller `_AA_cache/*.npz`
archives instead. See
[`notes/log/2026-08-27_rfull_vs_mmn_position_matrix.md`](../notes/log/2026-08-27_rfull_vs_mmn_position_matrix.md)
for which inputs are genuinely required and which were retired.

---

## Where things came from

This tree was assembled from scripts previously scattered across three
calculation directories. Old path -> new path:

| Was | Is |
| --- | --- |
| `tests/test_berry_curvature_covariance.py` | `unit/test_berry_curvature_covariance.py`, later folded into `unit/test_curvature.py` |
| `tests/test_wannier_io.py` | `unit/test_wannier_io.py` |
| `tests/test_refactored_pipeline.py` | `unit/test_pipeline_end_to_end.py` |
| `tests/test_transl_inv_full_validator.py` | `unit/test_position_matrix_compare.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_interpolation_inputs.py` | `inputs/check_interpolation_inputs.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_raw_mmn_wilson.py` | `inputs/check_raw_mmn_wilson.py` |
| `calculations/SrTiO3/compare_position_sources.py` | `inputs/compare_position_sources.py` |
| `calculations/SrTiO3/validate_transl_inv_full.py` | `inputs/validate_transl_inv_full.py` |
| `calculations/Y2Ir2O7/phonon/validate_external_terms.py` | `external_curvature/check_vs_pythtb.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_postw90_reference.py` | `external_curvature/check_vs_postw90_line.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_postw90_full_mesh.py` | `external_curvature/check_vs_postw90_full_mesh.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_interpolant_plaquette_flux.py` | `external_curvature/check_vs_wilson_flux.py` |
| `calculations/Y2Ir2O7/phonon/mode1_symmetry_checks.py` | `symmetry/check_magnetic_symmetry.py` |
| `calculations/SrTiO3/analyze_convention_symmetry.py` | `symmetry/check_convention_symmetry.py` |
| `calculations/SrTiO3/Ti_Q_2A/check_trial_curvature_geometry.py` | `symmetry/check_trial_curvature_geometry.py` |
| `calculations/SrTiO3/Ti_Q_2A/compare_trial_curvature.py` | `symmetry/compare_trial_curvature.py` |
| `calculations/SrTiO3/born_effective_charge.py` | `observable/born_effective_charge.py` |
| `calculations/SrTiO3/audit_mixed_born_response.py` | `observable/audit_mixed_born_response.py` |
| `calculations/SrTiO3/summarize_convention_results.py` | `observable/summarize_convention_results.py` |
| `calculations/SrTiO3/plot_berry_curvature_generic_line.py` | `observable/plot_curvature_generic_line.py` |

Result directories moved too: `interpolation_validation -> results/interpolation_inputs`,
`raw_mmn_wilson_comparison -> results/raw_mmn_wilson`,
`postw90_reference_comparison -> results/postw90_line`,
`postw90_full_mesh_comparison -> results/postw90_full_mesh`,
`plaquette_flux_interpolant_comparison -> results/wilson_flux`,
`curvature_trial_comparison -> results/trial_curvature`,
`trial_curvature_geometry -> results/trial_geometry`,
`convention_analysis -> results/convention_symmetry`,
`mixed_born_validation -> results/mixed_born`,
`mode1_symmetry -> results/magnetic_symmetry`
(and the `_nk12` / `_rfull` / `_full_A_support` variants alongside each).

`calculations/` later held only first-principles data and the production drivers. Since 2026-09-24 the
data is in `data/` (git-ignored) and `calculations/` holds the run registry, `run_axion.py` and the analysis
notebooks; see `calculations/README.md`.

Every CLI script above (everything except `unit/`) was later converted
one-for-one to a `.ipynb` of the same name — `check_interpolation_inputs.py`
became `check_interpolation_inputs.ipynb`, and so on. Where a script's
functions were imported by another check, the functions moved to a `*_lib.py`
module next to it and the notebook imports from that; see
[Running them](#running-them) above. No script's logic changed in the
conversion, only how it's invoked.

---

## Related notes

| Note | What it holds |
| --- | --- |
| [`external_terms.md`](../notes/theory/external_terms.md) | The math of the external terms as implemented, with the code boundaries. |
| [`theory/decomposition.md`](../notes/theory/decomposition.md) | Symbol map, the two things that fail silently, and the 13-point check list mapped to these scripts. |
| [`theory/derivation.md`](../notes/theory/derivation.md) | The full derivation. §10.4 covers what the k–beta blocks cannot be checked against. |
| [`code.md`](../notes/code.md) | Which module owns which stage of the pipeline. |
| `wannier/position_matrix.py` (docstring) | The four A(R) sources, what each reads, and how to add one. |
| [`2026-08-20_handoff_wannier_conventions.md`](../notes/log/2026-08-20_handoff_wannier_conventions.md) | Why the aliasing-class representative is a choice of *interpolant*, not a gauge, and why mixing two choices is fatal. |
| [`2026-08-27_rfull_vs_mmn_position_matrix.md`](../notes/log/2026-08-27_rfull_vs_mmn_position_matrix.md) | Where `<0i\|r\|Rj>` should come from; why `_r.dat` cannot carry it. |
| [`2026-08-27_Ti_Q_2A_full_validation.md`](../notes/log/2026-08-27_Ti_Q_2A_full_validation.md) | The written verdict of the `inputs/` and `external_curvature/` checks on Ti_Q_2A. |
