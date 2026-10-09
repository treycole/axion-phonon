# Code map

Three stages, one direction of dependence:

```text
wannier/                  STAGE 1  Wannier90 -> seed_hr.dat, seed_r.dat, seed_wsvec.dat
  README.md  run_plot_stage.sh  check_outputs.py
  wannier_io.py  position_matrix.py      the older in-house route to A(R) (legacy, see below)
        |   (PythTB, external branch, reads the files:  W90(dir, seed).model(orb_vecs=tau))
        v
modules/curvature.py      STAGE 2  the curvature of the occupied bands, note §3-§7
        |
        v
modules/axion.py          STAGE 3  finite difference in beta, c2, BZ integral, note §8-§9
        |
        v
calculations/run_axion.py    the whole calculation on one page: edit the settings block at the top, then run it
```

`run_axion.py` has no command-line options: the two run folders, their A(R) caches, `N_OCC`, `DBETA`, `NKS`,
`SOURCE` and `SHARED_WS` are variables at the top of the file, with the choices listed next to them.
Paths into `data/` are written out in full. `calculations/README.md` explains the layout of `calculations/`.

The dated notes up to 2026-09-23 give older forms of the same run (environment variables on
`calculations/Y2Ir2O7/phonon/run_axion.py`). They map one to one onto the settings: `AXION_POSITION_SOURCE=x` is
`SOURCE = "x"`, `AXION_SHARED_WS=1` is `SHARED_WS = True`, `AXION_NKS=4,6` is `NKS = [4, 6]` and
`AXION_OUTPUT_DIR=d` is `OUT = "d"`. The output folders are named as before.

Nothing in `modules/curvature.py` or `modules/axion.py` reads a `.chk`, `.mmn`,
`.amn` or a Wannier90 file; they see a PythTB model. Everything they compute is
a transcription of `theory/derivation.md` (the "note"), with the
note's equation numbers in the docstrings.

## Using it yourself

```python
from pythtb import W90
from modules.axion import dtheta

base, mode = W90("base_dir", "seed"), W90("mode_dir", "seed")
tau = base.lattice.orb_vecs                      # ONE fixed tau for both (note A3)
model_base = base.model(orb_vecs=tau, min_hopping_norm=1e-5)
model_mode = mode.model(orb_vecs=tau, min_hopping_norm=1e-5)
result = dtheta(model_base, model_mode, nk=12, dbeta=0.01, n_occ=156)
result["dtheta"]                                 # (2,) base-end and mode-end estimates
```

To look at the curvature itself for one structure, see the top docstring of
`modules/curvature.py`.

## Stage 2: `modules/curvature.py` (seven functions)

Direction index first: `dH`, `A` are `(3, nk, J, J)`, `Omega` is
`(D, D, nk, Nv, Nv)` with `D = 3` (kx, ky, kz) or `4` (adding beta).

| function | returns | note |
|---|---|---|
| `position_terms(model)` | `X(R)=<0s\|r-tau_s\|Rt>`, bond vectors `b=R+tau_t-tau_s`, curl weights; once per structure | (8), (9) |
| `fields(model, pos, k)` | Wannier gauge, k space: `H`, `dH/dk` (from PythTB), `A(k)`, plain curl of `A` | (7)-(9) |
| `connection(f, n_occ, beta=None)` | band gauge: internal `d_cv`; external `a_cv`, `A_vv`, `Omega_vv` | (10), (28), (31), 7.2 |
| `omega_internal(c)` | `i[d_mu^† d_nu - d_nu^† d_mu]` | (29) |
| `omega_cross(c)` | `(d_mu^† a_nu + h.c.) - (mu<->nu)`, no factor of i | (30) |
| `omega_external(c)` | `Omega_vv - i[A_vv_mu, A_vv_nu]` | (27) |
| `omega(c)` | internal + cross + external | (24) |

`beta` (a `BetaTerms`) adds the fourth direction: `dH/dbeta`, the mixed curl
`Omega_{l,beta}`, and optionally `A_beta` (the parametric connection). Without
`A_beta` the run is in the **frozen Wannier gauge**, `A_beta = 0`, recorded in
`connection(...).frozen_gauge`; know which one a run is using.

PythTB supplies `H(k)`, `dH/dk`, the model reading and the tau embedding. It also
has its own `berry_curvature(non_abelian=True, include_external=True)` for the
three k planes, which returns cross and external summed; that is the independent
test of `omega_internal` and `omega_cross + omega_external`
(`validation/external_curvature/check_vs_pythtb.ipynb`).

## Stage 3: `modules/axion.py`

| function | does | note |
|---|---|---|
| `dtheta(model_base, model_mode, nk, dbeta, n_occ, ...)` | one loop over k-batches: `fields` for both structures, `beta_terms`, `connection` and `omega` at each end, `c2`, integral | §8-§9 |
| `beta_terms(f_base, f_mode, dbeta, ...)` | `dH/dbeta = (H_mode-H_base)/dbeta`, `dA/dbeta`, mixed curl `-dA/dbeta` (frozen) or `d_l A_beta - dA_l/dbeta` | 8.3, (35) |
| `c2_density(Omega)` | `(1/16 pi) eps^{mu nu rho sigma} tr[B_mu nu B_rho sigma]` | (38); equals (39) |
| `integrate_c2(c2)` | mesh sum and Simpson over the unit cube of reduced k | (38) |

The two entries of `dtheta` share **one** secant `dH/dbeta`, `dA/dbeta`. They differ
only in where the spatial fields are evaluated (beta = 0 or beta = dbeta). There is
no right endpoint and their spread is not a validity check.

## Conventions that fail silently

* **tau is one fixed array for both structures** (note A3). `dtheta` refuses
  models with different `orb_vecs`. The Wannier centres move physically through
  `X(R)`, not through tau.
* **Frozen gauge is explicit** (above).
* **`H(R)` and `A(R)` must come from the same Wannier90 run** with
  `use_ws_distance`, `transl_inv_full` and `write_ndegen_applied` all set
  (`wannier/README.md`). `wannier/check_outputs.py` verifies it.
* **Base and mode choose their Wigner-Seitz images separately** (each from its own
  centres). Measured, effect on dtheta not yet known:
  `log/2026-09-20_jaemo_ndegen_yio_mode1.md`.

## Verification

* `validation/unit/test_curvature.py`: the analytic covariance checks (each piece
  and the total against exact moving frames, including the beta versions and the
  non-Abelian commutator), (38) = (39), the Bloch sums against a brute-force sum
  over R, and the parametric connection A_beta against its own k derivative.
* Independent codes on the k planes: PythTB's `berry_curvature` (real YIO base, 3.6e-16) and
  WannierBerri's non-Abelian `Omega.nn` block (SrTiO3 trial04, <= 3e-8, internal and
  cross + external separately): `validation/external_curvature/check_curvature_vs_*.py`.
* Data regression: `run_axion.py` on the YIO pair reproduces the saved results
  of the older driver (see the regression note in `log/2026-09-20_jaemo_ndegen_yio_mode1.md`).

## Kept, not used by the new path

**Stage 1, the in-house route to `A(R)`** (`wannier/wannier_io.py`,
`wannier/position_matrix.py`). Python readers for the `.chk`, `.mmn` and `.nnkp`,
the `Lambda(R)` loader, the base/mode pair loader, and the builder that forms
`A(R)` from `.chk` and `.mmn` (five sources, a provenance-checked cache, an
in-house Wigner-Seitz mapping). The current path does not need any of it, because
Wannier90 writes `_r.dat` directly. It stays for the `.mmn` route, the validation
notebooks, and `ws_pair_consistency` (used by `wannier/check_outputs.py`).
Import as `from wannier.wannier_io import ...` and
`from wannier.position_matrix import ...`; there are no re-exports left in `modules/`.

## Removed (2026-09-20)

The earlier pipeline (`modules/berry_curvature.py`, `modules/axion_angle.py`, the
old Y2Ir2O7 driver `axion_response_vs_nk.py`) is gone. The validation libs,
notebooks and unit tests were ported to `curvature.py` / `axion.py`; see
[`log/2026-09-20_legacy_removal.md`](log/2026-09-20_legacy_removal.md) for the old-versus-new
regression and for which notebooks could be rerun.

**The frozen MnBi2Te4 driver is gone too (2026-09-24).** `calculations/MnBi2Te4/phonon/axion_response_nk_sweep.py`
and its copy of the old module (`frozen_legacy/`) were deleted, together with their outputs
(`data/MnBi2Te4/phonon/T2-/mode*/Q_0p01/symmetrized/nk_sweep_shared_embedding*`). They existed only to
reproduce the paper's `d_Q theta ~ 10 rad/A`, which
[`log/2026-08-20_handoff_wannier_conventions.md`](log/2026-08-20_handoff_wannier_conventions.md) (3a) already flagged as wrong. **The
paper's numbers are in flux:** every MnBi2Te4 and Y2Ir2O7 value will be recomputed on the new pipeline, so
keeping code to reproduce an old number is no longer a reason to keep it. A copy (code and outputs) is in
`tmp/calculations_archive_2026-09-24/deleted_frozen_mbt_driver/` (git-ignored).

**Still on the old TensorFlow code path (2026-09-24 reorganisation).** The per-run copies of these notebooks
were merged into one notebook per task in `calculations/analysis/`, driven by run names. Where the task only
reads results (`c2_density`, `dtheta_vs_nk`, `bands`, `pdos`) it no longer computes anything on the old path.
Two still do, and their numbers are exploratory until ported to `modules/`:

| file | what it is for |
|---|---|
| `calculations/analysis/theta.ipynb` | θ of one structure (PythTB `Wannier`, 3-form, TensorFlow); warns that its SCDM gauge is near-singular |
| `calculations/analysis/z2_index.ipynb` | Z2 from Wilson-loop spectra (PythTB); on `mbt_base` and `mbt_base_sym` it gives ν = 0, whereas the saved output of the old notebook (a model at an old, unidentifiable path) showed ν = 1 |

`calculations/legacy/` keeps the old TensorFlow dθ drivers (`mbt_axion_response.py`/`.ipynb`,
`mbt_c2_density.py`, `yio_axion_response.ipynb`). `run_axion.py` does their job (dθ and the c2 density), but
the MnBi2Te4 runs have no `_r.dat` or `.mmn` cache, so it cannot run on them until they are regenerated. Delete
`legacy/` once that works. The old `nk_sweep/` and `c2_density/` folders in `data/MnBi2Te4` and the
`superseded` Y2Ir2O7 runs came from this path.
