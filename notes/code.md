# Code map

Since 2026-10-30 (restructure, `notes/restructure/design.md`) the production code is one package, `axion_phonon/`,
in the order the calculation runs, and one script per cluster step, `pipeline/`. The physics files never open a file:
`axion_phonon/io/` reads QE and Wannier90 files into arrays, and the physics takes and returns arrays. The old
`modules/`, `wannier/`, `run_axion.py` and the dated diagnostics are in `archive/` (its README; `axion_phonon`
reproduces them, `validation/restructure/`).

```text
axion_phonon/
  io/qe.py           wfcN.dat, save XML, UPF, pw.x outputs (eigenvalues, bands, Fermi level, k-path), projwfc PDOS
  io/wannier90.py    .chk, .eig, .win, .nnkp, .amn; read_run(folder, prefix) -> one dict of arrays per Wannier run
  io/results.py      save_dtheta_results (uniform-mesh dθ tables)
  paw.py             step 4: PAW metric S, overlaps O(k) = <psi_structure|S|psi_base> between two saves, A = O V_base
  wannier.py         step 6a: H(R), A(R) (transl_inv_full from the .chk), WS fold, pinning -> PythTB models
  curvature.py       step 6b: Omega and its three parts                                note §3-§7
  axion.py           step 6c: beta terms, c2, BZ integral (uniform mesh)                note §8-§9
  tube.py            step 6d: the MBT small-gap tube k list and its quadrature
  chern_simons.py    check: theta_CS in a smooth projection gauge; [theta_CS(Q) - theta_CS(0)]/Q vs c2
pipeline/            1_scf.sh ... 7_dtheta.sh (README there) and the Python the jobs run
```

## A dθ/dQ by hand

```python
from axion_phonon import axion, curvature, tube, wannier
from axion_phonon.io import wannier90

base = wannier90.read_run("data/.../base_wannier_run", "MnBi2Te4")      # V, energies, M(k,b), k, lattice, atoms
mode = wannier90.read_run("data/.../mode_wannier_run", "MnBi2Te4")
mb, mm = wannier.pair_models(base, mode)          # shared positions = base centres' nearest atoms = tau (A3)
pb, pm = curvature.position_terms(mb), curvature.position_terms(mm)
k = tube.klist(16, 4, 3, 8)["k"]                  # or a uniform mesh
fb, fm = curvature.fields(mb, pb, k), curvature.fields(mm, pm, k)
beta = axion.beta_terms(fb, fm, 0.01)             # Q = 0.01 A
c2_base = axion.c2_density(curvature.omega(curvature.connection(fb, 82, beta)))
```

`pipeline/c2.py` is exactly this on a slice of the k list; `pipeline/dtheta.py` integrates. For one run alone (bands,
parities): `wannier.run_model(wannier90.read_run(folder, prefix, overlaps=False), pin=True)`.

## `wannier.py` (seven functions; replaces most of the 3,000 lines of the old `wannier/`)

| function | does |
|---|---|
| `class_representatives(mp_grid, lattice)` | one R per aliasing class R mod N: the shortest Cartesian member (CLAUDE.md pitfall 1) |
| `hamiltonian_R(V, E, k, R)` | H(R) = (1/N_k) Σ_k e^{-ik·R} V† E V, full precision from the .chk + .eig, no hopping cut (bug #22) |
| `fold(matrices, lattice, positions, mp_grid)` | each (i, j) entry to its nearest image of R + τ_j − τ_i, ties (≤ 5e-3 Å) split equally |
| `position_matrix_R(M, ...)` | A(R) by Wannier90's transl_inv_full formula per b, folded on the same positions |
| `pin_to_atoms(centres, atoms, lattice)` | each centre → its nearest atom (bug #21) |
| `tb_model(lattice, tau, H_R, A_R)` | the PythTB model (hoppings, `_pos_r`); checks R ↔ −R closure and A(R) = A(−R)† |
| `pair_models(base, mode)` | both models on the base's pinned positions, one τ for both |

## `curvature.py` (seven functions)

Direction index first: `dH`, `A` are `(3, nk, J, J)`, `Omega` is `(D, D, nk, Nv, Nv)` with `D = 3` or `4` (with β).

| function | returns | note |
|---|---|---|
| `position_terms(model)` | `X(R)=<0s\|r-tau_s\|Rt>`, bond vectors, curl weights; once per structure (complex64 by default) | (8), (9) |
| `fields(model, pos, k)` | Wannier gauge: `H`, `dH/dk`, `A(k)`, curl of `A` | (7)-(9) |
| `connection(f, n_occ, beta=None)` | band gauge: internal `d_cv`; external `a_cv`, `A_vv`, `Omega_vv` | (10), (28), (31) |
| `omega_internal(c)` | `i[d_mu^† d_nu - d_nu^† d_mu]` | (29) |
| `omega_cross(c)` | `(d_mu^† a_nu + h.c.) - (mu<->nu)` | (30) |
| `omega_external(c)` | `Omega_vv - i[A_vv_mu, A_vv_nu]` | (27) |
| `omega(c)` | internal + cross + external | (24) |

`beta` (a `BetaTerms`) adds the fourth direction. Without `A_beta` the run is in the **frozen Wannier gauge**,
`A_beta = 0` (exact at O(Q) in the projection gauge with PAW-corrected `.amn`: `Λ = 0`, bug #23).

## `axion.py`

| function | does | note |
|---|---|---|
| `beta_terms(f_base, f_mode, dbeta)` | `dH/dbeta`, `dA/dbeta`, mixed curl (one secant for both ends) | 8.3, (35) |
| `c2_density(Omega)` | `(1/16 pi) eps tr[B B]` | (38) |
| `integrate_c2(c2)`, `dtheta(...)` | uniform-mesh integral; the one-call driver for a pair | §8-§9 |

The two ends (base, structure) share one secant and differ only in where the spatial fields are evaluated; their
spread is not a validity check.

## Verification

* `validation/unit/` (pytest, no data): curvature covariance and (38) = (39) (`test_curvature.py`); tube quadrature,
  fold, b-shell weights, `.amn` layout (`test_axion_phonon_core.py`); PAW pieces (`test_axion_phonon_paw.py`).
* `validation/restructure/`: the new package against the code that made the current numbers: PAW overlaps 0.0; H(k)
  6e-16, A(k) identical, c2 4e-14 (complex128); production precision vs the cluster's c2 1.7e-6; Chern-Simons θ;
  the pipeline's job drivers and its dry run.
* The other `validation/` notebooks (external codes, symmetry, Born charge) were written against the old code and run
  against `archive/code` (`validation/_paths.py`); they validate the same curvature code.
