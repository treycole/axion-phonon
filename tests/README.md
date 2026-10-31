# tests

Checks that the code is right and that the production results satisfy what the paper claims about them. Run from the
repository root with the `axion` env:

    python -m pytest tests                 # everything (~10 min, ~8 GB); data tests skip, naming the missing path,
                                           # where a run is not in data/
    python -m pytest tests -m "not data"   # unit tests only (seconds, no data)
    python -m pytest tests -s              # also print the measured residuals

`conftest.py` builds the production models the data tests share, once per session: the MnBi₂Te₄ base/pattern-3 pair
(window `f10.962_w12.810_k8`, ~2 min) and the Y₂Ir₂O₇ base (run 214249, ~2 min). Their paths are its settings block.

| test | what it checks | data | ported from (`archive/validation/`) |
|---|---|---|---|
| `unit/` | `axion_phonon` on analytic models and fixed-seed fakes: curvature pieces against exact moving frames, (38) = (39), PAW overlaps, Wannier fold, tube k list, `tube.integrate_at`; the symmetry pass/fail gate catches each kind of defect | none | `unit/` |
| `inputs/test_interpolation.py` | H(R) reproduces the `.chk` Hamiltonian on the mesh (eigenvalues to 1e-8 eV); off the mesh H, A, curl A are Hermitian and covariant under k → k + G; analytic dH/dk and curl A match central differences | MBT, YIO bases | `inputs/check_interpolation_inputs` |
| `curvature/test_vs_pythtb.py` | our internal and cross + external curvature equal PythTB's own (`include_external` True − False), full occupied matrices, to 1e-8; the mixed k–β block with a nonzero A_β equals a hand assembly of every term, and is Hermitian | MBT, YIO bases; MBT pair | `external_curvature/check_vs_pythtb` |
| `curvature/test_wilson_flux.py` | ∫ Tr Ω over mesh plaquettes (Gauss–Legendre) against the occupied Wilson-loop phase from the `.chk` overlaps: agree to 5% rms (measured 2.5–2.7%, the finite-link error of 1/8-zone links) | YIO base (MBT's trace is zero by PT) | `external_curvature/check_vs_wilson_flux` |
| `symmetry/test_mbt_symmetry.py` | MBT models obey their magnetic point groups, read from each structure's `scf.out` (QE's matrices and time-reversal flags): base 24 operations (QE's 12 plus their products with the anti-translation {1′\|½,½,½}, which QE does not list for the supercell), each displaced structure 12: DFT eigenvalues, interpolated bands, T_ab = Re Tr[F_a F_b] by part; controls (C4z, broken inversion) fail. Thresholds: `symmetry_lib.THRESHOLDS` | session pair (p3); other patterns from saved cluster results | `2026-10-22_mbt_mode_symmetry` (archive/diagnostics) |
| `symmetry/test_yio_symmetry.py` | YIO base obeys m-3m′ (47 operations): DFT eigenvalues, interpolated bands, traced curvature by part; a fourfold operation with flipped priming fails | YIO base | `symmetry/check_magnetic_symmetry`, `2026-10-26_yio_symmetry_conv1e-10` |
| `mbt/test_dtheta.py` | MBT production dθ/dQ: window edges in gaps with the base's band counts; dense k list reproduces the production values, which agree to 0.2%; integrand real, gap > 0.15 eV; pattern-3 Q → 0 limit within 0.1% | `windows/f10.962_w12.810_k8/` | (new) |

## The MnBi₂Te₄ symmetry table (paper `tab:symmetry-mbt`)

`symmetry/symmetry_lib.py` does the tests; `symmetry/run_mbt_symmetry.py` runs them for one base/structure pair and
writes the table (its log), `symmetry_{base,mode}.npz`, and PASS/FAIL against the thresholds. Locally: set `BASE_RUN`
and `MODE_RUN` in its settings block and run it with the repository on `PYTHONPATH`. On the cluster (the `.chk` of
patterns 1, 2, 4 are only there): copy `run_mbt_symmetry.{sh,py}` and `symmetry_lib.py` to one folder on the login
node, set `STRUCTURE` in the `.sh` and run it (`SUBMIT=false` first: a dry run); the job writes
`windows/<window>/symmetry_<structure>/`. Copy that folder to `data/` with the job folder renamed `out/`, and
`test_mbt_symmetry.py` checks it.

## Not here yet

- Y₂Ir₂O₇ displaced structures (−4′3m′): their `.chk` are on the cluster only.
- Ports that need new runs: the curvature against postw90 (`kpath_task = curv`) and WannierBerri on a production
  `.chk`, and the SrTiO₃ Born effective charge (the only check of the mixed k–u curvature against a measured number),
  which needs SrTiO₃ base + Ti displacement through `pipeline/` steps 1–6.
- What was not ported, and why: `archive/validation/` compared A(R) sources (`.mmn`, `_r.dat`, PR #702, Jae-Mo's
  fork, postw90 `rfull`) and SrTiO₃ trial sets; production uses one source (CLAUDE.md pitfall 2), so those are a record.
