# `rdat_ndegen_applied` closes on the production SrTiO3 system

Status 2026-09-16. Resolves the "What could not be tested" gap left by
[2026-09-13_jaemo_transl_inv_full_ws_distance.md](2026-09-13_jaemo_transl_inv_full_ws_distance.md):
that note validated Jae-Mo Lihm's `write_ndegen_applied` fix on two small
Wannier90 test-suite systems (diamond, `testw90_basic2`) because the
production `.mmn` had already been deleted after caching (see
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md) §0b). This
note reruns the decisive test — the one that actually mattered, since PR #702
was *also* bit-exact at the file level and still 318,000x off on it — on the
real production system.

**Bottom line: `rdat_ndegen_applied` is symmetry-exact on SrTiO3 Ti_Q_2A
trial04 (48 WF / 38 occupied), and cross-validated against an independently
built postw90 `rfull` export to ~1e-8.** It is implemented as
`source="rdat_ndegen_applied"` in `wannier/position_matrix.py`, needs no
`.mmn`, and is safe to use standalone — **provided the `.win` sets all three
flags together**: `use_ws_distance=T`, `transl_inv_full=T`,
`write_ndegen_applied=T`. Any one missing reproduces an earlier, worse
failure mode (see `rdat-position-matrix-unusable` / `jaemo-wannier90-fork`
project memory for the intermediate `transl_inv_full`-missing result, which
broke C4z/Mx covariance to ~4.3e-2 relative).

## What was run

Ti_Q_2A trial04 (`.../trial_04_Sr_sp_Ti_p_O_sp/proj_gauge/191583bc889`), all
three flags set together, confirmed via pythtb's own readers
(`num_wann=38`, all H(R) degeneracies already 1, `n_R=729` for both `ham_r`
and `pos_r`, `_wsvec_write_ndegen_applied()` reads `True` from the header).
`validation/symmetry/check_rdat_ndegen_applied_symmetry.ipynb` run against
it, borrowing trial03's Jae-Mo `.mmn` as the `mmn` reference (same NSCF,
`check_mmn_compatibility` passes; H(R)/gauge for both sides still comes from
trial04's own `.chk`, so `max_abs_energy_difference_eV = 0.0` exactly — a
clean A(R)-only comparison).

## Symmetry result

C4z and Mx covariance of Tr Ω on 24 generic off-grid k-points:

| route | C4z (relative) | Mx (relative) |
|---|---|---|
| `mmn_pair` | 5.9e-10 | 5.6e-10 |
| `rdat_ndegen_applied` | 5.2e-10 | 5.1e-10 |

Matches `mmn_pair` to within noise. A(R) Hermiticity
(`ndegen_applied_hermiticity_relative`): **1.06e-12** (was 7.2e-3 without
`transl_inv_full` — see the intermediate result in project memory). This is
7-8 orders of magnitude better than the `transl_inv_full`-missing
intermediate result, with `transl_inv_full` the only change between them.

## Cross-validation against an independently built `rfull`

`_r_full.dat` for the identical trial04 run, generated locally (no cluster
needed): our own `~/Repos/wannier90` fork (`compare/transl-inv-full`,
unrelated codebase to Jae-Mo's), `postw90.x` with `write_aa_r=true`,
symlinking trial04's `.chk`/`.eig`/`.nnkp`/`_hr.dat`/`_centres.xyz`/
`_wsvec.dat` plus trial03's `.mmn` (borrowed, same NSCF). 19s on one core.

All three share trial04's own `.chk`, so H(R)/gauge is identical across all
three (confirmed: `mmn` vs `rdat_ndegen_applied` gives
`max_abs_energy_difference_eV = 0.0` exactly):

| comparison | yz/zx relative_l2 | xy relative_l2 | max|energy diff| |
|---|---|---|---|
| `rdat_ndegen_applied` vs `mmn` | 5.8e-02 | 3.80e-01 | 0.0 eV |
| `rfull` (local postw90) vs `mmn` | 5.8e-02 | 3.80e-01 | 6.47e-05 eV |
| **`rdat_ndegen_applied` vs `rfull`** | **1.0e-08** | **2.1e-08** | 6.47e-05 eV |

Two completely independent Fortran implementations — Jae-Mo's
`wannier90.x`-side `write_ndegen_applied` + `transl_inv_full`, versus our own
`postw90.x`-side `write_aa_r` + `transl_inv_full` — converge on the identical
answer, while both differ from `mmn_pair` by the identical ~5.8%/38%. This
confirms that gap is the
[rfull-vs-mmn family difference](2026-08-27_rfull_vs_mmn_position_matrix.md), not a
defect in either implementation: `rdat_ndegen_applied` is now in the `rfull`
family, not the `mmn` family, and the two families disagree by this much at
finite mesh density regardless of code correctness (§1 of
`2026-08-20_handoff_wannier_conventions.md` — this is the aliasing-class-as-interpolant
point, not a bug). `rfull` checked independently against the same symmetry
test: C4z=3.3e-10, Mx=3.3e-10 — also clean. The 6.47e-5 eV energy differences
involving `rfull` are `_build_rfull` redundantly reapplying the Python
orbital-aware R mapping to an already-correctly-mapped H(R) — floating-point
noise, not physics.

## Conclusion

The "not rely on `.mmn`" goal from
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md) is met.
`rdat_ndegen_applied` is a validated, symmetry-exact, `.mmn`-free `A(R)`
source — but it is a *different interpolant* from `.mmn` (the `rfull`
family), not a bit-identical replacement, so don't expect its pointwise
curvature to match `.mmn`'s beyond the known family gap. `.mmn` remains the
production default (closer to DFPT per
`2026-08-27_rfull_vs_mmn_position_matrix.md` point 3); `rdat_ndegen_applied` is the
right choice when disk or reproducibility argues against keeping `.mmn`
around. No further action needed on this question.

## Reproducing

```bash
python3 validation/symmetry/check_rdat_ndegen_applied_symmetry.ipynb  # run in Jupyter, axion kernel
```

See `wannier/position_matrix.py`'s `POSITION_SOURCES["rdat_ndegen_applied"]`
for the builder and its guard clause (refuses to run unless
`_wsvec.dat`'s header confirms `write_ndegen_applied=.true.`, rather than
silently returning the old wrong `rdat_direct` numbers).
