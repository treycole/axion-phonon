# Wannier90 PR #702 writes our `A(R)` exactly - at the wrong Wigner-Seitz representative

Status 2026-09-10.  Settles the `transl_inv_full` control left open in
`validation/inputs/validate_transl_inv_full.py`, and revisits section 0 of
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md).

Revised the same day: the first version of this note blamed the `6F12.6`
print format alongside `use_ws_distance`.  The symmetry measurement below
shows the format contributes ~1% of the residual.  It is `use_ws_distance`,
and only `use_ws_distance`.

**Bottom line:** PR #702's `_r.dat` *is* postw90's `A(R)` at
`use_ws_distance=F` - agreement is 4.4e-16, machine epsilon, not a tolerance.
But that matrix is the wrong representative for us, and the symmetry test says
so: PR #702 improves the old `_r.dat` by 15x on the base P+T-forced zero and is
still **318,000x worse than the `.mmn` route**.  **The print format is not the
problem** - printing the same matrix at `ES24.16E3` changes the symmetry
residual by ~1%.  `_r.dat` remains unusable for this work, and the reason is
`use_ws_distance` alone.

## What the PR changes

[wannier-developers/wannier90#702](https://github.com/wannier-developers/wannier90/pull/702)
(`sjhong6230`, open, 3 files) teaches `wannier90.x`'s `plot_write_rmn` the
translation-equivariant formula, gated on `transl_inv_full=T`, which
`wannier90.x` now accepts as its own keyword rather than ignoring it.

Where postw90's `get_AA_R` applies `exp(i b.(r_n + r_m)/2)` before the Fourier
transform and `exp(-i b.R/2)` after it, PR #702 folds both into a single
`exp(i b.(r_n + r_m - R)/2)` inside the k-sum, using each k-point's own `b`
instead of postw90's `nninv` canonical neighbour ordering.  The `R=0` diagonal
is overwritten with the Wannier centres, as `get_AA_R` also does.  The two
routes are therefore *algebraically* the same and, as measured below,
numerically the same.

## The test

Bit-exactness after rounding, not a tolerance: if PR #702's `6F12.6` output is
anything other than the rendering of postw90's `ES24.16E3` value, the two
formulas differ.  `validation/inputs/compare_pr702_rdat.py`.

| system | WF | disent. | SOC | nrpts | components | bit-exact |
|---|---|---|---|---|---|---|
| diamond (`testw90_example05`) | 4 | no | no | 93 | 8,928 | **100.000000%** |
| `testw90_basic2` | 4 | yes | no | 89 | 8,544 | **100.000000%** |
| SrTiO3 Ti_Q_2A trial-03 | 48 | yes | yes | 729 | 10,077,696 | **100.000000%** |

Max residual after rounding on the production case: `4.4e-16`, i.e. zero.
Max raw difference `5.0e-07`, which is exactly half the `6F12.6` step.

Rebuilding `plot_write_rmn` at `ES24.16E3` removes the format from the question
entirely: PR #702 then agrees with postw90's `use_ws_distance=F` `A(R)` to
**4.4e-16** on diamond, machine epsilon, differing only by summation order.  The
two are the same computation, not two formulas that happen to agree.

At the file level against the `use_ws_distance=T` export (production SrTiO3
base, full precision, matched on `(R,n,m)` because the row order differs), 8.1%
of blocks disagree, median relative difference `8.1e-11` but max relative
difference exactly **0.500** - the factor-of-two signature of a block postw90
split across two Wigner-Seitz representatives at weight one half each.

Negative control - the *stock* `_r.dat` from the same checkpoint, against the
same reference:

| | stock `_r.dat` | PR #702 `_r.dat` |
|---|---|---|
| bit-exact vs postw90 | 74.3% | **100.0%** |
| max abs difference | 2.87e-02 | 5.00e-07 |
| mean rel. difference, `|A|>1e-5` | 0.419 | - |
| max rel. difference, `|A|>1e-5` | 29.1 | 4.8e-02 |
| Hermiticity residual | 1.12e-02 | 0.0 |

The `1.12e-02` reproduces the `7.2e-03` recorded for the other SrTiO3 run, so
this is the same defect measured before, and PR #702 removes it.

## The symmetry test -- the decisive measurement

File-level agreement cannot say whether the matrix is *usable*, so the three
`_r.dat` variants were pushed through the same in-house pipelines as the two
good routes.  `validation/symmetry/check_pr702_rdat_symmetry.py`; set (c),
48 WFs / 38 occupied; max `|Tr Omega_xy|` on symmetry-forced zeros, A^2.

| pipeline | base generic (P+T) | base kx=1/2 mirror | mode kx=1/2 mirror |
|---|---|---|---|
| `mmn_pair` | 8.570e-09 | 4.419e-09 | 8.072e-08 |
| `rfull_pair` (postw90, ws=T) | 8.720e-09 | 5.037e-09 | 8.158e-08 |
| `rdat_pair` [stock] | 3.971e-02 | 1.651e-02 | 1.653e-02 |
| **`rdat_pair` [PR #702, 6F12.6]** | **2.727e-03** | **7.050e-03** | **7.168e-03** |
| **`rdat_pair` [PR #702, ES24.16E3]** | **2.763e-03** | **7.087e-03** | **7.092e-03** |
| `rdat_direct` [stock] | 4.057e-02 | 1.911e-02 | 1.922e-02 |
| `rdat_direct` [PR #702, 6F12.6] | 4.243e-03 | 9.685e-03 | 9.865e-03 |
| `rdat_direct` [PR #702, ES24.16E3] | 4.395e-03 | 9.922e-03 | 9.904e-03 |

Read three things off it.

1. **PR #702 is a real improvement over the old `_r.dat`**: 14.6x on the base
   P+T zero, 2.3x on both mirrors.  The missing link phase was genuinely costing
   us that much.
2. **It is nowhere near good enough.**  Still 318,000x the `.mmn` residual on
   the base generic path, 1,600,000x on the base mirror.  A symmetry-forced zero
   coming out at 2.8e-03 is not a zero.
3. **The print format is a red herring.**  Rebuilding PR #702's writer at
   `ES24.16E3` and rerunning changes the residual by 0.5-1.3% - and in two of
   three columns makes it very slightly *worse*.  The `6F12.6` rounding is
   simply not what breaks the symmetry.

That last point corrects the first version of this note, which called the format
"a one-line format change upstream" and implied fixing it would matter.  It
would not.  The whole residual is the `use_ws_distance` representative.

The 2.763e-03 is also an independent confirmation of the mechanism: it lands on
the 2.8e-03 already recorded in
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md) for
"export at `use_ws_distance=F` and remap in house", measured there by a
completely different route.  Same number, same cause.

## What PR #702 does *not* fix

Of the three reasons `_r.dat` failed, only the second is addressed.

1. **`6F12.6` is unchanged** - but measured, it does not matter.  On the
   production matrix **66.4%** of postw90's values round to exactly zero, and
   among the 6.2% of components above `1e-5` the format costs up to **4.8%**
   relative error, so the file certainly looks damaged.  It is not what breaks
   the physics: the symmetry table above puts the format's share of the residual
   at ~1%.  The destroyed entries are overwhelmingly numerical noise on
   symmetry-zero elements, which is why killing them is nearly free.
2. **No link phase.**  Fixed.
3. **`use_ws_distance` still never enters `plot_write_rmn`** - and this is the
   one that matters most, because for `transl_inv_full` the Wigner-Seitz remap
   *cannot be deferred*.  `get_AA_R`'s second phase `exp(-i b.R/2)` is applied
   in R-space **after** `operator_wigner_setup`, so it depends on which
   representative `R` the remap chose.  A file written with `use_ws_distance=F`
   has that phase frozen at the wrong representative, and remapping afterwards
   cannot undo it.  We measured the cost of exactly this before: only `2.5e-04`
   in `A(R)`, but the base P+T-forced zero moves `8.7e-09 -> 2.8e-03`, i.e.
   entirely symmetry-breaking.

   PR #702's `_r.dat` *is* the `use_ws_distance=F` export.  Measured on
   `testw90_basic2`, postw90's `A(R)` at `use_ws_distance=T` differs from the
   `F` version by up to `2.7e-03` absolute and 100% relative on 536 of 1424
   rows, and emits the R rows in a different order; PR #702 is 85.4% bit-exact
   against the `T` file and 100% against the `F` one.

   **So even at infinite print precision, PR #702's `_r.dat` is still the
   wrong matrix for our symmetry work** - now measured, not inferred: the
   `ES24.16E3` rebuild scores 2.763e-03 on the base P+T zero.  The stored `SrTiO3_r_full.dat` under
   `base/.../185116bc889` was written with `use_ws_distance=T` and is not
   comparable to a PR #702 `_r.dat` row-for-row.

So the operational conclusion of the earlier note stands unchanged: keep the
`.mmn` + `_AA_cache` route, and keep `write_aa_r` for the `_r_full.dat`
cross-check.  Defect (3) is the whole story and it is structural -
`plot_write_rmn` runs before any Wigner-Seitz mapping exists, so no amount of
precision or downstream remapping can recover the frozen phase.  What PR #702
buys us is that `transl_inv_full` finally means the same thing in both binaries,
which makes it a genuinely useful independent check on `write_aa_r` - that is
how it is used here, and it is not a replacement for it.

## Caveat worth watching

`get_AA_R` phases with `wannier_centres_from_AA_R` (rebuilt from the `.mmn`);
PR #702 phases with `wannier_data%centres` (from the `.chk`).  `get_AA_R` only
checks that these agree as `sum((...)**2) < 1e-8`, which for 48 WF permits a
per-component drift of ~`8e-6` - *above* the `1e-6` print floor.  Measured here
it is `9.6e-09`, far below, so it never surfaced.  If `_r.dat` ever gains
precision, this becomes a real difference.

## Reproducing

Branches in the `wannier90` checkout, both on current `develop` (`e7132e80`):

- `feature/write-aa-r` - our `write_aa_r` / `_r_full.dat` export (`5c3bd255`)
- `compare/transl-inv-full` - the above plus PR #702 cherry-picked (`93535be6`)

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWANNIER90_MPI=OFF \
  -DBLAS_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib \
  -DLAPACK_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib
cmake --build build -j8
```

The production case reuses the converged checkpoint instead of re-wannierising:
copy `SrTiO3.win`, `.eig`, `.chk` and the unpacked `.mmn` into a scratch dir,
set `restart = plot`, `write_rmn = true`, `use_ws_distance = false`, turn off
`bands_plot`/`write_hr`/`write_xyz`, then

1. run `wannier90.x` as-is -> reproduces the archived `SrTiO3_r.dat` *numerically
   exactly* (only signed-zero text differs), which validates `restart = plot`;
2. add `transl_inv_full = true`, run `wannier90.x` -> the PR #702 `_r.dat` (1m54s);
3. add `write_aa_r = true`, run `postw90.x` -> `_r_full.dat` (21s);
4. `python3 validation/inputs/compare_pr702_rdat.py <_r.dat> <_r_full.dat>`.

Note that `use_ws_distance` must not be *appended* to these `.win` files - all of
them already set it, and Wannier90 rejects the duplicate as an unrecognised
keyword.
