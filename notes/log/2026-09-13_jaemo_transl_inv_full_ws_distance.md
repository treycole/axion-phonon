# Jae-Mo's `write_ndegen_applied` fixes `use_ws_distance` under `transl_inv_full`

Status 2026-09-13. Extends
[2026-09-10_pr702_transl_inv_full_rdat.md](2026-09-10_pr702_transl_inv_full_rdat.md),
which left one thing unresolved: PR #702 gives `wannier90.x` the
translation-equivariant `_r.dat` formula, but `use_ws_distance` never enters
`plot_write_rmn`, so the file is permanently frozen at the wrong
Wigner-Seitz representative and is unusable for symmetry work.

**Bottom line:** Jae-Mo Lihm's `plan5-write-ndegen-applied` branch
(`jaemolihm/wannier90`, HEAD `88134317`) fixes exactly that. With
`transl_inv_full=T, use_ws_distance=T, write_ndegen_applied=T`, its `_r.dat`
matches postw90's `write_aa_r` reference at `use_ws_distance=T` **bit-exact
after 6-decimal rounding, 100.000000% of components, on both test systems
checked**, with Hermiticity exactly `0.000e+00` and every row's `-R` partner
present. PR #702's `_r.dat` is confirmed **byte-for-byte identical** whether
`use_ws_distance` is `T` or `F` in the `.win` -- direct proof that the flag
is silently ignored, not just numerically inert.

**What this does not settle:** the production SrTiO3 symmetry-forced-zero
test (`validation/symmetry/check_pr702_rdat_symmetry.py`, set (c)) could not
be rerun. See "What could not be tested" below.

## What changed vs. PR #702

Jae-Mo's branch is PR #702's own branch (`rmn_file`, `sjhong6230`) plus his
own commits on top: `e9434700` "share the Wigner-Seitz expansion between
wannier90.x and postw90.x", `692ea6e0` "one owner for `<0m|r|Rn>`, shared by
`_r.dat` and `_tb.dat`", `6f9d3aa9` "feat: write_ndegen_applied, self-contained
real-space output files", plus review/test follow-ups (`008f5719`, `305a2484`,
`88134317`).

Practically: requesting `use_ws_distance=T` under `transl_inv_full=T` without
`write_ndegen_applied=T` is now a hard error --

```
Exiting.......
transl_inv_full=T with use_ws_distance=T needs write_ndegen_applied=T:
_r.dat/_tb.dat cannot hold <0m|r|Rn> on the folded R grid
```

-- rather than the silent wrong-representative file PR #702 produces. Setting
`write_ndegen_applied=T` changes the file format: a Wigner-Seitz-degenerate
bond can now appear at more than one `R` image, each row *already* divided by
its `ndegen` weight (self-contained, unlike the classic `_hr.dat` +
`_wsvec.dat` convention where the reader must divide separately). Measured on
diamond, this takes `nrpts` from 93 (`use_ws_distance=F`) to 141; on
`testw90_basic2` it stays at 89 (this system apparently has no WS-boundary
bonds that need splitting). postw90's own `write_aa_r` reference reorders --
and, on diamond, also lengthens -- its row list the same way at
`use_ws_distance=T`, so the two files must be matched by `(R, n, m)` key, not
by row position. `validation/inputs/compare_jaemo_ndegen_rdat.py` does this.

## The test

Two Wannier90 test-suite fixtures that still ship their `.mmn` (the production
system's was deliberately deleted -- see below): `testw90_example05`
(diamond, seedname `diamond`, 4 WF, no disentanglement) and `testw90_basic2`
(seedname `wannier`, 4 WF, disentangled). Each was wannierised fresh (not
`restart=plot`, to avoid any question of cross-branch `.chk` compatibility) by
three binaries built from three branches, all release 4.0.3:

- `wannier90.x` from `/Users/treycole/Repos/wannier90` `compare/transl-inv-full`
  (`feature/write-aa-r` + PR #702 cherry-picked) -- the PR #702 candidate and,
  via its `postw90.x`, the `write_aa_r` reference.
- `wannier90.x` from `/Users/treycole/Repos/wannier90-JaeMo`
  `plan5-write-ndegen-applied` -- the Jae-Mo candidate.

The reference for the Jae-Mo leg needed a `.chk` gauge-matched to Jae-Mo's own
wannierisation. Tried the simple route first: feed Jae-Mo's `.chk` directly to
the `compare/transl-inv-full` build's `postw90.x` (`write_aa_r=T`). It loaded
without complaint on both systems -- the two branches' `.chk` format is
compatible at 4.0.3, so no patch-porting was needed.

### `use_ws_distance=F` -- regression check

Reproduces the known PR #702 result (now also for a from-scratch
wannierisation, and now also for Jae-Mo's branch on its `ws=F` path, which
shares the same formula):

| system | candidate | bit-exact vs 6F12.6(ref) | max abs diff | Hermiticity |
|---|---|---|---|---|
| diamond | PR #702 | 100.000000% | 4.83e-07 | 0.00e+00 |
| diamond | Jae-Mo | 100.000000% | 4.83e-07 | 0.00e+00 |
| basic2 | PR #702 | 100.000000% | 4.99e-07 | 0.00e+00 |
| basic2 | Jae-Mo | 100.000000% | 4.99e-07 | 0.00e+00 |

(`validation/results/jaemo_rdat_file/file_level_summary.md`, via
`compare_pr702_rdat.py`.)

### PR #702 ignores `use_ws_distance` -- confirmed directly

`diff`ing each system's PR #702 `_r.dat` written with `use_ws_distance=F`
against the one written with `use_ws_distance=T` (same `.chk`, only that one
`.win` line changed): **identical apart from the run timestamp in the header**,
on both diamond and basic2. Not "numerically indistinguishable" -- literally
the same bytes. This is the sharpest possible confirmation of the mechanism
identified in the PR #702 note: `use_ws_distance` never reaches
`plot_write_rmn`.

### `use_ws_distance=T` -- the new result

Jae-Mo `_r.dat` (`write_ndegen_applied=T`, folded R grid) vs. the
`write_aa_r` reference at `use_ws_distance=T`, matched by `(R, n, m)` key
(`compare_jaemo_ndegen_rdat.py`):

| system | candidate/reference rows | keys only on one side | bit-exact vs 6F12.6(ref) | max abs diff | Hermiticity (keyed) |
|---|---|---|---|---|---|
| diamond | 2256 / 2256 | 0 / 0 | 100.000000% | 4.83e-07 | 0.00e+00 (0 unmatched) |
| basic2 | 1424 / 1424 | 0 / 0 | 100.000000% | 4.99e-07 | 0.00e+00 (0 unmatched) |

Every row in Jae-Mo's file has a matching key in the reference and vice versa,
every value rounds to the same 6-decimal number, and the file is exactly
Hermitian under its own folded convention. This is the same standard PR #702
only meets at `ws=F` (`max abs diff` here, ~4.8e-7 to ~5.0e-7, is the same
half-a-6F12.6-step floor seen throughout this series -- i.e. this is a
print-format-limited match, not an approximate one).

## What could not be tested

**Resolved 2026-09-16** — see
[2026-09-16_jaemo_rdat_production_closure.md](2026-09-16_jaemo_rdat_production_closure.md).
The production `.mmn` question below turned out not to block a production-scale
rerun; `rdat_ndegen_applied` is symmetry-exact on Ti_Q_2A trial04 and
cross-validated against an independently built `rfull` export to ~1e-8.

The production SrTiO3 symmetry-forced-zero pipeline
(`check_pr702_rdat_symmetry.py`, set (c), 48 WF / 38 occupied) needs its
`.mmn` to run `wannier90.x` from scratch. That `.mmn` was deliberately deleted
after caching into `_AA_cache/*.npz` --
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md) section 0b,
a verified-bit-identical 40.5 GB cleanup -- and there is no local QE
wavefunction data to regenerate it (would need rerunning the NSCF +
`pw2wannier90` step, not attempted this session). So the "is it usable"
question that made PR #702's file-level bit-exactness insufficient on its own
(PR #702 was *also* bit-exact at `ws=F`, and still 318,000x off on the
symmetry-forced zero) is not yet re-answered for Jae-Mo's fix on the
production system. The file-level result above is necessary but the same
prior note shows it is not automatically sufficient; the small-system result
is a genuine but partial validation. If the production `.mmn` becomes
available again (regenerated, or restored from an archive), rerun
`check_pr702_rdat_symmetry.py --variant JaeMo=<base _r.dat>,<mode _r.dat>`
(after rebuilding both structures' `_r.dat` with `write_ndegen_applied=T`)
alongside the already-saved `PR702`/`stock` curves in
`validation/results/pr702_rdat_symmetry/symmetry_curves.npz`.

## Reproducing

Branches:

- `/Users/treycole/Repos/wannier90` `compare/transl-inv-full` (`feature/write-aa-r`
  + PR #702, unchanged from the prior note)
- `/Users/treycole/Repos/wannier90-JaeMo` `plan5-write-ndegen-applied`, tracking
  `origin/plan5-write-ndegen-applied` (`jaemolihm/wannier90`, HEAD `88134317`)

Both built with the same recipe:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWANNIER90_MPI=OFF \
  -DBLAS_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib \
  -DLAPACK_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib
cmake --build build -j8
```

Per system (`testw90_example05`/`diamond`, `testw90_basic2`/`wannier`), per
`use_ws_distance` setting: copy `.win`/`.amn`/`.eig`/`.mmn` into a scratch dir,
set `use_ws_distance`, add `write_rmn = .true.`, `transl_inv_full = .true.`,
and -- only for the Jae-Mo build at `use_ws_distance=T` -- `write_ndegen_applied
= .true.`; run `wannier90.x <seed>` (full wannierisation, not `restart=plot`).
For the reference, copy the resulting `.chk`/`.eig`/`.mmn` into another
scratch dir with a minimal `.win` (`write_aa_r = .true.`, matching
`use_ws_distance`/`transl_inv_full`, the same `kpoints` block as the original
`.win`) and run the `compare/transl-inv-full` build's `postw90.x <seed>`.

```
python3 validation/inputs/compare_pr702_rdat.py --case LABEL=CANDIDATE,REFERENCE ...   # ws=F
python3 validation/inputs/compare_jaemo_ndegen_rdat.py CANDIDATE REFERENCE              # ws=T
```

## Upstream: both fixes are now on `wannier-developers/wannier90` develop (2026-09-22)

**PR #702 is merged** (`gh pr view 702 --repo wannier-developers/wannier90`: `state=MERGED`,
`mergedAt=2026-09-16T10:48:57Z`, merged by maintainer Jerome Jackson, `rmn_file` -> `develop`) —
superseding this note's and `2026-09-10_pr702_transl_inv_full_rdat.md`'s "open" status.

**Jae-Mo's `write_ndegen_applied` work landed the same day**, on top of it: commit
[`e345784c`](https://github.com/wannier-developers/wannier90/commit/e345784cd543948e8d6cf65bd534314d2a2707c8)
("Merge pull request #2 from jaemolihm/plan5-write-ndegen-applied — Add write_ndegen_applied:
self-contained _hr.dat/_r.dat/_tb.dat, and make transl_inv_full work with use_ws_distance"),
authored/merged by Seung-Ju Hong (`sjhong6230`, same author as PR #702) into the mainline repo
directly (not `jaemolihm/wannier90`). Confirmed by ancestry, not just the commit message:
`git merge-base --is-ancestor e345784c origin/develop` returns true against current develop HEAD
`c4a924a5` (28 commits further). Confirmed in the actual merged source
(`/Users/treycole/Repos/wannier90`, `git show origin/develop:src/hamiltonian.F90`): the
`write_ndegen_applied` path threads `use_ws_distance` into `ws_apply_ndegen` and gates on
`transl_inv_full .and. write_ndegen_applied` together — the same three-flag combination this
project's recipe uses, not merely the keyword name present with different wiring.

**Practical consequence:** production likely no longer needs `jaemolihm/wannier90` at all — stock
`wannier-developers/wannier90` on `develop` (or whatever release/tag first includes `c4a924a5` or
later) should carry the full fix. **Not yet done:** rerunning this note's own decisive test (the
production SrTiO3 symmetry-forced-zero check, `validation/symmetry/check_pr702_rdat_symmetry.py`
set (c)) against a `develop`-built binary before actually switching production over — "merged
upstream" confirms the code is there, not that it behaves identically to the already-validated
Jae-Mo-fork build on our specific test. `wannier/run_plot_stage.sh`'s default `WANNIER90_X` path
and `wannier/README.md`'s attribution to "Jae-Mo's fork" are now stale if this switch happens and
should be updated together with it, not left pointing at a fork that's no longer necessary.
