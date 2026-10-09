# A bug in QE 7.5's noncollinear DFT+U symmetriser (2026-09-21)

**Where:** `PW/src/new_ns.f90`, subroutine `new_ns_nc`. **What:** it symmetrises the Hubbard occupation matrix
`ns` with the rotation applied in the wrong orientation relative to the atom it pairs it with. The error
vanishes for any operation that is its own inverse (E, inversion, 2-fold rotations, mirrors) and is present for
every operation of order 3, 4 or 6 that permutes Hubbard atoms. **Effect on our runs:** every Y2Ir2O7 SCF (base
and displaced) with `symmetry_with_labels = .true.` and four Ir species. This is the source of the 4 to 14%
symmetry floor of `notes/log/2026-09-21_yio_symmetry.md`, and it also shrinks the Ir Hubbard moment.

## How the operations reach `new_ns_nc`

`symmetry_with_labels = .true.` in a noncollinear run sets `colin_mag = 0` (`setup.f90`). `sgam_at` and
`checksym` (`symm_base.f90`) then match atoms by chemical symbol, `chem_symb(atm(ityp(nb))) ==
chem_symb(atm(ityp(na)))`, instead of by `ityp`. So `Ir1` ... `Ir4` count as equivalent, `sgam_at_mag` keeps the
operations that also preserve the moments `m_loc` (with or without time reversal), and `irt(isym, na)` can now
point at an atom of a different species. That is how QE finds the 48 operations of the AIAO pyrochlore.
`irt(isym, na)` is the image of `na` under `sr(isym)`; checked on all 48 operations of our run.

## The bug

Collinear `new_ns` contracts the orbital rotation on its first index (`D^T nr D`):

    ns(m1,m2,is,na) += d(m0,m1,isym) * nr(m0,m00,is2,nb) * d(m00,m2,isym) / nsym        nb = irt(isym,na)

`new_ns_nc` contracts it on the second (`D nr D^T`) and puts `conj(d_spin_ldau)` on the left, with the same `nb`:

    nr1(m1,m2,is1,is2,na) += conj(d_spin_ldau(is1,is3,isym)) * d2(m1,m3,isym) * nr(m3,m4,is3,is4,nb)
                             * d_spin_ldau(is2,is4,isym) * d2(m2,m4,isym) / nsym          nb = irt(isym,na)

`nr` is the transpose of the density matrix (`init_ns_nc` and the projection loop both build it that way).
For an invariant state the forward relation is `nr(irt(g,a)) = U_g nr(a) U_g^dag` with `U_g = D_S (x)
conj(dS_g)` (and `nr(a)^T` for a time-reversed `g`), so the term for atom `a` must be taken from the atom that
`g` maps *into* `a`. The collinear routine does that through its transposed contraction; `new_ns_nc` does not.

## Evidence

Everything below uses QE's own data for the run in `data/Y2Ir2O7/base/soc/u_3.0/output/178382bc889`
(48 operations, `irt`, `t_rev`, `find_u`, `comp_dspinldau`, `d_matrix`), ported line for line to Python in
[`calculations/diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py`](../../calculations/diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py).

* The operation data are consistent. `find_u(sr)` implements `sr` for all 48; the AIAO moments obey
  `m(irt(a)) = +-det(sr) sr m(a)` for all 48 (`-` for time-reversed); QE's own starting ns is invariant under the
  forward action of all 48 to 6e-15, as is a random orbital-spin state built by group averaging.
* Feed QE's analytic AIAO ns to `new_ns_nc` and every one of the 48 terms must reproduce it. **28 do not:** the
  8 C3, 8 S6, 6 C4 T and 6 S4 T. The other 20 (E, i, all 2-fold rotations and mirrors) are exact. Site-fixing
  operations are fine even when they have order 3; only the atoms that an operation permutes go wrong.
* Full orbital x spin test, deviation from the fixed point:

  | | as written | `nb = irt(invs(isym), na)` |
  |---|---|---|
  | random invariant state | 0.44 | 1.6e-15 |
  | AIAO ns, |m| = 3 | 0.16 | 2.2e-15 |
  | idempotence `P(P(X)) - P(X)` | 0.87 | 1.7e-15 |
  | output invariant under the 48 forward actions | 1.7 | 4.8e-15 |

* Applied to the exact AIAO ns the average returns |m| = 1.0 instead of 3.0 on every Ir: the Hubbard moment
  is cut by a factor 3 per application. (An earlier figure of "15x" in the conversation was wrong: it compared
  a per-orbital moment with a per-site total.)
* The symptom in the run matches: iteration 1 (analytic guess) is exactly symmetric, the first computed ns is
  already asymmetric (|m| spread 4.4%, tilts up to 1.9 deg) and stays so (9%, up to 5 deg), while the
  density-based site moments from `report = 1` are symmetric to 7 digits because the density symmetriser is
  correct. `HUBBARD (atomic)` gives the same pattern, so the projector is not involved.

## The fix

    -  USE symm_base, ONLY : nsym, irt, time_reversal, t_rev
    +  USE symm_base, ONLY : nsym, irt, time_reversal, t_rev, invs
    ...
    -  nb = irt (isym, na)
    +  nb = irt (invs(isym), na)

in `new_ns_nc` only ([`qe_patches/new_ns_nc_preimage.patch`](../../calculations/qe_patches/README.md)).
`invs` is QE's index of the inverse operation; for all 48 operations `irt(invs(g), .)` is the inverse
permutation of `irt(g, .)`.

**Verification in a real run.** QE 7.5 built from a copy of the tree (gfortran 14, Open MPI, OpenBLAS, FFTW3) and
run on the actual YIO base input (22 atoms, 48 operations, 8 k-points, `HUBBARD (ortho-atomic)`,
`symmetry_with_labels`), unpatched and patched, same Davidson settings and seeds, so both receive **identical raw
occupations** and only the symmetriser differs. Hubbard-projected |m| on Ir1..Ir4
([`qe_patches/verification/`](../../data/Y2Ir2O7/qe_patches/verification), compared with `calculations/qe_patches/verification/compare_hubbard_moments.py`):

| block | unpatched | spread | tilt off the 3-fold axis | patched | spread | tilt |
|---|---|---|---|---|---|---|
| analytic start | 3.000 x4 | 0 | 0 | 3.000 x4 | 0 | 0 |
| first computed ns | 0.781, 0.782, 0.787, 0.767 | 2.5% | up to 0.7 deg | 1.069 x4 | **0.00%** | **0.00 deg** |
| second | 0.187, 0.187, 0.185, 0.184 | 2.0% | up to 0.8 deg | 1.447 x4 | **0.00%** | **0.00 deg** |
| unpatched, block 5 | 0.194, 0.194, 0.190, 0.170 | 12.6% | up to 3.6 deg | (run stopped at 2 iterations) | | |

The patched occupations are identical on the four Ir to the printed precision at every iteration. The unpatched map
is asymmetric and also shrinks the moment (0.78 against 1.07 on the same raw data at the first step, then 0.08 to 0.19).

**Caveat on that build.** The scratch build's SCF is numerically unhealthy: iteration 1 gives a total energy of
-4960 Ry against -4625.9 Ry in the cluster run, CG aborts with "too many bands are not converged", and later
energies overflow (probably the OpenBLAS/gfortran combination). So only the *comparison* is claimed, not any
physical value from it. Because the symmetriser is a projector acting on whatever raw occupation it receives, and
both runs receive the same data, the comparison isolates the patch. The cluster build converges these runs.

## Cluster confirmation (2026-09-22)

**Where:** `data/Y2Ir2O7/base/soc/u_3.0/output/hub_sym_fix/192271bc889`, the patched `pw.x` (built
alongside `pw_unpatched.x` in the same directory, not yet run) on the actual production base-state input: 22
atoms, `symmetry_with_labels = .true.`, `nosym = .false.`, `HUBBARD (ortho-atomic)`,
`mpirun -np 48 pw.x -input Y2Ir2O7.scf.in`. `48 Sym. Ops., with inversion, found (36 have fractional
translation)` — the order-3/4/6 operations that expose the bug are active, unlike the local build's own SCF
which never got numerically healthy enough to trust.

**Result:** the four Ir (atoms 5-8) are exactly equivalent at every iteration: equal |m|, each purely along its
own local 3-fold axis, to the printed precision. (An earlier version of this note quoted `|m| = 0.017916` for "all 8
Ir sites"; that is the sphere moment of atom 11, not Ir.) Ir moments by iteration, same input as the unpatched
`178382bc889` apart from `electron_maxstep`:

| iteration | patched sphere | patched Hubbard | patched scf acc. (Ry) | unpatched sphere | unpatched Hubbard |
|---|---|---|---|---|---|
| 1 | 0.439 | 1.070 | 6.0 | 0.439 | 0.348 |
| 10 | 0.420 | 0.547 | 1.8e-2 | 0.235 | 0.100 |
| 20 | 0.352 | 0.480 | 6.9e-3 | 0.133 | 0.064 |
| 47 (converged) | - | - | - | 0.1225 | 0.058 |

Sphere = `magnetization` in the `report` block (mu_B); Hubbard = |`Atomic magnetic moment`|. The patched values are
still drifting down at iteration 20 and are not a result.

**Caveat.** This run did not converge: `convergence NOT achieved after 20 iterations: stopping`, with
`electron_maxstep = 20`, `mixing_mode = 'plain'`, `mixing_beta = 0.03` from a cold atomic start (no
`startingpot`/`startingwfc = 'file'`); estimated scf accuracy was still oscillating in the 0.005-0.05 Ry range
against `conv_thr = 1.0d-7`. That reads as an under-provisioned run (too few steps and too conservative mixing
for noncollinear+SOC+DFT+U from a cold start), not a symmetriser problem — the moment pattern is already exact
at whatever iteration state QE printed at stop, and the symmetriser acts identically every iteration regardless
of how close the density is to self-consistency. No unpatched comparison was run here (unlike the local build);
`pw_unpatched.x` sits in the directory unused. Not yet done: rerun with `electron_maxstep` raised to ~60-100 (and
consider `mixing_beta`/`mixing_mode`) to get a converged patched base state, and the matching unpatched cluster
run for a paired comparison at production scale.

## What this changes for the project

* Every YIO DFT run so far (base, mode 1, ...) used `symmetry_with_labels` with four Ir species, so its Hubbard
  occupations, potential, bands and Wannier model carry the artefact. The 4 to 14% curvature floor, the 3 to 16
  meV band violations and the tilted Ir moments have this one cause.
* It is not just noise. On an exactly symmetric AIAO input the symmetriser cuts the Hubbard spin polarisation to
  a third; on the first raw data of the local run it returned 0.73 of the correct value, and over the next
  iterations the unpatched map drove |m| down (0.78, 0.19, 0.08) while the patched one grew (1.07, 1.45). So the
  U-driven exchange enhancement on Ir is presumably too weak in the SCF that converged (sphere moment 0.12
  mu_B, Hubbard-projected |m| 0.05). **Not measured**: it needs a *converged* patched (or `nosym`) SCF on the
  cluster — the 2026-09-22 cluster run above confirms the symmetriser but did not converge (see "Cluster
  confirmation" section).
* Not affected, from their SCF inputs: SrTiO3 has no Hubbard correction, and in MnBi2Te4 the `HUBBARD` card is
  commented out (checked `base/soc` and `T1+/mode1`); `new_ns_nc` only runs with DFT+U. Other MnBi2Te4 inputs
  were not all read.

## Workarounds without patching

* `nosym = .true.` for the SCF (64 k instead of 8; about 2 h against a 20 h NSCF). No symmetriser runs.
* Drop `symmetry_with_labels`: then only operations that fix every Ir species individually are kept, so no
  Hubbard atoms are permuted. More k-points than 8, fewer than 64.

## Upstream status (checked 2026-09-21, read-only)

* **Not fixed on `develop`.** `PW/src/new_ns.f90::new_ns_nc` on gitlab.com/QEF/q-e `develop` still has
  `nb = irt (isym, na)`, the same contraction, no `invs`. QE reports go to the GitLab tracker (`CONTRIBUTING.md`);
  the GitHub repository is a mirror.
* **No issue** matches `new_ns_nc` or `symmetry_with_labels` (GitLab API search, all states).
* **There is an open draft, MR !2512** "Draft: Update symmetrization of Hubbard occupations" (Ryota Masuki, opened
  2024-12-11, last updated the same day, no reviewer, one note; attachments `Symmetrizations.pdf`, an FeO input and
  two outputs). It edits only `new_ns_nc` and removes the same mismatch from the other side: it **keeps
  `nb = irt(isym,na)`** and transposes the index order of `d_spin_ldau`, `d1`, `d2`, `d3` (raw diff lines 50 to 73
  for l = 2), moving the `CONJG` to the other factor. Its test is a noncollinear run of collinear antiferromagnetic
  FeO without SOC, where it "gives the same result as develop". That test cannot tell the forms apart: there the spin
  part is diagonal and (my inference from the description, input not read) only operations that are their own
  inverse permute the two Fe.
* **Checked against our operations** (`qe_ns_symmetriser_check.py`, section 6; port only, not built into `pw.x`):

  | variant | fixed point, random invariant | fixed point, AIAO ns0 | idempotent | `\|m\|` from 3 | failing terms of 48 |
  |---|---|---|---|---|---|
  | develop as written | 0.44 | 0.16 | no (0.87) | 1.0 | 28 (order 3, 6, 4T) |
  | **MR !2512 as drafted** | 0.66 | 0.26 | **yes (1e-15)** | 1.0 | **36** (order 3, 6, and 8 of 12 each of 2T, 4T) |
  | transposed, derived CONJG, `nb = irt` | 1e-15 | 3e-15 | yes | 3.0 | 0 |
  | as written, `nb = irt(invs)` (our patch) | 1e-15 | 3e-15 | yes | 3.0 | 0 |

  The draft is a projector, but onto the wrong subspace: it maps QE's own symmetric AIAO state to one of a third of
  the moment, and it breaks 8 of the 12 time-reversed 2-fold terms that develop handles. Its `CONJG` placement is
  the opposite of what `X_b = Ut X_a Ut^dag` (`Ut = D (x) conj(dS)`) requires in both branches: unitary needs
  `dS(is3,is1) d(m3,m1) nr(m3,m4,is3,is4,nb) CONJG(dS(is4,is2)) d(m4,m2)`, time-reversed needs
  `CONJG(dS(is3,is1)) d(m3,m1) nr(m4,m3,is4,is3,nb) dS(is4,is2) d(m4,m2)`. With that placement the transposed form
  and our pre-image form are the same map (the group sum is re-indexed by `g -> g^-1`), so either fix is right.
* **Reported upstream (2026-09-22).** Filed as [`gitlab.com/QEF/q-e` issue #887](https://gitlab.com/QEF/q-e/-/issues/887),
  "new_ns_nc (noncollinear DFT+U Hubbard-occupation symmetriser) pairs the rotation with the wrong atom for
  order-3/4/6 symmetry operations" — the derivation, both evidence tables, the real-`pw.x` verification, and
  the MR !2512 comparison above. A linking comment on !2512 itself is the next step, not yet posted. **2026-09-22 evening:** the full write-up (summary, derivation, two-line fix, Fe4
  reproducer before/after, YIO paragraph without recipe details, !2512 note) was posted to #887 from
  `~/Repos/q-e/repro/new_ns_nc/ISSUE_887_draft.md`. The issue body itself is still the template; the text is in the
  comments (3 as of 21:56 UTC). Still
  open: a converged SCF and an unpatched cluster-scale comparison to cite alongside the issue (see "Cluster
  confirmation" section).
* **Minimal public reproducer (2026-09-22).** The issue template asks for input/pseudopotential/output, which
  the YIO production files satisfy but shouldn't be the vehicle for (unpublished U value, k-mesh, recipe).
  Built a generic, non-YIO stand-in instead: 4 Fe atoms on a C4 orbit in a simple-tetragonal cell (`ibrav=6`),
  noncollinear in-plane "radial-vortex" starting moments, DFT+U on Fe-3d (`U=2.0`, norm-conserving
  `pseudo/Fe.pz-n-nc.UPF` — already shipped in the q-e repo itself, no separate hosting needed),
  `symmetry_with_labels=.true.` over 4 same-pseudopotential species labels (same trigger mechanism as the YIO
  report). QE finds 16 Sym. Ops. including pure (no time reversal) 90° rotations about ±z that cyclically
  permute all 4 Fe atoms. No SOC needed — the bug lives in the SU(2) spin-rotation part of the symmetriser,
  present for any noncollinear run — and no convergence needed, 2 iterations reproduce the effect exactly as
  the earlier local-build table did. Built and run from `develop@61569eb480231b649c47f631df9c5c83537df461`
  (unpatched) and `fix/new_ns_nc-preimage@74bc68eb3` (patched), ~9 s wall time, 1 rank:

  | build | iter | Fe1 | Fe2 | Fe3 | Fe4 | spread |
  |---|---|---|---|---|---|---|
  | unpatched | 1 | 1.575360 | 1.455607 | 1.575360 | 1.455607 | 7.9% |
  | unpatched | 2 | 1.593306 | 1.448095 | 1.593306 | 1.448095 | 9.6% |
  | patched | 1 | 3.030967 | 3.030967 | 3.030967 | 3.030967 | 0.00% |
  | patched | 2 | 3.096393 | 3.096393 | 3.096393 | 3.096393 | 0.00% |

  Verified directly against both `.scf.out` files (not taken on faith from the build): the unpatched split
  (atoms 1,3 vs 2,4) is exactly the pairing the wrong `nb=irt(isym,na)` produces from this orbit; patched is
  equal to full printed precision at both iterations. Two earlier design attempts (uniform +z moment; a
  single shared species with orbital-anisotropic `starting_ns_eigenvalue`) showed **zero** asymmetry despite
  using the same detected operations — a moment exactly along the rotation axis, or a basis-independent
  trace/eigenvalue observable, can't distinguish a correctly- vs incorrectly-paired unitary conjugation of
  identical inputs. The bug only becomes visible in a scalar per-atom observable when the moment is genuinely
  non-collinear with the axis, which is what the radial-vortex design forces. Files (input,
  both outputs, the patch against `develop`, and a draft of the issue text) are kept
  outside this repo, in the QE fork: `~/Repos/q-e/repro/new_ns_nc/` (git-excluded there, not on the fix branch).
  Build note: this machine's `./configure` left `DFLAGS` without an FFT backend macro; needed
  `-D__FFTW3` and the FFTW include path added to `make.inc` by hand before either branch would build.
* **Still not established:** the same-species case *without* `symmetry_with_labels` — standard `pw.x` input
  assigns starting noncollinear direction per species type, so distinguishing per-atom moments on nominally
  "the same" species inherently requires separate type labels (i.e. `symmetry_with_labels`-style setup) in
  the first place; untested whether spontaneous (SCF-driven, not starting-condition-forced) non-collinear
  symmetry breaking on truly single-labeled equivalent atoms shows the same effect.

## Not established

* Whether the bug is present in other releases; the checkout used is a 7.5 tree (`version_number = '7.5'`).
* Whether the same species case (an order-3 operation permuting two atoms of one species, no `symmetry_with_labels`)
  shows it in a real run: the algebra does not use the species, but only the four-species YIO case was run.
* That !2512 fails in a real `pw.x` run: it was tested in the port only.
* `new_nsg.f90` (DFT+U+V) uses `irt` and `d2` in a similar pattern and was not tested.
* The size of the physical change (moment, gap, curvature) after the fix.
