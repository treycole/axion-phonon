# Where `<0i|r|Rj>` should come from, and how to stop storing 40 GB of `.mmn`

Status 2026-08-27.  Closes the open questions left by
[2026-08-20_handoff_wannier_conventions.md](2026-08-20_handoff_wannier_conventions.md).

**Bottom line, in the order that matters operationally:**

1. **`_r.dat` cannot carry the position matrix at all** - not with
   `use_ws_distance`, not with any embedding, not for any material.  Two
   independent structural reasons, both verifiable in the Wannier90 source.
   See section 0.  This retires 28 GB of `_r.dat` as dead weight for this
   purpose.
2. **The `.mmn` does not have to stay on disk.**  A warm `_AA_cache/*.npz`
   reproduces `A(R)` bit-for-bit with the `.mmn` deleted; the cache is ~15x
   smaller (2.6 GB of cache replaces 40.5 GB of `.mmn`).  Verified, and the
   code gate that used to require the file is gone.
3. Between the two routes that *do* work, `.mmn` and postw90's
   `_r_full.dat`, keep `.mmn`: it is closer to DFPT and its cache is
   4.7x smaller than the `_r_full.dat` export.  The 0.5% `Z*` difference
   between them is not by itself evidence that `.mmn` is more accurate - see
   the caveat on gate 6 - but nothing argues for paying more disk to switch.

## 0. `_r.dat` + `use_ws_distance` does not work

Tested because it is by far the cheapest option: every run already has a
`_r.dat`, so using it would need no `.mmn`, no new files and no reruns.
Implemented as `apply_rdat_pair_convention` in `wannier/wannier_io.py` -
strip Wannier90's scalar WS weights, re-distribute by `min |R + tau_j - tau_i|`
(which *is* the `use_ws_distance` rule), Hermitian-project, install.

### Result: fails on both systems

SrTiO3, set (c), `max |Tr Omega_xy|` on symmetry-forced zeros (A^2):

| diagnostic | `.mmn` | `_r.dat` + ws remap | `_r.dat` raw | `_r.dat`, tau=0 |
|---|---|---|---|---|
| base generic, P+T zero | 8.6e-09 | 3.97e-02 | 4.06e-02 | 1.65e+00 |
| base mirror, zero | 4.4e-09 | 1.65e-02 | 1.91e-02 | 9.51e-01 |
| mode kx=1/2 mirror, zero | 8.1e-08 | 1.65e-02 | 1.92e-02 | 9.51e-01 |

The remap buys essentially nothing (4.06e-02 -> 3.97e-02) and a tau=0 embedding
is 40x worse.  The physical curvature is wrong too: 1.62x the `.mmn` amplitude
with correlation 0.32, versus `_r_full.dat`'s 1.22x at correlation 0.9999.

Y2Ir2O7 mode-1, relative symmetry residual (this system has a ~5% DFT/Hubbard
symmetry floor already, see [[yio-symmetry-floor-hubbard-ns]], so only large
effects are visible):

| test | part | `.mmn` | `_r.dat` |
|---|---|---|---|
| base_inversion | total | 0.079 | **1.230** |
| mode1_C2z | total | 0.133 | **1.276** |
| mode1_C3_111_axis | total | 0.085 | 0.445 |
| median over all `total` rows | | 0.085 | 0.223 |
| any `internal` row | | *identical* | *identical* |

A relative residual above 1 means the symmetry violation is larger than the
curvature being tested.  The `internal` rows agreeing exactly is the control:
that part uses only `H(R)`, which both routes share.

### Why - and why no remapping can fix it

`plot_write_rmn` (`src/plot.F90:2693-2791` in the local Wannier90 checkout):

1. **The print format is `6F12.6`.**  Six decimals, fixed point.  Measured on
   the SrTiO3 trial-03 `_r.dat`: **51.9% of all 5,038,848 components are
   exactly zero** and 70.7% are at or below the 1e-6 A granularity, while the
   median off-diagonal element at `R=(1,0,0)` is 1.3e-05 A.  Over half the
   off-diagonal position matrix is destroyed by rounding before anything reads
   the file.
2. **No link phase.**  `m_matrix` is used raw - there is no
   `exp(i b.tau_j)` (our `.mmn` route) and no `exp(i b.(tau_i+tau_j)/2)`
   (postw90 `transl_inv_full`).  The link sum over `nn` happens *inside* the
   routine, so the per-link phase is not recoverable from the output.
3. `use_ws_distance` never enters this routine.  It does not change `_r.dat`;
   it only makes Wannier90 emit the mapping separately in `wsvec.dat`.

The direct measurement of (1) and (2) together: raw `_r.dat` violates
`A(R) = A(-R)^dag` at **7.2e-03** relative (SrTiO3) and **1.3e-02** (Y2Ir2O7).
The `.mmn` builder raises above 1e-10.  Hermiticity here is exact in the
continuum and is untouched by the Hubbard symmetry breaking, so this is a
property of the file, not of the physics.

Raising the format to `ES24.16E3` would fix (1) but not (2), and would make
`_r.dat` as large as `_r_full.dat`.  Not worth it.

## 0b. Dropping the `.mmn` without losing anything

A warm `_AA_cache/*.npz` already holds the finished `A(R)`; the raw inputs were
only being consulted for a size+mtime provenance stamp.
`build_position_matrix_from_mmn` now accepts a deleted `.mmn` (or `.chk`) when
the cache records a stamp for it, warns that the deleted file's staleness can no
longer be re-verified, and still compares every surviving input exactly.
Verified end-to-end: with the `.mmn` absent, `A(R)` comes back **bit-identical**
(`max |difference| = 0.000e+00` over all 729 R blocks).

`calculations/` is 116 GB.  What is in it:

| | size | keep? |
|---|---|---|
| `*_hr.dat` | 35.0 GB | yes - this is `H(R)` |
| `*.mmn` | 40.5 GB | **no** - all 5 are backed by caches carrying valid provenance |
| `*_r.dat` | 28.0 GB | **no** for the position matrix (section 0); only the legacy control uses it |
| `*.chk` | 6.8 GB | needed for `mp_grid`/centres unless those are recorded elsewhere |
| `*_wsvec.dat*` | 2.4 GB | no - we apply the rule in house |
| `AA_*.npz` | 2.6 GB | **yes - this is the compact correct format** |

Per structure, Y2Ir2O7: `.mmn` 18.57 GB vs cache 0.75 GB, a 25x saving.

Before deleting any `.mmn`, confirm its cache carries provenance - four legacy
SrTiO3 caches have no `cache_meta` and would be rejected, but each is
superseded by a `_wbnone_ws1em05` sibling that does.

## 0c. Why the WS remap cannot be deferred (transl_inv_full does not commute)

Natural question: `transl_inv_full` supplies the phase, our code supplies the
`min |R + tau_j - tau_i|` remapping -- so why not export with
`use_ws_distance = .false.` and remap in house?  Answer: **they do not commute**,
and deferring the remap costs six orders of magnitude of symmetry.

`get_AA_R` (`src/postw90/get_oper.F90`) applies **two** phases, not one:

| | what | where |
|---|---|---|
| `phase1` | `exp(i b.(tau_i + tau_j)/2)` | k-space, per link, **before** the Fourier transform (lines 258-262) |
| `phase2` | `exp(-i b.R/2)` | R-space, per link, **after** `operator_wigner_setup` (lines 297-304) |

`phase2` is a function of `R` (`crvec_pw90`), and the Wigner--Seitz remapping
*changes* `R`.  So the remap has to happen inside the `do nn` link loop, with the
phase evaluated at the **mapped** `R`.  Sum the links first and you have frozen
`phase2` at the pre-mapping representative; remapping afterwards moves `R`
without correcting it.

That is also why our `.mmn` route is safe doing it in the other order: its phase
`exp(i b.tau_j)` depends only on the orbital, never on `R`, so it factors out
before the Fourier transform and the remap is then a pure linear redistribution.
The symmetric phase cannot be factored that way -- it necessarily leaves an
R-dependent residue.

### Measured (SrTiO3 base, set (c), 48 WF)

Exported once with `transl_inv_full=T, use_ws_distance=F`, then remapped in house:

| comparison | rel | off-diagonal |
|---|---|---|
| our remap vs postw90 per-link remap | 2.485e-04 | 2.485e-04 |
| no remap at all vs postw90 per-link | 3.078e-04 | 3.078e-04 |
| our remap vs `.mmn` | 1.237e-03 | 1.237e-03 |
| postw90 per-link vs `.mmn` | 1.257e-03 | 1.257e-03 |

So the remap ordering is worth only 2.5e-04 -- five times *smaller* than the
1.25e-03 the phase itself contributes -- and our remap recovers only ~20% of the
3.1e-04 the mapping moves.  It looks negligible.  It is not:

| position matrix | base generic, P+T zero | base mirror, zero |
|---|---|---|
| export `use_ws_distance=F` + our remap | **2.763e-03** | **7.087e-03** |
| postw90 per-link (`use_ws_distance=T`) | 8.719e-09 | 5.037e-09 |
| `.mmn` | 8.570e-09 | 4.416e-09 |

The 2.5e-04 residual is **entirely symmetry-breaking**: freezing `phase2` at the
`min |R|` representative re-injects exactly that convention's non-covariance.
The ~11x amplification from A(R) error to curvature residual is the usual
conditioning of the k-k curvature (see [[berry-curvature-conditioning]]).

**Conclusion: always export with `use_ws_distance = .true.`.**  There is nothing
to gain by deferring the remap -- the file is the same size either way -- and
everything to lose.

### `use_ws_distance` is NOT a free aliasing choice under `transl_inv_full`

For an ordinary `A(R)`, moving a block from `R` to `R + T` (a supercell
translation) is invisible on the ab-initio mesh, because `exp(ik.T) = 1` there;
it changes only the interpolation between mesh points.  Measured on the
production export, moving every block by one supercell vector:

| | A(k) on the ab-initio mesh | A(k) off-mesh |
|---|---|---|
| pure aliasing shift, no R-dependent phase | **4.79e-15** | 1.39 |

That is the textbook statement, confirmed to machine precision.

It does **not** survive `transl_inv_full`, because `phase2 = exp(-i b.R/2)`
makes the block's *value* depend on which representative it sits at.  Comparing
the same structure exported with `use_ws_distance` T vs F:

| | A(R) | A(k) on the ab-initio mesh | A(k) off-mesh |
|---|---|---|---|
| `use_ws_distance` T vs F | 3.08e-04 | **1.71e-04** | 3.41e-04 |

1.7e-04, not 5e-15.  So under `transl_inv_full` the representative choice is
part of the *definition of the approximation*, not a choice of interpolant:
picking the nearest image is what keeps `b.(bond)` small, which is what makes
the finite difference translationally invariant.  This is the same fact as the
non-commutation above, seen from the k side.

### What is actually remapped in the production `rfull` path

Worth stating explicitly, because "our remap breaks it" is the wrong reading:

| operator | who remaps it | how |
|---|---|---|
| `A(R)` | **postw90, internally** | per link, inside the `do nn` loop, with `phase2` evaluated at the mapped `R` |
| `H(R)` | **us** | `apply_orbital_dependent_R_mapping` before `W90.model` |

`install_postw90_position_matrix` performs **no** remapping -- it validates
Hermiticity and the `R=0` centres and assigns the blocks verbatim.  It also
*refuses* to run unless `H(R)` has already been pair-mapped (every
`_ham_deg == 1`), because the two operators must share one convention.

Both halves are load-bearing (SrTiO3 base, P+T-forced zero):

| configuration | max abs Tr Omega_xy |
|---|---|
| `H` pair-remapped + `_r_full.dat` (production) | 8.719e-09 |
| `H` left in Wannier90's `min \|R\|` + `_r_full.dat` | **8.997e-04** |
| production + our remap applied to `A(R)` on top | 8.665e-09 |

So our remap applied *to `_r_full.dat`* is harmless and nearly redundant -- it
moves `A(R)` by 1.8e-04 (different tie-breaking tolerance, both choices valid and
covariant) and leaves the symmetry untouched.  The failure in the table above is
not the remap failing; it is that a `use_ws_distance=F` export is **already
damaged**.  Its `phase2` was evaluated at the `min |R|` representative and baked
into the *values*.  Remapping changes which `R` a block is filed under; it cannot
un-bake a wrong number inside the block.


## 1. The cluster is not needed for the export

The rebuilt Wannier90 does not have to live on Rutgers.  A serial build with
Homebrew gfortran and Accelerate reproduces the cluster export exactly:

```
cmake -S ~/Repos/wannier90 -B <build> -DCMAKE_BUILD_TYPE=Release \
      -DWANNIER90_MPI=OFF -DWANNIER90_TEST=OFF -DWANNIER90_INSTALL=OFF \
      -DCMAKE_Fortran_COMPILER=gfortran -DBLA_VENDOR=Apple
cmake --build <build> -j 10
```

`postw90.x` needs only `seedname.{win,chk,eig,mmn,nnkp}` and takes **25 s** per
structure for SrTiO3 (48 WF, 100 bands, 8x8x8, one core).  The `.mmn` may stay
where it is and be symlinked in, so the export costs no extra disk.

The `.win` needs three keywords; everything else can be turned off:

```
transl_inv_full = true
use_ws_distance = true
write_aa_r      = true
```

This matters for the deployment plan: the export is a *local, cheap*
post-processing step, not a cluster job.  What still has to happen on the
cluster is the thing that was always expensive - the SCF/NSCF and the `.mmn`
itself.

## 2. The local build reproduces the cluster export

Running the whole Born-charge workflow on the locally produced
`_r_full.dat` gives `Z*_zz = 7.33296` at nk=16 - **identical to the value the
cluster-built exporter produced** (see the handoff).  The exporter is
compiler- and platform-independent.

## 3. What the two constructions actually differ by

`validation/inputs/compare_position_sources.py` (new) answers this directly.

| | set (c) base | set (c) mode | set (b) base | set (b) mode |
|---|---|---|---|---|
| `A(R)` rel. difference | 1.257e-03 | 1.265e-03 | 2.144e-03 | 2.162e-03 |
| ... diagonal (centres) only | 4.63e-06 | 2.29e-05 | 1.12e-05 | 1.65e-04 |
| ... off-diagonal only | 1.257e-03 | 1.265e-03 | 2.144e-03 | 2.156e-03 |
| `A(k)` rel. diff **on the ab-initio mesh** | 1.249e-03 | 1.259e-03 | 2.143e-03 | 2.162e-03 |
| `A(k)` rel. diff off-mesh | 1.316e-03 | 1.324e-03 | 2.251e-03 | 2.279e-03 |

Two things follow, and they are the whole story:

1. **The difference is purely off-diagonal.**  Both constructions reproduce the
   Wannier centres (the `R=0` diagonal) to ~1e-5 relative.  They disagree only
   about `<0i|r|Rj>` for `i != j`.
2. **The difference is already there on the ab-initio mesh**, at the same size
   as off-mesh.  So this is *not* an aliasing-class/representative choice - if
   it were, `A(k)` would agree exactly at the source k-points.  It is a
   genuinely different finite-difference discretisation of the off-diagonal
   connection: postw90 gives each link the symmetric phase
   `exp(i b.(tau_i+tau_j)/2)` and keeps the links separate through the Fourier
   transform, while our route applies the right-centre phase `exp(i b.tau_j)`
   and sums the links first.

**Consequence for convergence:** a denser *interpolation* mesh cannot shrink
this.  Only a denser ab-initio `mp_grid` can.  Confirmed empirically - the
`Z*` gap between the two routes is 0.03196 / 0.03195 / 0.03195 at nk = 8 / 12 /
16, i.e. completely flat in the BZ mesh.

## 4. Physical observables

### 4a. Born effective charge, `Z*_{Ti,zz}` (DFPT reference 7.28829)

Same base and displaced directories for every row; nk = 16.

| position matrix | set (c), 48 WF | set (b), 30 WF, occupied-only |
|---|---|---|
| `.mmn`, right-centre phase | **7.30101  (+0.17%)** | 7.370126295  (+1.12%) |
| `.mmn` + MV diagonal terms | 7.30280  (+0.20%) | - |
| postw90 `transl_inv_full` | 7.33296  (+0.61%) | 7.370126307  (+1.12%) |

Term by term for set (c), `rfull - mmn`: `z_int` **+0.00000000** (H-only, so
identical by construction), `z_cross` +0.00358715, `z_ext` +0.02835863.  The
external term carries 89% of the difference, exactly where `A(R)` enters.

**The occupied-only set is immune.**  For set (b) the two routes agree to
**1.24e-08** in `Z*` although their `A(R)` differ by 2.2e-03.  With
`n_occ = num_wann` the occupied projector is the identity, the response
collapses onto `Tr A(R=0)` - the diagonal - and the entire off-diagonal
disagreement cancels identically.  This is the sharpest confirmation that the
difference lives only in the empty-state sector.

### 4b. Cubic-forbidden components of `Z*`

| | `z_ext,x` | `z_ext,y` |
|---|---|---|
| `.mmn` | 1.639e-05 | 1.639e-05 |
| postw90 `transl_inv_full` | 4.22e-11 | 6.10e-11 |

`rfull` is **six orders of magnitude cleaner** here, and the `.mmn` violation is
nk-independent, so it is an R-space construction defect rather than a BZ
integration error.  The symmetric `(tau_i+tau_j)/2` phase is point-group
covariant on the off-diagonal blocks; the right-centre phase is not.

### 4c. Berry-curvature symmetry, `max |Tr Omega_xy|` (A^2)

`validation/symmetry/check_convention_symmetry.py --conventions mmn_pair
rdat_direct rfull_pair`.  Every entry except the last row of each block is a
symmetry-enforced zero.

| set / diagnostic | `mmn_pair` | `rfull_pair` | `rdat_direct` |
|---|---|---|---|
| (b) base generic, P+T zero | 1.17e-09 | 1.26e-10 | 3.23e-16 |
| (b) base mirror, zero | 6.95e-10 | 2.99e-10 | 8.59e-16 |
| (b) mode kx=1/2 mirror, zero | 6.71e-10 | 9.17e-11 | 7.58e-16 |
| (b) mode generic, **physical** | 3.365e-03 | 2.955e-03 | 3.204e-03 |
| (c) base generic, P+T zero | 8.57e-09 | 8.72e-09 | 4.06e-02 |
| (c) base mirror, zero | 4.42e-09 | 5.04e-09 | 1.91e-02 |
| (c) mode kx=1/2 mirror, zero | 8.07e-08 | 8.16e-08 | 1.92e-02 |
| (c) mode generic, **physical** | 6.812e-03 | 8.335e-03 | 4.023e-02 |

`rfull_pair` passes every symmetry test at the same level as `mmn_pair`.
`rdat_direct` fails at 1e-2, reproducing the known convention bug.

### 4d. The pointwise curvature is *not* convention-independent

On the displaced generic path the two good conventions give:

* set (c): correlation **0.99992**, but a near-uniform rescaling by **1.2207**
  (rms difference 22.1% of the mmn rms);
* set (b): correlation **0.4927** - completely different shapes at the same
  order of magnitude (rms difference 92% of the mmn rms).

Yet every *integrated* quantity is robust: `Z*` moves by 0.4% for set (c) and by
1e-8 for set (b), and the mean of each curve is zero to 1e-10.

So: **BZ-integrated observables are ~50x less sensitive to the position-matrix
discretisation than the pointwise curvature is.**  `dtheta/dQ` is an integral,
so it should inherit the ~0.5% robustness.  Any figure that plots `Omega(k)`
along a line, however, is convention-dependent at tens of percent at
`mp_grid = 8^3` and should say which construction produced it.

## 5. Acceptance gate for retiring `.mmn`

| # | gate | threshold | measured | verdict |
|---|---|---|---|---|
| 1 | exporter vs independent reconstruction | rel < 1e-10 | 6e-15 [1] | pass |
| 2 | `A(R) = A(-R)^dag` | rel < 1e-10 | 7.2e-14 / 7.6e-14 | pass |
| 3 | `R=0` diagonal vs `.chk` centres | < 1e-6 A | 1.3e-09 / 1.2e-08 | pass |
| 4 | Berry symmetry residuals | <= existing tol (1e-07) | <= 8.2e-08 | pass |
| 5 | occupied-only cross-check vs `.mmn` | < 1e-06 | 1.24e-08 | pass |
| 6 | `Z*` vs DFPT, no worse than `.mmn` | <= +0.17% | **+0.61%** | **fail** |

Five of six pass.  Gate 6 is the one that matters and it fails, so `.mmn` stays.

[1] Carried over from the 2026-08-20 session, which compared the cluster-built
exporter against an independent Python reconstruction.  Not re-measured here.
The independent check this session is different and arguably stronger: the
locally built exporter reproduces the cluster export's `Z*` to all printed
digits (section 2).

**Caveat on gate 6, stated honestly.**  DFPT is not an exact benchmark for this
comparison.  Both routes are O(b^2) at `mp_grid = 8^3` and the Wannier `Z*` also
carries a forward-difference error in `Q` at `D_BETA = 0.01 A`.  That `.mmn`
lands closer could be genuine accuracy or partial error cancellation; a single
benchmark point cannot separate those.  The decisive test is the one thing not
doable here - **rerun the `.mmn`/Wannierization at a denser `mp_grid`** (10^3 or
12^3) and check which route's `Z*` moves less.  Until then "keep `.mmn`" is the
conservative choice, not a proven one.

## 6. What is still open

* Denser ab-initio `mp_grid` for SrTiO3 (needs new NSCF + `.mmn` on the cluster).
  This is the only measurement that can settle gate 6.
* A second displacement amplitude (0.005 A) to separate the `Q` finite-difference
  truncation error from the position-matrix error.
* Nothing has been checked for Y2Ir2O7 or MnBi2Te4; the exports there do not
  exist yet.  The recipe in section 1 applies unchanged.

## 7. Reproducing everything here

```
# exports (25 s each; ~296 MB for 48 WF, ~110 MB for 30 WF)
postw90.x SrTiO3          # in a dir with .win/.chk/.eig/.mmn/.nnkp

# what the two constructions differ by
python validation/inputs/compare_position_sources.py

# Born charges, both routes, identical directories
python validation/observable/born_effective_charge.py \
    --position-source rfull --nks 8 12 16 --trial trial_03_Sr_sp_Ti_pd_O_sp
python validation/observable/born_effective_charge.py \
    --position-source mmn   --nks 8 12 16 --trial trial_03_Sr_sp_Ti_pd_O_sp

# Berry-curvature symmetry, three conventions, two sets
python validation/symmetry/check_convention_symmetry.py \
    --trials b c --conventions mmn_pair rdat_direct rfull_pair \
    --output-dir calculations/SrTiO3/convention_analysis_rfull
```

`_r_full.dat` now exists for: base `185116bc889` and `188725bc889` and mode
`185119bc889` (trial_03); base `183813bc889` and mode `183965bc889` (trial_02).
All are gitignored.
