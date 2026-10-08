# `rdat_ndegen_applied` on Y₂Ir₂O₇ mode 1: the first real base/mode pair

Status 2026-09-20. Follows
[2026-09-16_jaemo_rdat_production_closure.md](2026-09-16_jaemo_rdat_production_closure.md),
which validated the route on **one** SrTiO₃ structure. This note runs it on a
finite-difference **pair** (YIO base + Γ₂⁻ mode 1, Q = 0.01 Å, trial_02 set,
176 WF, U = 3.0) and reports what that exposes. **The mechanics are validated;
the physics is not settled** (see "What is and is not established").

## What was run

No DFT, no `pw2wannier90`, no `.mmn`/`.amn`. Jae-Mo's `wannier90.x`
(`plan5-write-ndegen-applied`, `88134317`, serial build) with `restart = plot`
on the existing `.chk` of each structure (written by Wannier90 **3.1.0** on the
cluster; the 4.0.3-based fork read them without complaint). In a `jaemo/`
subdirectory beside each production run: `.chk`/`.eig` symlinked, `.win` copied
with `restart = plot`, `transl_inv_full = true`, `write_ndegen_applied = true`
appended and `bands_plot`/`write_tb` set false. `restart = plot` never reaches
`write_chkpt` (`wannier_prog.F90:317,326`, both gated on `ldsnt`/`lwann`), so
the production `.chk` is untouched.

| | base | mode 1 |
|---|---|---|
| wall time | 213 s | 213 s |
| peak RSS | 3.95 GB | 4.46 GB |
| `_hr.dat` / `_r.dat` / `_wsvec.dat` | 1.70 / 3.33 / 0.92 GB | same |
| R vectors | 1097 (617 before folding) | 1097 |

`_wsvec.dat` header carries `write_ndegen_applied=.true.`; the `_hr.dat` ndegen
block is all ones. Then, in `calculations/Y2Ir2O7/phonon/`:

```bash
AXION_POSITION_SOURCE=rdat_ndegen_applied AXION_NKS=4,6,8,10,12 AXION_WS_PAIR_TOL=none \
  python axion_response_vs_nk.py        # 25 min, 7.7 GB peak; -> .../nk_sweep_ext_ndegen
```

## Mechanical checks (all pass)

- `A(R)` Hermiticity 0.00e+00 (exact, as in the small-system tests).
- `R = 0` diagonal of `A(R)` vs Wannier centres: **5.0e-7 Å**, the six-decimal
  print floor. This is the check added to `_build_rdat_ndegen_applied`; a wrong
  b-vector ordering (the risk in re-reading a 3.1.0 `.chk` with a 4.0.3 fork)
  would show here as ~Å. It did not.
- Minimum gaps agree with the `.mmn` route to 1e-5 eV (0.71840 vs 0.71842), so
  `H(k)` is the same; **every difference below is in `A(R)`**.

## Result: dθ (rad/Å), real part; end0/end1 = the two endpoint estimates

| nk | `.mmn` end0 | `.mmn` end1 | ndegen end0 | ndegen end1 |
|---|---|---|---|---|
| 4 | +0.003768 | +0.000572 | +0.038744 | +0.041573 |
| 6 | −0.006597 | −0.007130 | −0.004509 | −0.004501 |
| 8 | −0.005554 | −0.005669 | −0.003518 | −0.001519 |
| 10 | −0.006237 | −0.006756 | −0.003764 | −0.003401 |
| 12 | −0.006224 | −0.006746 | −0.003774 | −0.003409 |

Both routes converge by nk = 10–12 (0.3%), same sign, endpoint spread 8–10% in
each. The ndegen value is **0.61× / 0.51×** the `.mmn` value. Presumably the
"different interpolant family" gap of
[2026-08-27_rfull_vs_mmn_position_matrix.md](2026-08-27_rfull_vs_mmn_position_matrix.md), now seen on
dθ for the first time (previously only on `A(R)` and Tr Ω). It is expected to
exist; that it is a factor of two on this small number is new.

## New: base and mode do not share a Wigner-Seitz rule

`write_ndegen_applied` makes Wannier90 choose the Wigner-Seitz images **per
run, from that run's own centres**, splitting weight equally among images tied
within 1e-5 Å. The symmetric base has exact ties. In mode 1 only the four Y atoms
move (0.0173 Å each); the 24 Y-p WFs follow them by 0.0152 Å (88%), 96 of 176
centres move by > 1 mÅ and **172 of 176 by more than the 1e-5 Å tolerance**
(largest single Cartesian component 8.9e-3 Å), so ties break and a bond spread (½, ½) in the base can sit
on one image in the mode. `(H_mode − H_base)/Δβ` then contains a change of
interpolant as well as the response. The production `.mmn` route has the *same*
property for `H(R)`: `apply_orbital_dependent_R_mapping(w90)` uses each
structure's own `lattice.orb_vecs` (`wannier_io.load_wannier_pair`); the
numbers below were measured on `H(R)`. `A(R)` was **not** measured.

`position_matrix.ws_pair_consistency` measures it **for `H(R)` only** (per aliasing class, compare
which images carry weight, only where both structures have a class-summed
hopping above 1e-4 eV; report `‖dH on changed entries‖/‖dH‖`). On this pair:

| case | entries changed / compared | contamination |
|---|---|---|
| each structure from its **own** centres — Jae-Mo files | 92,736 / 6,436,834 (1.44%) | **5.50%** |
| same, `.mmn`-route in-house mapping of the stock `_hr.dat` | 92,706 / 6,436,812 (1.44%) | **5.50%** |
| **base centres for both** (one shared rule) | 0 | 0.0000 |
| mixed rules (base mapped, mode raw min\|R\|) — the CLAUDE.md #1 failure | 1,178,234 | **15.74%** |
| both raw min\|R\| (centre-independent, not symmetric) | 0 | 0.0000 |

At `significance = 1e-3` eV the own-centres value drops to 1.85%: it is
dominated by many sub-meV far hoppings. Because the tie-broken part of `dH` is
independent of Q while the physical part scales with Q, this fraction should
grow like 1/Q at fixed tie set (reasoning, **not measured**).

`load_wannier_pair(max_ws_pair_contamination=0.10)` raises above that. 0.10 sits
between the 5.5% every current pipeline already carries and the 15.7% of mixed
rules; it is a **tripwire for gross rule mismatch, not an accuracy bound**.

## How the β derivative is formed (endpoint agreement is not a check)

`compute_dtheta` takes **one** secant, `dH/dβ = (H_mode − H_base)/Δβ` and
`dA/dβ = (A_mode − A_base)/Δβ`, and evaluates `H(k)`, `∂_k H`, `A`, `curl Ω` at
β = 0 (end0) and β = Δβ (end1) with that same secant. There is no right
endpoint: the two agree only in the linear limit, and their spread measures how
much `c2` varies across the step. An older debugging note used
`dtheta[0] ≈ dtheta[1]` as a sanity check on independent Wannierizations; it
flags gross gauge trouble, not correctness, and is not used as a criterion here.

Centre shifts *are* part of this derivative: `transl_inv_full` sets the `R = 0`
diagonal of `A(R)` to the run's Wannier centres and puts `exp(i b·(r_i+r_j)/2)`
(that run's centres) into every `A(R)` entry, so both are in `dA/dβ`. The Bloch
sum's embedding τ is deliberately held at the base's (CLAUDE.md #3). The tie
inconsistency is a different use of the centres: as input to the discrete
Wigner-Seitz rule, which as a function of β is a step (exactly at a tie for the
base, (½, ½) → (1, 0) for any β ≠ 0), so its finite difference is
jump/Δβ and grows as Δβ shrinks instead of converging. That is reasoning, not
a measurement; the Δβ test below decides it.

## Why `A(R)` cannot be refolded from class sums (a retracted experiment)

An earlier version of this note reported a "refold": sum each aliasing class of
`H(R)` and `A(R)`, redistribute with the in-house rule at the base centres, and
rerun. It changed dθ at nk = 10 from (−0.0038, −0.0034) to (−0.0019, +0.0006).
**That result is void.** It assumes the images of a class are copies of one
another. True for `H(R)`, false for `A(R)`:

`hamiltonian_get_rmn` (`transl_inv_full` + `write_ndegen_applied`) applies the
real-space half of the equivariant phase, `exp(−i b·R/2)`, *per b, at the final
expanded lattice vector* ("has to be evaluated at the final lattice vector").
Here `b = ±G_i/8` and images of a class differ by `L = 8n`, so the factor
differs by `(−1)^{n_i}` between images. Checked on the base
(`calculations/diagnostics/2026-09-20_ndegen_pair/image_structure_check.py`):

| | `H(R)` | `A(R)` |
|---|---|---|
| tied two-image entries (140,412, all with odd `n`) | images identical (max diff 0.0 eV) | images differ, up to a full sign flip (`max |A_i−A_j|/|A| = 2.0`) |
| refold onto *own* centres, entries differing > 5e-6 | 26 + 8 + 2 of 15.8 M | 26 of 15.7 M single-image entries, but **104,131 of 142,046** two-image and 2,992 of 3,922 three-image entries |

The 2–3e-4 identity error I first read as noise is the refold destroying the
per-b phase structure at exactly the tie entries. `H(R)` refolds exactly;
`A(R)` cannot be re-assigned from the b-summed file at all.

## What is and is not established

Established: the route runs on a real pair from existing `.chk`s; `H(k)` matches
production; dθ converges by nk = 10–12, to (−0.0038, −0.0034) vs (−0.0062,
−0.0067) for `.mmn`; the base/mode Wigner-Seitz tie inconsistency is real,
measurable in `H(R)`, and present in the production route as well.

**Not** established: what that inconsistency does to dθ. The one experiment
that tried to remove it was invalid (above). `A(R)` tie consistency has not been
measured. Treat the sign and order of magnitude of the YIO mode-1 dθ as the
result until these are done.

*Added 2026-09-21:* the symmetry test of this pair was run
([`2026-09-21_yio_symmetry.md`](2026-09-21_yio_symmetry.md)). The Jae-Mo route obeys the
base and mode-1 magnetic point group at the same 4–14% floor as the `.mmn` route, so it is
usable symmetry-wise. The two routes disagree most on the inversion-odd part of the
mode-1 curvature (suggestive, seven k-points), which is the part dθ depends on.

Open, in order of cost:

1. Refold **`H(R)` only** (valid: images are copies), leave `A(R)` as written
   (`refold_driver.py`, `REFOLD_A=0`, untested). Says whether the `H` tie
   inconsistency alone moves dθ.
2. Δβ dependence (0.005 / 0.02 Å): decisive for whether the tie jump is an
   artefact (∝ 1/Δβ) or converges. **Not cheap:** YIO has only Q = 0.01 Å, so
   this needs a new displaced structure (SCF + NSCF, `pw2wannier90`, rigid-shift
   `.amn`, Wannier90) on the cluster.
3. `A(R)` tie consistency needs a fix inside Wannier90: Wigner-Seitz rule from
   fixed (base) centres, phases and the `R = 0` diagonal from the run's actual
   centres. It cannot be done after the b-sum. In the fork, `plot.F90:198` hands
   `wannier_data%centres` to `ws_translate_dist` (its result is cached in
   `ws_distance` and used for both `H(R)` and `A(R)`), while
   `hamiltonian_get_rmn` (`plot.F90:~346`) gets the same array separately for the
   phases and the `R = 0` diagonal. So the change is to give the *first* call
   reference (base) centres only. There is no existing option (only
   `ws_search_size`, `ws_distance_tol`). Not written or tested.
4. `.mmn` route with the base's centres for the mode
   (`apply_orbital_dependent_R_mapping(..., centers_cartesian=base_centres)`)
   — `H` side only: the `.mmn` is deleted and the `_AA_cache` was built with
   own centres.

## Regression of the rewritten code (`curvature.py`, `axion.py`)

The same `jaemo/` pair, run through the new three-stage code
(`calculations/run_axion.py`, complex64, numpy in place of
TensorFlow, frozen gauge) against the saved results above:

| | new | old | relative difference |
|---|---|---|---|
| nk=4 dtheta (base end, mode end) | 0.03874383, 0.04157307 | 0.0387438, 0.04157265 | 7.6e-7, 1.0e-5 |
| nk=6 dtheta (base end, mode end) | -0.00450949, -0.00450055 | -0.0045091, -0.00450058 | 8.5e-5, 6.5e-6 |
| c2 density, rel. L2 (nk=4, 6) | | | 3.7e-5, 4.8e-5 |
| minimum gap and its k-point | identical | | |

That is complex64 round-off, amplified about ten times in the integral by the
cancellation between regions (sum of |c2| d3k is ten times |dtheta|). Load plus
nk=4 and 6 took 240 s and 7.0 GB.

Independent check of the k planes against PythTB's own
`berry_curvature(non_abelian=True)` on the real base structure (176 WF, 156
occupied, six generic k-points, `calculations/diagnostics/2026-09-20_ndegen_pair/check_curvature_vs_pythtb.py`):
`omega_internal` equals its internal curvature, and `omega_cross + omega_external`
equals its `include_external=True` minus `False`, on all three planes to 3.6e-16
relative, comparing band traces and Frobenius norms (the matrices sit in the
occupied-band basis, so elementwise comparison of two eigensolvers is not
meaningful). The beta sector has no PythTB counterpart; it is covered by the
analytic moving-frame tests and the regression above.

## Reproducing

Scripts are in `calculations/diagnostics/2026-09-20_ndegen_pair/`
(`ws_pair_calibration.py`, `refold_driver.py`, `image_structure_check.py`); run
with the `axion` env. Unit tests: `validation/unit/test_ndegen_pair_checks.py`.
The two `jaemo/` directories hold ~6 GB each and the disk was at 97% when they
were written.
