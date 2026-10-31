# 19. QE k-point reduction drops −k stars for the MBT T₂⁻ group (energy 0.3 mRy too low)

**Status: FIXED** (found 2026-10-09; cause 2026-10-10; validated 2026-10-11, job 212618). A regression in QE's
`PW/src/irreducible_BZ.f90` from commit 1ce38add3 ("Simplify irreducible_BZ", Sept 2025, in 7.5CN and 7.6). It drops
whole k-stars when every lattice rotation is in the magnetic group but some only with time reversal. One-line fix:
[`qe_patches/irreducible_BZ_unfold_trev.patch`](../../qe_patches/irreducible_BZ_unfold_trev.patch);
patched pw.x `~/q-e-kfix` (md5 `db4b99b4`) = `~/q-e-segnifix` + this patch. The original
evidence is in [`../log/2026-10-08_mbt_noncollinear_gga_lda.md`](../log/2026-10-08_mbt_noncollinear_gga_lda.md) §11.1.

## Symptom

All four MnBi₂Te₄ T₂⁻ mode SCFs (LDA+U, Q = 0.01 Å, jobs 210629-210632) end **0.18-0.27 mRy below the base**, where a
displacement from a zero-force base must cost +½kQ². They also share the same **\|m\|abs of 10.27 μB** (base 9.99),
whichever atoms move.
- **First iteration:** from the converged base density it is off by 0.2-3.9 Ry (the base run itself: 2.5e-7).
- **Q-dependence:** the effect does not shrink with Q. At Q = 0.002 Å, E − E_base = −0.311 mRy.

## Why it is an artifact

The same Q = 0.002 Å structure, started from the "broken" density with **`nosym = .true.`** (job 210656: no
symmetrization, 64 k):

| | E − E_base | \|m\|abs | −½(F_b + F_m)·Δr |
|---|---|---|---|
| symmetrized (12 ops, 210645) | −0.311 mRy | 10.27 | +0.0024 mRy |
| `nosym` (210656) | **+0.0028 mRy** | **9.99** | +0.0031 mRy |

- **The `nosym` run falls back within the first iteration and converges in 6.** Its energy equals the base-branch
  harmonic value.
- **The symmetrized run's forces are nearly right** (to ~1e-4 Ry/bohr). Only the total energy and \|m\|abs are
  wrong.
- **The starting density is not the problem:** it is invariant under every mode operation to 4e-7
  (`archive/diagnostics/2026-10-08_mbt_noncollinear_gga/density_symmetry_qe.py`). So the G-space density
  symmetrizer, applied to it, is a no-op.
- **The Hubbard ns stay correct** and mirrored.

## Root cause: `irreducible_BZ` unfolds with time reversal

`setup.f90` builds the k points in two steps:
1. **`kpoint_grid`** runs *before* `find_sym`, with the lattice group (`nrot`, plain rotations, `t_rev = 0`).
   For MBT's rhombohedral lattice that is D3d, which contains inversion, so k and −k end up as one point.
2. **`irreducible_BZ(nrot, s, nsym, ...)`** then unfolds each point with the `nrot` lattice operations and reduces
   again with the crystal's `nsym` magnetic operations. For those, k → S k, or −S k when the operation carries time
   reversal.

Since 1ce38add3, the unfold step also flips k → −S k for every operation with `t_rev(irot) = 1`
(`irreducible_BZ.f90:117`). The input set is irreducible under the *plain* lattice rotations, so the unfold must use
them. Before 1ce38add3, `irrek_nc` applied `t_rev` only in the reduction.

**When it bites:** `nsym == nrot`, with some rotations present only with T. Then the unfold actions are
{S : unitary} ∪ {−S : with T}. If those don't span the lattice group, part of each star is never generated.
- **MBT T₂⁻ group:** unitary E, C3±, C2′ (with f); with T: i, S6±, σd. The unfold actions are E, C3±, C2′, −i = E,
  −S6 = C3, −σd = C2′, i.e. only D3, so −k is never generated. Result: **13 k instead of 16**, three −k orbits
  missing and their partners double-weighted (12/64 where the orbit is 6/64). The orbits cover 52 of 64 grid points.
  The density is still symmetrized consistently, but it's built from an incomplete, re-weighted zone. That gives a
  self-consistent, non-variational state: E 0.3 mRy too low, |m|abs 10.27, forces nearly right.
- **MBT base:** inversion is unitary, so the unfold actions span D3d; 13 k is correct there.
- **YIO (all runs):** the base has unitary inversion. The modes have nsym = 24 < nrot = 48, and the remaining lattice
  operations carry `t_rev = 0`. All checked lists are exact (base YP, modes 1 and 2: 8 k, weights = orbit sizes,
  64/64 covered).

**Evidence:** [`archive/diagnostics/2026-10-10_qe_irreducible_bz/kpoint_weights_check.py`](../../archive/diagnostics/2026-10-10_qe_irreducible_bz/kpoint_weights_check.py)
(log next to it).
- **Orbits check:** reads QE's printed operations and k list and compares the weights to the true orbits of the
  magnetic group.
- **Python port:** of `kpoint_grid` + `irreducible_BZ`. The develop version reproduces QE's 13-point list exactly;
  the patched version gives the correct 16.
- **Patched pw.x:** prints `number of k points= 16`, with the −k partners of k5, k6 and k12 added (212618).

**Earlier suspects, cleared:**
- **`PAW_symmetrize` (becsum magnetization):** pairs the image atom `irt(S,ia)` with the inverse spin rotation
  `s(·,·,invs(S))`, a consistent pairing. It flips m for every improper (`inv…`) and every time-reversed op.
- **`sym_rho`:** the density check found it consistent.

## Consequences

- **Every MBT T₂⁻ run with symmetry and the unpatched binary is wrong**, LDA and the earlier PBE ones alike.
  Rerun them with the patched pw.x.
- **Not affected:** the YIO runs (k lists checked exact; modes 1 and 3 pass energy vs force with ratio 1.000).
- **Also hit:** any QE 7.5CN/7.6 noncollinear magnetic run whose magnetic group contains every lattice rotation,
  some only with time reversal.
- **No base run is affected:** no fractional translations.
- **The mode Wannier chain is held** (210633-210640).

## Validation (212618, 2026-10-11)

Mode 3 at Q = 0.002 Å, symmetry on, the patched pw.x, the same broken-state start as 210645. Output:
`data/MnBi2Te4/archive/qe_patches/irreducible_BZ/mbt_t2m_mode3_Q0p002_patched_212618.scf.out`.

| | k points | E (Ry) | E − E_base | \|m\|abs |
|---|---|---|---|---|
| unpatched, symmetry (210645) | 13 | −6947.39582107 | −0.311 mRy | 10.27 |
| `nosym` (210656) | 64 | −6947.39550713 | +0.0028 mRy | 9.99 |
| **patched, symmetry (212618)** | **16** | **−6947.39550724** | **+0.0027 mRy** | **9.99** |

- **Energy:** agrees with `nosym` to 1.1e-7 Ry (conv_thr 1e-7). It falls back from the broken state within 3
  iterations and converges in 10. Total magnetization is 0.
- **Forces:** on the displaced Te they agree with `nosym` to ≤ 5e-5 Ry/bohr; the unpatched run was off by up to
  1.6e-4. The Mn forces (±3e-4, near zero) differ at the level of the SCF correction (1.3e-3 total).

## Remaining

- **Production binary: done 2026-10-11** (user). `db4b99b4` is in CLAUDE.md; source branch `production/axion-phonon`
  (fa9b19c8f). Dry runs of both binaries on every current input: byte-identical k lists for the YIO base, YIO
  modes 1-3 and the MBT base; only the MBT T₂⁻ modes change (13 → 16).
- **MBT: reruns submitted 2026-10-11:** `phonon_lda/T2-/mode_{1..4}/Q_0p01_kfix_2026-10-11/`, jobs 212708-212711
  (SCF only; the NSCF waits for the Wannier mesh). **Converged 2026-10-12:** 16 k, |m|abs 9.99, E above the base,
  energy vs force 1.01-1.04 (Mn mode 1.24, within force precision); log 10-11 §5.1. The old runs 210629-210632 and their held Wannier chain
  (210633-210640, removed) are invalid.
- **Upstream:** MR draft [`../upstream/qe_irreducible_BZ_unfold_trev.md`](../upstream/qe_irreducible_BZ_unfold_trev.md),
  branch `fix/irreducible_BZ-unfold-trev` rebased on develop c5a3cf95b (308a2f8b2). 1ce38add3 is Jae-Mo Lihm's.
