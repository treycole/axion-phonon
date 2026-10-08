# #18 — Noncollinear PBE has a spurious lower-energy, layer-asymmetric MnBi2Te4 state

> **Status: WORKED AROUND for MnBi2Te4** (2026-10-08: base moved to LDA+U, `soc_lda/u_4.0`, job 209653, converged in
> 12 iterations). **YIO (#17), same family:** LDA converges the YIO base in 21-22 iterations where PBE wanders (2026-10-08).
> The functional switch for YIO is not decided.
> - **Not a code regression:** the QE 7.2 → 8.0dev audit found nothing that explains it. It's a limitation of
>   QE's noncollinear GGA.
> - **Full record:** [`../progress/2026-10-08_mbt_noncollinear_gga_lda.md`](../progress/2026-10-08_mbt_noncollinear_gga_lda.md).

**Found:** 2026-10-05 to 10-08. **Affects:** MnBi2Te4 base SCF with SOC (therefore noncollinear) and PBE, from a
random start (`startingwfc = 'atomic+random'`, the default). Every MBT PBE+SOC base run on the new QE build.

## Symptom
- **No convergence:** base SCFs stall at 1e-5 to 1e-3 Ry and never reach `conv_thr`.
- **|m|abs climbs:** from about 10.3 to 12.6-13.5 μB.
- **Layers stop mirroring:** the induced Bi/Te moments of the Mn1 layer and the Mn2 layer no longer mirror
  (Te7/Te8 +0.016/−0.032 μB against ±0.027), and Mn1/Mn2 ns differ.
- **Lower energy:** the state is **55 mRy below** the symmetric one, with the gain in the grid XC energy (−285 mRy).

## Root cause
1. **The symmetry isn't enforced.** QE drops MBT's anti-translation {1′|½,½,½} ("This is a supercell"), so nothing
   enforces the equivalence of the two Mn layers.
2. **The random start breaks it.** The random part of the starting wavefunctions breaks the equivalence by about
   4e-3 in iteration 1.
3. **Noncollinear PBE has a lower spurious state.** With `lsign = .true.` QE's GGA uses n↑,↓ = ½(n ± sign(m_z)|m|).
   A transverse texture (∫|m_⊥| = 3.46 μB, against 0.10 in the symmetric state) makes m_z cross zero at finite
   |m|. There the GGA spin density **jumps** (84,258 grid links; total variation +48.5% over m_z). PBE's gradient
   dependence rewards it, and the potential omits ∂s/∂m, so the SCF can't settle.

## Why it is a bug (for our purposes)
- **The same physics without SOC converges.** Collinear PBE (C1) converges symmetric in 13 iterations. Noncollinear
  PBE with moments along ±z (C2), the same physical problem, finds a lower, non-converging state.
- **LDA has no such state.** E1 (random start) and E2 (atomic start) converge in 12 iterations to the same energy
  (4e-8 Ry).
- **It breaks experiment.** The state breaks the experimentally established magnetic group R_I-3c.

## Fix / workaround
- **MnBi2Te4 → LDA+U** with `rel-pz` PAW: Mn and Bi official; Te generated with `ld1.x` from the official
  `rel-pbe` recipe (validated). Base 209653: converged in 12 iterations, symmetric to 3e-5 μB.
- **Avoiding it within PBE:** `startingwfc = 'atomic'` (test D converged in 18 iterations), but only for structures
  that keep the anti-translation. The T2⁻ modes break it physically, so PBE can't protect them.

## Remaining
- **The LDA U value:** U = 4.0 eV kept for now (user, 2026-10-08). The DMC benchmark optimum for LDA+U is 4.4 eV.
- **YIO:** LDA converges (see #17 update); the decision is open.
- **Atomic start in PBE is knife-edge:** the symmetric state is a saddle (progress note §4.1).
- **An optional QE patch** to enforce the anti-translation.
- **Cleanup of the PBE test runs.**
