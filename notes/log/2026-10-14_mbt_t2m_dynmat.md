# 2026-10-14 — MBT T₂⁻ dynamical matrix, LDA+U (first pass)

> Superseded by [2026-10-15_mbt_rigid_layer_mode](2026-10-15_mbt_rigid_layer_mode.md): reruns at Q = 0.03 Å, conv_thr 1e-10 give 17.8 / 82.7 / 128.7 / 143.3 cm⁻¹, all stable.

Script: `archive/diagnostics/2026-10-14_mbt_t2m_dynmat/t2m_dynmat.py`. It reads the run folders in
`data/MnBi2Te4/` directly. Same construction as YIO (log 10-12 §16).

## Inputs

- Base forces: 210628. Mode SCFs: 212708-212711 (Q = 0.01 Å along the raw SMODES vectors, entries ±1, `db4b99b4`).
- **The base has residual forces** (experimental cell, not relaxed in LDA+U): about 0.26 eV/Å on Bi and Te, in a
  Θ-even pattern. So Φ_ij = −v_i·[F(Q_j) − F(0)]/Q.
  - The T₂⁻ projection of F(0) is ≤ 3e-13 eV/Å, zero by symmetry, as it must be.
- SCF force corrections of the mode runs: 0.010 / 0.012 / 0.023 / 0.034 eV/Å (Mn / Bi / Te1 / Te2). The force
  changes are 0.05-0.3 eV/Å, so the noise is a few % to ~20% for the Mn mode.

## Result

- Φ (eV/Å²), basis Mn, Bi, Te1, Te2: diagonal 5.36 / 29.91 / 20.37 / 12.71. Largest off-diagonal Bi-Te1 ≈ 19.7.
  Φ − Φᵀ is 2.5% of Φ.
- Masses (12 m per basis mode, as for YIO's v normalization): 109.9 / 835.9 / 510.4 / 510.4 amu.

| ω (cm⁻¹) | meV | THz | Mn / Bi / Te1 / Te2 share of \|Δr\|² | spread from force noise (cm⁻¹) |
|---|---|---|---|---|
| −20.4 (imaginary) | — | — | 0.23 / 0.24 / 0.24 / 0.29 | 13 |
| 77.9 | 9.7 | 2.34 | 0.30 / 0.07 / 0.24 / 0.39 | 3 |
| 124.1 | 15.4 | 3.72 | 0.77 / 0.03 / 0.05 / 0.15 (Mn) | 2 |
| 139.9 | 17.3 | 4.19 | 0.39 / 0.23 / 0.37 / 0.01 | 1 |

- **The lowest mode is the rigid antiphase motion of the two septuple layers along c** (M-metric overlap 0.996).
  - In the basis it is (1, −1, 1, −1).
  - It is the longitudinal acoustic branch at the non-magnetic Z point, folded to Γ by the AFM doubling: the
    interlayer (van der Waals) breathing mode.
  - Its stiffness, c·Φ·c / c·M·c = −0.0012 eV Å⁻² amu⁻¹, is a ~3% remainder of entries of order 5-30. Φ's asymmetry
    and the force noise are that size, so **its frequency and sign are not resolved**. The imaginary value is not
    evidence of an instability.
  - The other three modes are stable to the noise.
- The spread column comes from Gaussian noise of SCF-correction size per run, a rough model.

## Consequences

- **The interlayer breathing mode is T₂⁻, so it is axion-active.**
  - Its frequency is low, so for a given energy its amplitude (zero-point or thermal) is large.
  - Its in-phase partner is the acoustic translation, not an optical Γ mode.
  - So the "Davydov partner" handle in `paper/pitch.md` applies to the three optical modes only.
- **To resolve it:** one SCF along (1, −1, 1, −1) at a larger amplitude (0.03-0.05 Å, both signs). It gives the
  curvature directly from E and F without the cancellation. The base residual forces matter most for this mode.
  A check against the c-axis relaxation would also tell whether the experimental cell is under interlayer tension in
  LDA+U.
- **dθ per physical mode:** set `DTHETA_DQ` in the script once the four basis-mode dθ/dQ exist. Until the soft mode
  is resolved, its eigenvector, and so the mixing of the other three, is uncertain at the few-% level.
- **Not compared with experiment.** IR, Raman or neutron data on the zone-boundary or folded modes of MnBi₂Te₄:
  to look up.

## Rigid-layer SCFs submitted (2026-10-14 22:53)

Jobs **213463** (Q = 0.03 Å, n136) and **213464** (Q = 0.05 Å, n116), started 23:00. They replace 213461/213462, which
landed on n247 and n138 and were cancelled after 6 min; the stubs were removed. The job scripts now carry
`-l h_vmem=3.6G,hostname=!(n138|n126|n246|n247)`. Copy-back from n138 did work for the cancelled job. Cluster
`MnBi2Te4/phonon_lda/T2-/rigid_layer/Q_{0p03,0p05}_2026-10-14/` (pw.x `db4b99b4`; NOTE in each folder).
- Each atom moves by Q along z, SL1 up and SL2 down: the pattern (1, −1, 1, −1).
- Θ maps +Q to −Q, so a second amplitude replaces the ±Q pair and also tests whether E is quadratic.
- Input: the mode-1 SCF with these positions, and conv_thr tightened from 1e-7 to 1e-9 to cut the force noise. The
  mode-1 SCF took 4.8 h.
- **Reading:** E(Q) − E_base (base 210628: −6947.39550994 Ry) > 0 and ∝ Q² means stable.
  - If the mode is a pure eigenvector, ω² = 2ΔE / (Q² Σm), with Σm = 1966.6 amu (all 14 atoms).
  - The force projection Σ r·[F(Q) − F(0)] gives the same curvature, as a check.
  - Expected scale: ΔE ≈ 0.8 mRy at 0.05 Å for 35 cm⁻¹, against a conv_thr of 1e-9 Ry.
