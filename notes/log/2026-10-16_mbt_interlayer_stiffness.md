# 2026-10-16 — MBT interlayer stiffness: corrected gap counting and experiment

**Result.** The conversion from the rigid-layer generalized stiffness to a spring per van der Waals gap in
[log 10-15](2026-10-15_mbt_rigid_layer_mode.md) was too large by two. The corrected effective spring density is
about **2.9–3.0 × 10¹⁹ N/m³**. The computed normal-mode frequencies are unchanged. This is roughly half the
experimental chain-model stiffness, and close to the earlier DFT value quoted by the experimental paper.
This audit uses existing structures, energies, forces, and literature; no new electronic-structure runs.

## Coordinate and gap counting

The 14-atom magnetic primitive cell contains two septuple layers. Its basal area is 16.2640 Å², volume
443.7479 Å³, and effective height 27.284 Å. Each layer has mass 983.310 amu. The rigid displacement moves seven
atoms by +Q along z and seven by −Q. There are two van der Waals gaps per periodic cell, including the gap across
the boundary; their widths change by +2Q and −2Q.

For identical nearest-neighbor gap springs K (energy per squared relative layer displacement),

\[
\Delta E=\tfrac12 K(2Q)^2+\tfrac12 K(-2Q)^2=4KQ^2=\tfrac12 k_Q Q^2,
\qquad K=k_Q/8.
\]

The kinetic energy is \(T=\tfrac12(2M_{\rm SL})\dot Q^2\), so the rigid-layer Rayleigh frequency obeys
\(\omega^2=k_Q/(2M_{\rm SL})=4K/M_{\rm SL}\). The existing phonon calculation already uses the full mass
matrix correctly; the factor-of-two error was only in the written interpretation of K.

| Existing calculation | k_Q (eV/Å²) | K = k_Q/8 (eV/Å²) | K/area (10¹⁹ N/m³) |
|---|---:|---:|---:|
| Rigid Q = 0.03 Å, energy | 2.39853 | 0.29982 | 2.954 |
| Rigid Q = 0.05 Å, energy | 2.40113 | 0.30014 | 2.957 |
| Rigid runs, forces | 2.40764–2.40972 | 0.30096–0.30122 | 2.965–2.967 |
| Production basis force matrix, rigid projection | 2.34570 | 0.29321 | 2.888 |

The rigid energy curvature gives 18.21–18.22 cm⁻¹ (about 0.546 THz). Diagonalizing the production force matrix
gives 17.7588 cm⁻¹ (0.532396 THz), as before. The normal mode permits internal adjustment within the layers,
and the independent force constructions also have the previously documented small numerical difference.

Raw inputs: `data/MnBi2Te4/phonons/base_forces_2026-10-09/210628bc889/`,
`data/MnBi2Te4/phonons/rigid_layer/Q_0p03_2026-10-14/213463bc889/`, and
`data/MnBi2Te4/phonons/rigid_layer/Q_0p05_2026-10-14/213464bc889/`.
Production basis runs and their validation are in the final section of log 10-15 and
`analysis/phonon_frequencies.ipynb`.

## Experimental comparison

[Bartram et al., npj Quantum Materials 7, 84 (2022)](https://www.nature.com/articles/s41535-022-00495-x)
report an interlayer breathing stiffness per area of (5.7 ± 0.1) × 10¹⁹ N/m³ and a fitted chain frequency
scale f₀ = 760 ± 5 GHz (25.35 cm⁻¹). Their Discussion also quotes a previous bilayer DFT stiffness of
2.85 × 10¹⁹ N/m³ from [Rodriguez-Vega et al. (2022)](https://doi.org/10.1021/acs.jpclett.2c00070).

The experimental modes are finite-flake standing waves. Comparing their fitted f₀ with our alternating-layer
bulk mode is a nearest-neighbor-chain extrapolation, not a direct observation or optical assignment of our
T₂⁻ mode. Our stiffness is about 51–52% of the experimental fit; our frequencies are about 70–72% of its
zone-edge scale. The earlier DFT value provides context, not independent validation of our calculation.

With longer-range interlayer interactions, K inferred from this single alternating displacement is an
effective spring: it does not independently determine the complete dispersion or long-wavelength elastic
response. Substrate, structure, and electronic-structure approximations are possible sources of differences,
but this comparison does not identify their cause.

## Consequence for this paper

Keep the calculated frequency and state the experimental comparison and model assumptions. Do not claim
percent-level agreement with experiment or use the numerical force-convergence precision as an experimental
accuracy estimate. The stiffness conversion does not change dθ/dQ. A sensitivity estimate using the two
frequency scales can show how assumed phonon stiffness affects a response estimate, without presenting it as
a new self-consistent prediction. Resolving the stiffness discrepancy is not necessary to report the present
microscopic coupling with these limits stated.
