# Paper figures

One notebook per figure of `paper/main2.tex`, `fig<n>.ipynb` for Fig. n (numbering as of 2026-10-21; renumber when
figures move). Each writes `paper/figures/fig<n>.pdf` (+ `.png`). Settings at the top, code inline.

| notebook | label | status |
|---|---|---|
| fig1 | `fig:overview` | stub: structures, patterns, schematic θ(Q) |
| fig2 | `fig:selection` | stub: schematic |
| fig3 | `fig:bands` | **done**: DFT vs Wannier bands (copy of `analysis/bands.ipynb`; run all cells) |
| fig4 | `fig:workflow` | stub: schematic |
| fig5 | `fig:harmonic` | **done**: YIO ΔE(Q), θ(Q) from the 3-form vs the four-curvature slope, θ/(Q ∂θ/∂Q) |
| fig6 | `fig:yio` | stub: needs an nk sweep and a symmetry rerun at production settings; (c), (d) data exist |
| fig7 | `fig:mbt` | stub: waits for the MBT dθ chain |
| fig8 | `fig:detection` | stub: (b) YIO data exist |
| c2_distribution | (draft, not placed) | c₂ in Cartesian k: MBT map + radial concentration (mode 3, 1e-10 tube 214627, since 10-22), YIO (110) plane (10³, interpolated), concentration curves; every mode: `../analysis/c2_distribution_modes.ipynb` |
