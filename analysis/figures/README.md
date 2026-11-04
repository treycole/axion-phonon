# Paper figures

One notebook per figure of `paper/draft/main_draft.tex`, `fig<n>.ipynb` for Fig. n (numbering as of 2026-10-21; renumber when
figures move). Each writes `paper/figures/fig<n>.pdf` (+ `.png`). Settings at the top, code inline.

| notebook | label | status |
|---|---|---|
| fig1 | `fig:overview` | **done**: (a) YIO Ir₄ tetrahedron + O(48f) + 8 Y seen along [001], (b) MBT two septuple layers in the Fig. 7 style, with moments (black) and axion charges Z^θ_κ (red; directions from the Q = 0.01 Å pattern inputs, z_κ from the paper tables), (c) schematic θ(Q) |
| fig2 | `fig:selection` | stub (schematic); waits for the selection-rule section to settle |
| fig3 | `fig:bands` | **done**: DFT vs Wannier bands with a PDOS column per material (projwfc.x 220778 YIO, 219776 MBT) |
| fig4 | `fig:workflow` | **done** (schematic): each stage, its condition, and what violating it did |
| fig5 | `fig:thetaQ` | **done** (one column, 11-03): (a) YIO θ(Q) from the 3-form against the four-curvature slope, (b) each estimate over its Q → 0 limit: YIO 3-form chord and undistorted-end four-curvature, MBT pattern 3 undistorted end |
| fig6 | `fig:yio` | **done**: (a) nk sweep (2026-11-03), (b) symmetry residuals per operation (pattern 1), (c) the three modes on the tetrahedron, (d) c₂ on the plane ⊥ [110] (12³) |
| fig7 | `fig:mbt` | **done**: (a) the four T₂⁻ modes with couplings, (b) k-density convergence (dense tube), (c) c₂ maps of two modes on the k_z = 0 and Z planes |
| fig8 | `fig:detection` | stub (schematic + estimates); waits for the detection section |
| fig9 | `fig:harmonic` | **done** (appendix, 11-03): ΔE(Q) against −½Qv·ΔF along every basis pattern, YIO (a) and MBT (b) (the panels a, d of the earlier Fig. 5) |
| c2_distribution | (draft, not placed) | c₂ in Cartesian k on the densest data: MBT pattern 3 on the dense tube (16, 16, 5, 32), YIO 12³; concentration curves. Every mode: `../mbt/c2_distribution_modes.ipynb`, `../yio/c2_distribution_modes.ipynb` |

**Colours** (shared): magnetic ion (Ir, Mn) `#7b3294`, other cation (Y, Bi) `#2b6cb0`, anion (O, Te1) `#d9821e`, Te2
`#f7dc6f` (lines `#d4ac0d`); arrows `#c0392b`. Fig. 3's PDOS uses the element colour, solid for the main valence
orbital and dashed for the second.
