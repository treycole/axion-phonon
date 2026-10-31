Figures for `../2026-10-27.md` (made 2026-10-27).

Made by `make_figures.py` in this folder (axion env), from the saved results:
- `mbt_dtheta_convergence.png`: MBT tube k convergence (in-plane level, k_z planes, refined box) for all four T₂⁻
  patterns, and pattern 3's two estimates vs Q with a + bQ² fits. Data: `data/MnBi2Te4/archive/phonon/T2-/mode*/Q_*/lda/dtheta_tube_conv1e-10_*/{box_H5,convergence_M8_N16}/dtheta_convergence.npz`.
- `mbt_window_study.png`: MBT pattern 3 dθ/dQ vs frozen and outer window tops. Data: `data/MnBi2Te4/archive/phonon/T2-/mode3/Q_0.01A/lda/window_study_2026-10-22/*/tube_dtheta_tube.log`.
- `symmetry_pinned.png`: curvature symmetry per operation; YIO base + modes with Wannier-centre vs atom-pinned
  positions (`archive/diagnostics/2026-10-26_yio_symmetry_conv1e-10/*.npz`), MBT modes pinned
  (`archive/diagnostics/2026-10-22_mbt_mode_symmetry/mode*_symmetry.npz`, generic k).

Copies (the script copies them; rerun the notebook to regenerate the source):
- `c2_distribution.png`: `analysis/figures/c2_distribution.ipynb` (paper draft; MBT now the 1e-10 pattern 3 tube 214627)
- `mbt_c2_patterns.png`, `yio_c2_patterns.png`: `analysis/c2_distribution_modes.ipynb`
- `pbe_vs_lda_bands.png`: `analysis/bands_pbe_vs_lda.ipynb`
