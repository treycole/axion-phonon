# calculations/

Code that runs on the first-principles data. The data itself (QE and Wannier90 runs, and everything
computed from them) lives in `data/` at the repository root. That folder is not in git.

```
calculations/
  run_axion.py    dθ/dQ versus k-mesh for one base/mode pair
  analysis/       one notebook per task: bands, pdos, c2_density, dtheta_vs_nk, theta, z2_index,
                  wannier_centres, sto_validation_summary
  setup/          building the structures and QE inputs (SMODES, displacements, positions)
  diagnostics/    dated one-off investigations, each next to the note it supports
  qe_patches/     the QE 7.5 new_ns_nc fix and its check
  legacy/         the old TensorFlow dθ drivers; delete once MnBi2Te4 has been rerun through run_axion.py
```

## Convention: a settings block at the top, full paths, no flags

Every script and notebook starts with its settings as plain variables: the run folders written out in full
under `data/`, `PREFIX`, `N_OCC`, the amplitude, and each option with its choices listed in a comment next to
it. To run on something else, edit those lines. There are no command-line options and nothing searches `data/`
for runs.

A path is derived instead of written out only when it is always the same relative to one that is:
`<run>/jaemo/` for Jae-Mo's files, `<scf folder>/<prefix>.scf.out`, `<run>/symmetrized/`, and the output
folder `<mode run>/nk_sweep_<source>_<own|shared>_ws/`.

Folder names in `data/` are cluster job IDs (`178365bc889`). A plain-text `NOTE` file in a run folder (or a
folder above it) records what no file can tell, e.g. "superseded" or "QE bug"; read it with `cat`.

## Computing dθ/dQ

Edit the settings block at the top of `run_axion.py`, then:

```bash
conda activate axion
python calculations/run_axion.py
```

```python
# ============================================================================== settings
PREFIX = "Y2Ir2O7"
BASE = DATA / "Y2Ir2O7/base/soc/u_3.0/output/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181538bc889"
MODE = DATA / "Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/178365bc889/trial_02_Y_p_Ir_d_O_sp/no_displacement/181565bc889"
CACHE_BASE = DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_base_181538bc889.npz"
CACHE_MODE = DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_mode_181565bc889.npz"
N_OCC = 156
DBETA = 0.01
NKS = [4, 6, 8, 10, 12]
SOURCE = "transl_inv_full_from_chk"   # "rdat_ndegen_applied" | "mmn" | "transl_inv_full_from_chk" | "right_centre_from_chk"
SHARED_WS = True
OUT = None
```

The result is saved next to the data, in `<MODE>/nk_sweep_<source>_<own|shared>_ws/` (`nk_*.npz` holds dθ,
the c2 density and the gap; there is also `summary.csv`). Plot it with `analysis/dtheta_vs_nk.ipynb` and
`analysis/c2_density.ipynb`.

**Production A(R) (2026-10-07): `transl_inv_full_from_chk`**, Jae-Mo's translationally invariant formula rebuilt
from the `.chk` overlap matrices, with `SHARED_WS = True`. Unlike Jae-Mo's own `_r.dat` files it can be folded on
the base's Wigner–Seitz images, and it needs no `.mmn`. `right_centre_from_chk` (the `.mmn` formula from the
`.chk`) was the setting from 2026-09-30 to 2026-10-06.

The `.mmn` route reads the `.mmn` from the Wannier run folder. For SrTiO3 the `.mmn` was written to the NSCF
folder instead, so use `SOURCE = "right_centre_from_chk"` there (the same formula, built from the `.chk`).

## Status of the analysis notebooks

`bands`, `pdos`, `c2_density`, `dtheta_vs_nk`, `wannier_centres` and `sto_validation_summary` read either
`modules/` results or QE files directly. `theta` and `z2_index` still compute with PythTB and TensorFlow, not
`modules/`: their numbers are exploratory (the θ notebook warns that its gauge is near-singular). The
notebooks they replaced are archived in `tmp/calculations_archive_2026-09-24/` (git-ignored).
