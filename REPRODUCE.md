# Reproducing the paper's data

A guide for someone new to this repository. It says where every number and figure of the paper comes from, which
code made it, and how to make it again. The results themselves (settings, every pattern's response, every convergence test) are on one page in
[`RESULTS.md`](RESULTS.md). The paper draft is `paper/draft/main_draft.tex` (see `paper/README.md`);
`README.md` at the root is David Vanderbilt's 2024 project proposal (the motivation), not a software readme.

## What is computed

For two magnetic insulators, MnBi₂Te₄ (θ = π) and Y₂Ir₂O₇ (θ = 0), the derivative ∂θ/∂Q of the Chern–Simons axion
angle along each symmetry-adapted displacement pattern that can change θ at linear order (MnBi₂Te₄: four zone-boundary
T₂⁻ patterns; Y₂Ir₂O₇: three zone-centre Γ₂⁻ patterns). Each derivative is a Brillouin-zone integral of a trace of
products of Berry curvatures (paper Eq. 4curv-reduced), evaluated on Wannier-interpolated tight-binding models built
from DFT for the undistorted structure and for the structure displaced by Q. From the pattern derivatives follow the
axion charges of the atoms and, with the frozen-phonon force constants, the coupling of each normal mode.

## Where things are

| folder | what | in git |
|---|---|---|
| `axion_phonon/` | the production Python package: readers (`io/`), Wannier models H(R), A(R) (`wannier.py`), PAW overlaps (`paw/`: `qe_port.py` from QE, `cross_structure.py` new; `notes/theory/paw_overlap.md`), curvature (`curvature.py`), β terms and c₂ (`axion.py`), k lists (`tube.py`), Chern–Simons 3-form (`chern_simons.py`), symmetry tests (`symmetry.py`). Code ↔ equation map: `notes/code.md` | yes |
| `analysis/pipeline/` | one script per step of the calculation, run on the cluster login node (steps 1–7 below; `README.md` there) | yes |
| `analysis/figures/` | one notebook per paper figure, `fig<n>.ipynb` → `paper/figures/fig<n>.pdf` (`README.md` there) | yes |
| `analysis/yio/`, `analysis/mbt/` | per-material notebooks: bands, PDOS, phonons, symmetry, c₂ maps, nk sweep | yes |
| `analysis/validation/` | dated study scripts (k-density and k sweeps, c₂ maps, PDOS) and the SrTiO₃ Born-charge check | no (ignored) |
| `tests/` | pytest suite: unit tests and checks of the production results (`tests/README.md`); also the symmetry-table runners `tests/symmetry/run_{mbt,yio}_symmetry.{sh,py}` | yes |
| `data/` | the DFT and Wannier runs and everything computed from them, per material (`data/<material>/README.md`) | no (ignored) |
| `paper/` | `main.tex` (Overleaf copy), `draft/main_draft.tex` (working draft), `figures/` | `draft/` ignored |
| `notes/` | status (`notes/README.md`), theory (`theory/derivation.md`), bugs (one file each), dated investigation log, how-tos | yes |
| `archive/` | superseded code and dated diagnostics, kept as a record | no (ignored) |

**Data locations.** The full data (wavefunction saves, `.mmn`, `.chk`) are on the Rutgers cluster under
`/mnt/wkdk/bc889/axion/<material>/`. `data/<material>/` here has the same tree with the small outputs (inputs, `.out`,
`.win`, `.eig`, logs, the c₂ chunks, the production `.chk` of the bases), which is everything the notebooks read. The
two trees share one layout, described in `data/MnBi2Te4/README.md`: `dft/<structure>/` holds what depends only on a
structure (SCF, NSCF, pw2wannier90, PAW overlaps), `windows/<window>_k<mesh>/` what depends also on the Wannier
energy window (wannier90 runs, `dtheta_<structure>/` c₂ results, `symmetry_<structure>/`).

## Software

| what | version | where |
|---|---|---|
| Quantum ESPRESSO | 8.0dev, branch `production/axion-phonon` of the `treycole/q-e` GitLab fork, commit `fa9b19c8f`, with four corrections (paper App. "Corrections to the DFT code"; patches and reports in `notes/upstream/`) | cluster `/mnt/wkdk/bc889/axion/bin/qe_db4b99b4/` (pw.x md5 `db4b99b4`, pw2wannier90.x `8f5368ba`, projwfc.x `9255d8a7`, ph.x `554fabc7`) |
| Wannier90 | 4.0.3 (develop v4.0.1-160-gc4a924a5), projection only (`num_iter = 0`, `dis_num_iter = 0`), `mpirun -np 1` | cluster `bin/`, wannier90.x md5 `26c8b506` |
| pseudopotentials | pslibrary 1.0.0 / 0.3.1 `rel-pz` PAW; Te generated with `ld1.x` (`data/MnBi2Te4/pseudo_gen/`) | cluster `/mnt/wkdk/bc889/axion/pseudo/rel-pz/` with `MD5SUMS` |
| Python | `axion` conda env (numpy, scipy, matplotlib, nbformat/nbconvert, pytest); PythTB on its `external` branch, commit `5776792` (needed for the position-matrix model; other branches fail with `no '_pos_r'`); WannierBerri 1.8.0 (validation only) | local `conda activate axion`; cluster `/mnt/wkdk/bc889/envs/axion/bin/python` |
| production code on the cluster | `/mnt/wkdk/bc889/axion/code/2026-10-29/` (`axion_phonon/`, `pipeline/`, `pythtb/`, `MD5SUMS`) | the repository's `axion_phonon/` is this code plus additions used only by notebooks (`.win` window parsing, `tube.integrate_at`, `symmetry.py`), with `paw.py` split into the `paw/` package (2026-11-03; the overlaps agree to 2e-16); the computation is identical |

Every job script checks the md5 of its binaries, pseudopotentials and code before it submits (`analysis/pipeline/common.sh`).

## The calculation, step by step

`analysis/pipeline/README.md` documents each step. For one material and one Wannier window:

| step | script | makes | for |
|---|---|---|---|
| 1 | `1_scf.sh` | SCF of a structure (conv_thr 1e-10) | base and every displaced structure |
| 2 | `2_nscf.sh` | NSCF on the 8³ Wannier mesh | base and every displaced structure |
| 3 | `3_pw2wannier.sh` | `.mmn`, `.eig` (base: also the atomic-orbital `.amn`) | base and every displaced structure |
| 4 | `4_overlap.sh` | PAW cross-structure overlaps with the base (metric S, first-order augmentation) | displaced structures |
| 5 | `5_wannier.sh` | wannier90: base from atomic orbitals; displaced structures projected onto the base's Wannier functions (`mode_amn.py`), the frozen Wannier gauge | base, then each displaced structure |
| 6 | `6_c2.sh` | c₂(k) of the base/displaced pair on a k list (both endpoint estimates), as N tasks | each displaced structure |
| 7 | `7_dtheta.sh` | ∂θ/∂Q = ∫ c₂ (`dtheta.log`, `dtheta.npz`, `c2.npz`) | each displaced structure |
| – | `cs_theta.sh`, `cs_theta_sum.sh` | θ from the Chern–Simons 3-form of each structure (θ = 0 materials only) | Y₂Ir₂O₇ check |

Run a step on the login node, wait for `OK` in its folder, then run the next; nothing is queued behind a step that has
not succeeded. Inputs written by hand: each structure's `scf.in`, and once per mesh the base's `nscf.in` and `.win`
(displaced structures' NSCF inputs and `.win` are generated from the base's).

**Production settings.**
- MnBi₂Te₄: LDA+U (U = 4 eV), window `f10.962_w12.810_k8`, 116 Wannier functions, 82 occupied, tube k list
  (16, 4, 3, 8) (`axion_phonon.tube`: c₂ is concentrated in a narrow tube around the k_z axis).
- Y₂Ir₂O₇: LDA+U (U = 3 eV), window `f12.422_w22.644_k8`, 176 Wannier functions, 156 occupied, uniform 10³ grid.
- Both: H(R) from the `.chk` + `.eig` with no hopping cut; A(R) from Lihm's translationally equivariant formula
  rebuilt from the `.chk` overlaps; Wannier positions pinned to the atoms, one embedding τ for both structures of a
  pair.

The eight ways to get a silently wrong answer are listed in `CLAUDE.md` ("The eight things..."); read them before
changing any of these settings.

## Paper item → code → data

Paths of data are relative to `data/<material>/`. Job numbers are SGE job ids on the cluster; each step folder's
`NOTE` and `OK` say which job made it.

### Figures

| figure | notebook | reads |
|---|---|---|
| Fig. 1 overview | `analysis/figures/fig1.ipynb` | base SCF inputs (positions, moments), `dft/<pattern>_Q0.01/scf/` inputs (axion-charge directions), z_κ from Tables yio-charges / mbt-charges |
| Fig. 3 bands + PDOS | `analysis/figures/fig3.ipynb` | `dft/base/bands/` (pw.x bands), production base `.chk` (Wannier bands), `dft/base/projwfc_k8/` (PDOS) |
| Fig. 5 θ(Q) | `analysis/figures/fig5.ipynb` | YIO `windows/f12.422_w22.644_k8/dtheta_*/dtheta.npz` and `cs/theta_cs.npz`; MBT `windows/f10.962_w12.810_k8/dtheta_T2-_p*_Q0.0{1,3}/` |
| Fig. 9 harmonic check (appendix) | `analysis/figures/fig9.ipynb` | YIO `dft/GM2-_p*_Q*/scf/`; MBT `phonons/amplitude_series/`, `phonons/force_constants/` |
| Fig. 6 Y₂Ir₂O₇ | `analysis/figures/fig6.ipynb` | `dtheta_nk{4,6,8,12}_*` (nk sweep), `symmetry_GM2-_p2_Q0.01/`, `mode_couplings.log`, base and pattern SCF inputs |
| Fig. 7 MnBi₂Te₄ | `analysis/figures/fig7.ipynb` | `phonons/force_constants/`, `windows/f10.962_w12.810_k8/dtheta_*`, `dtheta_dense_*` (k density), `c2_maps_*` |
| c₂ distribution (draft) | `analysis/figures/c2_distribution.ipynb` | MBT `dtheta_dense_T2-_p3_Q0.01/`, YIO `dtheta_nk12_*` |
| Fig. 4 workflow | `analysis/figures/fig4.ipynb` | none (schematic; numbers quoted from the text) |
| Figs. 2, 8 | schematics; placeholders in the draft | — |

### Tables and key numbers

| item | made by | data |
|---|---|---|
| Table patterns (YIO ∂θ/∂Q, both ends, 3-form, Q → 0) | pipeline steps 1–7 and `cs_theta.sh` at Q = 0.01–0.05 Å; fits a + bQ² + cQ⁴ in `analysis/figures/fig5.ipynb` | `Y2Ir2O7/windows/f12.422_w22.644_k8/dtheta_GM2-_p*_Q*/`, summary in `data/Y2Ir2O7/README.md` |
| Table yio-charges, yio-modes | `analysis/yio/phonon/mode_couplings.py` | `Y2Ir2O7/windows/f12.422_w22.644_k8/mode_couplings.log` |
| Table patterns-mbt, mbt-charges | pipeline steps 1–7 (tube list) | `MnBi2Te4/windows/f10.962_w12.810_k8/dtheta_T2-_p*_Q*/dtheta.log`, summary in `data/MnBi2Te4/README.md` |
| Table mbt-modes, phonon checks | `analysis/mbt/phonon/phonon_frequencies.ipynb`, `analysis/figures/fig7.ipynb` | `MnBi2Te4/phonons/` |
| Table symmetry-yio | `tests/symmetry/run_yio_symmetry.{sh,py}` (cluster; jobs 220828–30) | `Y2Ir2O7/windows/f12.422_w22.644_k8/symmetry_GM2-_p*_Q0.01/out/symmetry.log` |
| Table symmetry-mbt | `tests/symmetry/run_mbt_symmetry.{sh,py}` (cluster; jobs 220823–26) and `analysis/mbt/symmetry/mbt_symmetry.ipynb` (pattern 3) | `MnBi2Te4/windows/f10.962_w12.810_k8/symmetry_*/` and `symmetry/` |
| YIO k convergence (nk 4–12) | `analysis/validation/2026-11-03_k_sweeps/c2_variant.sh` (jobs 220779–220802) | `Y2Ir2O7/windows/f12.422_w22.644_k8/dtheta_nk*_*/` |
| MBT k convergence (dense tube), uniform meshes | `analysis/validation/2026-10-30_k_density/c2_dense.sh`; `c2_variant.sh` (jobs 220803–13) | `dtheta_dense_*`, `dtheta_uniform{4,8,16}_T2-_p3_Q0.01/` |
| window studies | `analysis/validation/2026-10-30_frozen_top/`, `..._window_parts/` (MBT); log 2026-10-17 (YIO) | `MnBi2Te4/studies/window_tests/`, `data/MnBi2Te4/README.md` |
| SrTiO₃ Born charge (validation of the mixed curvature) | `analysis/validation/sto_born_charge.ipynb`, `analysis/pipeline/dfpt.sh` | `data/SrTiO3/` |

## Re-running from scratch

1. Install the software above; check the md5s against `bin/MD5SUMS` and `pseudo/rel-pz/MD5SUMS`.
2. For each material, write the base `scf.in` (from `data/<material>/dft/base/scf/`), run steps 1–3 for the base and
   step 5 for the base in the production window.
3. For each displaced structure (`dft/<pattern>_Q<Q>/scf/<prefix>.scf.in`: the base with positions moved by Q along
   the SMODES vector, `data/<material>/phonons/SMODES*`), run steps 1–7.
4. Run the figure notebooks (`jupyter nbconvert --execute`, `axion` kernel); each starts with a settings block of
   full paths into `data/`.
5. `python -m pytest tests` checks the production results against the paper's claims (symmetry, k convergence,
   interpolation, the SrTiO₃ Born charge).

The current status of every calculation, the open questions, and the history of every decision are in
`notes/README.md` and `notes/log/`.
