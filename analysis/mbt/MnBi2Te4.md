# MnBi₂Te₄

θ = π (inversion + a composite antiunitary symmetry). Axion-active phonons: the zone-boundary **T₂⁻** modes.
Where the runs stand is in [`../README.md`](../README.md); this file is the material's setup reference.

## Current setup (2026-10-11)

- **Functional: LDA+U**, U = 4 eV on Mn 3d (`HUBBARD (ortho-atomic)`), `rel-pz` PAW, noncollinear + SOC.
  Never noncollinear PBE ([bug #18](../bugs/18_qe_noncollinear_gga_spurious_state.md)).
- **Base:** `MnBi2Te4/base/soc_lda/u_4.0/`, job 209653. Pseudopotentials and how the Te `rel-pz` file was generated:
  `soc_lda/pseudo_gen/` and [`../log/2026-10-08_mbt_noncollinear_gga_lda.md`](../log/2026-10-08_mbt_noncollinear_gga_lda.md) §11.
- **Every run:** the exact hex-CIF cell (a = 4.33360, c = 81.85200 Å, atoms at u·c), `startingwfc = 'atomic'`,
  `diago_thr_init = 1.0d-5`. Check layer mirroring and ∫|m⊥| with
  `archive/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py`.
- **Wannier trial set (2026-10-09): Bi s;p, Te s;p, Mn1 d, Mn2 d — 116 spinor WFs, bands 57-172 of the LDA base,
  n_occ = 82** (bands 57-138).
  - **Rule:** keep every occupied band unless it is truly flat core.
  - **Left out as flat core:** Mn 3s/3p (≤ 0.6 meV wide) and Bi 5d (≤ 7.3 meV).
  - **Included:** the dispersive Te 5s / Bi 6s bands (up to 1.75 eV wide), which the old 92-WF p-d set (below) left
    out.
  - **Windows:** outer −7..20 eV, frozen −7..10 eV, projection only.
  - **Run:** `base/lda/w90_sp_d_k4_210651` (sites with z > c/2 moved down by
    a₁+a₂+a₃, [#20](../bugs/20_w90_centre_wrap_long_cell.md)); 4×4×4 is too coarse near Γ, denser mesh pending.
  - **WS fold and τ: atom-pinned positions** (`PIN_TO_ATOMS = True` in `run_axion.py`). The p_x/p_y centres sit off-axis, so per-WF centres break C3
    ([#21](../bugs/21_mbt_wf_centres_not_c3_covariant.md)). The base's traced curvature is PT-forced zero, so check
    symmetry on the interpolated eigenvalues, not on Tr Ω.
- **T2⁻ modes on LDA:** `MnBi2Te4/phonon_lda/T2-/mode_{1..4}/Q_0p01_kfix_2026-10-11/`, jobs 212708-212711 (pw.x
  `db4b99b4`, SCF only, 16 k). The first runs (`Q_0p01/`, 210629-210632, 2026-10-09) had bug
  [#19](../bugs/19_qe_sym_ftrans_trev_energy.md) and are invalid.
  - **Displacement:** the SMODES basis patterns (`setup/MnBi2Te4/mbt_qe_inputs.ipynb` cell 28), Q = 0.01 Å
    along z, on the exact cell. Partners in the two Mn layers move opposite.
    - mode 1: Mn1 +, Mn2 −.
    - mode 2: Bi 17.39 +, 58.32 −, 23.54 −, 64.46 +.
    - mode 3: Te 5.46 +, 46.38 −, 35.47 −, 76.39 +.
    - mode 4: Te 12.03 +, 52.96 −, 28.89 −, 69.82 +.
  - **Each run:** SCF with `tprnfor` from the base state (`startingpot = 'file'`, `tmp.tar.gz`) → NSCF (64 k,
    `-nk 8`).
  - **Symmetry:** each mode keeps 12 operations, 6 with fractional translations. The base forces run is job 210628.
  - **Wannier chain:** the held 210633-210640 was removed 2026-10-11 (built on the invalid runs); the folders
    `Q_0p01/w90_trial_sp_d_no-displacement/` are templates.
    - Each untars its save on the node, makes the Julia `--no-displacement` `.amn` from the base
      (`w90_trial_sp_d/210622bc889`), then runs pw2wannier90 and wannier90.
    - The `.win` is the base's with the mode's displaced `atoms_cart`, YIO's convention.
- **Runs** (local `data/MnBi2Te4/`, reorganized 2026-10-13; cluster paths unchanged under `MnBi2Te4/`):

  | | local | cluster |
  |---|---|---|
  | base SCF + 4×4×4 NSCF | `base/lda/scf_nscf_209653` | `base/soc_lda/u_4.0/209653bc889/` |
  | base forces | `base/lda/scf_forces_210628` | `base/soc_lda/u_4.0/forces_2026-10-09/` |
  | base 8×8×8 NSCF (running) | `base/lda/nscf_k8_212997` | `base/soc_lda/u_4.0/209653bc889/nscf_k8_2026-10-12/` |
  | base Wannier, 4×4×4 | `base/lda/w90_sp_d_k4_210651` (windows −7..20/..10), `w90_sp_d_k4_win15_212985` (−7..15.047/..10.009) | `209653bc889/w90_trial_sp_d_recentred/` (`.mmn`/`.eig` from `w90_trial_sp_d/210622`) |
  | T₂⁻ mode SCFs 1-4 | `phonon/T2-/mode{1..4}/Q_0.01A/lda/scf_21270{8,9,10,11}` | `phonon_lda/T2-/mode_{1..4}/Q_0p01_kfix_2026-10-11/` |
  | T₂⁻ 8×8×8 NSCFs (running) | `phonon/T2-/mode{1..4}/Q_0.01A/lda/nscf_k8_{212998..213001}` | `.../Q_0p01_kfix_2026-10-11/nscf_k8_2026-10-12/` |
  | T₂⁻ force-constant SCFs (Q 0.03 Å, conv_thr 1e-10; frequencies only) | `phonon/T2-/mode{1..4}/Q_0.03A/lda/scf_2135{55..58}` | `phonon_lda/T2-/mode_{1..4}/Q_0p03_2026-10-15/` |
  | rigid-layer SCFs (Q 0.03 / 0.05 Å; check) | `phonon/T2-/rigid_layer/Q_{0.03A,0.05A}/lda/scf_{213463,213464}` | `phonon_lda/T2-/rigid_layer/Q_{0p03,0p05}_2026-10-14/` |

  The SCF/NSCF folders hold text only (inputs, outputs, job script, NOTE). Every PBE run (the paper's old pipeline,
  the U scan, QE 7.6 debugging, T₁⁺) is in `data/MnBi2Te4/archive/`, whose README lists each one. Templates and the
  cluster workflow scripts: `archive/setup_templates/MnBi2Te4/`.
- **Cluster:** compute node **n138 has no `/mnt/wkdk` mount** (and no `/scratch`). Exclude it (`-l hostname='!n138'`)
  from any job that reads the shared filesystem. n126/n246/n247 were excluded with it until 10-17 without a reason:
  they mount `/mnt/wkdk` and have NVMe `/scratch` (checked 10-17), so they are usable.
- Past MBT-specific bugs: [#9](../bugs/09_mbt_t2_rounded_cell_and_labels.md), [#14](../bugs/14_mbt_stale_binary.md),
  [#15](../bugs/15_mbt_u3_base_rounded_cell.md), [#18](../bugs/18_qe_noncollinear_gga_spurious_state.md).

---

## Setup notes (early 2026)

The working notes from the initial setup (SMODES, mode selection, Wannier windows). The procedure still holds; the
specific windows were for PBE without SOC. The files the text calls `SMODES.in.txt`/`SMODES.out.txt` are now `data/MnBi2Te4/phonons/SMODES.in`/`.out`.

### Phonon mode calculation

### Packages

### SMODES
- SMODES is a symmetry analysis package for phonon modes
- https://iso.byu.edu/smodes.php
- Can find irreps of phonon modes at different k-points
- Input: crystal structure (lattice vectors, atomic positions, space group)
- Output: list of irreps, degeneracies, and displacement patterns for each mode

#### Input
See SMODES.in.txt for input file. The (a, b, c, alpha, beta, gamma) parameters are obtained from MnBI2Te4_hex.cif. This file is from springer https://materials.springer.com/isp/crystallographic/docs/sd_1940632. The Wyckoff positions are from PhysRevMaterials.3.064202.
For some reason, using Wyckoff from cif file gives different atomic positions than those in nscf output. 

#### Output
See SMODES.out.txt for output file. This contains the list of irreps, degeneracies, and displacement patterns for each mode.

### Axion active phonons

1. Run SMODES with primitive non-magnetic cell
- choose Gamma and Z modes
- TR odd modes will be at Z, since atoms with opposite spin have opposite displacement

2. Find axion active irreps 
- Odd under TR
- Odd under improper rotations
- Check literature
- Check character tables for irreps
- Focus on 1d irreps at first (easier to see) 
    - T2- is our irrep of interest

3. Compute the phonon modes with Quantum Espresso or Phonopy
- Give you list of modes (set of displacements of atoms in Cartesian coordinates) and frequency
- Tells you symmetry of each mode (eigenvectors of force constant matrix)
- Filter axion active irreps

4. Apply distortion to unit cell multiplying by small constant so that we are in linear response regime
- Forward and back (positive and negative amplitude)

### Steps
--------
1. 5e-3 -> 0.02 Angstrom amplitude
2. displace smodes output in cartesian
3. use angstrom cartesian coords in QE input file
4. Set tprnfor = :true:
5. Get total forces on each atom
6. Project forces onto mode

- Dynamical matrix is 4x4 matrix in basis of modes of T2-.
- General DM entry i,j is slope of F_j due to displacing atom i. Index i/j 
  is composite for atom i and direction x/y/z.
- Now in 4 dim space entry i,j mode di (basis mode of T2-)
  slope of the force on atom n comp m wrt di on each QE calc 
  dotted into mode dj -> 4x4 matrix
- Diagonalize, 4 eigenvectors give you linear comb of d's, these
  physical modes, eigenvalues are frequencies 
- Displace atoms again based on physical modes
- Wannier Hamiltonian
---------------------------------
- Can't use non-collinear and Hubbard U with tprnfor = .true..
- Use collinear with Hubbard U, then do non-collinear without U
- Forces should be similar with and without U


### Wannier windows

- Distentanglement: -7-20 eV
- Frozen: 0-10 eV

### Trial wavefunctions
Number of bands from -7-20 eV = 75 (without SOC)
Number of bands from 0-10 eV ~ 46 (without SOC)

- Mn1 d: 5 (l orbs) x 1 (atoms) x 2 (spins) = 10
- Mn2 d: 5 (l orbs) x 1 (atoms) x 2 (spins) = 10
- Bi p: 3 (l orbs) x 4 (atoms) x 2 (spins) = 24
- Te p: 3 (l orbs) x 8 (atoms) x 2 (spins) = 48
Total number of trial states = 92 (/2 = 46)

Remaining states:
- Mn1 s2: 1 (l orbs) x 1 (atoms) x 2 (spins) = 2
- Mn2 s2: 1 (l orbs) x 1 (atoms) x 2 (spins) = 2
- Bi s: 1 (l orbs) x 4 (atoms) x 2 (spins) = 8
- Te s: 1 (l orbs) x 8 (atoms) x 2 (spins) = 16


- The gap is smaller when SOC is included, as expected for MnBi2Te4. VdW + U leads to a larger gap.
### Notes from 2026-01-21

1. PT symmetry enforces Kramers degeneracy at each k-point, even with distortion. Quantum Espresso doesn't detect this symmetry and the bands are split at each k-point. In part due to this, the convergence is slow. We need to move origin to halfway between Mn atoms to enforce PT symmetry. New flag should be added in QE calculation.

2. Check if undistorted structure is Z2 odd. Should be, but need to confirm.

3. We need small smearing for intermediate metallic phases.


### TODO (early 2026)


- We want to plot the $C_2$ density for each symmetry mode to show that it is strongly peaked along the z-axis (111 direction in reduced coordinates). This is because of the small band gap and large Berry curvature along this direction. 
- We need to use consistent gauge as in YIO and figure out why the symmetry detection isn't working as expected.