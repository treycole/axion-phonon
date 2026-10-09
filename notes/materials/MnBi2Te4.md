# MnBi₂Te₄

θ = π (inversion + a composite antiunitary symmetry). Axion-active phonons: the zone-boundary **T₂⁻** modes.
Where the runs stand is in [`../README.md`](../README.md); this file is the material's setup reference.

## Current setup (2026-10-08)

- **Functional: LDA+U**, U = 4 eV on Mn 3d (`HUBBARD (ortho-atomic)`), `rel-pz` PAW, noncollinear + SOC.
  Never noncollinear PBE ([bug #18](../bugs/18_qe_noncollinear_gga_spurious_state.md)).
- **Base:** `MnBi2Te4/base/soc_lda/u_4.0/`, job 209653. Pseudopotentials and how the Te `rel-pz` file was generated:
  `soc_lda/pseudo_gen/` and [`../log/2026-10-08_mbt_noncollinear_gga_lda.md`](../log/2026-10-08_mbt_noncollinear_gga_lda.md) §11.
- **Every run:** the exact hex-CIF cell (a = 4.33360, c = 81.85200 Å, atoms at u·c), `startingwfc = 'atomic'`,
  `diago_thr_init = 1.0d-5`. Check layer mirroring and ∫|m⊥| with
  `calculations/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py`.
- **Wannier trial set:** Bi p, Te p, Mn1 d, Mn2 d — 92 spinor WFs (below). n_occ = 58 within the p-d manifold
  (bands 81-138 of the LDA base).
- Past MBT-specific bugs: [#9](../bugs/09_mbt_t2_rounded_cell_and_labels.md), [#14](../bugs/14_mbt_stale_binary.md),
  [#15](../bugs/15_mbt_u3_base_rounded_cell.md), [#18](../bugs/18_qe_noncollinear_gga_spurious_state.md).

---

## Setup notes (early 2026)

The working notes from the initial setup (SMODES, mode selection, Wannier windows). The procedure still holds; the
specific windows were for PBE without SOC. The files the text calls `SMODES.in.txt`/`SMODES.out.txt` are now `data/MnBi2Te4/phonon/SMODES.in`/`.out`.

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