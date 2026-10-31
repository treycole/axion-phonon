# Y₂Ir₂O₇

θ = 0 (all-in-all-out order, correlation-driven): the non-topological counterpoint to MnBi₂Te₄. Axion-active
phonon: the zone-centre **Γ₂⁻ (A₂ᵤ)** mode, IR- and Raman-silent. Where the runs stand is in
[`../README.md`](../README.md); this file is the material's setup reference.

## Current setup (2026-10-12)

- **Functional: LDA+U** (since 2026-10-12, user), U = 3 eV on Ir 5d, `rel-pz` PAW (Y `spn` 1.0.0, Ir `n` 0.2.3, O `n`
  0.1), noncollinear + SOC. Why: curvature symmetry exact (PBE 0.013%, [bug #22](../bugs/22_hopping_cut_and_hr_rounding.md)),
  band energies ~100× more symmetric, SCFs converge (PBE mode 2 stalled like [#17](../bugs/17_segni_fix_base_scf_stall.md)).
- **Binary:** pw.x `db4b99b4` (cluster `~/q-e-kfix`) for new runs. Runs on the segni-fixed `e270b05b` stay valid (same k
  lists). Every YIO run before the segni fix is superseded ([bug #16](../bugs/16_qe_paw_segni_stale_ux.md); earlier
  [#7](../bugs/07_qe_new_ns_nc_symmetrizer.md)).
- **Runs** (local `data/Y2Ir2O7/`, reorganized 2026-10-13; cluster paths unchanged under `Y2Ir2O7/phonon/`):

  | | local | cluster |
  |---|---|---|
  | LDA base | `base/lda/{scf_210527, nscf_210590, w90_trial02_212615}` | `base/u_3.0/lda_test_2026-10-08/YL2_lda_relpz/` |
  | LDA modes 1-3 | `phonon/GM2-/mode{1,2,3}/Q_0.01A/lda/{scf_nscf_21275{5,6,7}, w90_trial02_{212982,212994,212995}}` | `mode_{1,2,3}/fromYL2_2026-10-12/` |
  | LDA amplitude series (10-17) | not copied yet | `mode_{1,2,3}/Q_0p0{1,2,3,5}_conv1e-10_2026-10-17/`, base `.../YL2_lda_relpz/conv1e-10_2026-10-17/` ([log 10-17](../log/2026-10-17_yio_amplitude_series.md)) |
  | PBE base (comparison) | `base/pbe/{scf_210525, nscf_210589, w90_trial02_212616}` | `.../YP_pbe/` (W90 in `w90_trial02_base_frozfix_2026-10-10/`) |
  | PBE modes 1, 3 | `phonon/GM2-/mode{1,3}/Q_0.01A/pbe/{scf_nscf_21059{6,8}, w90_trial02_{212996,213002}}` | `mode_{1,3}/fromYP_2026-10-08/` (W90 `..._frozfix`) |

  Each `w90_*` folder has the Wannier90 output, the dθ sweeps (`nk_sweep_*`, modes) and `job_scripts/` (the
  cluster `.win`, `pw2wan.in`, job scripts and NOTE). The SCF/NSCF folders hold text only (inputs, outputs, NOTE).
  PBE mode 2 never converged (212707, [#17](../bugs/17_segni_fix_base_scf_stall.md)). Superseded runs:
  `data/Y2Ir2O7/archive/` (its README lists each one and why). Templates and the cluster workflow scripts:
  `archive/setup_templates/Y2Ir2O7/`.
- **Wannier windows:** every `dis_*` edge in an eigenvalue gap of that NSCF (tens of meV); 15.0 eV sat 0.16 meV from a
  band and broke the threefold operations ([log 10-08 §12.8](../log/2026-10-08_mbt_noncollinear_gga_lda.md)).
- **WS tie tolerance 5e-3 Å** (default; with base Wannier centres YIO's near-ties end at 1.7e-3, genuine differences start at 2.1e-2). **Atom-pinned positions since 2026-10-26:** ties exact (≤ 1e-12 Å), next distinct distance 0.39 Å, so the tolerance no longer matters (`archive/diagnostics/2026-10-26_yio_pinned/tie_window_pinned.log`).
- **Structures:** base `Y2Ir2O7/phonon/base/u_3.0/` and the three Γ₂⁻ mode structures `mode_1..3`
  (SMODES files: `data/Y2Ir2O7/phonon/SMODES.Y2Ir2O7.*`).
- **Magnetic point groups:** base m-3m′, displaced −4′3m′.
- **Wannier set: "Option A"** — `Y: p`, `Ir: d`, `O: s;p`, 176 WFs, projection only, **n_occ = 156**. The diffuse
  empty shells (Y 4d, Ir 6s/6p) close the gap and must stay out ([bug #1](../bugs/01_wannier_diffuse_projections_gap_collapse.md)).
- **Mode `.amn` from the base Wannier functions, unshifted:** `generate_amn.jl --no-displacement`, on the login
  node (the rigid-shift variant, [`../theory/rigid_shift.md`](../theory/rigid_shift.md), is no longer used).
- Spreads of every trial set tried: [`Y2Ir2O7_trial_orbital_spreads.txt`](Y2Ir2O7_trial_orbital_spreads.txt).

---

## Structure

Y₂Ir₂O₇ crystallizes in the cubic Fd̅3m space group. Y³⁺ is bonded in a body-centered cubic geometry to eight O²⁻ atoms. There are two shorter (2.21 Å) and six longer (2.44 Å) Y-O bond lengths. Ir⁴⁺ is bonded to six equivalent O²⁻ atoms to form corner-sharing IrO₆ octahedra. The corner-sharing octahedral tilt angles are 53°. All Ir-O bond lengths are 2.01 Å. There are two inequivalent O²⁻ sites. In the first O²⁻ site, O²⁻ is bonded to two equivalent Y³⁺ and two equivalent Ir⁴⁺ atoms to form a mixture of distorted corner and edge-sharing OY₂Ir₂ tetrahedra. In the second O²⁻ site, O²⁻ is bonded to four equivalent Y³⁺ atoms to form a mixture of corner and edge-sharing OY₄ tetrahedra.

- Space group: Fd-3m (227)
- Ferromagnetic (Materials Project's label; the order used here is all-in-all-out)
- Total magnetization: 2.04 µB/f.u.

- Conventional cell
a = 10.19 Å
b = 10.19 Å
c = 10.19 Å
α = 90.00 º
β = 90.00 º
ɣ = 90.00 º
Volume = 1056.73 Å³

- Primitive cell
a = 7.20 Å
b = 7.20 Å
c = 7.20 Å
α = 60.00 º
β = 60.00 º
ɣ = 60.00 º
Volume = 264.18 Å³

- Atomic positions
Wyckoff	Element	x	y	z
8b	O	1/4	3/4	1/4
16c	Ir	5/8	7/8	3/8
16d	Y	1/8	7/8	3/8
48f	O	0	0	0.212932

## k-path

```text
0.000000  0.000000  0.000000  30   ! G
0.500000  0.000000  0.500000  30   ! X
0.250000  0.500000  0.750000  30   ! W
0.375000  0.375000  0.750000  30   ! K
0.000000  0.000000  0.000000  30   ! G
0.500000  0.500000  0.500000  30   ! L
0.250000  0.625000  0.625000  30   ! U
0.250000  0.500000  0.750000  30   ! W
0.500000  0.500000  0.500000  30   ! L
0.375000  0.375000  0.750000   1   ! K
```

## Wannierization notes (early 2026)

> Superseded in part: the 224-WF set below (with Ir s/p and Y d) is the one that closed the gap (bug #1).
> Production is Option A (above). The band-character notes still hold.


### Atoms in unit cell

There are 

- 4 Y atoms
- 4 Ir atoms
- 14 O atoms

### Character of bands w/o SOC

- Y s: form deeper valence states around -30 eV
- Y p: flat manifold around -10 eV
- Y d: predominantly conduction states above E_F
- O p: predominantly valence states around [0, 15] eV
- O s: predominantly valence states around [-12, -5] eV
- Ir s, p, d: span valence and conduction from [-12, 15] eV

### Trial wavefunctions

We want all dispersive occupied bands and some conduction bands to be captured by the trial wavefunctions. The following are the trial wavefunctions we will use:

- O p: large DOS in primary valence bands from [0, E_F] eV
- O s: large DOS in semi-core valence states from [-12, 0] eV
- Ir p: span valence (semi-core and primary) and conduction states. Important for effective spin-1/2 around Fermi level.
- Ir s: similar DOS spread as Ir p.
- Ir d: largest DOS that spans primary valence and conduction bands.
- Y d: predominantly conduction states above E_F. 

Total number of trial wavefunctions when spin-orbit coupling is included: 

- O p = 14 atoms * 3 orbitals * 2 spin = 84
- O s = 14 atoms * 1 orbital * 2 spin = 28
- Ir p = 4 atoms * 3 orbitals * 2 spin = 24
- Ir s = 4 atoms * 1 orbital * 2 spin = 8
- Ir d = 4 atoms * 5 orbitals * 2 spin = 40
- Y d = 4 atoms * 5 orbitals * 2 spin = 40

Total = 224 trial wavefunctions. 

This is a large number of trial wavefunctions, but we want to capture all the relevant physics in the valence and conduction bands.

### Disentanglement and frozen windows

We want enough conduction bands (`nbnds`) in the nscf calculation so that the frozen window can be large enough to include most of the conduction bands, minus a few Wannier bands, and the remaining disentanglement window is large enough so that the iterative subspace selection can find localized Wannier functions. 

We have been using `nbnds = 300` but the remaining disentanglement window is somewhat small. We are now increasing to `nbnds = 350` so that the disentanglement window is larger.

### Frozen window

The frozen window is chosen to include all occupied bands and most of the conduction bands so that most of the Wannier bands exactly match the DFT bands.

Previously, we chose the inner window to be from [0, 23] eV, which makes most of the Wannier bands exactly match the DFT bands. The disentanglement window was chosen to be [0, 28] eV. This is small enough so that the unconverged upper most conduction bands are excluded. However, now that we are including Os trial orbitals, the windows need to be adjusted to incorporate the semi-core valence states. 

Looking at the the bands of the undistorted structure (with SOC), the semi-core valence states are within [-14, 0] eV and there are 52 of them. If we extend the bottom of the inner window to there, and leave the upper end of the inner window as is, then there are at most 244 bands across the BZ which is too many. Changing the upper end of frozen window to be 19 eV means there are at most 216 states between [-14, 19] eV across the BZ, which is less than the number of trial wavefunctions (224). This is a reasonable choice for the frozen window.

**Frozen window** = [-14, 19] eV

### Disentanglement window

The remaining states above the frozen window are included in the disentanglement window. The disentanglement window is chosen to be so that the unconverged upper most conduction bands are excluded. The upper end of the disentanglement window was previously chosen to be 28 eV. Now that we are including more bands, we can increase the upper end of the disentanglement window to 30 eV for a first try. This is small enough so that the unconverged upper most conduction bands are excluded.

**Disentanglement window** = [-14, 30] eV




## Early roadblocks

- The Hubbard U value is chosen to be 3.0 eV, which is a reasonable value for Ir 5d orbitals. We want it large enough so that the gap is sizable, and so that we are in the trivial insulator phase. 

- There are several local minima in the Wannierization procedure that is due to the DFT+U procedure. Different combinations of Hubbard occupations lead to different ground state energies. Sometimes, depending on the starting magnetization, the scf will converge to a different AIAO configuration with quenched moments. We need to make sure that the undistorted and distorted scf ground states are in the same AIAO configuration and in the same local minimum branch. We need small step sizes as well for this same reason.
