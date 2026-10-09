# YIO symmetry floor: a QE PAW noncollinear-GGA bug that singles out Ir1 (2026-10-05)

Follows [`2026-10-05_yio_ws_tolerance_full_group.md`](2026-10-05_yio_ws_tolerance_full_group.md),
which showed that the curvature residual left at the physical WS tolerance (~1.5% total / ~26%
external, rotations) isn't WS tie-breaking. This note traces it down the pipeline to the DFT
eigenvalues. **Corrects the 09-29 note's "the floor moved into A(R)"**: A(R) isn't a separate source.

Base run throughout: `data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889` (SCF + NSCF, QE 7.6 + `new_ns_nc`
patch, **without** the PAW `v_rad`/`g_rad` fix), Wannierization `.../trial_02_Y_p_Ir_d_O_sp/192473bc889`.
Scripts: `calculations/diagnostics/2026-10-05_ws_tolerance_window/`.

## 1. The DFT eigenvalues break the symmetry (`eig_symmetry.py`)

Comparing `.eig` (NSCF, full 8³ mesh, `nosym = .true.`) at `k` and `g k` (`−g k` for primed `g`) over
the 48 ops of m-3m', the mesh maps to itself exactly. Bands 1-156 (occupied):

| class | n | max \|Δε\| (meV) |
|---|---|---|
| inversion | 1 | 0.87 |
| C3 / S6 about [111] | 4 | 1.7 |
| C2' / m' with axis ⟂ [111] (⟨110⟩ type) | 6 | 1.0-1.9 |
| all other 36 ops | 36 | 7.1-7.6 |

The 12 ops that roughly survive are exactly **−3m' about [111] through Ir1** (atom 5, at the origin),
i.e. Ir1's site group, and the 12 ops without a fractional translation. The other 36 are broken by
~7.5 meV. This doesn't come from diagonalization noise: it's in occupied bands, and it has a clean
group structure.

**Ruled out at this layer:**
- **SCF symmetry detection.** The SCF found all 48 ops ("48 Sym. Ops., with inversion, found (36 have
  fractional translation)").
- **The SCF state.** Its moments are exactly AIAO (all four Ir at 0.172951 μB along ⟨111⟩). The final
  Hubbard ns of the four Ir have identical traces (up/down swapped as the moment directions require)
  and identical eigenvalues.
- **FFT commensurability.** All fractional translations are 0 or ½ in crystal coordinates; the grids
  are 100³ (dense) and 64³ (smooth).
- **The b-vector shell.** Wannier90 used the 8 b-vectors (±1,±1,±1) · 0.0771 Å⁻¹, a set closed under Oh.

So the breaking enters at the **NSCF**, which uses the SCF's density and ns without symmetrizing
anything, and whose eigenvectors are what Wannier90 is built on.

## 2. It's Ir1 (`onsite_spectra.py`)

Spectra of the on-site H(R=0) blocks (Wannier `_hr.dat`) of symmetry-equivalent atoms. These are
gauge-free, since g maps one atom's Wannier set onto the other's:

| atoms | spread |
|---|---|
| Y ×4 | ≤ 0.1 meV |
| **Ir1 vs Ir2, Ir3, Ir4** | **11.6-11.9 meV**. Ir2, Ir3, Ir4 agree with each other to ~0.3 meV |
| O(8b) ×2 | 0.04 meV |
| O(48f) ×12 | two groups of 6 about 1.6 meV apart. This is the D3d split (the 6 O bonded to Ir1 vs the other 6) |

The Hamiltonian treats Ir1 differently from the three other Ir, which still look equivalent to each
other. That's exactly what reduces the symmetry to Ir1's site group.

## 3. It accounts for the whole curvature floor (`dft_vs_curvature.py`)

Per operation, DFT eigenvalue asymmetry vs the 5e-3 curvature residual (47 ops):

| part | corr. with DFT rms Δε |
|---|---|
| internal | +0.89 |
| cross | +0.92 |
| external | +0.93 |
| total | +0.93 |

Total curvature residual ≈ 2.9% per meV of rms Δε, with a relative spread of 17% across ops. Inversion,
at 0.098 meV rms, gives 0.18%; the broken ops, at ~0.55 meV rms, give ~1.5%. In absolute units the
internal, cross and external residuals are of the same size (e.g. z-line: internal 0.0048, external
0.0050). External only looked dominant as a *relative* residual because its scale (~0.024) is 25×
smaller than internal's (~0.59). **One source, the DFT Hamiltonian, drives all three parts.**

Older run for comparison (`eig_symmetry_any.py`): the July base (178382, QE 7.5, before the
`new_ns_nc` patch) breaks every operation by 3.6-6.7 meV, with no surviving subgroup and inversion at
1.2 meV. The `new_ns_nc` patch changed the *pattern* (everything now consistent with Ir1's site
symmetry) but not the size.

## 4. Cause: unguarded `segni_rad` in QE's PAW one-centre noncollinear GGA (found 2026-10-05)

**Mechanism.**
1. **`ux` is left at Ir1's moment.** `PW/src/compute_ux.f90` sets `ux` to the starting moment of the
   *first* magnetic atom. That's atom 5, Ir1, since the Y atoms 1-4 have none, so `ux ∝ (−1,−1,−1)`.
   It then sets `lsign = .false.` because the other moments aren't parallel to it, but **it never
   resets `ux`**.
2. **The smooth-grid GGA is guarded.** `compute_rho.f90:38` (`gradcorr`) uses `segni = SIGN(1, m·ux)`
   only `IF (lsign)`, and `segni = 1` otherwise.
3. **The PAW one-centre GGA is not.** `paw_onecenter.f90::compute_rho_spin_lm` (line 2340) and its
   linear-response twin `compute_drho_spin_lm` (2554) **always** use
   `segni_rad = SIGN(1, m·ux)`. In every PAW sphere, the up/down channels fed to the GGA are swapped
   wherever the local m has a negative projection on Ir1's axis. The GGA then differentiates up/down
   densities that jump across the surface m·ux = 0.

A global flip of `segni` is harmless, because the functional is ↑↔↓ symmetric. So only operations
that map the line of `ux` to itself survive: −3m′ about Ir1's [111], as observed.

**Tests (cluster, `base/u_3.0/test_nscf_sym_2026-10-05/`; outputs copied to
`data/Y2Ir2O7/base/soc/u_3.0/output/test_nscf_sym_2026-10-05/`; `nscf_test_analyse.py`,
`nscf_L_analyse.py`).** Each is an NSCF only, with 48 procs and `nbnd = 180`.

| run | setup | breaking |
|---|---|---|
| A | 192344 SCF save, **fixed** pw.x (`c41c409`), 60 k | 7.1-7.6 meV, identical to the original |
| B | 192344 SCF save, old pw.x (`1541022`), 60 k | identical to A and to the 64-proc original (0.08 meV) |
| C | 193145 SCF save (fixed SCF), fixed pw.x, 60 k | 4.6-4.8 meV, same Ir1-axis pattern |
| L0 | C at the 4 L points only | 4.6 / 4.8 / 4.8 meV (reproduces C) |
| L1 | L0 without the HUBBARD card | 5.1 / 5.6 / 5.2 |
| L2 | L0 with `atomic` projectors | 4.3 / 4.5 / 4.4 |
| **L3** | L1 with `input_dft = 'pz'` (LDA) | **0.0 / 0.0 / 0.0** |
| **L5** | L0 with a 0.01 starting moment on Y along +z, so `ux = +z` | **0.2 / 0.1 / 0.1** |

The L columns are max |Δε| (bands 1-156) of L(001), L(010) and L(100) against L(½,½,½), which lies on
Ir1's axis. In L5, `lsign` also became true because `is_parallel` uses an absolute threshold
(|a×b|² < 1e-6) and the seed moment was tiny. That isn't a clean single-variable change, but it shows
the breaking follows `ux` and not the physics.

**Ruled out, with the evidence:**
- **The PAW `v_rad`/`g_rad` bug.** A = B; that fix only reduces the size via the SCF, as C shows.
- **Heap or MPI layout.** B at 48 procs = the original at 64.
- **The Hubbard term and the ortho-atomic construction.** L1 and L2 are still broken.
- **The stored SCF inputs.** All exactly covariant under the 48 ops:
  - ρ and m(r) in G-space to 1e-15 (`density_symmetry.py`);
  - Hubbard ns including phases, to 3e-10, once QE's Condon-Shortley signs on xz/yz are used
    (`hubbard_ns_covariance.py`);
  - PAW becsum magnitudes to 1e-14 (`paw_becsum_magnitudes.py`).
- **The signed-zero `SIGN(1, −0.0)` route.** I suspected it first, but ifort 2021.11 at the build
  flags (`-O2`, default `-assume nominus0`) returns +1 for −0.0, scalar and vectorized. Only
  `-assume minus0` gives −1. That route would matter only if `ux` were 0, and it isn't.

**MnBi2Te4 is not affected the same way.** Its Mn starting moments are collinear (±z), so `lsign` is
true and both GGA paths use `segni` consistently, as designed. Bi and Te start at zero.

## 5. Fix, explanation and reruns: see the next day's note

The patch was validated on 2026-10-06 (worst eigenvalue asymmetry 4.83 meV → 46 μeV over all 47 operations), and the
production reruns were submitted that day. Both, plus a plain-language explanation of the bug with worked GGA examples
and why MnBi2Te4 is unaffected, are in
[`2026-10-06_segni_fix_validated_and_reruns.md`](2026-10-06_segni_fix_validated_and_reruns.md).
