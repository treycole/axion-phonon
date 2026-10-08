# #16 — QE PAW noncollinear GGA uses a stale `ux` without the `lsign` guard (root cause of the YIO symmetry floor)

> **Status: FIXED** (patch validated 2026-10-06). `paw_segni_lsign.patch` makes `segni_rad = SIGN(1, m·ux)` apply
> only `IF (lsign)`, and 1 otherwise.
> - **Build:** cluster `~/q-e-segnifix` (branch `fix/paw-segni-lsign`, `4b480ac68`), pw.x md5 **`e270b05b`**.
> - **Validated:** worst DFT eigenvalue asymmetry 4.83 meV → 46 μeV.
> - **Production YIO reruns in progress** (jobs 207546-207549). Not yet reported upstream.

**Found:** 2026-10-05. **Affects:** QE noncollinear magnetic + GGA + PAW calculations whose starting moments are
not all (anti)parallel. In this project: **every Y2Ir2O7 DFT run**, base and modes 1-3, with any binary.
MnBi2Te4 is **not** affected (see the end).

## Symptom
- **DFT eigenvalues** (`.eig`, NSCF on the full 8³ mesh) broke YIO's m-3m′ symmetry by up to 7.6 meV (4.8 meV
  after the #12 fix). Only 12 operations survived (to 1-1.8 meV): −3m′ about [111] through **Ir1**, the atom at
  the origin.
- **On-site H(0):** Ir1's spectrum was 12 meV off Ir2, Ir3 and Ir4, which agreed with each other to 0.3 meV. The 12
  O(48f) split into the 6 bonded to Ir1 and the other 6.
- **Curvature:** a symmetry floor of ~1.5% total (≈ 26% external, relative), tracking the DFT asymmetry operation
  by operation (corr ≈ 0.9, equally for internal, cross and external).
- **The SCF looked perfect:** all 48 operations found, all four Ir moments exactly 0.172951 μB along ⟨111⟩, and the
  stored ρ, m, ns and becsum exactly covariant.

## Root cause
GGA needs n↑ and n↓, which a noncollinear code builds from n and the vector m (explained in §"Background" below).
1. **`ux` isn't reset.** `PW/src/compute_ux.f90` sets the GGA reference axis `ux` to the starting moment of the
   **first magnetic atom**, Ir1, `∝ (−1,−1,−1)` (atoms 1-4, the Y atoms, have none). It then finds the other Ir
   moments not parallel and sets `lsign = .false.`, but **never resets `ux`**.
2. **The smooth-grid GGA is fine.** `PW/src/compute_rho.f90` uses `segni = SIGN(1, m·ux)` only `IF (lsign)`, and
   `segni = 1` (the local frame) otherwise.
3. **The PAW one-centre GGA is not.** `PW/src/paw_onecenter.f90::compute_rho_spin_lm` (~line 2340) and its
   linear-response twin `compute_drho_spin_lm` (~2554) have **no `lsign` guard**. They always use
   `segni_rad = SIGN(1, m·ux)`, i.e. n↑ = (n + s|m|)/2 with s = sign(m · Ir1-axis).

## How the code works, step by step (QE `934f8cd`, line numbers in `PW/src/`)

**1. Starting moments: `setup.f90:270-283`.** For noncollinear runs, each atom gets
`m_loc = starting_magnetization × (sin angle1 cos angle2, sin angle1 sin angle2, cos angle1)` from its species. YIO:
Y and O get 0; Ir1 gets 0.1·(−0.577, −0.577, −0.577), Ir2 0.1·(−0.577, 0.577, 0.577), and so on. Then
`ux = 0` and, for any GGA, `CALL compute_ux(m_loc, ux, nat)`.

**2. Choosing the GGA axis: `compute_ux.f90`.**
- **First magnetic atom:** sets `lsign = .FALSE.`, then loops over the atoms; the first with |m_loc|² > eps12 (atom 5,
  Ir1) gives `ux = m_loc(:, Ir1)` and `lsign = .TRUE.`.
- **Parallel test on the rest:** `lsign = lsign .AND. is_parallel(ux, m_loc(:, na))`. `is_parallel`
  (`divide_class.f90`) returns |a×b|² < 1e-6; for Ir2 this is 8.9e-5, so `lsign = .FALSE.`.
- **`ux` not cleared:** `IF (lsign)` normalizes `ux` and prints "Fixed quantization axis". There is no `ELSE`, so
  `ux` keeps Ir1's (unnormalized) moment. This is the first half of the bug.

**3. Smooth-grid GGA (correct): `v_of_rho.f90:553` → `gradcorr.f90:86` → `compute_rho.f90`.**
- **Densities:** `IF (lsign)` uses `segni = SIGN(1, m·ux)`, and `ELSE segni = 1`. Then
  `rhoout(↑,↓) = ½(n ± segni·|m|)`.
- **Potential:** `gradcorr` evaluates the gradient correction on these and maps the result back to a magnetic field
  `v(:,2:4) += segni·½(v↑ − v↓)·m/|m|` (`gradcorr.f90:233-235`).
- **Result:** for YIO, `segni = 1` everywhere, which is the local frame.

**4. PAW one-centre potential: `paw_onecenter.f90:61` `PAW_potential(becsum, ddd_paw, …)`.**
For each atom, it builds the sphere's density in spherical harmonics (`rho_lm`) from becsum, computes Hartree, then
calls `PAW_xc_potential` (line 199). That does the LDA part on a radial × angular grid with no `segni`, and for a
GGA calls `PAW_gcxc_potential` (lines 619 → 635):

1. **`compute_rho_spin_lm` (717 → 2267).** It maps `rho_lm` to the radial × angular grid (`PAW_lm2rad`). At every
   grid point, `mag = |m|`. If `mag < eps12`, `segni_rad = 1`; otherwise
   **`segni_rad = SIGN(1, m̂·ux)` (line 2340), with no `lsign` check.** This is the second half of the bug. Then
   `rhoout_rad(↑,↓) = ½(n ± segni_rad·mag)`, mapped back to a *truncated* spherical-harmonic expansion
   (`PAW_rad2lm`).
2. **Gradients (761-762).** `PAW_lm2rad` + `PAW_gradient` take the gradients of n↑ and n↓ from that expansion. A
   jump in n↑ becomes a large, smeared gradient.
3. **PBE (790).** `xc_gcx` evaluates the spin-polarized PBE correction on (n↑, n↓, ∇n↑, ∇n↓), giving the up/down
   potentials, plus a divergence term (`PAW_rad2lm3`).
4. **Back to a scalar plus a field: `compute_pot_nonc` (892 → 2369).**
   `v(1) = ½(v↑ + v↓)`, `vs = ½(v↑ − v↓)`, and **`v(2:4) = vs · segni_rad · m/|m|`** (line 2439).
   `segni_rad` enters twice. In a region of uniform `segni_rad = −1`, n↑ ↔ n↓ swap, so v↑ ↔ v↓ swap, vs → −vs,
   and the factor −1 restores the field direction. That region is fine. Only where `segni_rad` *changes sign*
   inside the sphere does step 2 produce a spurious gradient term.
5. **Back in `PAW_potential`:** each atom's potential becomes its nonlocal coefficients,
   `D_ij = ∫ v_Hxc p_ij − ∫ ṽ_Hxc (p̃_ij + aug_ij)`, stored as `ddd_paw(ij, atom, spin)`.

**5. Into the Hamiltonian.** `newd` adds `ddd_paw` to `deeq` (`electrons.f90:973`), and H|ψ⟩ includes
`Σ_ij |β_i⟩ deeq_ij ⟨β_j|ψ⟩` (`add_vuspsi`). Every sphere gets some spurious term, but the terms of
symmetry-equivalent atoms are no longer images of each other (see "Why it was a bug"). So `ddd_paw` of Ir1 differs
from that of Ir2-4, and the O atoms bonded to Ir1 differ from the rest.

**6. SCF (`electrons.f90:883-884`, `928-929`): the error is averaged, not removed.** Each iteration calls
`PAW_potential` **and then `PAW_symmetrize_ddd(ddd_paw)`**, which averages D_ij over the 48 operations. ρ and m are
symmetrized by `sym_rho`, becsum by `PAW_symmetrize` (`sum_band.f90:237`), and ns by `new_ns_nc`.
- **The SCF Hamiltonian is symmetric,** but it contains the group average of the spurious term.
- **The converged state is symmetric but shifted.** Every printed and stored quantity is exactly symmetric, which is
  why the bug was invisible there.

**7. NSCF (`run_pwscf` → `init_run.f90:173` → `potinit.f90:283`): the raw error is used.** `PAW_potential` runs once
on the stored becsum ("`PAW_pot … 1 calls`" in `nscf.out`), with **no `PAW_symmetrize_ddd`**; with `nosym` there
is only the identity anyway. The NSCF Hamiltonian therefore contains each atom's own spurious D_ij, unaveraged. It
is solved at all 512 k-points, giving eigenvalues and wavefunctions that break m-3m′ down to −3m′ about Ir1's axis.
Everything downstream (pw2wannier90, Wannier90, curvature, dθ) inherits it.

**What the patch changes.** In step 4.1 (and in `compute_drho_spin_lm`), `segni_rad = 1` unless `lsign`, the same
rule as step 3. No sphere has sign changes, so steps 4.2-4.5 produce no spurious D_ij, and step 7's Hamiltonian is
symmetric.

## Why it was a bug
- **Spurious jumps:** inside each PAW sphere, n↑ and n↓ are swapped wherever m points away from Ir1's axis. Because
  m varies in direction within a sphere (spin-orbit texture, induced O moments), s flips across the surface
  m ⊥ ux, and n↑, n↓ **jump** there. GGA differentiates them, so the jumps produce a spurious E_xc and a spurious
  potential D_ij.
- **Not covariant:** `ux` is one fixed direction that doesn't rotate with the crystal. Every sphere has flip
  surfaces, because m's direction varies inside it: on the Ir1 → Ir4 line, m turns more than 90° from Ir1's own axis
  within 0.6-0.9 Å of Ir1 while |m| is still 13-32% of its peak; see "What m(r) looks like" below. What breaks
  symmetry is *where* the flip plane sits relative to each atom's own texture:
  - **Ir1:** the plane is perpendicular to Ir1's own moment, so its flip pattern is symmetric under Ir1's site group.
  - **Ir2-4 and O:** the plane is still perpendicular to *Ir1's* axis, not their own, so their patterns are not the
    images of Ir1's under the operations that map Ir1 onto them.
  - **Result:** equivalent atoms get inequivalent spurious potentials. Ir2-4 stay equivalent to each other, because
    C3 about Ir1's axis maps them onto each other and leaves the plane fixed.
- **Which operations survive:** a global flip of s is harmless (the functional is ↑↔↓ symmetric). Only operations
  that map the *line* of `ux` to itself survive, which is −3m′ about [111].
- **Why it was missed:** the SCF symmetrizes ρ, m, ns, becsum *and* D_ij (`PAW_symmetrize_ddd`) every iteration,
  so its Hamiltonian and outputs are symmetric, though shifted by the averaged spurious term. Only the NSCF, which
  builds D_ij once and never symmetrizes it, exposes the bug.

## How it was isolated (cluster NSCF tests, `base/u_3.0/test_nscf_sym_2026-10-05`)

| test | setup | breaking |
|---|---|---|
| A | old SCF + #12-fixed pw.x, 60 k | 7.1-7.6 meV (= original) |
| B | old SCF + old pw.x, 48 procs | identical to A and to the original 64-proc run |
| C | #12-fixed SCF + fixed pw.x | 4.6-4.8 meV, same pattern |
| L1 / L2 | no Hubbard / atomic projectors | 5.1-5.6 / 4.3-4.5 meV |
| **L3** | **LDA instead of PBE** | **0.0 meV** |
| **L5** | `ux` moved to +z (tiny Y seed moment) | **0.1-0.2 meV** at L |

**Ruled out:**
- **#12:** A = B.
- **Heap / process count:** B at 48 procs = the original at 64.
- **The Hubbard term and the ortho-atomic construction:** L1, L2.
- **The stored SCF data:** ρ and m covariant to 1e-15, ns with phases to 3e-10, becsum to 1e-14.
- **Signed zero:** `SIGN(1, −0.0)` is +1 with ifort's default `nominus0` at the build flags.

## What m(r) looks like (`magnetization_profiles.py`)

m(r) evaluated from the stored G-space n, m (smooth + augmentation density, not all-electron; the PAW code builds
its own in-sphere density from becsum, so exact flip locations differ).

**YIO** (193145), Ir1 → bridging O → Ir4. |m| relative to Ir1's on-site value, and the angle to Ir1's axis:

| from Ir1 (Å) | \|m\| | angle to Ir1 axis | sign(m·ux) |
|---|---|---|---|
| 0 (Ir1) | 1.00 | 0° | + |
| 0.34 | 0.20 | 22° | + |
| 1.01 | 0.009 | 80° | + |
| 1.34 | 0.026 | 162° | − |
| 2.01 (O) | 0.048 | 55° (along −z) | + |
| 3.02 | 0.009 | 171° | − |
| 4.02 (Ir4) | 1.00 | 109.5° | − |

The size falls to a few percent between the atoms, while the direction rotates continuously. The O moment points
along −z, the sum of its two Ir neighbours' moments. The sign relative to `ux` flips where |m| is finite, which is
the jump case.

**MnBi2Te4** (202833, iteration 115, not converged), Mn1 → Mn2 along the 3-fold axis: m_z = +1.21 at Mn1, drops
~200× within 1.4 Å, small induced moments of either sign on Te and Bi (up to ~0.1 near Bi), −1.21 at Mn2.
m_⊥ = 0 (1e-16) on this line by symmetry. Sign changes happen through m_z = 0, so s·|m| = m_z is smooth: the
collinear case.

**During the calculation:** starting moments only seed the SCF. Ir1's moment went from 0.43 μB (iteration 1) to
0.30 μB (converged), with the direction held on ⟨111⟩ by symmetry. `ux` is computed once from the input starting
moments and never follows m. The NSCF keeps the SCF's m(r) fixed.

## Fix
`calculations/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/paw_segni_lsign.patch`. In both routines,
`use_ux = lsign` (passed to the OpenACC loop as `firstprivate`), and `IF (use_ux) segni_rad = SIGN(1, m·ux)`,
`ELSE segni_rad = 1`, mirroring `compute_rho.f90`.

Not sufficient on its own: clearing `ux` in `compute_ux` when `lsign` is false. With a compiler that honours
signed zero (gfortran does by default), `SIGN(1, m·0)` is −1 in the (−,−,−) octant.

## Verification
Test C2 (job 205711) is identical to C except for the patched binary. Full-precision XML eigenvalues, bands
1-156, all 47 non-identity operations:

| run | max | median | Ir1 site group | other 36 ops |
|---|---|---|---|---|
| C (unpatched) | 4.83 meV | 4.80 meV | 0.74-1.5 meV | 4.6-4.8 meV |
| **C2 (patched)** | **0.046 meV** | 0.045 meV | 19-46 μeV | 0.08-46 μeV |

That's 100× smaller, and the Ir1 pattern is gone. The ≤ 46 μeV left is a numerical floor of unverified origin.

## Background: GGA, the local frame and the fixed axis
- **GGA** (PBE): `E_xc = ∫ f(n↑, n↓, ∇n↑, ∇n↓)`. It needs **gradients** of n↑ and n↓; LDA (`∫ f(n↑, n↓)`)
  doesn't.
- **Local frame (QE's intended rule for noncollinear orders):** n↑ = (n + |m|)/2. Example: n = 1.0 and
  m = (0.3, 0, 0.4) give |m| = 0.5, n↑ = 0.75, n↓ = 0.25, for *any* direction of m.
  - *Why only |m| although we use SOC:* SOC is in the Hamiltonian (the fully relativistic pseudopotentials), not in
    the XC functional. PBE is built for spin only and is invariant under rotating all spins together, so it can
    depend only on n and |m|, the eigenvalues of ½(n + m·σ).
  - *What it claims:* this is the standard (Kübler) noncollinear generalization. It is exact for LDA and
    rotation-covariant, but an approximation for GGA: it ignores gradients of the direction m̂.
- **Fixed axis (`lsign = .true.`, for collinear *orders*, i.e. all moments ∥ or anti-∥ one axis like MBT's Mn ±z,
  even when run with `noncolin = .true.`):** n↑ = (n + s|m|)/2, s = sign(m·ux). It avoids the
  kink of |m| where m passes through zero. In a 1D antiferromagnet with m_z = sin x, |m| = |sin x| is V-shaped at
  x = 0, while s|m| = sin x is smooth.
- **Wrong axis (this bug):** take m(θ) = 0.4 (cos θ, 0, sin θ), rotating smoothly with constant length, and
  `ux = x̂`.
  - s = sign(cos θ) flips at θ = 90°, and n↑ jumps from 0.7 to 0.3 although nothing physical happens.
  - In the local frame, n↑ = 0.7 everywhere.
  - A *uniform* swap would be harmless; only the sign *changes* hurt.
  - LDA has no gradients, so it doesn't care (test L3).

**Kink vs jump.** `SIGN(1, x)` is Fortran for ±1 with the sign of x. The local-frame n↑ = ½(n + |m|) (the
eigenvalues of ½(n + m·σ)) has a **kink** where m passes through zero. n↑ stays continuous and only its slope
changes, so the gradient PBE sees stays finite. That's a known, mild, symmetric weakness of the local-frame GGA. The
fixed axis removes that kink for a *collinear* order, because there m·ux changes sign only through m = 0, so s|m| =
m·ux is smooth. For a *noncollinear* order, m·ux changes sign where |m| is finite, so n↑ **jumps** by |m|, and its
gradient spikes. That's the bug.

| rule | where n↑ misbehaves | how | breaks symmetry? |
|---|---|---|---|
| local frame ½(n ± \|m\|) | where \|m\| = 0 | kink (finite gradient) | no |
| fixed axis, collinear order (MBT) | nowhere | none | no |
| fixed axis, noncollinear order (YIO bug) | where m ⊥ ux, \|m\| finite | jump (gradient spike) | yes |

The patch removes the jumps and leaves the standard kink; YIO's GGA is not made exact.

**What the sign means, and why not each atom's own axis.**
- **The sign:** in the fixed-axis mode, "up" means spin along +ux everywhere. s = +1 where the local majority spin is
  on the +ux side (m·ux > 0) and −1 where it's on the −ux side. For a collinear order, s|m| = m·ûx exactly, a signed
  projection that passes smoothly through zero.
- **The local frame uses no axis at all:** at each point, n↑,↓ are the eigenvalues of ½(n + m·σ), with "up" along
  m(r) at that point, in spheres and between atoms alike.
- **A per-atom axis is not the correct alternative.** Each atom's moment as the axis for its own sphere would be
  symmetric, since each axis rotates with its atom, but it still creates jumps wherever m turns more than 90° from
  the atom's axis. That happens within 0.6-0.9 Å of Ir1 with |m| still 13-32% of the peak. It also can't be defined
  between atoms. QE has only the two modes, one global axis or none, and the bug is the PAW part using one global
  axis (Ir1's) when QE had chosen none.

| scheme | where n↑ misbehaves | symmetric? | used by QE when |
|---|---|---|---|
| local frame (s = +1) | kink only where \|m\| = 0 | yes | `lsign = .false.`: smooth grid; spheres after the patch |
| one global axis | nowhere if m ∥ ±ux | yes, if symmetry preserves the axis | `lsign = .true.` (MBT) |
| global axis on a noncollinear order (**bug**) | jumps where m ⊥ ux, \|m\| finite | **no** | unpatched PAW spheres, `lsign = .false.` |
| each atom's own axis (hypothetical) | jumps where m turns > 90° from the atom's axis | yes | never |

## When it triggers
The trigger is a **noncollinear magnetic order**. QE correctly picks the local frame (`lsign = .false.`), the
smooth-grid GGA follows it, and the PAW GGA applies the collinear-only fixed-axis trick anyway with a leftover axis.
- **All of these:** noncollinear and magnetic, a GGA functional, PAW, and starting moments **not all parallel or
  antiparallel**: all-in-all-out pyrochlores, canted and spiral magnets.
- **Which axis breaks:** it's set by whichever magnetic atom is first in the input.
- **The `is_parallel` threshold:** `is_parallel` (`PW/src/divide_class.f90`) uses an absolute threshold,
  |a×b|² < 1e-6, so tiny starting moments can be judged parallel.

## MnBi2Te4 is not affected
- **lsign is true:** Mn1/Mn2 start at `angle1` = 0° / 180° (antiparallel along z), so `lsign = .true.` and
  `ux = z`. Every MBT spin-orbit SCF prints `Fixed quantization axis for GGA: 0 0 1`.
- **Consistent:** both GGA paths use s = sign(m_z), as designed.
- **Symmetric:** every MBT magnetic operation maps z to ±z, so s stays put or flips globally, which is harmless.
- **No-SOC MBT runs** are collinear (`nspin = 2`).
- **Future starts:** a canted MBT starting configuration would hit this bug. The patched build is safe either way.

## Remaining
1. **Finish the production reruns** (base 207546, modes 207547-207549) and the Wannier stage, then rerun the
   `.eig` and curvature symmetry checks and dθ. Every earlier YIO number is superseded.
2. **Replace `~/q-e` or use `~/q-e-segnifix` as the production QE.** Add the fix to CLAUDE.md's "correct versions"
   list and update pitfall 7.
3. **Report to QEF:** the code is unchanged on develop as merged 2026-09-23.
4. **Delete the old buggy YIO runs** after validation (~3 TB), with the user's OK.

**Sources:** [`../progress/2026-10-05_yio_symmetry_floor_dft_origin.md`](../progress/2026-10-05_yio_symmetry_floor_dft_origin.md),
[`../progress/2026-10-06_segni_fix_validated_and_reruns.md`](../progress/2026-10-06_segni_fix_validated_and_reruns.md);
scripts in `calculations/diagnostics/2026-10-05_ws_tolerance_window/`.
