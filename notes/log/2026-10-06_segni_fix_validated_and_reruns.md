# The PAW `segni` fix, validated; how the bug works; production reruns (2026-10-06)

Follows [`2026-10-05_yio_symmetry_floor_dft_origin.md`](2026-10-05_yio_symmetry_floor_dft_origin.md), which found
bug #16 (numbering of [`../bugs/README.md`](../bugs/README.md)): QE's PAW one-centre noncollinear GGA uses a stale `ux` without
the `lsign` guard. This note covers four things:
- the fix and its validation (§1);
- a plain-language explanation of the bug, with worked examples (§2);
- why MnBi2Te4 is not affected (§3);
- the production reruns started today (§4).

## 1. Fix, validated

**Patch:** `paw_segni_lsign.patch`. In `compute_rho_spin_lm` and `compute_drho_spin_lm` (`PW/src/paw_onecenter.f90`),
`use_ux = lsign` is passed to the OpenACC loop as `firstprivate`, and `segni_rad = SIGN(1, m·ux)` only
`IF (use_ux)`, 1 otherwise, mirroring `compute_rho.f90`.
- **Built** on the cluster in **`~/q-e-segnifix`**: a copy of `~/q-e` 934f8cd with `make.inc` TOPDIR/BUILDDIR
  repointed, committed as `4b480ac68` on branch `fix/paw-segni-lsign`.
- **Binary:** `pw.x` md5 **`e270b05b9424b3dca8145a4b84f26be4`**. The production `~/q-e` (`c41c409`) is
  unchanged.
- **Files:** the patch and the build/submit script are in
  `calculations/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/`.

**Validation: run C2** (job 205711). It's identical to test C (193145 fixed SCF save, the same 60 k,
byte-identical `nscf.in`) except for the patched binary. Full-precision eigenvalues come from each run's
`data-file-schema.xml` (`nscf_test_analyse_xml.py`), bands 1-156, max |Δε| over the 47 non-identity ops:

| run | max | median | Ir1 site group (C3, C2′⟨110⟩, I) | other 36 ops |
|---|---|---|---|---|
| C (unpatched) | 4.83 meV | 4.80 meV | 0.74-1.5 meV | 4.6-4.8 meV |
| **C2 (patched)** | **0.046 meV** | 0.045 meV | 19-46 μeV | 0.08-46 μeV |

That's a ~100× reduction, and the Ir1 pattern is gone: C3 about [111] is now no better than the other
operations. The remaining ≤ 46 μeV is a numerical floor of unverified origin (NSCF diagonalization tolerance or
the PAW angular quadrature are the likely candidates). With the ~2.9%-per-meV(rms) scaling of the 10-05 note §3,
it should put the curvature residual well below 0.1%, though that hasn't been measured.

**Consequences:**
- **Every YIO DFT run (base, modes 1-3) must be redone with the patched pw.x, SCF included.** The SCF's ρ is
  symmetrized, but the state it converged to was driven by the asymmetric PAW GGA potential.
- **The Wannier, WS and A(R) layers were never the floor.**
- **dθ will move.** All earlier YIO numbers predate the fix.
- **The bug is probably upstream.** The code is unchanged on QEF develop as merged 2026-09-23. It should be
  reported to QEF; the user decides.
- **`CLAUDE.md`**'s "correct versions" list needs the segni fix added for YIO.
- **Cluster cleanup:** `test_nscf_sym_2026-10-05/` went from 39 GB to 658 MB. `nscf.out` for all runs, plus the
  XMLs for A, B, C and C2, are in `data/Y2Ir2O7/base/soc/u_3.0/output/test_nscf_sym_2026-10-05/`.

## 2. How the bug works

### 2.1 What GGA is and what it needs

DFT's one unknown piece is the exchange-correlation energy E_xc. Every practical functional writes it as a sum
over points r of an energy density that depends on the electrons near that point:

- **LDA** (local density): `E_xc = ∫ f(n(r)) dr`. It uses only the density at that point.
- **LSDA** (local spin density): `∫ f(n↑(r), n↓(r)) dr`. It uses the up and down densities separately.
- **GGA** (generalized gradient, e.g. PBE, used here): `∫ f(n↑, n↓, ∇n↑, ∇n↓) dr`. It also uses how fast
  those densities change in space.

The key point: **GGA needs spatial derivatives of n↑ and n↓**, so it compares each point with its neighbours.
LDA doesn't. These functionals were built for *collinear* spin and take two numbers, n↑ and n↓.

### 2.2 Noncollinear magnetism: where n↑ and n↓ come from (local frame)

A noncollinear run stores the charge n(r) and the magnetization **m(r)**, a 3D vector that can point anywhere
at each point. To feed the collinear functional, QE takes "up" along the local m:

    n↑ = (n + |m|)/2,   n↓ = (n − |m|)/2

**Example:** n = 1.0 and m = (0.3, 0, 0.4) give |m| = 0.5, so n↑ = 0.75 and n↓ = 0.25. For
m = (0.4, −0.3, 0), |m| is still 0.5, with the same n↑ and n↓. **The direction of m doesn't enter, only its
length.**

*Why only |m|, even though we use SOC:* SOC is in the Hamiltonian, through the fully relativistic pseudopotentials
(V_ion) and PAW's small-component terms. That's what locks YIO's moments to ⟨111⟩ and twists m(r) inside each
atom. The XC functional itself (PBE) contains no SOC: it's built for spin only and is unchanged if all spins rotate
together. So E_xc can depend only on frame-independent quantities, n and |m|, the two eigenvalues of the local
spin-density matrix ½(n + m·σ).

*What "local frame" does and doesn't claim:* this construction (Kübler et al.; used by QE, VASP and others) is the
standard noncollinear generalization.
- **Exact for LDA:** the energy density depends only on the local eigenvalues.
- **Rotation-covariant:** so it respects the crystal symmetry.
- **An approximation for GGA:** it uses gradients of n↑ and n↓ (the magnitudes) but ignores the energy of the
  *direction* m̂ varying in space. More complete noncollinear GGAs exist (Scalmani-Frisch, Peralta et al.) but
  aren't in QE.

It's what QE intends for a noncollinear order like YIO's (`lsign = .false.`).

### 2.3 The alternative: a fixed axis, and why QE has it

|m| is never negative, so it has a kink wherever m passes through zero. **1D example:** an antiferromagnet with
m_z(x) = sin x, the moment flipping from +z on one atom to −z on the next:

    x:              −π/2    0    +π/2
    m_z = sin x:    −1      0    +1      <- smooth
    |m| = |sin x|:   1      0     1      <- V-shaped kink at x = 0

GGA would see an artificial slope jump between the atoms.

"Collinear" here means the *magnetic order*: all moments parallel or antiparallel to one axis, like MnBi2Te4's
antiferromagnet (Mn ±z). It does not mean the calculation mode. MBT is run with `noncolin = .true.` (two-component
spinors, vector m), because SOC requires it, but its order is collinear. YIO's all-in-all-out order is
noncollinear: four moments along four different ⟨111⟩ axes. For collinear orders, QE chooses one global axis
`ux` and attaches a sign:

    n↑ = (n + s·|m|)/2,   s = sign(m · ux)

Here s·|m| = sin x, which is smooth again. In effect "up" means +ux everywhere, as in a collinear calculation.
QE uses this **only when all starting moments are parallel or antiparallel** (`lsign = .true.`), and then
prints `Fixed quantization axis for GGA`.

### 2.4 What a wrong axis does: a 2D example

Take a moment that rotates smoothly with constant length, as happens inside an atom with spin-orbit:
m(θ) = 0.4 (cos θ, 0, sin θ), with θ = 0° → 180° across a region and n = 1.0.

- **Local frame (correct):** |m| = 0.4, so n↑ = 0.7 everywhere. Smooth, no spurious gradient.
- **Fixed axis with a bad choice, `ux = x̂`:** s = sign(cos θ).

      θ:    0°    60°   89° | 91°   120°  180°
      s:    +1    +1    +1  | −1    −1    −1
      n↑:   0.7   0.7   0.7 | 0.3   0.3   0.3

  n↑ **jumps from 0.7 to 0.3** at θ = 90°, the surface where m ⊥ `ux`, although nothing physical happens
  there. GGA differentiates n↑, and a jump has an enormous gradient. In a PAW sphere, the finite
  spherical-harmonic expansion smears it into a ripple, which still gives a spurious energy and potential.

A *uniform* swap (s = −1 throughout a region) is harmless, because GGA is invariant under relabelling all of ↑
as ↓. **Only the places where s changes sign do damage.** LDA uses no gradients, so even the jumps are
invisible to it.

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

### 2.5 The bug in those terms

The trigger is a **noncollinear order**. QE correctly decides there's no common axis (`lsign = .false.`, local
frame), and the smooth-grid GGA follows that decision. The PAW part applies the collinear-only fixed-axis trick
anyway, with a leftover axis that means nothing for this order. For a collinear order (MBT), `lsign = .true.`, using
`ux` is intended, and there's no bug.

1. **`ux` isn't reset.** `compute_ux` (`PW/src/compute_ux.f90`) sets `ux` to the starting moment of the **first
   magnetic atom**. In YIO that's Ir1 (atoms 1-4 are Y with no starting moment), so `ux ∝ (−1,−1,−1)`. It then
   finds Ir2 not parallel and sets `lsign = .false.`, but **never resets `ux`**. It only normalizes and prints
   it when `lsign` is true.
2. **The smooth-grid GGA is fine.** `compute_rho.f90` uses `segni = SIGN(1, m·ux)` only `IF (lsign)`, and
   `segni = 1` otherwise: the local frame.
3. **The PAW one-centre GGA is not.** `paw_onecenter.f90::compute_rho_spin_lm` (line ~2340) and its
   linear-response twin `compute_drho_spin_lm` (~2554, used by phonon/DFPT) have **no `lsign` guard**. They
   always set `segni_rad = SIGN(1, m·ux)`, which is §2.4's wrong-axis case, with the axis fixed to Ir1's moment.

**Why that singles out Ir1.** `ux` is one fixed direction in space, the same for every atom. It doesn't rotate
with the crystal, so symmetry-equivalent atoms sit at different angles to it:
- **Ir1:** the flip plane m ⊥ ux is perpendicular to Ir1's *own* moment. Ir1's sphere still has flip surfaces (m
  turns more than 90° from Ir1's axis within 0.6-0.9 Å while |m| is 13-32% of its peak), but they're arranged
  symmetrically under Ir1's site group.
- **Ir2-4:** their moments sit at 109.5° to `ux` (cos = −1/3), so the same plane cuts their textures at an angle
  unrelated to their own axes.
- **O atoms:** small induced moments of varying direction (along −z at the Ir1-O-Ir4 bridge) cross that plane too.
- *(Corrected 2026-10-06: an earlier version said Ir1's sphere has no flips; the m(r) profiles in bug #16 show it
  does.)*

Equivalent atoms get different spurious GGA potentials. Under a symmetry operation g, m rotates, so the result
stays covariant only if g maps the *line* of `ux` to itself (a global flip of s is harmless). Those operations
are exactly −3m′ about [111]: C3, the primed C2′ ⟂ [111], and inversion. Everything else breaks.

**Effects in YIO:**
- DFT eigenvalues break m-3m′ by 4.8-7.6 meV.
- Ir1's on-site levels are 12 meV off Ir2-4.
- The 6 O(48f) bonded to Ir1 split from the other 6.
- The curvature floor is ~1.5% total (≈ 26% external), because the Wannier pipeline faithfully reproduces an
  asymmetric Hamiltonian.

**Why it hid:**
- **The SCF covers it up.** Every iteration it symmetrizes ρ, m, ns, becsum *and* the PAW D_ij
  (`PAW_symmetrize_ddd`, `electrons.f90:884, 929`). So the SCF Hamiltonian and every stored quantity are exactly
  symmetric, though the state is shifted by the group-averaged spurious term.
- **Only the NSCF shows it.** The NSCF builds D_ij once (`potinit`) and never symmetrizes it. With `nosym` it applies
  that Hamiltonian at all k. The line-by-line code walk is in `../bugs/16_qe_paw_segni_stale_ux.md`.
- **The earlier clearance couldn't catch it.** "The PAW fix doesn't matter for YIO" compared symmetrized SCF
  quantities.

**How it was pinned down** (10-05 note §4):
- **Not these:** the PAW `v_rad` fix, the process count, the Hubbard term and the ortho-atomic projectors.
  Each left the breaking unchanged.
- **LDA removes it.**
- **Moving `ux` moves it.** A tiny seed moment on Y moved `ux`, and the breaking moved with it.
- **The patch removes it:** 100× smaller (§1).

**When it triggers:** all of these must hold:
- noncollinear and magnetic;
- a GGA functional;
- PAW pseudopotentials (the faulty line exists only in the PAW one-centre code);
- starting moments not all parallel or antiparallel.

That covers all-in-all-out pyrochlores, canted and spiral magnets, and any noncollinear starting configuration.
The breaking axis is set by whichever magnetic atom comes first in the input. One wrinkle: `is_parallel`
(`PW/src/divide_class.f90`) uses an absolute threshold, |a×b|² < 1e-6, so very small starting moments can be
judged parallel. That happened in test L5.

## 3. MnBi2Te4 is not affected

- **Starting moments:** Mn1 has `angle1 = 0°` and Mn2 has `angle1 = 180°`, antiparallel along z. Bi and Te
  start at zero, which passes the parallel test.
- **Result:** `lsign = .true.` with `ux = z`. Every MBT spin-orbit SCF prints
  `Fixed quantization axis for GGA: 0 0 1`, checked on:
  - local: 152188, the qe76 debug runs, T2⁻ and T1⁺ modes;
  - cluster: 195245, 202833, 193147, the u_0.0 runs.
- **Consequence:** both GGA paths use s = sign(m_z) consistently, as QE designed (§2.3).
- **No symmetry breaking:** every MBT magnetic operation maps the z line to itself, so s either stays put or flips
  everywhere, and both are harmless.
- **No-SOC MBT runs** (`base/wo_soc`) are collinear (`nspin = 2`) and never enter this code.

**Caveats:**
1. **The fixed-axis mode is itself an approximation.** s jumps where m_z changes sign, e.g. in small induced
   moments on Bi or Te. That is standard QE behaviour and symmetric, and it is identical in the base and mode
   runs, so it cancels consistently in the finite difference.
2. **A canted or noncollinear MBT start would hit this bug.** The patched build is safe either way: it changes
   nothing when `lsign` is true.

## 4. Production reruns (submitted 11:15)

The patched pw.x was copied over the old `pw.x` (md5 `1541022`) in each structure folder, with a `NOTE` written
there. Inputs and `send_job.sh` are unchanged (SCF → NSCF → bands):

| structure | job | queue | previous run (buggy) | previous NSCF wall time |
|---|---|---|---|---|
| base/u_3.0 | 207546 | dko64m4, 64 | 192344 | 10.4 h |
| mode_1 | 207547 (started 11:34 on n109) | dko64m4, 64 | 192345 | 10.8 h |
| mode_2 | 207548 | dko48m4, 48 | 192346 | 23.7 h |
| mode_3 | 207549 | dko48m4, 48 | 192347 | 23.8 h |

**Freeing a node for mode_1:** MBT u_4.0 job **205075 was killed** for stagnation: 82 iterations in 17.5 h, SCF
accuracy flat at 1.4-2.8e-3 (conv_thr 1e-7). Its `scf.out` and `scf.in` are saved in `MnBi2Te4/base/soc/u_4.0/`
as `MnBi2Te4.scf.out.205075_killed_2026-10-06` and `MnBi2Te4.scf.in.205075`. MBT u_3.0 (202971) was left running:
115 iterations, accuracy 0.6 → 1.9e-4 → 5.5e-4, slow and noisy but still about 10× lower over the last 50
iterations.

**Status at 18:30.**
- **mode_1 and mode_3:** SCFs converged with the patched pw.x (68 and 93 iterations); NSCFs running.
- **base (207546):** the SCF **did not converge**. It hit `electron_maxstep = 500`, stagnating at 2-9e-5 Ry from
  iteration ~150 on (target 1e-7; plain mixing, β 0.03). The job script went on to the NSCF anyway. That NSCF is
  kept for the DFT symmetry check (the SCF density is symmetrized regardless), but **not for dθ**.
- **mode_2 (207548):** the SCF stagnated the same way (4-9e-6 Ry over iterations 100-307) and was killed. Its
  scf.out is saved as `mode_2/Y2Ir2O7.scf.out.207548_killed_2026-10-06`. Resubmitted as **209009** with
  `mixing_mode = 'local-TF'`, `mixing_beta = 0.1`; the old input is kept as `Y2Ir2O7.scf.in.plain0.03.bak`.
- **Before the fix:** the same inputs converged with the buggy binaries (192344: 252 iterations; 193145: 411). The
  patch changes the convergence behaviour of the high-symmetry structures; the cause is unknown.
- **mode_2 cleanup:** deleted the stale 178380bc889 (330 GB, QE 7.5), 192346 (554 GB, segni bug) and the killed
  207548bc889. mode_2 went from 883 GB to 135 MB. Their small text outputs are in
  `data/Y2Ir2O7/phonon/GM2-/mode2/Q_0.01A/output/`, and the Wannier-stage templates (rigid-shift script, `.win`,
  `pw2wan.in`, job scripts) are in `mode_2/wannier_templates_from_192346/`.
- **Base symmetry check queued:** `eig_symmetry_xml_full.py` runs on base's NSCF XML when 207546 finishes. It was
  validated on the old 192344 NSCF: 7.586 meV, matching the `.eig` test.

**mode_2, second attempt (21:00).**
- **The run:** job 209009 (`local-TF`, β = 0.1, from scratch) was 3× faster early on and reached 6e-6 to 1e-5 Ry at
  iterations 41-47.
- **Then a period-2 oscillation:** accuracy alternated 5e-5 / 1.9e-4 Ry, and the total energy flipped by 1.2e-4 Ry
  every step.
- **Diagnosis** (`calculations/diagnostics/2026-10-05_ws_tolerance_window/scf_oscillation.py`): the oscillation is
  in the **Ir Hubbard occupations**. Tr ns↑ and ns↓ swing by ~2e-4 and the Hubbard moment by 4e-4, about 7× the
  noise. The charges and integrated moments of all atoms (Ir, Y, both O types) do not oscillate.
- **Secondary signal:** CG "eigenvalues not converged" warnings went from about 0.5 to 3 per iteration in the
  oscillating stretch.
- **Clean stop:** the run was stopped at iteration 75 with a `Y2Ir2O7.EXIT` file (the job's NSCF/bands inputs were
  renamed first, so the job ended after the SCF). QE wrote the density, ns and becsum.
- **Restart:** job **209022** uses `startingpot = 'file'` from that state, `local-TF`, β = 0.05.

**Davidson fixes the SCF stall (23:30).**
- **The hint:** mode_2 restart 209022 (local-TF, β 0.05, CG) still sat at 1e-5 to 1.4e-4 Ry, with ~6.6 CG
  "eigenvalues not converged" warnings per iteration.
- **Restarted with Davidson** (job **209040**): `diagonalization = 'david'` from its saved state, everything else
  unchanged. Warnings dropped ~10× (0.7 per iteration), iterations got faster (~60 s vs ~100 s), and the SCF
  **converged in 39 iterations** (7e-8 Ry).
- **The final step was abrupt** (3e-5 → 7e-8), but sound: the Hubbard ns landed exactly between the two values it had
  been alternating between, and the total energy trend was continuous (−4624.80627 Ry, consistent with earlier
  plateaus to ~2e-5).
- **Base:** rerun 209024 (CG) was still at 2.5-6.7e-5 Ry after 88 iterations. It was stopped cleanly and restarted
  the same way with Davidson as **209042**.
- **Working conclusion (not yet confirmed for base):** the stall in the patched runs came from CG not converging some
  bands in each iteration, which gave the SCF a noise floor of ~1e-5 Ry, not from the mixing.

**Still to do:**
1. **Wannier stage.** Base trial_02: `-pp` → pw2wannier90 → wannier90 (`num_iter = 0`) → `restart = plot` job.
   Each mode: the `generate_rigid_shift_amn.sh` `.amn` → pw2wannier90 → wannier90 → plot.
2. **Checks:** `.eig` symmetry on the new base (`eig_symmetry_any.py`), curvature symmetry, then dθ.
3. **Apply the 5e-3 WS tolerance** (#13).
4. **Delete the old buggy runs** once the new ones validate (YIO is ~3 TB on the cluster); this needs the user's
   OK.
5. **Update CLAUDE.md** pitfall 7 and its "correct versions" list.
6. **Decide whether to report the `segni` bug to QEF.**
