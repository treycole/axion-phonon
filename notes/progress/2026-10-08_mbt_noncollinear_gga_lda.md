# MnBi2Te4: noncollinear PBE has a spurious lower-energy state; the base moves to LDA+U (2026-10-05 → 10-08)

The MBT base SCF would not converge on the fixed QE build. This note records the investigation:
- the settings audit (§1);
- the code audit against QE 7.2 (§2);
- the symmetry the SCF can't enforce (§3);
- the starting-wavefunction finding (§4);
- the lower-energy asymmetric state and its energy breakdown (§5);
- the LDA, collinear and noncollinear tests that pin it on the noncollinear GGA (§6, §7);
- the real-space magnetization analysis: transverse m, jumps rather than kinks (§8);
- the full run table (§9);
- why noncollinear PBE is fragile and the LDA trade-off (§10);
- the decision and the new LDA+U base (§11);
- the YIO connection and the YIO LDA tests (§12);
- open items (§13).

Bug number: **#18** ([`../bugs/18_qe_noncollinear_gga_spurious_state.md`](../bugs/18_qe_noncollinear_gga_spurious_state.md)).
All MBT runs are U = 4 eV (`HUBBARD (ortho-atomic)`, Mn 3d) unless stated. Cluster paths are under
`/mnt/wkdk/bc889/axion/MnBi2Te4/base/`.

## Bottom line

- **With SOC, MBT must run noncollinear.** Noncollinear PBE then has a spurious self-consistent state:
  - **Lower energy:** 55 mRy below the physical one.
  - **Asymmetric:** the two Mn layers stop being time-reversed copies.
  - **Transverse magnetization:** 3.5 μB of transverse m in the interstitial.
  - **Never converges:** the SCF stalls around 1e-5 Ry at best.
- **A random start falls into it, and stalls.** QE does not enforce the symmetry that forbids it. An exactly
  symmetric start (`startingwfc = 'atomic'`) avoids it.
- **Collinear PBE (no SOC) and LDA (with SOC) have no such state.** Both converge from random starts in 12-13
  iterations.
- **The cause:** QE's noncollinear GGA builds n↑/n↓ = (n ± sign(m·z)|m|)/2. Once m has a transverse part, that
  spin density **jumps** where m_z changes sign at finite |m|. The gain is entirely in the grid XC energy.
  This is a jump, not the |m| kink (§8).
- **Decision:** the MBT base is now **LDA+U** with `rel-pz` PAW. Base 209653 converged in 12 iterations, symmetric
  to 1e-5, NSCF done (§11).
- **U:** stays at 4.0 eV (decided 2026-10-08).
- **Why the atomic start converges:** it starts exactly in the symmetric subspace, which the SCF preserves. In PBE
  the symmetric state is a saddle and the start only works because convergence beats the growth of numerical
  noise (§4.1).
- **YIO:** the same pattern. LDA converges the base in 21-22 iterations with an exact AIAO order and a 0.68 eV gap
  on the SCF mesh; PBE wanders with |m|abs creeping up (§12).

## 1. Starting point: the settings audit

Every MBT base SCF on the fixed binaries stalled at 1e-4 to 1e-3 Ry:
- 195245 (u_4.0);
- 202833 (u_4.0, plain β 0.2);
- 205075 (u_4.0, local-TF β 0.1 restart);
- 202971 (u_3.0, plain β 0.2).

The only converged base with SOC+U is **152188** (QE 7.2, u_4.0, 13 iterations).

| setting | 152188 (converged) | stalled runs |
|---|---|---|
| QE | 7.2 | 8.0dev (`~/q-e` c41c409; later `~/q-e-segnifix` e270b05) |
| diagonalization | **david** | **cg** |
| mixing | local-TF β 0.4, default ndim 8 | plain β 0.2 or local-TF β 0.1, ndim 20 |
| degauss (m-v) | 0.001 | 0.005 or 0.02 |
| conv_thr | 1e-6 | 1e-7 |
| `symmetry_with_labels` | off | on |
| iteration-1 `ethr` | **1e-5** (see below) | 1e-2 |

**Two details mattered:**
- **CG, ndim 20 and the wide smearing were workarounds** for the old QE 7.5/7.6 iteration-1 blow-up (bug #12), and
  were never reverted. No converged MBT run (≈60, all QE 7.2) ever used CG.
- **152188 had `restart_mode = 'restart'` with no restart files.** It printed "restart disabled" and started from
  atomic, but `input.f90` had already set `startingpot = 'file'`. So `setup.f90` chose the strict first-iteration
  threshold `ethr = 1e-5` (24.4 Davidson passes per k). Every later cold start got 1e-2 (3.5 passes).

**Restoring the 152188 recipe did not fix it.** 209025 (segni-fixed binary, david, local-TF 0.4, degauss 0.001)
still separated from 152188 at iteration 4 and stalled. Its absolute magnetization climbed from 10.35 to 12.6 μB,
where 152188 stays at 10.0-10.26 μB.

## 2. Code audit, QE 7.2 vs the build (934f8cd)

Routines on the MBT path, compared function by function:

| routine | finding for MBT |
|---|---|
| `new_ns_nc` (Hubbard ns symmetrizer) | identical to 7.2 except our `invs` patch (#7), which is inert for MBT: all 12 ops map each Mn to itself |
| `sym_rho` noncollinear branch | unchanged; edits only touch the collinear `colin_mag == 2` path |
| `PAW_symmetrize` | only the collinear branch changed; `symmetry_with_labels` sets `colin_mag = 0`, which only matches atoms by element |
| `mix_rho` | restructured; default path equivalent to 7.2 (`simple_magn_mix` default `.FALSE.`, spin ranges default `1:nspin`, mix-file record layout correct) |
| `electrons.f90` ns copies (`scf_ns_copy`) | same direction as 7.2 |
| PAW one-centre noncollinear GGA | OpenACC rewrite (2024); equivalent apart from the two known fixes. `v_rad` (#12) is fixed. `segni` (#16) never affects MBT, since `lsign = .true.` ("Fixed quantization axis for GGA: 0 0 1") |

**One real path change.** Since `ec2b117bc` (QE 7.3), noncollinear +U without J runs as `lda_plus_u_kind = 0`, the
Binci-Marzari noncollinear Dudarev implementation, which the release notes call "experimental". 7.2 ran it as
`kind = 1` (Liechtenstein, averaged j radial WFs). At J = 0 the two are the same functional. Matched 7.2/7.6 test
runs (qe_ab, 2026-09-26) give identical first iterations. **Not the cause.**

## 3. The symmetry QE cannot enforce

- **Same group everywhere.** All runs (QE 7.2 and the new build, any cell, labels on or off) find the same 12
  operations: E, C3, three C2′ and their inversion partners, six of them with time reversal. The full blocks
  (integer and Cartesian matrices, the time-reversal flags) diff clean.
- **The missing operation.** MBT's magnetic space group (R_I-3c, BNS 167.108, in the CIF as centring
  `x+1/2, y+1/2, z+1/2, −1`) also contains the anti-translation {1′ | ½,½,½}. It maps the Mn1 layer onto the Mn2
  layer with spins reversed.
- **Why QE drops it:**
  - QE stores at most one operation per rotation matrix, and this one has the identity rotation.
  - With `symmetry_with_labels`, `sgam_at` (`symm_base.f90:499-521`) finds "identity + (½,½,½)", prints
    "This is a supercell, fractional translations are disabled", and drops it. `inverse_s` (`symm_base.f90:101`)
    pairs inverses by rotation matrix alone, so the operation can't simply be forced in.
- **Consequence:** nothing makes the two Mn layers equivalent. Only the SCF itself can keep them mirrored.

**The cell is not the issue.** The u_4.0 production cell is the exact hex CIF geometry: a = 4.33360,
c = 81.85200 Å, every atom at u·c, a1 = (2.1668, −1.25100256328, 27.284). It has exact C3, inversion and c/2
half-translation. 152188's cell went through the rhombohedral CIF's 7-digit a0 and α and is off by up to 1e-6 Å.
Both give 12 ops. Test B (209034: 152188's cell, no labels) matched 209025 to 4 digits at iterations 1-2. The
current u_4.0 cell is kept (`mbt-symmetry-with-labels` memory corrected).

## 4. The random start breaks the Mn1 ↔ Mn2 equivalence

`startingwfc` defaults to `'atomic+random'`. After iteration 1:

| run | code, procs | Mn1 / Mn2 minority Tr ns | split |
|---|---|---|---|
| 152188 | 7.2, 64 | 0.30058 / 0.30058 | 0 |
| 192889 (qe_ab test) | 7.2, 48 | 0.28290 / 0.28651 | 3.6e-3 |
| 209025 | new, 64 | 0.29982 / 0.30402 | 4.2e-3 |
| B2 (209047) | new, 48 | 0.29924 / 0.30314 | 3.9e-3 |

152188's random start happened to be symmetric; that's why it converged. **Test D** (209052: B2's settings +
`startingwfc = 'atomic'`) mirrors exactly after iteration 1 (0.30060 / 0.30060) and converges in 18 iterations to
−7140.60174 Ry. That's 0.16 mRy from 152188, with every layer pair mirrored to 3e-5 μB.

### 4.1 Why the atomic start converges

1. **The symmetric densities form an invariant subspace.** The exact-CIF structure is symmetric under the
   anti-translation {1′|½,½,½}, so the Kohn-Sham map (density in → density out) commutes with it, even though QE
   doesn't impose it. A density that is exactly anti-translation-symmetric produces an exactly symmetric output:
   start inside the subspace and you stay inside, up to numerical noise.
2. **`startingwfc = 'atomic'` starts inside it.**
   - **The orbitals:** the superposed atomic orbitals have Mn2's spinors rotated by `angle1 = 180°`, so they are
     exactly the time-reversed, half-shifted images of Mn1's.
   - **No random part:** without the random component, the starting ρ, m, ns and becsum are exactly symmetric.
   - **Evidence:** D's iteration-1 Mn ns are 0.30060 / 0.30060.
3. **The spurious state lies outside it.** It breaks the layer mirror, and its transverse texture grows only
   together with the asymmetry (∫|m_⊥| 0.10 in D against 3.46). From inside the subspace the SCF can't reach it.
4. **In PBE the symmetric state is a saddle, and it converges before the noise grows.**
   - **PBE:** random starts seed the antisymmetric direction at about 4e-3, and it grows.
   - **LDA:** the same seed decays (E1: 4e-3 → 3e-5). So in LDA the symmetric state is a real minimum; in PBE it is
     unstable in that direction.
   - **The noise in D:**
     - the 375-point FFT grid is odd, so the half-cell shift is not a grid translation, and the real-space XC
       evaluation breaks the symmetry slightly;
     - D also used 152188's cell, off by up to 1e-6 Å.
   - **Result:** D's final layer mismatch is 2.6e-5 μB. It converged in 18 iterations, before that noise could be
     amplified into the spurious state.
5. **So in PBE the atomic start is knife-edge.**
   - **The base:** it works because the base keeps the anti-translation.
   - **The T2⁻ modes:** they break it physically, so the antisymmetric amplitude is there from the start and nothing
     protects them.
   - **Why that justifies the switch:** in LDA the symmetric state is stable from any start (E1).
   - **Hardening the base further:** `nr1 = nr2 = nr3 = 384`, an even grid, would make the half-cell shift exact on
     the grid and remove the main noise source.
6. **The leftover restart flag** gave 152188 and D a well-converged first diagonalization (`ethr = 1e-5`). That
   helps the first density but is not what selects the state. 152188 converged from a random start because its
   random part happened to come out symmetric (Mn1 = Mn2 ns at iteration 1). Why QE 7.2's generator did that on that
   run is not known.

## 5. The asymmetric state is lower in energy

Left to run, the random-start PBE runs settle into a definite state:
- **Layers not mirrored:** Te7/Te8 +0.016/−0.032 μB, where the symmetric state has ±0.027 μB.
- **|m|abs ≈ 13.3 μB**, against 10.2 μB in the symmetric state.
- **Accuracy plateau:** they never get below about 5e-6 to 1e-5 Ry.

209025 was stopped at iteration 149 with EXIT. **209380** then evaluated that state at once
(`startingpot = 'file'`, `conv_thr = 1e-3`, accuracy 8.6e-6 Ry), which prints the energy breakdown:

| term (Ry) | asymmetric (209380) | symmetric (D) | Δ (mRy) |
|---|---|---|---|
| **total** | −7140.65681 | −7140.60174 | **−55.1** |
| one-electron | −651.24918 | −651.55694 | +307.8 |
| Hartree | 403.04045 | 403.16008 | −119.6 |
| **XC (grid)** | −318.17895 | −317.89433 | **−284.6** |
| one-centre PAW | −5901.85611 | −5901.89759 | +41.5 |
| Hubbard | 0.08339 | 0.08349 | −0.1 |
| Ewald | −672.20712 | −672.20716 | +0.04 (cells) |
| D3 | −0.28929 | −0.28929 | 0 |

B2 reached the same state independently (−7140.6569 Ry, |m|abs 13.5 μB).

**The gain is in the grid XC term.** The PAW one-centre term rises, so this isn't another one-centre bug.

## 6. LDA has no such state (E1, E2)

Diagnostic runs with `input_dft = 'pz'` on the `rel-pbe` PAW set, no D3, the exact cell with labels, and
`diago_thr_init = 1e-5`:

| run | start | converged | E (Ry) | \|m\|abs | Mn1/Mn2 minority ns (iteration 1 → final) |
|---|---|---|---|---|---|
| E1 (209378) | random | 12 iterations | −6949.61020821 | 9.99 | 0.32319/0.32743 → 0.31990/0.31993 |
| E2 (209379) | atomic | 12 iterations | −6949.61020825 | 9.99 | 0.32454/0.32454 → 0.31993/0.31993 |

The random start's 4e-3 split decays to 3e-5, and both runs reach the same energy to 4e-8 Ry.

## 7. Collinear vs noncollinear PBE, without SOC (C1, C2)

`test_C_collinear_vs_nc_2026-10-07`:
- PBE+U, scalar-relativistic pslibrary PAW (Mn `pbe-spn` 0.3.1, Bi `pbe-dn` 1.0.0, Te `pbe-n` 1.0.0);
- random start, no labels, otherwise as above.

With every moment along ±z these are **the same physical problem**:

| run | spin treatment | converged | E (Ry) | \|m\|abs | layer pairs |
|---|---|---|---|---|---|
| C1 (209651) | collinear `nspin = 2` | 13 iterations | −7099.74790 | 10.16 | mirrored (Te7/Te8 ±0.0240) |
| C2 (209652) | noncollinear, moments ±z | **no** (59 iterations, 1e-2 to 2e-2) | ≈ −7099.79 (±7 mRy) | 12.3 | broken (+0.0266/−0.0234) |

Only the noncollinear treatment finds the lower, asymmetric, non-converging state, and it does so without SOC.

## 8. What the state is: transverse magnetization, and jumps (not kinks)

`charge-density.dat` holds n, m_x, m_y, m_z in G-space. Transformed to the 375³ grid:

| | asymmetric PBE (209025) | symmetric PBE (D) | LDA base (209653) |
|---|---|---|---|
| ∫n | 138.000 | 138.000 | 138.000 |
| ∫\|m\| (= QE's \|m\|abs) | 13.21 | 10.17 | 9.99 |
| ∫\|m_z\| | 10.32 | 10.15 | 9.98 |
| **∫\|m_⊥\|** | **3.46** | 0.10 | 0.08 |
| **max \|m_⊥\|** (μB/bohr³) | **0.31** | 0.0018 | 0.0023 |
| grid fraction with \|m_⊥\| > \|m_z\| | 66% | 10% | 4% |

The extra ∫|m| (3.03 μB) is almost all transverse: ∫(|m| − |m_z|) = 2.88. Net ∫m_x = ∫m_y = 0. **The low state is
a noncollinear texture between the atoms** that collinear PBE cannot represent.

**Which spin density the GGA sees.** QE's noncollinear GGA needs n↑/n↓. MBT has `lsign = .true.`, so both the smooth
grid (`gradcorr`) and the PAW sphere use the fixed axis, n↑,↓ = ½(n ± s|m|) with s = sign(m_z). Total variation
(Σ|Δf| over grid links, a proxy for ∫|∇f|):

| | symmetric (D) | asymmetric (209025) |
|---|---|---|
| m_z (what a collinear code uses) | 88,037 | 95,237 |
| \|m\| (local frame; kinked at m = 0) | 87,636 | 113,721 |
| **s\|m\| (what QE's GGA uses)** | 88,120 (**+0.1%**) | 141,459 (**+48.5%**) |
| links where s flips with \|m\| > 1% of max on both ends | **0** | **84,258** |
| share of the s\|m\| excess on those links | — | 60% |

**Reading the table:**
- **Symmetric state:** s|m| = m_z. The sign flips only where |m| → 0, so the GGA spin density is smooth: no jump, and
  no kink either.
- **Asymmetric state:** the transverse m makes m_z cross zero at finite |m|, and s|m| **jumps** by up to 2|m_⊥|
  there. GGA differentiates n↑/n↓, and a jump is a gradient spike.
- **Why that lowers the energy (inference, not separately verified):** PBE's exchange enhancement F_x grows with
  the reduced gradient, so spurious gradients make the XC energy more negative. That fits the −285 mRy.
- **Why it doesn't converge:** QE's potential omits the derivative of s, so energy and potential are inconsistent
  and the SCF never settles.

**Is it the |m| kink?** Not for MBT. The kink of the local-frame |m| at m = 0 is the hypothesis for YIO's stall
(#17). With `lsign = .true.`, MBT never uses |m| alone. What hurts MBT is the jump surface of s|m|, which exists
only once transverse m has formed. Both are ways the noncollinear GGA's spin density stops being smooth; LDA uses no
gradients and is immune to both. For YIO, smoothing the kink (R1-R3, δ up to 1e-3 e/bohr³) did **not** remove the
stall (§12).

### 8.1 What a "jump" is, with numbers

GGA needs n↑,↓ = ½(n ± X). For MBT QE uses X = sign(m_z)·|m|: the length of m, with the sign of its z component.
- **A kink** is a continuous function whose slope changes abruptly, like |x| at 0.
- **A jump** is a function whose value changes abruptly between neighbouring points.

**Walking through the interstitial of the asymmetric state.** Grid points are about 0.14 bohr apart along the
rhombohedral vectors (|a1| = 51.8 bohr / 375). Take a transverse part m_x = 0.30 μB/bohr³ (the asymmetric state's
max |m_⊥| is 0.31) while m_z drifts slowly through zero:

| point | m_x | m_z | \|m\| | sign(m_z) | X = sign(m_z)·\|m\| |
|---|---|---|---|---|---|
| 1 | 0.30 | +0.02 | 0.3007 | + | +0.3007 |
| 2 | 0.30 | +0.01 | 0.3002 | + | **+0.3002** |
| 3 | 0.30 | −0.01 | 0.3002 | − | **−0.3002** |
| 4 | 0.30 | −0.02 | 0.3007 | − | −0.3007 |

**What happens between points 2 and 3:**
- **Physically almost nothing.** m turns by 4° and |m| is unchanged.
- **X changes by 0.60 in one grid step.** That's a slope of about 4.3 μB/bohr⁴; m_z itself changes by only 0.02.
  So n↑ and n↓ each step by about 0.3 between neighbours.

**Why that matters:**
- **Energy:** PBE's exchange enhancement grows with the reduced gradient |∇n_σ|/n_σ^{4/3}, so an artificial step
  lowers E_xc.
- **Potential:** QE's XC field omits ∂sign/∂m, so energy and potential disagree at the step, and the SCF can't
  settle.

**When it occurs:**
- **The symmetric state can't produce it.** m_⊥ ≈ 0, so m_z crosses zero only where |m| → 0, and X = m_z is
  smooth: 0 jump links.
- **The transverse texture produces it.** It makes m_z cross zero at finite |m|: 84,258 jump links.

**For contrast, the kink.** YIO's local frame uses X = |m|. Crossing m = 0 gives 0.02 → 0.01 → 0.01 → 0.02: a V,
continuous, with a finite gradient. That's the #17 hypothesis, and smoothing it (R1-R3) didn't remove the YIO stall.

## 9. All runs

| run | functional | SOC | spin | start | converged | E (Ry) | ΔE vs symmetric (same functional, SOC) | \|m\|abs | layer mirror (worst \|m₁+m₂\|, μB) |
|---|---|---|---|---|---|---|---|---|---|
| 152188 (QE 7.2) | PBE | yes | nc | random* | 13 | −7140.60190 | −0.2 mRy† | 10.21 | 4e-5 ✓ |
| D 209052 | PBE | yes | nc | atomic | 18 | −7140.60174 | reference | 10.17 | 3e-5 ✓ |
| 209025 (eval 209380) | PBE | yes | nc | random | **no** (149, ~1e-5) | −7140.65681 | **−55.1 mRy** | 13.26 | 1.7e-2 ✗ |
| B2 209047 | PBE | yes | nc | random | **no** (127, ~1e-5) | −7140.6569 | **−55.2 mRy** | 13.48 | 7.6e-3 ✗ |
| C1 209651 | PBE | no | **collinear** | random | 13 | −7099.74790 | reference | 10.16 | 1e-4 ✓‡ |
| C2 209652 | PBE | no | nc | random | **no** (59, 1e-2) | ≈ −7099.79 | **≈ −44 mRy** | 12.28 | 3.1e-3 ✗ |
| E1 209378 | LDA (PBE PAW) | yes | nc | random | 12 | −6949.61021 | +4e-5 mRy | 9.99 | 3e-5 ✓ |
| E2 209379 | LDA (PBE PAW) | yes | nc | atomic | 12 | −6949.61021 | reference | 9.99 | 3e-5 ✓ |
| **LDA+U base 209653** | LDA (`rel-pz`) | yes | nc | atomic | 12 | −6947.39551 | (production) | 9.99 | 3e-5 ✓ |
| 202833, 205075, 202971 | PBE | yes | nc | random | no (1e-4 to 1e-3) | — | — | 12-13 | ✗ |

Footnotes:
- \* 152188's leftover restart flag gave it `ethr = 1e-5`, and its random start happened to be symmetric.
- † Different QE version, and a cell off by 1e-6 Å.
- ‡ The collinear output prints moments to 4 decimals.

**Comparing energies:**
- **Within a block only.** Rows with different functionals or pseudopotentials differ by hundreds of Ry for
  method reasons.
- **Symmetry:** all runs find the 12 operations. The anti-translation is never enforced; the mirror column is what
  each run converged to.

## 10. Why noncollinear PBE is fragile, and the LDA trade-off

**The spin densities, per functional:**
- **Collinear:** n↑/n↓ = ½(n ± m_z) are well defined, for LDA and GGA alike.
- **Noncollinear LDA:** it needs only n and |m| (the local spin-frame eigenvalues). That is an exact, rotation-invariant
  extension (Kübler et al.).
- **Noncollinear GGA:** it needs gradients of n↑/n↓, which a rotating spin axis does not define uniquely. QE uses
  either:
  - the local frame ½(n ± |m|), kinked at m = 0, for noncollinear orders (`lsign = .false.`); or
  - a fixed axis ½(n ± sign(m·u)|m|), for collinear orders (`lsign = .true.`). This one jumps wherever m·u = 0 at
    finite |m|.
- **Better noncollinear GGAs** exist (Scalmani-Frisch, Peralta et al.); QE doesn't have them.

**When it bites:**
- **All of these hold:** a noncollinear calculation (forced by SOC), a GGA, and a magnetization free to grow
  transverse parts or with large regions of small m·u. An antiferromagnet's interlayer nodal regions are such
  regions.
- **And nothing holds the magnetization collinear.** Here the anti-translation isn't enforced, and the random
  start seeds it.
- **What's safe:** collinear runs without SOC (C1) are fine.

**Unique ground state.** The exact functional's ground state is the lowest-energy state. An approximate functional
that is not well defined for a class of densities can have a lower minimum that exploits the defect. Several
self-consistent solutions are common in DFT+U anyway. Evidence that the low PBE state is unphysical:
- C1 = C2's physics, but only C2 finds it.
- LDA, whose noncollinear form is exact, has none.
- It never converges, because energy and potential are inconsistent.
- 0.6-0.75 eV per cell from induced Te/Bi moments is implausible.
- It breaks the experimentally established magnetic group R_I-3c.

**Why the MBT literature uses PBE.** No paper found justifies it. PBE+U (U ≈ 3-5.34 eV, VASP PAW, DFT-D3) is
convention: VASP's default PAW sets, and D3 parameterized for GGAs to relax the vdW gap. Neither applies at our fixed
experimental geometry in QE.

**The one direct comparison.** A DMC benchmark of MBT (arXiv:2408.03248) finds that LDA+U vs PBE+U differences are
mostly a 1-1.5 eV shift in U. The DMC-optimal bulk values are **PBE+U ≈ 3.5-4.0 eV and LDA+U ≈ 4.4 eV**.

**LDA trade-offs:**

| | LDA | PBE |
|---|---|---|
| noncollinear + SOC + PAW in QE | exact extension, no axis or sign trick | fragile (this note, #16, #17) |
| equilibrium structure | overbinds 1-2% | slightly underbinds, needs D3 |
| relevant here? | no: fixed CIF geometry and prescribed displacements; D3 is energy-only | — |
| magnetism, bands | U means about 1 eV less; small gap and band inversion shift; θ = π protected | literature standard |
| consistency | differs from YIO (PBE) and most MBT literature | matches |

## 11. Decision and the LDA+U base

**MBT moves to LDA+U.** Base: `soc_lda/u_4.0/`, job **209653** (dko64m4, 64 procs).

**Pseudopotentials:**
- **Mn and Bi:** official pslibrary `Mn.rel-pz-spn-kjpaw_psl.0.3.1` and `Bi.rel-pz-dn-kjpaw_psl.1.0.0`.
- **Te:** no `rel-pz` Te of any kind is on the QE site, so `Te.rel-pz-n-kjpaw_psl.1.0.0.UPF` (md5 37a9d751) was
  generated in `soc_lda/pseudo_gen/`:
  - **How:** with `ld1.x` (built in `~/q-e-segnifix`; needs `module load intel/2024 intel/ompi`), from the
    `PP_INPUTFILE` embedded in the official `Te.rel-pbe-n-kjpaw_psl.1.0.0`, changing only `dft='PZ'`.
  - **Check:** the same recipe with PBE reproduces the official file. Header identical, total_psenergy to 1e-7 Ry,
    arrays to 8e-10.

**Input:** the u_4.0 PBE input with:
- no `vdw_corr`;
- `diago_thr_init = 1e-5` and `startingwfc = 'atomic'` (SCF).

Unchanged: the exact cell, labels, david, local-TF 0.4, degauss 0.001, nbnd 166, conv_thr 1e-7. Binary: segnifix
pw.x e270b05.

**Result:**
- **SCF:** converged in 12 iterations. E = −6947.39550993 Ry, |m|abs 9.99, Mn ±4.5536 μB.
- **Symmetry:** every layer pair mirrored to 3e-5; Mn1/Mn2 ns identical.
- **NSCF:** done. `tmp/` is 29 GB; `tmp.tar.gz` is not yet made.

## 12. YIO

**The two YIO bugs in these terms:**
- **#16:** the segni bug was the fixed-axis jump applied with a meaningless axis in the PAW sphere.
- **#17:** the base SCF stall after the fix. The leading hypothesis was the local-frame |m| kink in the PAW sphere.

**The kink hypothesis is not supported.** `stall_test_2026-10-07` R1-R3 (|m| regularization, δ = 1e-5, 1e-4,
1e-3 e/bohr³) each ran 500 iterations at ~3e-5 Ry, like T1c. Only the unpatched binary (T1) converged, in 153
iterations. The patched runs also carry more |m|abs (3.35-3.37 vs 3.32 μB), the same direction as MBT.

**LDA tests (2026-10-08)** in `Y2Ir2O7/phonon/base/u_3.0/lda_test_2026-10-08/`:
- **Setup:** base SCF from scratch, U = 3, segni-fixed binary, production settings (Davidson, plain β 0.1, labels).
- **YP_pbe** (210525): PBE control.
- **YL1_lda_pbepp** (210526): `input_dft = 'pz'` on the `rel-pbe` PAW (functional only).
- **YL2_lda_relpz** (210527): official `rel-pz` PAW (Y `spn` 1.0.0, Ir `n` 0.2.3, O `n` 0.1).

**Result: LDA converges and PBE does not.** All three runs find 48 operations (36 with fractional translations).

| run | converged | E (Ry) | Ir moments (all four) | off ⟨111⟩ | O \|m\| max | \|m\|abs | gap, 8-k SCF mesh |
|---|---|---|---|---|---|---|---|
| YP PBE | **no**: 72+ iterations, wandering 2e-7 to 1.5e-5 Ry around conv_thr 1e-7 | — | 0.3099 μB, equal | 0.000° | 0.0132 | 3.32 → **3.40**, creeping up | — |
| YL1 LDA (`rel-pbe` PAW) | **21 iterations** | −4501.58867 | 0.2830 μB, equal | 0.000° | 0.0127 | 3.21 | 0.683 eV |
| YL2 LDA (`rel-pz` PAW) | **22 iterations** | −4499.85691 | 0.2829 μB, equal | 0.000° | 0.0127 | 3.21 | 0.685 eV |

**Reading it:**
- **Same family as MBT.** YIO's base stall (#17) follows the MBT pattern. PBE wanders with |m|abs creeping up:
  YP 3.40, the patched restarts T1c/R1-R3 at 3.35-3.37, against the converged unpatched T1 at 3.32. LDA converges
  in about 20 iterations.
- **The mechanism here is different.** YIO keeps all 48 operations, enforced every iteration, so there is no layer
  asymmetry to fall into. The drift has to happen inside the symmetric subspace: a symmetry-allowed noncollinear
  texture rewarded by the GGA.
- **What isn't established:** which part of the GGA does it. The PAW kink regularization (R1-R3) didn't help; the
  smooth-grid `gradcorr` term was not touched.
- **LDA keeps the AIAO insulator at U = 3.** Four equal moments exactly along ⟨111⟩, 9% smaller than PBE, and a
  0.68 eV gap on the SCF mesh. For reference, PBE's 512-k NSCF gap was about 0.42 eV; the meshes differ.
- **θ stays 0.** A trivial AIAO phase is what the project needs for YIO (memory `y2ir2o7-axion-phonon-goal`).

**Caveat:** LDA+U at U = 3 is not PBE+U at U = 3. The YIO AIAO moment is non-monotonic in U, so if YIO moves to
LDA, check the gap on the full NSCF mesh first. **Not decided:** whether YIO switches. YP was still running at
15:25 on 10-08.

## 13. Open

1. **LDA U for MBT: decided 2026-10-08, U = 4.0 for now** (user).
   - **Literature:** published DFT+U values for MBT span about 3-5.34 eV, almost all PBE+U.
   - **Benchmark:** the DMC-optimal LDA+U value is 4.4 eV.
   - **Revisit** if a benchmark comparison is needed.
2. **MBT T2⁻ modes on LDA:** exact cell + Q·mode, same settings. These modes break the anti-translation physically,
   so in PBE no symmetric start could have protected them. That's the main reason for the switch.
3. **YIO functional:** LDA converges (§12). Decide whether YIO switches, and check the full-mesh gap at LDA U = 3
   first.
4. **A QE patch to enforce the anti-translation** (project charge onto even h+k+l, magnetization onto odd h+k+l, and
   pair-average ns/becsum). Not needed for LDA; would make PBE usable for the base only.
5. **Cleanup:** C2 and B2 were killed 2026-10-08; their final outputs are kept as
   `MnBi2Te4.scf.out.<job>_killed_2026-10-08` in their test folders. Still to do: delete the PBE test folders and superseded PBE runs (`u_4.0/2028xx`, `205075`,
   `209025`, tests B/B2/D/E/C, `eval_209025_*`); make the LDA base's `tmp.tar.gz`.
6. **Report upstream?** The noncollinear GGA spurious state is a QE limitation more than a bug; the user decides.

**Script:** `calculations/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py` (settings block;
run on the cluster login node, ~4 GB, ~1.5 min for three runs). It reproduces every number in §8.
- **Reads:** QE's `charge-density.dat` (Fortran records: `gamma_only, ngm, nspin`; b; Miller indices; nspin complex
  G-space components).
- **Grid:** inverse-FFT to the 375³ grid.
- **Validation:** ∫n = 138.000, and ∫|m| reproduces QE's |m|abs.
- **LDA base for reference:** ∫|m_⊥| 0.077 μB, total-variation excess 0.0%, 0 jump links.
