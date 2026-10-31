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
- **YIO is not the same failure.** A fresh PBE start (segni-fixed binary) converges, but slowly: 80 iterations, about
  30 of them hovering near 1e-6. LDA converges in 21-22, with an exact AIAO order and a 0.68 eV gap on the SCF
  mesh. YIO's #17 stall is in restarts from the stalled 207546 state (§12).

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

**Eigenvalue symmetry of the LDA base NSCF (2026-10-08).** Checked with
`archive/diagnostics/2026-10-08_mbt_noncollinear_gga/eig_symmetry_qe.py`: the 12 operations QE found, 4×4×4
`nosym` NSCF (64 k), the 138 occupied bands, full-precision xml eigenvalues.
- **The 12 operations: worst |Δε| = 0.043 μeV** (median 0.043). The rotations move 48-62 of the 64 k, so the test
  is not vacuous.
  - Caveat: inversion makes ε(k) = ε(−k), so eigenvalues test the rotations but cannot test the time-reversal
    flags. Flipping a flag only composes the operation with inversion.
- **The anti-translation (layer equivalence)**, which QE can't impose, is tested via its consequence: inversion ×
  anti-translation forces every band to be twofold degenerate at every k. **Pair splitting: max 20.2 μeV, median
  0.59 μeV** over bands 1-138.
  - The residual is consistent with the 3e-5 μB layer-moment mismatch and the odd 375 FFT grid
    (§4.1, `nr1 = nr2 = nr3 = 384` would remove it).
- **For comparison:** YIO with the segni fix had a 46 μeV floor (4.8 meV before). The MBT LDA base is symmetric for
  every practical purpose.
- **Gap on the 4×4×4 NSCF mesh: 198 meV** (band 138 → 139, direct). That is in line with MBT's bulk gap of about
  0.2 eV.

**Wannierization of the LDA base (2026-10-09 01:00):** `209653bc889/w90_trial_sp_d/`, prep 210621, Wannier job 210622.

**Trial set: Bi s;p, Te s;p, Mn1 d, Mn2 d = 116 spinor WFs = bands 57-172, n_occ = 82.** The user's rule: keep every
occupied band unless it is truly flat core.

| bands | eV | character | widest band | treatment |
|---|---|---|---|---|
| 1-16 | −74 to −38 | Mn 3s/3p | 0.6 meV | flat core: excluded |
| 17-56 | −17 to −14 | Bi 5d | 7.3 meV | flat semicore: excluded |
| 57-80 | −5.9 to −0.96 | Te 5s + Bi 6s | 1752 meV | dispersive valence: **included** |
| 81-138 | 1.62 to 7.27 | Te 5p + occupied Mn 3d | 2139 meV | included |
| 139-172 | 7.47 to ~11 | Bi 6p, empty Mn 3d | | included (entangled at the top) |

Every earlier MBT Wannier run (u_0.0, u_3.0, u_4.0, T1⁺, T2⁻) used the 92-WF p-d set, which leaves out the dispersive
occupied s bands. That's the same kind of omission as YIO's old n_occ = 104 (bug #1).

**Windows:**
- **Outer −7..20 eV:** −7 sits in the 8.2 eV gap above Bi 5d.
- **Frozen −7..10 eV:** at most 112 bands at any k; band 173 starts at 10.97 eV.
- `dis_num_iter = num_iter = 0` (projection only).

**`.win` fixes against the u_4.0 template:**
- `num_bands` 240 → 300;
- the `write_ndgen_applied` typo;
- the empty `unit_cell_cart` and `kpoints` blocks, filled with the exact cell and the NSCF's 64-k list in order;
- `ang` units for the atoms.

**Binaries:** wannier90.x 26c8b506 and pw2wannier90.x 8f5368ba.

**First attempt (92-WF set, prep 210618 / Wannier 210619) failed.** Its prep job ran on **n138, which doesn't mount
`/mnt/wkdk`**. All jobs that read the shared filesystem now exclude n126/n246/n247/n138, including the held YIO prep,
Wannier and forces jobs (changed with qalter).
- **That qalter broke those 9 YIO jobs.** `qalter -l hostname=...` replaces the whole `-l` list and dropped
  `h_vmem`/`h_rt`. Every such job then failed at launch with "can't set additional group id: Cannot allocate
  memory" (Eqw, qacct `failed 38`). Fixed 2026-10-09 02:40 by restating all three resources.

**Result of the base Wannierization (210622, finished 04:59): not usable as it stands.**
- **Gap:** interpolated on the band path, −346 meV indirect and 11 meV smallest direct, against DFT 198 meV on the
  4×4×4 mesh. The model matches DFT on the mesh, but swings ±0.3 eV between mesh points.
- **Spreads:** total 1978 Å².
  - Mn d: 0.6-5.6 Å².
  - Te s/p: medians 1.5 / 3.6 Å², but some reach 54.
  - **Bi s/p: medians 26-30 Å², up to 56.**
- **Likely cause:** the outer window −7..20 eV holds 210-232 bands for 116 functions. Projection-only (no
  disentanglement iterations, by design) then draws the Bi 6p functions from a very large space, and 4×4×4 is coarse
  for this cell.
- **Window scan (wannier90 only, same .amn/.mmn/.eig),** `w90_trial_sp_d/window_scan_2026-10-09/`, jobs
  210646-210649:

  | variant | outer | frozen | bands in outer per k |
  |---|---|---|---|
  | o11.5_f10 | −7..11.5 | −7..10 | 116-122 |
  | o12_f10 | −7..12 | −7..10 | 118-124 |
  | o13_f10 | −7..13 | −7..10 | 128-136 |
  | o12_f9 | −7..12 | −7..9 | 118-124 |

  **Result: the windows don't matter.** All four variants give Bi median spreads of 28 Å² and path gaps of −412 to
  −433 meV.
- **Cause: wannier90's centres wrap along the long axis.**
  - **Mechanism:** with `num_iter = 0`, a centre comes from Im ln M_nn(k,b), which is defined only modulo 2π. On the
    4×4×4 mesh, the shortest axial b-vector, (b₁+b₂+b₃)/4, has a wrap period of 4c/3 = **109.1 Å**. Any function
    at crystal x > 2/3 (z > 54.6 Å: Bi 58.32/64.46, Te 69.82/76.39) has that shell's phase wrap.
  - **Consequences:**
    - The reported centre becomes an inconsistent average: Bi at 58.32 is reported at −50.4, against
      58.32 − 109.1 = −50.8.
    - Ω_D balloons to 1597 of 1978 Å².
    - Bi1 at 17.39 Å comes out fine (1.5 Å²).
  - **Why it breaks the interpolation:** the WS-image choice uses those centres, both wannier90's
    `use_ws_distance` and our production folding on the `.chk` centres (`wannier/wannier_io.py`). Wrong images
    give exact bands on the mesh and garbage between.
- **Fix under test:** job 210651, `w90_trial_sp_d_recentred/`.
  - The `.win` places the six atoms with z > c/2 at z − 81.852 (lattice vector a₁+a₂+a₃), so every projection site
    has \|x\| ≤ ½.
  - New `.amn` only; `.mmn`/`.eig` reused from 210622. QE positions are untouched; only the Wannier home cell
    moves.
  - Cross-check: 210650, the production windows with `use_ws_distance = false`.
- **If the fix holds:** base and all four modes must use the recentred `.win`, since the Julia `.amn` reads the
  base `.chk`/`.nnkp`. The mode chain 210633-210640 (held) was built from the old `.win` and must be rebuilt
  before release.
- **The fix holds (210651, 05:36).**

  | | 210622 (original) | 210651 (recentred) |
  |---|---|---|
  | total spread | 1978 Å² | **382 Å²** |
  | Ω_D | 1597 Å² | **0.46 Å²** |
  | Bi spreads, median/max | 28.8 / 56.3 Å² | **5.2 / 6.2 Å²** |
  | Te spreads, median/max | 3.4 / 53.9 Å² | 3.3 / 3.9 Å² |
  | path gap, indirect / smallest direct | −346 / 11 meV | **+104 / 185 meV** |

  - **Centres:** every centre sits on its atom; inversion partners at exactly ±z.
  - **Cross-check:** `use_ws_distance = false` on the original gives +83 meV, so the WS images chosen from wrong
    centres caused most of the damage.
  - **Open question:** the path VBM sits at (0.06,0,0) on Γ-L, 108 meV above Γ (VB 7.363 vs Γ 7.255; mesh VBM
    7.268 at T). Camel-back or interpolation wiggle? DFT bands at 9 k on Γ-L: job 210652
    (`soc_lda/u_4.0/bands_GL_check_2026-10-09/`).
  - **Held mode chain updated:** each `.win` has the base's image choice per atom (old kept as
    `MnBi2Te4.win.pre_recentre`), and `send_job.sh` reads the base `.chk`/`.nnkp` from `210651bc889`.
- **DFT check on Γ-L (job 210652, 9 k from the base state):**

  | x | DFT VB / CB | DFT gap | W90 VB / CB | W90 gap |
  |---|---|---|---|---|
  | 0 (Γ, mesh) | 7.2547 / 7.4730 | 218 | identical | 218 |
  | 0.04 | 7.3064 / 7.5402 | 234 | 7.3419 / 7.5305 | 189 |
  | 0.06 | 7.2957 / 7.5730 | 277 | 7.3625 / 7.5533 | 191 |
  | 0.10 | 7.2356 / 7.5978 | 362 | 7.3308 / 7.5813 | 250 |
  | 0.25 (mesh) | 7.0137 / 7.9829 | 969 | identical | 969 |
  | T (mesh) | 7.2681 / 7.4662 | 198 | identical | 198 |

  - **DFT has a real camel-back VB:** +52 meV above Γ near x = 0.04. The true indirect gap is about 160 meV, not
    the mesh value of 198.
  - **The recentred model overshoots it:** +56 meV at the peak and up to 95 meV error at x = 0.06-0.10. Its
    smallest direct gap near Γ is 30-45 meV too small.
  - **Cause:** the 4×4×4 mesh has nothing between Γ and 0.25, exactly where the axion curvature concentrates.
  - **Options for the user:** a denser NSCF/Wannier mesh (6×6×6 = 216 k, or 8×8×8 = 512 k), at roughly 3.4-8× the
    NSCF and Wannier cost.
- **General lesson:** in a long cell on a coarse mesh, put every projection site at |crystal x| ≤ ½. In this
  geometry, wannier90's num_iter = 0 centres wrap for x > 2/3. Nothing warns; only the spreads (Ω_D) and the
  interpolated gap show it.

### 11.1 T2⁻ modes on LDA: the first SCF iterations (2026-10-09 03:00)

**Setup:** `MnBi2Te4/phonon_lda/T2-/mode_{1..4}/Q_0p01/`, jobs 210629-210632.
- Each mode starts from the base state (`startingpot = 'file'`, base `tmp.tar.gz`) with `tprnfor`.
- Base forces run: 210628. It converged in 3 iterations to the base energy (to 1e-8 Ry); the Mn forces are 0.

**The mode SCFs don't start where they should:**

| | iteration 1 accuracy | after 6 iterations | \|m\|abs |
|---|---|---|---|
| base forces run | 2.5e-7 Ry | (converged at 3) | 9.99 |
| mode 1 (Mn) | 3.9 Ry | 0.20 Ry | 10.26 |
| mode 2 (Bi) | | 0.15 Ry | 10.27 |
| mode 3 (Te) | 0.23 Ry | 1.2e-3 Ry | 10.27 |
| mode 4 (Te) | | 1.3e-3 Ry | 10.27 |

- **All four settle at the same \|m\|abs, 10.27,** whichever atoms move.
  - That change is quadratic in any symmetry-breaking component, so a 0.01 Å displacement can't produce +0.28 μB
    linearly.
- **Mode 3's energy is already 0.16 mRy below the base,** where a harmonic displacement costs +½kQ². Its SCF is not
  yet converged (1.2e-3 Ry).
- **The Hubbard occupations stay mirrored** (Tr ns 4.9707/0.3182 on both Mn, as in the base).

**The symmetry groups differ.** In QE's labelling (T = time reversal):

| | operations without a translation | with f = (½,½,½) |
|---|---|---|
| base (12 ops) | E, C3±, i, S6±; C2′ and σd with T | none: QE finds the pure translation (½,½,½), calls the cell a supercell and drops every fractional translation |
| mode (12 ops) | E, C3±, σd with T | C2′ (unitary), i and S6± with T |

- Each mode operation is a product of a base operation and the anti-translation {E|½,½,½}′.
- The mode group shares 6 operations with the base run's group. It **lacks plain inversion**, so it allows density
  components that are odd under i but even under {i|½,½,½}′, the "T₂⁻ channel". The base run forbade them by
  symmetry.
  - Every LDA base run (atomic or random start) had them forbidden. The base state's stability against them has
    never been tested.

**Check: is the base density already mode-symmetric?** Yes.
`archive/diagnostics/2026-10-08_mbt_noncollinear_gga/density_symmetry_qe.py` applies each mode operation to
the base `charge-density.dat` in G-space:
- **Mode operations with fractional translations:** |Δn|/|n| = 6e-8 and |Δm|/|m| = 4e-7.
- **The rest:** ≤ 7e-9.
- **Validation:** the base's own operations give ≤ 7e-9.

So QE's G-space symmetrization of the starting density is a no-op. The bad first iteration is not a density
symmetrizer bug.

**Hypotheses, to be decided by energy vs force once modes 3/4 converge:**
- **(a) LDA+U has a soft electronic mode in the T₂⁻ channel.** The symmetric base would then be a saddle there, and
  every mode falls into the same broken-symmetry state. Prediction: E_mode < E_base, and ΔE ≠ −½F·Δr.
- **(b) A code path that only the Mn-swapping operations exercise** (PAW becsum, or the off-diagonal Hubbard ns with
  fractional translations and T). Prediction: erratic convergence.
- **(c) Benign startup transient.** Prediction: ΔE ≈ −½F·Δr > 0, and |m|abs returns toward 9.99.

**Mode 3 converged (04:11, 11 iterations): hypothesis (a).**

| | E (Ry) | E − E_base | \|m\|abs |
|---|---|---|---|
| base forces run 210628 | −6947.39550994 | 0 | 9.99 |
| mode 3 (converged) | −6947.39575631 | **−0.246 mRy** | 10.26 |
| mode 4 (4.7e-7 Ry) | −6947.39576768 | −0.258 | 10.27 |
| mode 1 (1.8e-4 Ry) | −6947.39577634 | −0.267 | 10.26 |
| mode 2 (1.3e-5 Ry) | −6947.39569372 | −0.184 | 10.27 |

- **Energy vs force for mode 3:** ΔE = **−0.246 mRy**, but −½(F_b + F_m)·Δr = **+0.070 mRy**. F_base·Δr = 6e-14.
  It is not on the base branch.
- **The forces break the base's inversion.** On the displaced Te (Ry/bohr):
  - Te 5.46 / 35.47: ±0.00589. The 5.46 force points *along* its displacement.
  - Te 46.38 / 76.39: ±0.00962.
  - Base-branch linear response would make all four equal in magnitude, since base inversion maps 35.47 ↔ 46.38.
- **The density (`density_symmetry_qe.py` on the mode 3 state):**
  - Invariant under its own 12 operations to 7e-9.
  - Breaks the base's plain inversion and C2′T by |Δn|/|n| = 5.2e-3 and **|Δm|/|m| = 1.65e-2**.
- **The state is saved** at `phonon_lda/T2-/mode_3/Q_0p01/scf_state_2026-10-09/` (from the node scratch, before the
  NSCF).

**Decisive test running: does the state exist as Q → 0?**
- **Run:** job 210645, mode 3 at Q = 0.002 Å, from that state.
- **Broken-state prediction:** E ≈ base − 0.25 mRy, \|m\|abs ≈ 10.27, inversion-odd \|Δm\|/\|m\| ≈ 1.6e-2.
- **Linear-response prediction:** \|m\|abs ≈ 9.99, and the odd part falls about 5×.
- **If it survives:** LDA+U MBT also has a lower state that breaks i and the anti-translation (keeping their
  product). The symmetric base, imposed via plain inversion in every base run, is then not the LDA+U ground state,
  which matters for θ = π.

**Result (210645, 07:00): it survives, and is lower at smaller Q.**

| | E − E_base | \|m\|abs |
|---|---|---|
| Q = 0.01 (mode 3) | −0.246 mRy | 10.26 |
| **Q = 0.002** | **−0.311 mRy** | 10.27 |

- **Extrapolation:** E(Q) = E₀ + ½k′Q² gives **E₀ ≈ −0.31 mRy (−4.3 meV/cell) at Q = 0**. The 10-iteration SCF
  converged cleanly to 5e-8 Ry.
- **What it would mean:** in QE LDA+U, MBT has an electronic state below the symmetric base that breaks i and
  {E|½,½,½}′ but keeps their product. That is spontaneous magnetoelectric order with θ unpinned from π, which is
  physically very surprising.
- **What it can't be:** a Hubbard-ns effect. The Mn d occupations are inversion-even by parity, and they stay
  mirrored.

**Remaining doubt: a QE symmetrization artifact** from the operations that combine a fractional translation with
time reversal, which no base run exercised.
- **Diagnostic:** `nosym` (no symmetrization at all) from the broken state at Q = 0.002. Job 210656 (dko64m4, 64 k,
  -nk 8, `Q_0p002_nosym_test_2026-10-09/`), ~1 h per iteration.
  - Stays at −0.31 mRy and \|m\|abs 10.27 ⇒ a genuine stationary LDA+U solution.
  - Drifts back ⇒ an artifact.
- **Input lesson:** the first submission (210655) had a second `nosym = .false.` further down the input; the
  namelist's last value wins. It was deleted; 210656 is verified "No symmetry found", 64 k.

**Preliminary result (08:50, 3 iterations of 210656): the broken state is a symmetrization artifact.**
- **Run:** started from the converged broken-state density (\|m\|abs 10.27).
- **First printed iterations:** \|m\|abs **10.00** (base 9.99), E = −6947.39536 → −6947.39544 Ry (base
  −6947.39551), SCF accuracy 2.3e-4 → 1.3e-4.
- **Hubbard:** Mn Tr ns 4.97067 / 0.31818, identical to the base.

Without symmetrization the state disappears at once. So the −0.25 to −0.31 mRy "state" exists only under QE's
symmetrization with the T₂⁻ group.
- **Where the bug must sit:** the mode-only operations, i.e. the fractional translation (½,½,½) combined with time
  reversal (i, S6), or with a Mn1↔Mn2 swap under a unitary C2. Candidates are sym_rho (nc magnetization), PAW
  becsum symmetrization, and the IBZ weights with these operations.
- **Not this bug:**
  - The starting density is invariant under them (`density_symmetry_qe.py`, 4e-7).
  - The YIO modes (ratios 1.000) are unaffected.
- **Consequence:** the T₂⁻ mode runs 210629-210632 (and the earlier PBE T₂⁻ runs) are not trustworthy as run. Next
  steps to decide: run the modes with `nosym` (~5× cost) or find and patch the QE routine. Final confirmation:
  210656 converging to ≈ base + 0.003 mRy.

**Confirmed (09:38, 210656 converged in 6 iterations):**

| | E − E_base | \|m\|abs | −½(F_b + F_m)·Δr | Te z-forces |
|---|---|---|---|---|
| symmetrized (210645) | −0.311 mRy | 10.27 | +0.0024 | ±0.00742, ±0.00804 |
| `nosym` (210656) | **+0.0028 mRy** | **9.99** | +0.0031 | ±0.00733, ±0.00815 |

- **`nosym` lands exactly on the base branch.** Its +0.0028 mRy matches the harmonic prediction.
- **The symmetrized run's forces are nearly right;** only its energy and \|m\|abs are off. Filed as
  [bug #19](../bugs/19_qe_sym_ftrans_trev_energy.md).
- **Correction to the mode 3 analysis above:** the unequal Te force pairs (0.0059 vs 0.0096) are **not** evidence
  of broken inversion. The `nosym` run has unequal pairs too, because the unrelaxed CIF base already has nonzero
  forces on Te (F_base ≈ 0.0077 Ry/bohr). The evidence was the energy (−0.25 mRy against +0.07 from the forces) and
  \|m\|abs.

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

**Result: both converge; LDA about 4× faster.** (Corrected 19:08: YP had looked stalled at 72 iterations.) All three runs find 48 operations (36 with fractional translations).

| run | converged | E (Ry) | Ir moments (all four) | off ⟨111⟩ | O \|m\| max | \|m\|abs | gap, 8-k SCF mesh |
|---|---|---|---|---|---|---|---|
| YP PBE | **80 iterations** (hovering 2e-7 to 1.5e-5 from about 45 to 79) | −4624.80821 | 0.3102 μB, equal | 0.000° | 0.0132 | 3.39 | — |
| YL1 LDA (`rel-pbe` PAW) | **21 iterations** | −4501.58867 | 0.2830 μB, equal | 0.000° | 0.0127 | 3.21 | 0.683 eV |
| YL2 LDA (`rel-pz` PAW) | **22 iterations** | −4499.85691 | 0.2829 μB, equal | 0.000° | 0.0127 | 3.21 | 0.685 eV |

**Reading it:**
- **No spurious state in YIO.** YIO keeps all 48 operations, enforced every iteration, so there is no layer
  asymmetry to fall into, and a fresh PBE start converges (YP, 80 iterations).
- **The |m|abs trend was the patched functional's own value, not a drift.** The patched restarts (T1c, R1-R3,
  3.35-3.37) were approaching YP's converged 3.39. The unpatched T1's 3.32 belongs to the buggy functional.
- **So #17 is narrower than thought:** slow PBE convergence, and stalls when restarting from 207546's state. YP is
  the first converged base on the segni-fixed binary and could replace 207546 for production (its NSCF has not been
  run).
- **LDA converges about 4× faster** (21-22 iterations), consistent with its noncollinear form being exact.
- **LDA keeps the AIAO insulator at U = 3.** Four equal moments exactly along ⟨111⟩, 9% smaller than PBE, and a
  0.68 eV gap on the SCF mesh. For reference, PBE's 512-k NSCF gap was about 0.42 eV; the meshes differ.
- **θ stays 0.** A trivial AIAO phase is what the project needs for YIO (memory `y2ir2o7-axion-phonon-goal`).

**Decision test (started 2026-10-08 21:45).** NSCFs from identical production inputs (512 k, nosym, 300 bands)
compare the eigenvalue symmetry floor and the full-mesh gap. Each folder has an `eig_symmetry_qe.py` with its
settings (N_BANDS 164); results go here when done.
- **PBE:** YP → `YP_pbe/nscf_from_210525/` (job 210587).
- **LDA:** YL2 → `YL2_lda_relpz/nscf_from_210527/` (job 210588).

**Literature:** the founding Y2Ir2O7 calculations (Wan, Turner, Vishwanath, Savrasov, PRB 83, 205101 (2011)) and later
U scans use LSDA+U+SO. There, AIAO is a trivial insulator above U ≈ 1.8-2 eV. So LDA is the YIO literature's own
method, and U = 3 sits in the trivial phase.

**Caveat:** LDA+U at U = 3 is not PBE+U at U = 3. The YIO AIAO moment is non-monotonic in U, so if YIO moves to
LDA, check the gap on the full NSCF mesh first. **Not decided:** whether YIO switches. Since PBE does converge, there is no forcing reason; LDA would only buy
speed and robustness.

### 12.1 YIO PBE production set rebuilt on base YP (2026-10-08 22:05)

**Nothing had been Wannierized from a converged base.**
- **Base:** the only base Wannierization was 207546, the stalled SCF, about 1 mRy above YP's branch.
- **Modes 1 and 3:** stage 1 only (`.mmn`/`.eig`, in place).
- **Mode 2:** its stage 1 never ran; job 209381 died on a node without `/mnt/wkdk`.
- **No mode stage 2** (rigid-shift `.amn` + `wannier90`) had ever run.

**Branch check of the existing mode SCFs against YP** (all segni-fixed). Calibrated with the old same-binary pair
(192344 base and its modes):

| | old pair, E − base | current, E − YP | Ir \|m\| (YP 0.31021) |
|---|---|---|---|
| mode_1 (207547) | +1.674 mRy (192345) | +1.668 mRy | 0.31019 ✓ |
| mode_2 (209040, restart chain) | +1.807 mRy (192346) | +1.947 mRy | 0.30653 ✗ |
| mode_3 (207549) | unusable (192347 landed elsewhere) | +0.571 mRy | 0.30843 ? |

**Rebuild:** every mode SCF was rerun **from YP's converged state**, so base and modes share one branch.
- **Start:** `startingpot = 'file'`, with YP's density, ns and becsum in `tmp.tar.gz`.
- **Settings:** YP's (Davidson, plain β 0.1). The NSCF runs with `-nk 8` in the same job.
- **Folders and jobs:** `mode_X/fromYP_2026-10-08/`, jobs 210596-210598, all pw.x e270b05.
- **Read before Wannierizing:** each mode's E − E(YP) should match the calibration column, and Ir |m| should be about
  0.310.

**Held Wannier chain** (`phonon/setup_w90_chain_YP.sh`):
- **Rule:** every job that runs pw2wannier90 stages its NSCF save as `tmp.tar.gz` in its own submit folder and untars
  it on the node. The Julia `.amn` code reads the base save untarred on the shared filesystem (CLAUDE.md cluster
  practice).

| step | jobs | waits on |
|---|---|---|
| base YP: tar the save into `w90_trial02_base/`, then gunzip the save in place | 210599 | NSCF 210589 |
| base YP: `-pp` → pw2wannier90 `-nk 8` → wannier90 `-np 1` | 210600 | 210599 |
| modes: tar each NSCF save into `w90_trial02_no-displacement/` | 210601 / 210603 / 210605 | mode reruns |
| modes: `-pp` → Julia `generate_amn.jl --no-displacement` from YP's `.chk`/`.nnkp` + save → pw2wannier90 (mmn, eig) → wannier90 | 210602 / 210604 / 210606 | own prep and base 210600 |

The jobs delete the unpacked `tmp` before copy-back. Superseded: the base Wannierization of 207546 and the in-place
mode stage-1 folders.

### 12.2 Branch check by forces (2026-10-08 23:00)

**The check.** E_mode − E_base ≈ −½ (F_base + F_mode)·Δr, exact through third order, with Δr the Γ₂⁻ displacement.
- **The base term drops out.** The base keeps inversion and Γ₂⁻ is odd, so F_base·Δr = 0, and the check is
  ΔE ≈ −½ F_mode·Δr.
- **No reference needed.** It doesn't rely on the old buggy-binary pair. A mismatch means the mode is in a different
  self-consistent state from the base.

**Hubbard forces work for this setup.** `force_hub.f90` refuses only these: projectors other than atomic or
ortho-atomic, Hund's J, noncollinear with background states, and gamma_only with ortho-atomic. None apply. The run
printed a "Hubbard contrib. to forces" term.

**Base forces:** `YP_pbe/forces_2026-10-08`, job 210609. A short SCF from YP's state: 9 iterations to 9.4e-10 Ry,
E = −4624.80821431, 0.8 μRy from YP.
- **Pattern:** only the 12 O(48f) carry force, 0.01074 Ry/bohr (0.28 eV/Å) each along the free 48f coordinate; Y, Ir
  and O(8b) are zero to 1e-8. The base is the unrelaxed experimental structure.
- **Hubbard part:** −0.00176 Ry/bohr of each. SCF correction 1e-6.

**Mode forces:** jobs 210611/210613/210615, each held on its rerun via a 1-slot prep job that tars the converged
state.

### 12.3 Mode reruns from YP: results (2026-10-09 01:05)

| mode | job | SCF from YP | E − E(YP) | calibration | Ir \|m\| (YP 0.31021) |
|---|---|---|---|---|---|
| 1 | 210596 | 13 iterations | +1.666 mRy | +1.674 | 0.31032 ✓ |
| 3 | 210598 | 12 iterations | +0.563 mRy | — (force check pending) | 0.30928 |
| 2 | 210597 | **not converged**: 173 iterations, 4e-5 to 2.7e-4 Ry | — | +1.807 | 0.30689 |

**Mode 2:** plain β 0.1 again gave the period-2 Hubbard-ns oscillation this mode showed on 2026-10-06.
- **Stopped cleanly** with EXIT, with the NSCF input held back.
- **Restarted** as **210623** (`fromYP_2026-10-08/restart_localTF/`) from its own final state, still descended from YP:
  local-TF β 0.05, Davidson, the settings that converged it before.
- **Chain rebuilt:** Wannier prep 210624 → Wannier 210625 (held on base 210600); forces prep 210626 → forces 210627.
  The old 210603/210604/210612/210613 were removed.

### 12.4 YIO PBE vs LDA: eigenvalue symmetry and gap (2026-10-09 02:20)

NSCFs from identical production inputs (512 k, nosym, 300 bands, `-nk 8`): YP (PBE) → 210589, YL2 (LDA, `rel-pz`) →
210590. Checked with `eig_symmetry_qe.py` over all 48 operations and the 164 occupied bands; the operations move
448-510 of the 512 k.

| | worst \|Δε\| | median | gap on the full mesh |
|---|---|---|---|
| **PBE (YP)** | **137 μeV** | 112 μeV | 774 meV indirect (1054 direct) |
| **LDA (YL2)** | **1.18 μeV** | 1.12 μeV | 685 meV indirect (957 direct) |

**Reading it:**
- **LDA is ~100× more symmetric.** PBE's floor is also above the 46 μeV of the segni-fix validation run C2: YP is a
  different, fully converged SCF.
- **Both stay AIAO insulators.**
- **This is the strongest evidence so far for moving YIO to LDA as well;** the user decides.

### 12.5 YIO modes 1 and 3: energy vs force (2026-10-09 02:55)

Forces at each converged mode state (SCF from the mode's NSCF state, `tprnfor`, `conv_thr` 1e-9; jobs 210611,
210615), checked against the base YP forces (210609) with ΔE ≈ −½(F_base + F_mode)·Δr:

| mode | E − E_base | −½(F_b + F_m)·Δr | ratio | Ir \|m\| |
|---|---|---|---|---|
| 1 | 1.6667 mRy | 1.6674 mRy | 1.000 | 0.3103 |
| 3 | 0.5641 mRy | 0.5642 mRy | 1.000 | 0.3093 |

- Both are **on the base's branch.** The midpoint rule leaves an O(Δr³) error, so agreement to 4×10⁻⁴ is all that
  can be asked.
- F_base·Δr = 0 to 1e-13, as inversion requires for an odd mode.
- Each forces SCF converged in 9-10 iterations from the stored state.

### 12.6 YIO mode 2 (Ir along local [111]) is not cleanly on the base branch (2026-10-09 04:15)

**209040 was never converged.** Its "39 iterations" was one lucky dip, from ~3e-5 to 7e-8 Ry, inside an
oscillation.
- **Forces from its state with conv_thr 1e-9 (210641):** oscillated at 5e-5 to 1.7e-4 Ry for 126 iterations.
  Killed.
- **Forces on the oscillating state (210643):** `electron_maxstep = 40`, `scf_must_converge = .false.`. QE's SCF
  correction to the forces is only 1.7e-5 Ry/bohr.

| | E − E_base | −½(F_b + F_m)·Δr | ratio | Ir \|m\| |
|---|---|---|---|---|
| mode 1 | 1.6667 | 1.6674 mRy | 1.000 | 0.3103 |
| mode 3 | 0.5641 | 0.5642 mRy | 1.000 | 0.3093 |
| **mode 2** | **2.12-2.22** (oscillating) | **1.948 mRy** | **1.09-1.14** | **0.3065-0.3072** |

**Why this is suspect:**
- **No cubic term.** E(Q) is even for an odd mode, so the midpoint formula has no cubic error; modes 1 and 3 hit
  1.000.
- **The Ir moment is too low.** Each Ir sits on an inversion center, so \|m_Ir\| is even in Q. A 1% drop at
  0.017 Å Ir displacement would need an enormous quadratic coefficient.
- **The 209040 coincidence.** Its lucky-dip energy (+1.947) matches the force estimate. That may be a coincidence.

**Test:** half amplitude, Q = 0.005 Å, job 210644, from the YP state. Same branch ⇒ Ir moment drop ~4× smaller,
ΔE ≈ ¼ × 1.95 mRy, ratio ≈ 1.

**Result (210644, 05:43): mode 2 is on the base branch; the Q = 0.01 SCF just doesn't converge.**

| | Q = 0.01 (oscillating, ~1e-4 Ry) | Q = 0.005 (1e-8 to 3e-7 Ry) |
|---|---|---|
| E − E_base | 2.12-2.22 mRy | 0.508 ± 0.006 mRy |
| −½(F_b + F_m)·Δr | 1.948 mRy | 0.501 mRy |
| ratio | 1.09-1.14 | **1.015** |
| Ir \|m\| drop | 1.0% | **0.24%** |

- **The Ir moment drop is quadratic:** 0.24% × 4 = 0.96%. The Ir moves along its own moment axis, so the
  quadratic coefficient is large but physical.
- **The energies are consistent:** 4 × 0.508 = 2.03 mRy agrees with the full-amplitude force estimate (1.95) to
  within quartic terms. The ~0.2 mRy excess at Q = 0.01 is the energy error of the oscillating SCF; its forces were
  fine (SCF correction 1.7e-5).
- **Production candidate:** mode 2 at Q = 0.005, SCF (conv_thr 1e-7, local-TF 0.05, from YP) + NSCF, job 210654,
  `mode_2/Q_0p005_fromYP_2026-10-09/`. No Wannier chain attached. Using a different amplitude for mode 2 is the
  user's decision.
- **`tmp_bands/` duplicates:** the template's `cp -a tmp/Y2Ir2O7.save tmp_bands/` left about 48 GB copies of the
  save in 210589, 210590, 210596, 210597 and 210598. They await deletion with the user's OK; the line is commented
  out in the new job.

**The Q = 0.005 candidate is done (210654, 10:24).**
- **SCF:** converged for real in 52 iterations (conv_thr 1e-7, no forced stop). E − E_base = 0.5078 mRy;
  −½(F_b + F_m)·Δr = 0.4999 mRy (ratio 1.016). Ir \|m\| 0.30953.
- **NSCF:** JOB DONE; `tmp/` is 105 GB.
- **Ready:** usable as production mode 2 if the user accepts Q = 0.005. It would need its own Wannier chain:
  prep (tar), then w90 with the Julia `.amn` from the base.

**Hold:** the mode 2 chain stays on hold (210624-210627, qhold) until this is settled. Restart 210623 stopped at
120 iterations unconverged; its NSCF is not production.

### 12.7 YIO base Wannierization on YP (2026-10-09 10:39)

Job 210600 (`YP_pbe/nscf_from_210525/w90_trial02_base/`): trial_02, 176 WFs, n_occ 156, projection only.
- **Spreads:** total 188.6 Å², Ω_D 0.015 Å²; per function 0.55-1.88 Å² (median 1.11).
- **Gap:** interpolated on the path, 774.4 meV indirect and 1054.3 meV smallest direct, identical to the DFT 8×8×8
  NSCF (the extrema lie on mesh points).
- **Chain:** the mode 1 and 3 Wannier jobs (210602, 210606) started on its completion. Mode 2 is held (Q-amplitude
  decision).

### 12.8 YIO base curvature symmetry: the external term breaks the threefold axes (2026-10-10)

`base_full_group_run.py` (48 operations of m-3m′, 7 generic k, `transl_inv_full_from_chk`, WS tolerance 1e-3 Å)
on the YP base Wannier run 210600. Copied to
`data/Y2Ir2O7/base/soc/u_3.0/output/YP_210525/w90_trial02_base/210600bc889/`.
- **PythTB:** needs the `external` branch. The local checkout is on `dev`, so it was run from a worktree of
  `external`.

| medians (max) | old base 192344 (pre-segni) | YP 210600 |
|---|---|---|
| internal | 1.06% (1.41) | 0.62% (0.76) |
| cross | 18.7% (25.4) | 8.9% (13.3) |
| external | 44% (57) | **86% (104)** |
| total | 2.0% (2.5) | **4.0% (5.0)** |

**Control:** the same script on the old data (plain `_hr.dat`/`_r.dat`) reproduces the October table to every
digit. So the code is unchanged, and the new data really differs.

**Pattern:** the new external term is exactly symmetric (≤0.7%) under mmm: inversion, C2 about x, y, z, and the
three mirrors.
- About 14% off under C4/S4 about y and C2 about [101].
- **76-104% off under every threefold operation**, and C4 about x and z.
- An axis-respecting defect in the position matrix; the old run broke all operations about equally (~45%).

**Ruled out:** identical `.win`, cell, atoms and projections; centres agree to 2e-4 Å; identical b-shells (8 along
⟨111⟩); same wannier90.x (26c8b506); same `_hr.dat` degeneracies.

**Difference found:**
- **Old:** pw2wannier90 **v7.6** (4581573f), no pools.
- **New:** pw2wannier90 **v8.0dev** (8f5368ba, built 09-28: the copy in ~/q-e, ~/q-e-segnifix and ~/q-e-segnireg)
  with **`-nk 8` pools**.
- The MBT LDA base Wannierization (210622/210651) used the same 8.0dev binary with `-nk 8`.

**Isolation runs:** both are the same tarball, `.win` and wannier90 without pools.
- 212512: v7.6.
- 212513: v8.0dev.

**Result (2026-10-10 21:50): pw2wannier90 is not the cause.** Compared the raw outputs directly; no need to wait for
the wannier90 step.
- `.amn`, 212512 (v7.6, no pools) vs 212513 (v8.0dev, no pools): max |Δ| 1e-12, relative 7e-14. That is print
  rounding.
- `.amn`, 212512 vs 210600 (v8.0dev, `-nk 8`): max |Δ| 1e-12, at all 512 k, so pools on other k are fine too.
- `.mmn`, 212512 vs 210600, first 3.0M elements (k = 1-5, including neighbours that live on other pools):
  max |Δ| 1e-12.
- 210525 found all 48 operations (36 with fractional translations), so it isn't a reduced SCF group either.
- 212512/212513 were killed at 22:13 once this was settled, and their test folders deleted (the result is above).

**Cause (2026-10-10 22:10): the frozen-window edge cuts through a k-star.** In the YP NSCF, band 178 lies within
0.16 meV of `dis_froz_max = 15.0` on a whole star of symmetry-equivalent k-points (eigenvalues in eV relative to 15.0):

| k (cart, 2π/a) | e − 15.0 | frozen? |
|---|---|---|
| (0.18, 0.71, ±0.18), (0.53, 0, 0.53), (−0.18, 0.71, 0.18) | −0.01 to 0.00 meV | yes |
| (0.71, 0.18, ±0.18), (0.18, ±0.18, 0.71), (0.53, 0.53, 0), (0, 0.53, 0.53) | +0.01 to +0.16 meV | no |

- **Why it matches the pattern:** sign flips (mmm) keep each k in its row, so the frozen subspace is mmm-symmetric. Axis
  permutations (the 3-fold operations, C4 about x and z) move k between rows, which breaks them; C4 about y keeps the
  rows, hence only ~14%.
- **Why the old base was fine:** same `.win`, but its eigenvalues didn't sit on the edge.
- **The DFT isn't at fault:** the 137 μeV symmetry floor only decides which side of the edge each star member lands on.
- **Outer edge:** `dis_win_max = 25` had bands 265/266 at +0.44 meV and others at −4.0 meV. No star was split, but the
  clearance was thin.

**Rule:** put every disentanglement edge in an eigenvalue gap of the actual NSCF, with clearance far above its
symmetry floor. Check this for every new NSCF, since the edges are absolute energies.

**Runs (2026-10-10 22:15):**
- **212616:** wannier90-only rerun of 210600 in `YP_pbe/nscf_from_210525/w90_trial02_base_frozfix_2026-10-10/`.
  - `dis_froz_max = 14.858`: gap 14.838-14.879, 41 meV wide, 168-170 bands frozen.
  - `dis_win_max = 25.155`: gap 53 meV wide.
  - The windows only enter wannier90.x, so the `.amn`/`.mmn` are hard links to 210600's.
- **LDA base YL2:** prep 212614, then Wannier job 212615, in `YL2_lda_relpz/nscf_from_210527/w90_trial02_base/`.
  - `dis_froz_max = 14.396`: gap 14.374-14.417, 168-170 bands.
  - `dis_win_max = 24.522`: gap 40 meV wide.
  - LDA needed new windows anyway: 15.0 would freeze up to 180 > 176 bands.

**Result of 212616 (2026-10-11 00:30): the window edge was the 86%.** Same `.amn`/`.mmn`, edges moved into gaps:

| medians (max) | 210600, edge on a band | 212616, edges in gaps | old base 192344 |
|---|---|---|---|
| internal | 0.62% (0.76) | 0.81% (0.97) | 1.06% (1.41) |
| cross | 8.9% (13.3) | 3.5% (3.6) | 18.7% (25.4) |
| external | 86% (104) | **16.4% (17.2)** | 44% (57) |
| total | 4.0% (5.0) | **0.67% (0.84)** | 2.0% (2.5) |

- **Best YIO total yet.** Log: `archive/diagnostics/2026-10-05_ws_tolerance_window/base_full_group_YP_frozfix_212616.log`.
- **What is left is axial:** the 16 operations of a 4/mmm subgroup stay ≤ 0.6%; the other 32 sit at 16-17%.

**What is left is not the Wannier step.**
- With `dis_num_iter = 0`, the 6 non-frozen states at each k are the intersection of the trial subspace (176) with
  the non-frozen block, so they are fixed by dimension counting: eigenvalue exactly 1, then a gap of ≥ 0.73 to the
  next.
- The smallest singular value of A is 0.60. No near-tie, so the gauge is fully determined by the data.

**The PBE NSCF Hamiltonian itself has an axis.** `eig_symmetry_qe.py` on 210589, per operation:
- **Exact (0.000 μeV):** only E, C2 about x, and the mirrors normal to y and z.
- **Broken:** inversion by 45 μeV, and the others by 48-137 μeV.
- So cartesian x is special in PBE; LDA (210590) has 1.2 μeV.
- **Not in the smooth real-space GGA:** `compute_rho` with `lsign = .false.` is isotropic (segni = 1).
- **Prime suspect (not checked):** the PAW one-centre GGA, integrated on a θ-φ product grid that is not cubic.

Practical test: the LDA base Wannierization (212615) with the same windows logic.

### 12.9 YIO base curvature symmetry: LDA, WS tolerance, pinned centres (2026-10-11)

**Runs:**
- **LDA base Wannierization 212615** (`YL2_lda_relpz/nscf_from_210527/w90_trial02_base/`): Ω_D 0.016 Å², total
  spread 194 Å². Local copy: `data/Y2Ir2O7/base/soc/u_3.0/output/YL2_210527/w90_trial02_base/212615bc889/`.
- **Tables:** `transl_inv_full_from_chk`, 48 operations, 7 generic k. Logs `base_full_group_YL2_212615*.log` and
  `base_full_group_YP_frozfix_212616.log` in `archive/diagnostics/2026-10-05_ws_tolerance_window/`.

| medians (max), % | internal | cross | external | total |
|---|---|---|---|---|
| PBE 212616, tol 1e-3 | 0.81 (0.97) | 3.5 (3.6) | 16.4 (17.2) | 0.67 (0.84) |
| LDA 212615, tol 1e-3 | 0.58 (0.82) | 2.4 (2.8) | 17.8 (18.1) | 0.68 (0.84) |
| **LDA, tol 5e-3** | 0.54 (0.79) | 1.06 (1.20) | **0.13 (0.18)** | **0.52 (0.82)** |
| LDA, atom-pinned positions, tol 1e-3 | 0.56 (0.84) | 0.84 (1.28) | 0.14 (0.17) | 0.53 (0.85) |
| LDA, pinned + full-precision H(R) from `.chk` | 0.70 (0.90) | 1.04 (1.19) | 0.16 (0.22) | 0.67 (0.92) |

**The pattern (PBE and LDA alike):** the 16 operations of 4/mmm about **z** (C4z, C2 about x, y, [110] and [1−10],
and their mirrors) are exact to ≤ 0.015% in every part. The other 32 break.
- Those 16 are exactly the operations that keep spin quantized along z, the spin axis of the spinor projections.
- The other 32 map each Wannier function into a combination of its same-site spin and orbital partners.
- **Individual centres:** they can't be covariant under the spin-mixing operations. Spin partners on one atom
  differ by up to 8.6e-4 Å (median 1.2e-4). This is the "5.7e-4 Å centre noise" of bug #13.

**The external term (18%) was WS tie-breaking on those centres, and is fixed.**
- **Two fixes, same result:** tolerance 5e-3 (inside the physical window of bug #13) or atom-pinned positions in
  the H/A folds and τ. Both give external 0.13-0.14%.
- **Correction:** the 2026-10-05 conclusion that the rotations' floor "is not tie-breaking" was made while the DFT
  bugs (#7, #16) dominated. On clean DFT it *is* tie-breaking.
- **Use 5e-3.**

**What remains is internal (H only), 0.5-0.9%, on the same 32 operations.**
- **Not the fold or τ:** pinned positions don't change it.
- **Not `_hr.dat` print rounding:** H(R) rebuilt at full precision from `.chk` + `.eig` (class sums match `_hr.dat`
  to 3.4e-6 eV; this `_hr.dat` is `write_ndegen_applied`: each R holds its pair-wise WS fraction, all deg = 1)
  gives 0.70%.
- **Not the disentanglement:** the selection is unique (§12.8).
- **Hypothesis (untested):** NSCF eigenvector accuracy.
  - QE's atomic starting spinors are quantized along z, so Davidson errors would be identical at k-points related
    by z-preserving operations and different otherwise.
  - The LDA eigenvalues are symmetric to 1.2 μeV ≈ 9e-8 Ry, at the level of ethr 1e-8 Ry with `diago_full_acc`.
    Eigenvector errors are first order where eigenvalue errors are second order, so ~1e-4 residuals are possible.
  - **Test:** NSCF with `diago_thr_init` 1e-11 to 1e-12, then rerun pw2wannier90 + wannier90 (~9 h).

## 13. Open

1. **LDA U for MBT: decided 2026-10-08, U = 4.0 for now** (user).
   - **Literature:** published DFT+U values for MBT span about 3-5.34 eV, almost all PBE+U.
   - **Benchmark:** the DMC-optimal LDA+U value is 4.4 eV.
   - **Revisit** if a benchmark comparison is needed.
2. **MBT T2⁻ modes on LDA:** exact cell + Q·mode, same settings. These modes break the anti-translation physically,
   so in PBE no symmetric start could have protected them. That's the main reason for the switch.
3. **YIO:** PBE converges from scratch (YP, 210525), and LDA converges 4× faster. Choose: (a) make YP the PBE
   production base (run its NSCF, then rerun the modes from fresh starts as needed), or (b) switch YIO to LDA.
   Check the full-mesh gap first if switching.
4. **A QE patch to enforce the anti-translation** (project charge onto even h+k+l, magnetization onto odd h+k+l, and
   pair-average ns/becsum). Not needed for LDA; would make PBE usable for the base only.
5. **Cleanup:** C2 and B2 were killed 2026-10-08; their final outputs are kept as
   `MnBi2Te4.scf.out.<job>_killed_2026-10-08` in their test folders. Still to do: delete the PBE test folders and superseded PBE runs (`u_4.0/2028xx`, `205075`,
   `209025`, tests B/B2/D/E/C, `eval_209025_*`); make the LDA base's `tmp.tar.gz`.
6. **Report upstream?** The noncollinear GGA spurious state is a QE limitation more than a bug; the user decides.

**Script:** `archive/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py` (settings block;
run on the cluster login node, ~4 GB, ~1.5 min for three runs). It reproduces every number in §8.
- **Reads:** QE's `charge-density.dat` (Fortran records: `gamma_only, ngm, nspin`; b; Miller indices; nspin complex
  G-space components).
- **Grid:** inverse-FFT to the 375³ grid.
- **Validation:** ∫n = 138.000, and ∫|m| reproduces QE's |m|abs.
- **LDA base for reference:** ∫|m_⊥| 0.077 μB, total-variation excess 0.0%, 0 jump links.

## Cluster cleanup: PBE U = 0 runs deleted (2026-10-19)

`MnBi2Te4/base/soc/u_0.0/` job folders 127470, 130655 (with Wannier 136914, 137698), 150301, 160618 and 160635 deleted
by the user (no local copy; 150300's inputs/outputs are in `data/MnBi2Te4/archive/base/pbe_d3_U_scan/U_0.0/`). All were
QE 7.2 PBE + D3, U = 0 (HUBBARD card commented out), rounded cell. Last observation from their outputs, relevant to #18:
- 127470, 130655 (rounded-bohr cell, 4 operations): converged symmetric, layer pairs mirror to 1e-4 μB, |m|abs 9.90/9.88.
- 150300 (rounded-Å cell, 12 operations, the paper's old dθ base): Bi/Te pairs mirror only to 2-3% (Te 0.0488 vs
  −0.0495 μB), |m|abs 10.37: partly drifted toward the layer-asymmetric state.
- 150301 (T₂⁻, Mn1 +0.01 Å, conv_thr 1e-9): unconverged after 400 iterations (4.5 days); NSCF + .amn/.mmn made anyway.
- 160618 killed after 3 iterations (|m|abs 7.43); 160635 the QE 7.5 crash after 2 (|m|abs 14.71; bug #12).
