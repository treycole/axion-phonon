# 2026-10-15 — MBT soft T₂⁻ mode resolved: rigid interlayer mode at 18.0 cm⁻¹

> Superseded for the per-gap stiffness conversion only by [log 10-16](2026-10-16_mbt_interlayer_stiffness.md): K_gap = k/8, not k/4; the computed phonon frequencies are unchanged.

Follows [2026-10-14_mbt_t2m_dynmat](2026-10-14_mbt_t2m_dynmat.md), where the lowest T₂⁻ mode came out at −20 ± 13 cm⁻¹
(unresolved). Script: `archive/diagnostics/2026-10-14_mbt_t2m_dynmat/t2m_dynmat.py`, now with the `RIGID`
setting. The same numbers, and each method explained, are in `analysis/phonon_frequencies.ipynb`.

## Runs

| job | Q (Å per atom) | iterations | E (Ry) | ΔE vs base 210628 | k points | sym ops | \|m\|abs | Mn1 / Mn2 m_z (μB) |
|---|---|---|---|---|---|---|---|---|
| 213463 | 0.03 | 18 | −6947.39543061 | +0.0793 mRy (1.079 meV) | 16 | 12 | 9.99 | +4.6530 / −4.6530 |
| 213464 | 0.05 | 23 | −6947.39528934 | +0.2206 mRy (3.001 meV) | 16 | 12 | 9.99 | +4.6527 / −4.6527 |

- Displacement: SL1 +Q along z, SL2 −Q, every atom. That is exactly (1, −1, 1, −1) in the SMODES basis (least-squares
  residual 1e-24).
- pw.x md5 `db4b99b4` (checked in both job folders). conv_thr 1e-9; SCF correction 1.2e-4 / 6.8e-5 Ry/bohr.
- Total magnetization 0. The Mn layer moments mirror exactly.
- Data: `data/MnBi2Te4/archive/phonon/T2-/rigid_layer/Q_{0.03A,0.05A}/lda/scf_{213463,213464}/` (scf.in/out, job log,
  NOTE, send_job.sh). Cluster: `MnBi2Te4/phonon_lda/T2-/rigid_layer/Q_{0p03,0p05}_2026-10-14/`.

## Result

- **The mode is stable.** ΔE > 0, and ΔE(0.05)/ΔE(0.03) = 2.782 against (0.05/0.03)² = 2.778, so E is quadratic
  to 0.1%.
- The base's residual forces have zero projection on this pattern (v·F(0) = −6e-14 eV/Å), as symmetry requires, so
  there is no linear term.
- Stiffness along the rigid vector c (eV/Å² per unit Q):

| from | Q = 0.03 | Q = 0.05 |
|---|---|---|
| energy, 2ΔE/Q² | 2.3985 | 2.4011 |
| forces, −c·[F(Q) − F(0)]/Q | 2.4097 | 2.4076 |
| Q = 0.01 basis runs (log 10-14) | −2.4275 | |

- The Rayleigh quotient on c alone gives 18.2 cm⁻¹.
- With the measured row Φ·c substituted into the 4×4 matrix (basis c plus three directions M-orthogonal to it,
  the c row/column taken from the mean of the two rigid runs):

| ω (cm⁻¹) | meV | THz | Mn / Bi / Te1 / Te2 share of \|Δr\|² | change vs 10-14 |
|---|---|---|---|---|
| **18.0** | 2.23 | 0.539 | 0.14 / 0.27 / 0.28 / 0.31 | was −20 ± 13 |
| 77.7 | 9.64 | 2.330 | 0.32 / 0.06 / 0.22 / 0.41 | −0.2 |
| 124.0 | 15.37 | 3.717 | 0.79 / 0.04 / 0.04 / 0.13 (Mn) | −0.1 |
| 139.8 | 17.33 | 4.191 | 0.39 / 0.21 / 0.39 / 0.01 | −0.1 |

- Overlap of the 18.0 cm⁻¹ eigenvector with the rigid vector: 0.9995. The other three barely move.
- From the two amplitudes and the E-vs-F agreement, the soft-mode frequency is good to ~0.1 cm⁻¹ at this cell and
  functional.

## Interpretation and caveats

- **Interlayer force constant.** Each cell has two vdW gaps, and the relative displacement across each is 2Q. So
  K_gap = k/4 = 0.60 eV/Å² per in-plane cell (area 16.26 Å²), about 5.9 × 10¹⁹ N/m³ per unit area. That is the
  usual vdW range, somewhat below graphite or MoS₂. *Not yet compared with measured MnBi₂Te₄ breathing modes.*
- **The Q = 0.01 runs had the wrong sign along c.** Their error in c·Φ·c was 4.8 eV/Å², about 2.5% of Σ|Φ_ij|.
  The 10-14 noise model predicted only 13 cm⁻¹ of spread.
  - Φ_old·c − Φ_meas·c = (−3.2, +1.5, −0.3, −0.2) eV/Å² with the raw Φ, and (−2.1, +1.6, −1.0, +0.1) with its
    symmetric part. The largest error sits in the Mn row; the Mn SCF was the one whose E vs F ratio was 1.24.
  - Φ_MnMn is 6.67 eV/Å² from the energy against 5.36 from the forces. The energy has no linear or cubic term
    here, so it is the better number. Raising the Mn diagonal of the method-4 matrix by that difference moves the
    modes to 21.7 / 80.4 / 130.5 / 143.4 cm⁻¹, so the Mn-mode error is ~6 cm⁻¹, not 2. That is a sensitivity
    test, not a result.
  - So the error bars on the three optical frequencies (1-3 cm⁻¹ in log 10-14) are probably too small. The
    124 cm⁻¹ Mn mode is the most exposed: its Mn diagonal is only 5.4 eV/Å².
  - Cheapest fix: rerun the four basis modes at Q = 0.03 Å with conv_thr 1e-9, ~5 h each. Not submitted.
- The cell is the experimental one, unrelaxed in LDA+U. The soft mode is the one most sensitive to that, since
  the interlayer distance sets K_gap.
- **For the paper.** This axion-active mode sits at 0.54 THz, the lowest of the four. For a given energy its
  amplitude is the largest: the zero-point amplitude scales as ω^(−1/2), and the thermal one at 10 K as 1/ω.
  Its dθ comes from the same four basis-mode dθ/dQ, via `DTHETA_DQ`.

## Basis modes at Q = 0.03 Å submitted (2026-10-15 12:40)

The paper protocol (user, 10-15) has one frozen-phonon method for both materials:
- Φ from the forces of one SCF per symmetry-adapted basis vector.
- Q is chosen so that three criteria pass: E vs F diagonal within 1%, Φ symmetric within 1%, SCF force error
  below 1% of ΔF.
- YIO already passes at 0.01 Å. MBT gets the four T₂⁻ basis vectors at 0.03 Å with conv_thr 1e-10.
- The rigid-layer runs above become the independent check: the new Φ must give c·Φ·c = 2.40 eV/Å² without them.
  The row splice (method 4 in the notebook) is not part of the paper.

Jobs, all started 12:39, `db4b99b4`, 16 k, 12 operations, conv_thr 1e-10 confirmed in the outputs:

| job | basis mode | node |
|---|---|---|
| 213555 | 1 (Mn) | n135 |
| 213556 | 2 (Bi) | n137 |
| 213557 | 3 (Te1) | n134 |
| 213558 | 4 (Te2) | n117 |

- Cluster folders: `MnBi2Te4/phonon_lda/T2-/mode_k/Q_0p03_2026-10-15/` (NOTE in each). pw.x, the UPFs and
  tmp.tar.gz are hard links.
- Input: the Q_0p01_kfix scf.in with positions base 210628 + 0.03 v_k, and conv_thr 1e-7 → 1e-10.
- These runs give forces only. dθ stays on the 0.01 Å runs and their 8×8×8 NSCFs.
- When done: copy to `data/MnBi2Te4/archive/phonon/T2-/modeN/Q_0.03A/lda/scf_<job>/`. In
  `analysis/phonon_frequencies.ipynb`, set `BASIS_RUNS` to them and `Q_BASIS = 0.03`, keeping the rigid runs as
  `EXTRA_RUNS`.

## Does conv_thr 1e-7 also limit dθ? (2026-10-15)

Every dθ SCF (MBT and YIO, base and modes) used conv_thr 1e-7. Script:
`scratchpad force_parts.py`, which uses QE's per-contribution force blocks (verbosity high).

**The SCF error relative to the electronic response.** Units are eV/Å² per unit Q, projected on each basis
vector. "err" is QE's SCF-correction term.

| run | total | ionic | electronic | err / total | err / electronic |
|---|---|---|---|---|---|
| MBT Mn | 5.36 | 1176 | −1171 | 4.9% | 0.02% |
| MBT Bi | 29.9 | 2693 | −2663 | 2.8% | 0.03% |
| MBT Te1 | 20.4 | 855 | −834 | 12% | 0.30% |
| MBT Te2 | 12.7 | 651 | −639 | 5.5% | 0.11% |
| YIO Y | 439 | 9160 | −8721 | 0.18% | 0.009% |
| YIO Ir | 492 | 4468 | −3975 | 0.12% | 0.015% |
| YIO O | 149 | 3418 | −3269 | 1.9% | 0.09% |

- The SCF error has the same absolute size in both materials.
- The MBT force constant is a 0.2-2.5% remainder of ionic and electronic terms of order 10³. That cancellation
  amplifies the error 40-200×. YIO's remainders are 5-12%.
- dθ is purely electronic, with no ionic term to cancel against.
- The YIO dθ BZ sum has little cancellation of its own: Σ|c2 density| / |Σ c2 density| = 1.2 / 2.6 / 1.1 for
  modes 1-3 at nk 12.
- Expected SCF-threshold error in dθ: ≲ 0.2% for YIO, assuming dθ sees the density error the way the electronic
  force does. That is an estimate, not a test.
- The direct test would be one YIO mode at conv_thr 1e-10, then NSCF, Wannier and dθ, about 5 h in total (SCF
  12 min, NSCF 3.5 h). If that change is small, the base needs no rerun.
- For MBT, check the same c2 cancellation ratio once the MBT dθ exists.

**The Mn E/F = 1.24 is not explained by the SCF error.**
- ΔE_Mn = 333 μeV. The mismatch is ~80 μeV, but the SCF accuracy is 8e-8 Ry (~1 μeV).
- QE's force correction is 0.26 against a 1.31 eV/Å² gap.
- The Hubbard occupations match the base to 1e-5.
- The Hubbard force term along the Mn mode is −0.61 eV/Å², the largest of all runs (≤ 0.13 elsewhere). The rigid
  runs, where the octahedra move rigidly, have −0.008 and E/F 0.995-0.997.
- So a systematic error in the noncollinear ortho-atomic Hubbard force is a candidate, *not established*.
- The 0.03 Å Mn run at 1e-10 (213555) decides it. E/F → 1 means SCF noise. E/F staying near 1.24 means a
  systematic force error. The frequency method would then use energies for the Mn entries.
- This does not affect dθ, which uses no forces.

## Basis modes at Q = 0.03 Å: results (2026-10-15 21:30)

All four finished (17-23 iterations, ~5.5-7 h). Data: `data/MnBi2Te4/archive/phonon/T2-/mode{1..4}/Q_0.03A/lda/scf_2135{55..58}/`.
Notebook: `analysis/phonon_frequencies.ipynb`, which now reads these runs and uses method 1 as the
production result.

**Checks.** pw.x `db4b99b4` (md5 in each job folder), 16 k, 12 operations, conv_thr 1e-10 reached, total m 0,
|m|abs 9.99, Hubbard Mn moments ±4.6524 (mirrored). SCF force correction 0.0004-0.0009 eV/Å, against 0.010-0.034 at
Q = 0.01.

**Acceptance criteria (paper protocol), all pass:**

| | Mn | Bi | Te1 | Te2 |
|---|---|---|---|---|
| ΔE (meV) | 2.986 | 13.930 | 9.312 | 5.735 |
| Φ_ii from E (eV/Å²) | 6.635 | 30.956 | 20.692 | 12.745 |
| Φ_ii from F | 6.626 | 30.959 | 20.691 | 12.732 |
| E / F | 1.001 | 1.000 | 1.000 | 1.001 |

- Φ asymmetry max|Φ − Φᵀ| / max|Φ| = 0.07% (2.5% at Q = 0.01).
- SCF correction < 0.5% of the force change (Φ_ii Q ≥ 0.2 eV/Å).

**The Mn E/F of 1.24 was SCF noise, not a Hubbard-force error.** At conv_thr 1e-10 the Mn ratio is 1.001. Φ_MnMn
is now 6.63 eV/Å², which matches the Q = 0.01 *energy* value (6.67), not its force value (5.36). The Hubbard-force
candidate in the section above is ruled out.

**Φ (eV/Å², basis Mn, Bi, Te1, Te2):**

```
[[ 6.626  3.12   0.824  4.074]
 [ 3.119 30.959 19.073 -6.788]
 [ 0.808 19.096 20.691  2.078]
 [ 4.062 -6.782  2.101 12.732]]
```

**Frequencies (method 1):**

| ω (cm⁻¹) | meV | THz | Mn / Bi / Te1 / Te2 share of \|Δr\|² | Q = 0.01 value |
|---|---|---|---|---|
| **17.8** | 2.20 | 0.532 | 0.14 / 0.27 / 0.28 / 0.31 | 18.0 (spliced) |
| **82.7** | 10.26 | 2.480 | 0.21 / 0.06 / 0.24 / 0.50 | 77.7 |
| **128.7** | 15.96 | 3.859 | 0.67 / 0.09 / 0.15 / 0.10 | 124.0 |
| **143.3** | 17.77 | 4.297 | 0.68 / 0.12 / 0.19 / 0.02 | 139.8 |

**Rigid-layer check (independent runs, not used in Φ).**
- c·Φ·c = 2.346 eV/Å² from method 1, against 2.398-2.410 measured: 2% low.
  - The soft mode is 17.8 cm⁻¹, against the Rayleigh value of 18.2.
  - c·Φ·c is a ~2% remainder of entries of order 5-30, so a 0.05 eV/Å² error is 0.2% of Φ's entries.
- Φ·c from method 1 is (0.254, −1.970, 0.333, 0.211). The two rigid runs give (0.291, −1.952, 0.368, 0.201) and
  (0.235, −1.986, 0.390, 0.204), so they scatter about as much among themselves as method 1 differs from them.
- Splicing the rigid row in (method 4) moves only the soft mode, to 18.0. The other three don't change at 0.1 cm⁻¹.
- The soft mode's overlap with c is 0.9995.

**What changed against Q = 0.01:**
- The three optical modes moved up by 3.5-5 cm⁻¹, more than the 1-3 cm⁻¹ error bars of log 10-14. The section above
  had already flagged those error bars as too small.
- The eigenvectors changed:
  - The 129 cm⁻¹ mode's Mn share went from 0.79 to 0.67.
  - The 143 cm⁻¹ mode's Mn share went from 0.39 to 0.68.
  - So the two upper modes now share the Mn motion, which matters for dθ per normal mode.
- The soft mode is unchanged.

**Consequences.**
- The paper protocol holds with method 1 alone. The rigid-layer runs serve as an independent validation, accurate to
  2% in the softest stiffness and 1% in its frequency.
- dθ/dQ still comes from the Q = 0.01 runs and their 8×8×8 NSCFs. Only Φ and the eigenvectors come from these runs.
- The SCF threshold question for dθ is closed without the YIO 1e-10 test (user, 10-15).
- These runs' `tmp/` folders on the cluster (1.9 GB each, force-only runs) were deleted 10-15.
