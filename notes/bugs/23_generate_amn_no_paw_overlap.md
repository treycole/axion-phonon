# 23. `generate_amn.jl` overlaps PAW pseudo-wavefunctions without the augmentation term

**Status: OPEN, fix in progress** (found 2026-10-28; Python fix `wannier/paw_overlap.py`, user's choice over a QE patch). Not fixed in any production run. Evidence:
`archive/diagnostics/2026-10-28_selftest_paw_overlap/` `check_selftest_amn.py` and its `.log` (inputs in
`data/MnBi2Te4/archive/base/lda/selftest_paw_overlap_2026-10-28/`).

## Symptom

The Q = 0 self-test doesn't reproduce the base it is projected from (log `2026-10-17_yio_amplitude_series.md`,
"The mode pipeline does not reproduce the base at Q = 0").
- **MBT:** Ω_I is 344.46 Å² for the self-test (214301) against 325.53 for the atomic-trial base (214300). The extra
  spread sits in the Mn d WFs (+30%).
- **YIO:** the gap is smaller, 116.42 against 117.10 Å² (−0.6%).

The cause was not found at the time.

## Cause

`QEWavefunctions.jl/scripts/generate_amn.jl` computes
A_mn(k) = Σ_m′ ⟨ψ̃_m(k)|ψ̃′_m′(k)⟩ V_m′n(k)
- ⟨ψ̃|ψ̃′⟩ is the plain plane-wave dot product of QE's pseudo-wavefunctions on shared G vectors.
- V = U_dis U comes from the reference `.chk`.

With PAW, the pseudo-wavefunctions are orthonormal only in the metric S = 1 + Σ_ij |p̃_i⟩ q_ij ⟨p̃_j|. pw2wannier90
includes S; the Julia package has no augmentation term anywhere. So for the self-test, A = O V with O = ⟨ψ̃|ψ̃⟩ ≠ 1,
not A = V.

## Evidence (MBT 8³: base 214300 `.chk`, self-test 214301 `.amn`)

| quantity | value |
|---|---|
| \|A − V\| / \|V\| | 0.22 |
| anti-Hermitian part of V†A, relative | 6e-12 (so A = O V with O Hermitian, as the cause predicts) |
| eigenvalues of V†A | 0.61 to 1.36 (median 1.08) |
| part of A outside span(V) | 7.8% |

| WFs | diagonal ⟨w\|w′⟩ | outside span(V) |
|---|---|---|
| Bi s,p | 1.10 | 5.7% |
| Te s,p | 1.15 | 6.9% |
| Mn d | 0.66 | 15.4% |

The Mn 3d norm deficit (large augmentation charge) is where the self-test's extra spread appears.

**YIO 8³** (base 214162 `.chk` = the self-test's `base_ref`, self-test 214249 `.amn`;
`check_selftest_amn_yio.log`). The defect is as large as in MBT:

| quantity | value |
|---|---|
| \|A − V\| / \|V\| | 0.17 |
| anti-Hermitian part of V†A | 8e-12 |
| eigenvalues of V†A | 0.78 to 1.18 (median 0.89) |
| part of A outside span(V) | 7.6% |

| WFs | diagonal ⟨w\|w′⟩ | outside span(V) |
|---|---|---|
| Y p | 1.16 | 5.2% |
| Ir d | 0.80 | 4.0% |
| O s,p | 0.95 | 9.0% |

But YIO's Ω_I moves only 0.6%, against 6% for MBT. So the size of the `.amn` error doesn't predict how much the
disentangled subspace moves, and Ω_I alone can't show which material's dθ is affected more. The YIO self-test changed
dθ by 0.5-1% from the atomic-base values (log 10-17).

## Consequences

1. **The self-test is not a fixed point.** Projection-only disentanglement on A = O V picks a subspace different from
   span(V), which is the 6% Ω_I shift. Modes are paired with the self-test, so both ends of dθ use the same
   construction; the dθ values are internally consistent.
2. **The "Λ = 0 exactly" argument does not hold as stated.** It needs A(λ) = ⟨ψ_λ|g⟩ with orthonormal ψ_λ and a fixed
   trial g in the base subspace. Here the effective trial in the true metric depends on λ, through the augmentation
   spheres that move with the atoms, and it is not in the self-test subspace. So A_β = 0 is again an approximation
   for production runs. Its size is unknown.
3. **YIO has the same defect** (7.6% of A outside span(V)), although its Ω_I moves only 0.6%.

## Whose bug

The Julia tool's (QEWavefunctions.jl, `scripts/generate_amn.jl`): it assumes norm-conserving pseudopotentials, and
nothing in the package or README says so. QE is correct: PAW pseudo-wavefunctions are meant to be orthonormal in S,
and pw2wannier90's `compute_amn` applies S (`calbec` + `s_psi`) when `okvan`.

## Fix in progress: `wannier/paw_overlap.py` (2026-10-28)

The user chose Python over a pw2wannier90 patch. The module computes the cross-structure overlap
O_mn(k) = ⟨ψ_Q,mk|ψ_0,nk⟩ from the two NSCF saves, with the PAW metric rebuilt from QE's own conventions:
- projectors as `init_us_2`, with the β table and interpolation of `beta_mod`, and `ylmr2`;
- q_ij from `compute_qqr`;
- the spin-orbit `fcoef`/`qq_so` of `init_us_1`;
- S applied as `s_psi_nc`.

From O(k) come:
- the corrected mode `.amn`, A = O·V₀;
- Λ(k) = antiHerm(V₀†O†U_Q)/Q (derivation §8.4), for any existing pair.

**Variants** (where the projectors sit; all are equal at Q = 0):
- `current`: both sides on the displaced atoms.
- `reference`: both sides on the base atoms.
- `mixed`: each side on its own atoms.
- `exact`: `mixed` plus the partial-wave term d·D_ij, with D_ij = ⟨φ_i|∇φ_j⟩ − ⟨φ̃_i|∇φ̃_j⟩ from the UPF partial waves
  (large components). It is exact to first order in the displacement d, which is all dθ/dQ needs.

**Validation** (`archive/diagnostics/2026-10-28_selftest_paw_overlap/check_paw_overlap.{py,log}`):
- **Test system:** local pw.x 7.4.1, MBT geometry, 25 Ry, 3 general k, 150 bands. Also pattern 3 displaced by 0.03 Å.
- **S matches QE exactly:** ⟨ψ|S|ψ⟩ = 1 to 1e-14 at every k for both structures, spin-orbit included. The plain
  overlap is off by 0.60.
- **Cross-structure overlap:** two orthonormal sets cannot overlap with a singular value above 1. The overshoot
  is:

  | variant | max singular value − 1 |
  |---|---|
  | current / reference | 2.5-3.1e-3 |
  | mixed | 0.5-1.1e-5 |
  | exact | 4-8e-8 |

  This validates D independently, at O(d).
- **Unit tests** (`validation/unit/test_paw_overlap.py`): Y_lm orthonormality and QE order, simpson, G matching, and
  D against a direct 3D quadrature.

**Production pair (running 2026-10-28):** MBT pattern 3, Q = 0.01 Å, 8³.
- `setup_paw_overlap_pair.sh` → jobs 218716-218723, `MnBi2Te4/paw_overlap_2026-10-28/mode3_Q0p01_k8/`.
- `assemble_paw_overlap.py` then:
  - checks the k order;
  - checks that O_plain·V reproduces the Julia `.amn` (conventions);
  - runs the Λ = 0 zero test for the `.amn` wannier90 used;
  - computes Λ for the production pair (self-test base) and against the atomic base;
  - writes the corrected mode `.amn` for a new wannier90 run. The base of the corrected pair is the atomic base
    214300 itself, since A = V exactly at Q = 0.

## Result: MBT pattern 3, Q = 0.01 Å, 8³ (2026-10-28)

Logs: `assemble_paw_overlap.log`, `lambda_corrected.log`; the tube log is in
`data/MnBi2Te4/paw_overlap_2026-10-28/mode3_Q0p01_k8/`.

**Checks:**

| check | result |
|---|---|
| k order (wfc vs `.chk`) | 5.6e-16 |
| our plain overlap × V vs the Julia `.amn` | 5.9e-12 |
| zero test (Λ for the `.amn` wannier90 used) | 1.2e-10 Å⁻¹ |
| S-orthonormality of QE's wavefunctions, all tasks | ≤ 1.7e-14 |

**Λ and spreads:**

| quantity | old (Julia, self-test base) | corrected (atomic base) |
|---|---|---|
| self-test Ω_I vs base 325.534103937 Å² | 344.46 | 325.534103932 (218726) |
| mode Ω_I | 344.50 | 325.556 (218725) |
| Λ^II(k), rms / max (Å⁻¹) | 0.029 / 0.037 (0.071 / 0.19 without D) | 3e-10 / 1.2e-10 (218760) |

**dθ/dQ** (rad/Å, tube m 4, 8 planes; base end / mode end):

| pair | base end | mode end |
|---|---|---|
| old | 6.43441 | 6.60273 |
| corrected | 6.43227 | 6.60007 |
| change | −0.033% | −0.040% |

**Reading:**
- The bug changed the Wannier subspace (Ω_I +6%) and made Λ ≈ 0.03 Å⁻¹, but it moved pattern 3's dθ/dQ by only
  0.03-0.04%. The self-test pairing kept both ends consistent, and the Λ term contributes at that level.
- The end gap (0.168 rad/Å) is unchanged, so it is not from the bug. That fits the O(Q²) gap-closing explanation.
- **Open:**
  - the window study at 8³ with corrected `.amn` (jobs 218804-218888, `paw_overlap_2026-10-28/window_study_k8/`);
  - the other patterns;
  - YIO.

## Earlier scope: a pw2wannier90 patch (not chosen)

Compute the mode/self-test `.amn` in pw2wannier90 (QE fork), with the reference WFs as the trial functions:

1. **New input:**
   - `ref_outdir` / `ref_prefix`: the base save;
   - `ref_chk`: the base `.chk`;
   - `ref_tau`: the base positions;
   - `ref_overlap = 'mixed' | 'current' | 'reference'`.
2. **In `compute_amn`, replacing `generate_guiding_functions(ik)`:**
   - `read_collected_wfc(ref_dir, ik, evc_ref)` reads the base wavefunctions straight into this run's G
     distribution. Same cell and cutoff, so it is the same G list.
   - g̃ = evc_ref(selected bands) · V(ik), with V = U_dis U from the `.chk` and lwindow/ndimwin as in
     `wannier_io.read_wannier_checkpoint`. Assert that the `.chk` k-points match xk.
3. **Apply S through the existing USPP branch** (`calbec` → `s_psi` → ZGEMV, spinor path):
   - **`current`:** vkb of the displaced structure. A = ⟨ψ̃|S_Q|g̃⟩.
   - **`mixed`:** `becp = calbec(vkb(base τ), g̃)`, then vkb of the current τ and `s_psi`. That gives
     A = ⟨ψ̃|g̃⟩ + Σ ⟨ψ̃|p̃^Q⟩ q ⟨p̃^0|g̃⟩, each projector on its own structure. `s_psi_nc` reads the global
     becp/vkb, so this is a τ swap: `struct_fact` + `init_us_2` with the base τ, then restore.
   - **`reference`:** vkb of the base structure on both sides.
   - All three agree exactly at Q = 0. They differ at O(Q) through the partial-wave dipole term
     ⟨φ_i|∇φ_j⟩ − ⟨φ̃_i|∇φ̃_j⟩, which none of them includes. Their spread in dθ bounds that term, and so the size
     of Λ.
4. **Read the `.chk` in Fortran** (unformatted; `exclude_bands`, `lwindow`, `ndimwin`, `u_matrix_opt`, `u_matrix`).
   Set n_proj = num_wann and skip the projection parsing of the `.nnkp`.

**Effort:** ~300-400 lines on `production/axion-phonon`, a new `pw2wannier90.x` md5 (the binaries list changes), about
2-3 days with the tests.

**Run cost:** one `.amn`-only pw2wannier90 run per structure, reading two saves. That is much cheaper than the `.mmn`
(7.5 h at MBT 8³).

**Validation:**
- **(a) Self-test:** A = V to ~1e-10, and Ω_I = 325.53 Å² (MBT) exactly, with any `ref_overlap`.
- **(b) MBT pattern 3, 8³:** run all three `ref_overlap` choices. Compare dθ with 6.434 / 6.603 rad/Å and the end gap.
- **(c) YIO:** the same, if its check shows an effect.

Not chosen:
- Porting β projectors and q_ij (with SOC qq_so) to Julia would duplicate `init_us_2`.
- Atomic trials for the modes would break the projection-gauge Λ = 0 argument outright.
