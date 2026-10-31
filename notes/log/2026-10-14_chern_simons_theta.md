# 2026-10-14: θ itself: MBT parity indicator, and the Chern–Simons 3-form as a check of dθ/dQ

**Bottom line.**
- **MBT θ = π.** The LDA base Wannier model (210651, pinned) has η₄ = 2 from the inversion parities at the 8 TRIMs.
  `analysis/theta.ipynb` now computes this; its old SCDM 3-form was singular at Γ.
- **The Chern–Simons 3-form reproduces the YIO dθ/dQ.** With an exact connection in a smooth projection gauge,
  [θ_CS(Q) − θ_CS(0)]/Q matches the production c₂ integral for all three Γ₂⁻ LDA modes at nk = 10 (0.09%, 0.25%,
  0.003%), each value lying between the two c₂ estimates. The frozen gauge (A_β = 0) is good to ≲ 0.3%.
- **For θ = π the 3-form is not practical:** the gauge must break I and S, and the mesh sum needs nk ≳ 16.

Scripts and logs: `archive/diagnostics/2026-10-14_cs_theta/`. The construction, with equations, is in
`notes/meetings/2026-10-12.md` §6.

## 1. MBT: why the old `theta.ipynb` failed, and the parity count

- **Parities** (`parity_and_scdm.py`): n₊/n₋ = 48/34 at Γ, 46/36 at the three (0,½,½)-type TRIMs, 41/41 at the other
  four. η₄ = 22 mod 4 = 2. Inversion on the Wannier basis from the `.nnkp` projections:
  |D H(k) Dᵀ − H(−k)| = 4.5e-7 eV.
- **The old gauge:** SCDM trials (columns of the mesh-averaged projector) are within 0.2% of I- and S-symmetric, so
  σ_min of ⟨u_occ|g⟩ is 4e-5 at Γ and θ comes out ~1e12. Multiplying spin-down by i is a spin rotation, which
  commutes with Θ: it changes nothing (3e-5).
- **Smooth gauges for θ = π** (`theta_cs.py`, `designed_trials.py`; tight-binding positions only, which is enough for a
  quantization check):

  | trials | σ_min | θ/π at nk 6 / 8 / 10 / 12 / 14 / 16 |
  |---|---|---|
  | designed: ionic set, two Te p_z → (Te p_z + Bi p_z)/√2, one per spin (best of 2048 site/phase choices) | 0.13 | 3.11 / 1.76 / 1.13 / 0.76 / – / – |
  | optimized: SCDM + random, ascent of Σ log det O, then −Σ Tr O⁻ᵖ (p = 1, 2, 4) on 6³ | 0.39 | 2.21 / 1.40 / 1.16 / 1.02 / 0.966 / 0.965 |

  The ionic set has parities 46/36 at Γ, one doublet short of the occupied 48/34; its σ_min is 3e-8. The nk = 20 run
  was killed for memory (all k at once, ~7 GB); `theta_cs` now batches k. The 0.965π plateau is not explained.

## 2. YIO: 3-form finite difference vs c₂ (`cs_full.py`, `yio_fd_dtheta.py`)

- **Setup:** production pair loader (H and A(R) from the `.chk`, shared WS fold, base centres as τ), positions in
  complex128. One SCDM trial set from the base: σ_min = 0.552 on base and every mode.
- **Checks:** Ω̃ = ∂Ã − ∂Ã − i[Ã, Ã] with A(R) included: 2e-9 of the scale. Im θ ≤ 1e-12. θ_CS(base): −1.6e-3,
  +1.2e-5, −4.1e-6, −9e-8 at nk 4 / 6 / 8 / 10 (must be 0 mod 2π).

| mode | nk | 3-form difference (rad/Å) | c₂ base / mode end | c₂ mean | difference − mean |
|---|---|---|---|---|---|
| 1 | 6 | −0.018277 | −0.018284 / −0.018409 | −0.018346 | 0.4% |
| 1 | 8 | −0.018761 | −0.017608 / −0.017780 | −0.017694 | 6.0% |
| 1 | 10 | −0.018035 | −0.017992 / −0.018110 | −0.018051 | 0.09% |
| 2 | 6 | +0.053989 | +0.053570 / +0.054191 | +0.053881 | 0.2% |
| 2 | 8 | +0.053901 | +0.052880 / +0.053510 | +0.053195 | 1.3% |
| 2 | 10 | +0.054300 | +0.053829 / +0.054502 | +0.054165 | 0.25% |
| 3 | 6 | +0.067788 | +0.065529 / +0.066411 | +0.065970 | 2.8% |
| 3 | 8 | +0.065333 | +0.065630 / +0.066557 | +0.066093 | 1.1% |
| 3 | 10 | +0.065486 | +0.065017 / +0.065959 | +0.065488 | 0.003% |

- **Why nk = 8 is off:** the difference θ(Q) − θ(0) is ~2e-4 (mode 1) to ~7e-4 (mode 3), so each θ needs a mesh error
  well below 1e-5; the base θ error is 4e-6 at nk = 8 and 9e-8 at nk = 10.
- **What this does and doesn't test:** both routes use the same H(R) and A(R), so it tests the c₂ code, the frozen
  gauge and the k sums, not A(R) (a finite difference in k on the DFT mesh).
- **Cost:** ~25 min per structure pair at nk = 10 (176 WF), after a 7-minute load.
