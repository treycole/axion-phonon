# #1 — YIO Wannier model closed the gap (diffuse empty projections); wrong n_occ

> **Status: FIXED** (2026-07-27). Projections changed to "Option A" (`Y: p`, `Ir: d`, `O: s;p`, 176 WF) and
> `n_occ = 156`. Every YIO run since uses them.

**Found:** 2026-07-27 (date of the memory record). **Affects:** Y2Ir2O7 Wannier models, so every YIO dθ before
this date.

## Symptom
- dθ/dQ came out at ±thousands and jumped erratically with nk.
- In the `nk_sweep` output, `min_gap` collapsed to 0.006-0.017 eV at k = (0, 0, 0.75), against ~0.72 eV for a
  healthy model, and the Chern² density reached ~1e5.
- The DFT itself was a clean insulator (0.42 eV gap).

## Root cause
1. **Diffuse empty projections.** The projections included empty, diffuse shells (Y 4d, Ir 6s/6p). These can't
   be Wannierized in one shot (`dis_num_iter = 0`, `num_iter = 0`) and came out delocalized: Final State spreads
   of 9-21 Å² (Ir 6s/6p) and ~24 Å² (Y 4d), against ~0.7-1.1 Å² for real WFs. They injected spurious
   long-range hoppings into `_hr.dat`, which nearly closed the interpolated gap.
2. **Wrong n_occ.** The old value, 104, was carried over from a model without Y 4p / O 2s as WFs. Option A
   includes them, so there are 164 − 8 (Y 4s, excluded) = 156 occupied Wannier bands: Y 4p (24) + O 2s (28) +
   O 2p (84) + Ir 5d below E_F (20). Using 104 puts the chemical potential inside the O 2p manifold.

## Why it was a bug
The Kubo-type curvature goes as 1/gap², so a spurious near-closing of the *interpolated* gap blows up the integral.
The number tracked the Wannier model's defects, not the physics.

## Fix
- **"Option A"** (num_wann = 176, projection-only): `Y: p`, `Ir1..4: d`, `O: s;p`. All spreads ≲ 1.1 Å²,
  Ω_total 148.5. The rejected options were B (adds Y:d, 23 Å² junk, Ω 877) and C (adds Ir:s;p, 6-9 Å², Ω 555).
- **n_occ = 156.**

## Verification
Final State spreads all ≲ 1.1 Å². The minimum gap is comparable to the DFT gap. The Ir 5d conduction manifold
(20 states) is enough for the Kubo sum.

## Remaining / guidance
- **Check the model:** compare `min_gap` in the npz with the DFT direct gap, and check `.wout` spreads (any
  WF ≫ 3 Å² is a diffuse empty projection).
- **Don't rescue diffuse shells with `dis_num_iter`:** it breaks the projection gauge the finite difference
  relies on.

**Sources:** memory `wannier-gap-collapse-Os`; `notes/materials/Y2Ir2O7.md`.
