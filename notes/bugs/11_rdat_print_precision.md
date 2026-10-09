# #11 — `_r.dat`'s 6-decimal print costs ~3% of dθ

> **Status: FIXED (mitigated)** (2026-09-23). Production builds A(R) at full precision from the `.chk`
> (`transl_inv_full_from_chk`, `right_centre_from_chk`), not from the printed `_r.dat`. Building H(R) at full
> precision from `.chk` + `.eig` is a possible follow-up (not done).

**Found:** 2026-09-23. **Affects:** dθ from Jae-Mo's `rdat_ndegen_applied` route (#5), which reads `_r.dat`.

## Symptom
On the YIO mode-1 pair, dθ from Jae-Mo's files (−0.0045095 / −0.0045006) differed by ~3% from the same formula
evaluated at full precision (−0.0046468 / −0.0046070).

## Root cause
- **Print precision:** `_r.dat` prints A(R) with 6 decimals (`F12.6`, rounding errors ≤ 5e-7 Å).
- **Independent errors:** the errors are independent in the base and mode files, so they don't cancel in the
  finite difference.
- **Amplified:** they are divided by Δβ = 0.01 Å.

## Why it was a bug
The finite difference amplifies print-level noise by 1/Δβ. Rounding our own full-precision A(R) to 6 decimals
reproduces Jae-Mo's numbers to 2e-4, which confirms the cause.

## Fix
Use the full-precision A(R) built from the `.chk` (the same formula as Jae-Mo's, `transl_inv_full_from_chk`).

## Remaining
`_hr.dat` is also printed at 6 decimals (eV) and feeds every route. Building H(R) at full precision from `.chk`
+ `.eig` would remove the analogous error; not done.

**Sources:** [`../log/2026-09-23_ws_ties_in_beta_derivative.md`](../log/2026-09-23_ws_ties_in_beta_derivative.md).
