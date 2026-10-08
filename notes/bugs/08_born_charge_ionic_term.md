# #8 — Born effective charge: wrong ionic term

> **Status: FIXED** (2026-09-21, date of the record). Z^ion,eff = z_valence − (excluded electrons on the
> displaced atom). SrTiO3 trial_02 Z* = 7.370 stands; trial_03 moved to 6.665.

**Found:** recorded 2026-09-21 (from the SrTiO3 Z* validation, late August). **Affects:** comparing the
Wannier-based Born effective charge with `ph.x`.

## Symptom
Ti's Z* came out at 15.1 instead of ~7.1 when the ionic term was taken as `z_valence`.

## Root cause
- **Missing bands:** bands below `dis_win_min` aren't in the Wannier space.
- **What they do instead:** those localized on the **displaced** atom translate rigidly with it, contributing
  exactly −1 each to dP/du.
- **Where they belong:** they have to go into the ionic term, not be dropped.

## Why it was a bug
- **The electronic part:** Z* = ionic + electronic, and the Wannier electronic part covers only the Wannier
  bands.
- **The fix-up:** semicore electrons excluded from the Wannier space still move with their atom, so they act
  like extra core charge.
- **Getting it wrong:** with the bare `z_valence`, those electrons get the ion's +1 each and miss their own −1,
  double-counting their charge.

## Fix
    Z^ion,eff = z_valence − (excluded electrons on the displaced atom)

- **`dis_win_min` = −10 or −20 eV:** Ti 3s (−46 eV) and 3p (−22.5 eV) are both excluded, so Z^ion,eff = 12 − 8 =
  **+4** (also the chemical Ti⁴⁺ core charge), and n_occ = 30, not 40.
- **`dis_win_min` = −40 (trial_03):** only Ti 3s is excluded, so Z^ion,eff = +10 and n_occ = 38.

## Verification
trial_02 Z* = 7.370 (+1.1% vs DFPT) stands. trial_03 moved to 6.665 (−8.6%).

**Sources:** memory `born-charge-ionic-accounting`; `validation/observable/born_effective_charge.ipynb`.
