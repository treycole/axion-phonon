# #9 — MnBi2Te4 T2⁻ inputs: a rounded bohr cell loses C3, and `symmetry_with_labels` triggers #7

> **Status: FIXED** (checked 2026-10-13). The production T₂⁻ runs (`phonon_lda/T2-/mode_*/Q_0p01_kfix_2026-10-11`,
> 212708-212711, and their 8×8×8 NSCFs) use the exact `CELL_PARAMETERS angstrom` cell (a = 4.33360 Å, c/3 = 27.284 Å)
> and QE finds 12 operations (6 with fractional translations). `symmetry_with_labels` is still on; the two Mn-swap
> operations it adds were the #7 trigger, and #7 is patched in the production pw.x `db4b99b4`.

**Found:** 2026-09-22. **Affects:** MnBi2Te4 zone-boundary T2⁻ mode inputs (`phonon/T2-/mode*/...`). The T1⁺
inputs already use `CELL_PARAMETERS angstrom`.

## Symptom
- **Lost symmetry:** QE finds fewer symmetry operations for the distorted T2⁻ cells than the structure has.
  With `symmetry_with_labels` off, 6 ops; with it on, 12 ops, two of which are the #7 trigger.
- **Scale drift:** static atoms carry a uniform 7.2e-8 z-scaling drift relative to base run 151625.

## Root cause
1. **Rounded cell.** The T2⁻ inputs give `CELL_PARAMETERS bohr` rounded to 8-9 digits. a1 + a2 + a3 then has a
   ~5e-5 bohr in-plane component, so on-axis atoms get x ≠ y in crystal coordinates at the 1e-5 level, above
   QE's symmetry tolerance. C3 is lost.
2. **Labels.** With `symmetry_with_labels = .true.`, T2⁻ gains two order-6 operations (S6, with ft = ½) that
   swap Mn1 ↔ Mn2, which is exactly #7's trigger.
   - In the base cell, the flag is harmless: Mn1 and Mn2 share a chemical symbol, so QE finds identity + (½,½,½)
     ("this is a supercell") and disables all fractional translations, which gives the same 12 ops as with the
     flag off.

## Why it was a bug
- **Fewer operations:** they mean weaker symmetrization and a lower-symmetry SCF state than the structure has,
  just from rounding.
- **Wrong operations:** with an unpatched pw.x, the label flag added operations that #7 mishandles.

## Fix (known)
- **Regenerate the T2⁻ inputs on the precise ångström cell** (as in 151625), with positions = base + Q·mode.
- **Flag on only with a #7-patched pw.x.** All current builds are patched.

## Remaining
**Apply it**, and check `Sym. Ops.` in the next T2⁻ SCF output. #15 is the same rounded-cell defect, found in
the MBT base u_3.0 input.

**Sources:** memory `mbt-symmetry-with-labels` (from a Python port of QE's `symm_base.f90`, validated op-by-op on
26 MBT and YIO runs).
