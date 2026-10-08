# #14 — MnBi2Te4 u_3.0 ran a stale, buggy pw.x

> **Status: FIXED** (2026-10-04). The fixed pw.x (`c41c409`) was copied into `MnBi2Te4/base/soc/u_3.0/` and the run
> resubmitted. Process fix: `md5sum` pw.x before any source-level debugging, and monitor the total magnetization on
> every check.

**Found:** 2026-10-04. **Affects:** MnBi2Te4 base u_3.0 runs made between 2026-09-28 and 2026-10-04.

## Symptom
- **u_3.0 drifted:** a u_3.0 base SCF (job 201275) drifted from ~0 total magnetization at iteration 1 to
  (0.00, −0.45, −8.33) μB/cell by iteration ~25, still growing.
- **u_4.0 didn't:** a concurrent u_4.0 run (202833) with the same structure stayed at (0, 0, 0).

## Root cause
- **What the update missed:** the 2026-09-28 binary update (to the #12-fixed `~/q-e` build, md5 `c41c409`)
  touched only `u_4.0/`.
- **What u_3.0 had:** `u_3.0/pw.x` was still the old `q-e-nsfix` build (md5 `1541022…`), which has bug #12.
- **Confirmed by md5sum.**

## Why it was a bug
An operational error, not a code bug, but it produced wrong physics: a net moment in an antiferromagnet.
- **Time lost:** about 1.5 hours went into hand-auditing QE's symmetrizer (`sgam_at_mag`, `sym_rho_serial`) against
  the Mn geometry. The symmetrizer was correct.
- **The answer was already recorded:** the bad md5 had been logged on 09-28.

## Fix
- **Binary:** copied `~/q-e/bin/pw.x` over `u_3.0/pw.x`.
- **Rerun:** killed 201275 and resubmitted (202840, then 202971 after #15).

## Lessons (memory `monitor-magnetization-and-check-binaries`)
1. **Check the total magnetization on every check:** grep `total magnetization` alongside `estimated scf
   accuracy`.
2. **md5sum first:** compare pw.x against the known-good/known-bad hashes before auditing source.
3. **Fix root causes:** never propose disabling a mechanism (`nosym`) as the fix.

**Sources:** memories `qe-paw-nc-uninit-vrad`, `monitor-magnetization-and-check-binaries`, `mbt-scf-restart-hang`.
