# #15 — MnBi2Te4 base u_3.0 input used a rounded bohr cell (same defect as #9)

> **Status: FIXED** (2026-10-05). u_3.0 was rebuilt as an exact copy of u_4.0 (inputs, `send_job.sh`, binaries,
> UPFs) with only `U Mn{1,2}-3d 3.0` changed, and resubmitted as job 202971.

**Found:** 2026-10-05. **Affects:** MnBi2Te4 base u_3.0 SCFs before 2026-10-05 (e.g. 202840).

## Symptom
QE found **4** symmetry operations / 24 k-points for u_3.0, against u_4.0's 12 operations / 13 k-points for the
same structure.

## Root cause
- **Rounded cell:** `MnBi2Te4/base/soc/u_3.0/MnBi2Te4.scf.in` used a rounded `CELL_PARAMETERS bohr` cell. The
  in-plane a1/a2 lengths differ by 3e-6 bohr, and C3 maps the lattice onto itself only to 4.5e-6 in crystal
  coordinates, above QE's 1e-6 lattice tolerance. C3 was lost.
- **Different cutoff:** the input also had `ecutrho = 240`, against u_4.0's 500.

## Why it was a bug
- **Different calculations:** u_3.0 and u_4.0 were meant to differ only in U, but they used different
  symmetrization (4 vs 12 operations) and a different density cutoff.
- **Not comparable:** the u_3.0 results couldn't be compared with u_4.0, and its SCF was symmetrized over too small
  a group.

## Fix
u_3.0 rebuilt as an exact copy of u_4.0 (the precise ångström cell, ecutrho 500), with only U changed. 202840 was
killed, the old 201275bc889/202840bc889 deleted (old inputs backed up), and the run resubmitted as 202971.

## Remaining
- **Before using any u_3.0 base run,** check `Sym. Ops.` = 12.
- **The T2⁻ inputs still have this defect** (#9).

**Sources:** memory `mbt-symmetry-with-labels`.
