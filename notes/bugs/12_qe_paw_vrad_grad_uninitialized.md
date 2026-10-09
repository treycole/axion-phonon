# #12 — QE PAW `v_rad` / `g_rad` used uninitialized in noncollinear GGA

> **Status: FIXED** (2026-09-28). Fixed upstream (QEF issue #863, commit `87ac195cf`, 2026-07-10, after the 7.6
> tag). The cluster `~/q-e` (pw.x md5 `c41c409`) and `~/q-e-segnifix` (`e270b05b`) include it. MBT u_3.0 picked it
> up only on 2026-10-04 (#14).

**Found:** 2026-09-27. **Affects:** QE 7.4, 7.5 and 7.6 with noncollinear + spin-orbit + PAW + GGA: every MBT and
YIO run made with those releases.

## Symptom
- **MnBi2Te4:** the base SCF blew up at iteration 1 (energy ~1e6 Ry) on 64 processes. On 48 processes it ran but
  converged to a net moment of about −2 μB in an antiferromagnet.
- **QE 7.2 was fine.**
- **Depended on process count and problem size.**

## Root cause
`PW/src/paw_onecenter.f90::compute_pot_nonc` (the noncollinear-GGA PAW one-centre XC potential) ALLOCATEs `v_rad`
and `g_rad` and never zeroes them:
- **`v_rad(:,:,2:4)`** is written only where |m| > eps12; elsewhere it keeps whatever was in memory.
- **`compute_g`** *accumulates* into `g_rad` (`g_rad = g_rad − …`) whenever `with_small_so` (fully relativistic
  PAW, i.e. every rel-*.UPF used here).

In QE 7.2 the routine had `v_rad = 0` and `IF (has_so .AND. ae == 1) g_rad = 0`. Commit 4ce366b33 (2024-03-15,
"Vpaw_acc - cosmetics") removed both during the OpenACC port.

## Why it was a bug
- **Garbage in the potential:** the PAW one-centre B-field term picks up whatever the heap held, so results depend
  on heap state, i.e. process count and problem size.
- **MBT is hit hardest:** Bi and Te start at exactly m = 0, so at iteration 1 their whole B-field term is garbage.

## Fix
Zero `v_rad` and `g_rad` after allocation (upstream `87ac195cf`; also an upstream test-suite case for
PAW + SO + magnetic + GGA).

## Verification
Patched build on MBT input 192903: 48 and 64 processes agree (iteration 1 −7140.62251 / −7140.62245 Ry; total
magnetization 0.00, absolute 10.43). The new `~/q-e` build reproduces that to 2e-8 Ry (193147).

## YIO
- **Little effect on the SCF:** the 192344 vs 193145 comparison found a small change (iteration-1 energy to
  2e-5 Ry, Ir moment +1.5%).
- **It is not the cause of the YIO symmetry breaking (#16).** NSCF tests A vs B on 2026-10-05 were identical.
- **It does change the converged SCF state:** a fixed SCF gave 4.8 meV of #16 breaking, against 7.6 meV
  unfixed.

## Lesson
`md5sum` the pw.x first when a run misbehaves: known-good `c41c409…`, known-bad `1541022…` (#14).

**Sources:** memory `qe-paw-nc-uninit-vrad`;
[`../log/2026-10-05_yio_symmetry_floor_dft_origin.md`](../log/2026-10-05_yio_symmetry_floor_dft_origin.md) §4.
