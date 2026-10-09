# #4 — Stock Wannier90 `_r.dat` cannot carry the position matrix

> **Status: FIXED (worked around)** (2026-08-27). Production never takes A(R) from a stock `_r.dat`. A(R) came
> from `.mmn` (with `_AA_cache`), and now comes from `.chk` or Jae-Mo's `_r.dat` (see #5).

**Found:** 2026-08-27 (tested on SrTiO3 and Y2Ir2O7). **Affects:** any attempt to get A(R) from the
`_r.dat` that `wannier90.x` writes.

## Symptom
A(R) read from `_r.dat` gave wrong curvatures and failed symmetry-forced-zero tests, with or without remapping
by `use_ws_distance` (min |R + τ_j − τ_i|), at any embedding.

## Root cause
All three are verifiable in `plot_write_rmn` (`src/plot.F90`):
1. **Print format `6F12.6`.** 51.9% of all components are exactly zero, and 70.7% are ≤ the 1e-6 Å granularity,
   while the median off-diagonal element at R = (1,0,0) is 1.3e-5 Å. The off-diagonal block is destroyed.
2. **No link phase.** The raw `m_matrix` is used, and the sum over neighbours happens inside the routine, so
   `exp(i b·τ)` can't be applied afterwards. This is unrecoverable, not just imprecise.
3. **`use_ws_distance` is ignored by this writer.** It only emits `wsvec.dat`.

## Why it was a bug
The file looks like the position matrix but is not usable as one for curvature work: wrong in content (2), wrong
representative (3), and with most of the signal rounded away (1).

## Fix
- **A(R) from `.mmn`:** built from the overlap matrices, with the link phase applied.
- **No `.mmn` kept on disk:** a warm `_AA_cache/*.npz` reproduces A(R) bit-for-bit, 2.6 GB of cache replacing
  40.5 GB of `.mmn`.
- **Later:** A(R) from the `.chk`'s stored M matrices (`right_centre_from_chk`, `transl_inv_full_from_chk`), or
  from Jae-Mo's fixed `_r.dat` (#5).

## Remaining
CLAUDE.md pitfall 2 summarizes the three tiers: stock (unusable), PR 702 (still wrong, #5), and Jae-Mo's branch
(works).

**Sources:** [`../log/2026-08-27_rfull_vs_mmn_position_matrix.md`](../log/2026-08-27_rfull_vs_mmn_position_matrix.md)
§0, §0b; memory `rdat-position-matrix-unusable`.
