# #5 — Wannier90 PR 702's `_r.dat` silently ignores `use_ws_distance`

> **Status: FIXED** (2026-09-13, production closure 2026-09-16) by Jae-Mo Lihm's `write_ndegen_applied` branch,
> used with `transl_inv_full = T, use_ws_distance = T, write_ndegen_applied = T`. That is
> `source="rdat_ndegen_applied"` in `wannier/position_matrix.py`. Stock PR 702 alone is still unusable.

**Found:** 2026-09-10. **Affects:** `_r.dat` written by Wannier90 with PR 702 (`transl_inv_full = T`).

## Symptom
- **Matches postw90 exactly:** PR 702's `_r.dat` equals postw90's A(R) to 4.4e-16 (machine epsilon).
- **But fails the symmetry test badly:** on the base P+T-forced zero it improves on stock `_r.dat` by 15× yet is
  still **318,000× worse than the `.mmn` route**.
- **The flag does nothing:** the file is byte-for-byte identical whether `use_ws_distance` is T or F in the
  `.win`.

## Root cause
PR 702 adds the translation-equivariant formula (fixing #4's missing link phase). But `use_ws_distance` never
reaches `plot_write_rmn`, so the file is permanently frozen at `use_ws_distance = F`, the wrong Wigner–Seitz
representative.

## Why it was a bug
For symmetry work, A(R) must be at the same (minimum-distance) representative as H(R) (#3). The print format is
*not* the problem: printing the same matrix at `ES24.16E3` changes the symmetry residual by only ~1%. The flag is
silently ignored, so a user setting `use_ws_distance = T` gets the wrong matrix with no warning.

## Fix
- **Jae-Mo's branch** (`jaemolihm/wannier90`, plan5-write-ndegen-applied, HEAD `88134317`) applies the WS
  distance and degeneracy in the writer. With all three flags its `_r.dat` matches postw90's `write_aa_r` at
  `use_ws_distance = T`, bit-exactly after 6-decimal rounding, on two test-suite systems. Hermiticity is exactly
  0, and every row's −R partner is present.
- **Production closure (2026-09-16):** symmetry-exact on SrTiO3. Cross-checked against an independently built
  postw90 `rfull` export to ~1e-8, two separate Fortran implementations agreeing.

## Remaining / caveats
- **Interpolant-family gap:** `rdat_ndegen_applied` still differs from `.mmn` by ~6-38% pointwise. That's the
  known `rfull`-family vs `mmn`-family interpolant gap, not an error.
- **Tie consistency inside Wannier90 is open:** under `transl_inv_full`, images of a class differ by (−1)ⁿ, so
  the base/mode tie fix (#6) can't be applied after the fact. It would need a change inside Wannier90
  (`plot.F90:198`).
- **6-decimal print still costs ~3% of dθ** (#11).

**Sources:** [`../log/2026-09-10_pr702_transl_inv_full_rdat.md`](../log/2026-09-10_pr702_transl_inv_full_rdat.md),
[`../log/2026-09-13_jaemo_transl_inv_full_ws_distance.md`](../log/2026-09-13_jaemo_transl_inv_full_ws_distance.md),
[`../log/2026-09-16_jaemo_rdat_production_closure.md`](../log/2026-09-16_jaemo_rdat_production_closure.md);
CLAUDE.md pitfall 2.
