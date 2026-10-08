# #10 — H(R) and A(R) folded on slightly different centre arrays

> **Status: FIXED** (2026-09-23). `load_wannier_pair` folds H(R) on the checkpoint centres (from the A(R) cache
> for `.mmn`), the same array the A(R) fold uses: one centre array per structure.

**Found:** 2026-09-23. **Affects:** the in-house pipeline before 2026-09-23. The effect is tiny.

## Symptom
20 near-tie H(R) entries flipped image (max 2.4e-4 eV) depending on which centres were used.

## Root cause
- **Two sources:** H(R) was folded on PythTB's centres, read from the text output files, while A(R) was folded on
  the `.chk` checkpoint's centres.
- **The difference:** the two arrays differ by ~5e-9 Å (print precision), which is enough to tip entries that sit
  within the 1e-5 Å tie tolerance (#13).

## Why it was a bug
It's a small instance of #3. H(R) and A(R) must be folded by the *same* rule on the *same* centres, otherwise
their combination isn't the curvature of a single interpolated model (CLAUDE.md pitfall 1).

## Fix
`load_wannier_pair` folds H(R) on the checkpoint centres.

## Verification
The effect on dθ at nk = 6 was 2e-7, negligible, but the rule is now consistent. The earlier `.mmn`
shared-centre sweep had used one array for all four folds and was already consistent.

**Sources:** [`../progress/2026-09-23_ws_ties_in_beta_derivative.md`](../progress/2026-09-23_ws_ties_in_beta_derivative.md)
("Also fixed").
