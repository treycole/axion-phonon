# #3 — H(R) and A(R) used different Wigner–Seitz conventions

> **Status: FIXED** (2026-08-20). H(R) and A(R) are folded by the same rule, today
> `wannier.position_matrix.apply_orbital_dependent_R_mapping` (min |R + τ_t − τ_s|). The original fix was
> `ws_reassign()`.

**Found:** 2026-08-20. **Affects:** curvatures built with A(R) from WannierBerri and H(R) from `_hr.dat`
(SrTiO3 trials at the time).

## Symptom
Curvature symmetry violations, and disagreement with WannierBerri. Swapping one ingredient at a time on the
displaced trial_03 line located it:

| what | result |
|---|---|
| our fields + our assembly | corr 0.929 with WB |
| WB formula fed our fields | identical: the assembly was fine |
| our A(k), Ω(k) vs WB | rel 0.00%: the fields were fine |
| PythTB H(R) vs WB H(R) | rel 0.0105%, max 3.3e-4 eV: **the difference** |
| all-WB inputs + WB formula | corr 1.000000, max diff 1.8e-10 |

## Root cause
A finite k-mesh fixes H(R) and A(R) only up to the aliasing class [R] = {R + Σ nᵢNᵢaᵢ}. Choosing which image in
a class carries the value is a *choice of interpolant*. WannierBerri's A(R) used one representative rule and
`_hr.dat`'s H(R) another.

## Why it was a bug
With different rules, H(k) and A(k) are interpolations of different underlying models, so their combination is
not the curvature of any Hamiltonian. It's invisible at the ab-initio k-points, where every rule agrees, and
shows up only on the dense interpolation mesh. This is CLAUDE.md pitfall 1.

## Fix
- **Fold H(R) with A(R)'s rule** before building the model, entirely in-house. Originally
  `ws_reassign()` / `apply_ws_reassignment(w90, mp_grid)` rewrote `ham_r` with `deg = 1`.
- **Prune zero images**, or the R-set balloons from 729 to 15625.
- **The current code** uses `apply_orbital_dependent_R_mapping`, with the same fold applied to A(R).

## Verification
Against WannierBerri: rel 0.0004%, max abs 7.06e-7 (the residual is `_hr.dat`'s 6-decimal printing). The
symmetry violations disappeared: they were this convention bug, not the Wannier sets.

## Related
#10 is a later, much smaller instance: the two folds used slightly different *centre arrays*.

**Sources:** [`../log/2026-08-20_handoff_wannier_conventions.md`](../log/2026-08-20_handoff_wannier_conventions.md)
§1, §3b.
