# #6 — Base and mode didn't share a Wigner–Seitz rule (assumption A3′)

> **Status: FIXED for the in-house folds** (2026-09-23): `load_wannier_pair(shared_ws_centres=True)` folds both
> structures on base's centres, and `calculations/run_axion.py` uses `SHARED_WS = True`. **Open** for the
> Fortran `rdat_ndegen_applied` route, where the fix would have to go inside Wannier90.

**Found:** 2026-09-20 (diagnosed and fixed 2026-09-23). **Affects:** every dθ from a base/mode pair folded on
each structure's own centres.

## Symptom
- **5.5% of |ΔH|:** on the first YIO base/mode pair through Jae-Mo's route, the new `ws_pair_consistency`
  check found that 5.5% of |ΔH| came from aliasing classes split differently in base and mode. The `.mmn` route
  carried the same inconsistency.
- **Mesh dependence:** dθ depended on nk in a way that tracked the mesh.

## Root cause
Each run's fold used its *own* Wannier centres to choose and weight images (`apply_orbital_dependent_R_mapping`,
and `checkpoint["wannier_centres"]` in the A(R) construction), while the Bloch phase used the shared base τ.
- **In the base,** symmetry makes many bonds exactly tied, so the weight is split (½, ½).
- **In mode 1,** atoms move by up to 0.015 Å, so ties break and the same bond becomes (1, 0).

## Why it was a bug
The derivation's β derivative (§8.2) differentiates at fixed bond vectors b_st(R). That needs τ fixed (A3)
**and** the image weights fixed:

> **A3′.** τ_s and the image weights w_st(R) of every aliasing class are the same for both structures of the
> finite difference.

- **The extra term:** a tie split (½, ½) in base and (1, 0) in mode adds (h / 2Δβ)(e^{ik·b1} − e^{ik·b2}) to
  ∂_β H(k), and the same to ∂_β A(k).
- **Where it vanishes:** this is exactly zero on the ab-initio mesh and O(1/Δβ) off it. With `mp_grid` = 8,
  nk = 4 and 8 are tie-free, while 6, 10 and 12 are not.
- **What it means:** a change of interpolant divided by Δβ, not physics.

## Fix
- **`load_wannier_pair(shared_ws_centres=True)`** folds H(R) and A(R) of **both** structures on the base centres
  (one array for all four folds, which is also the embedding τ).
- **A(R) re-fold:** `remap_A_R_to_centres` re-folds A(R). This is exact for `.mmn`-family A(R), where the images
  of a class are copies.
- **Check:** `ws_pair_consistency` must then report zero changed entries (asserted).
- **Valid when:** one lattice, one-to-one WFs (same projections, order, rigid-shift gauge), small centre motion
  (YIO mode 1: max |Δr̄| ≈ 0.015 Å), and an `.mmn`-family A(R).
- **Used in:** the `.chk` routes (`right_centre_from_chk`, `transl_inv_full_from_chk`) with `SHARED_WS`
  (2026-09-29 onward).

## Remaining
- **Refused for `rfull` and `rdat_ndegen_applied`:** under `transl_inv_full`, images differ by (−1)ⁿ, so class
  sums can't be redistributed afterwards. That fix would have to be in Wannier90 (`plot.F90:198`: reference
  centres for the WS rule, actual centres for phases and R = 0). Not written.
- **Not a fix for the family gap:** it doesn't settle the `.mmn` vs `rfull` interpolant-family gap (factor ~2 in
  dθ). It only makes each family consistent across β.

**Sources:** [`../progress/2026-09-20_jaemo_ndegen_yio_mode1.md`](../progress/2026-09-20_jaemo_ndegen_yio_mode1.md),
[`../progress/2026-09-23_ws_ties_in_beta_derivative.md`](../progress/2026-09-23_ws_ties_in_beta_derivative.md);
unit test `validation/unit/test_shared_ws_centres.py`.
