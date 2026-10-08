# #2 — PythTB embedding was a silent no-op (base and mode didn't share τ)

> **Status: FIXED** (2026-08-19). Models are built with `W90.model(orb_vecs=tau)`, and `modules/axion.py:269`
> raises `ValueError` if base and mode τ differ.

**Found:** 2026-08-19 (recorded in `HANDOFF` 2026-08-20). **Affects:** every dθ/dQ computed before the fix (YIO
and the old MBT driver), and the Born-charge validation.

## Symptom
- **Spurious τ dependence:** the internal/cross/external decomposition depended on τ in a way earlier notes
  blamed on "Fourier-truncation τ artefacts".
- **Endpoints disagreed:** for YIO dθ the two endpoint estimates (base end, mode end) had opposite signs, e.g.
  173199: [−0.32, +0.30]; 174901: [−0.13, +0.21].

## Root cause
The scripts did `model.lattice.orb_vecs = tau`. `TBModel.lattice` returns `copy.copy(self._lattice)`
(`pythtb/tbmodel.py:233`), so that line mutated a throwaway copy. The external-term code received `tau`, while H
and v kept each run's own Wannier centres.

## Why it was a bug
- **Within one model:** the curvature pieces mixed two different embeddings, which gives an error linear in τ.
- **In the β finite difference:** base and mode used different embeddings, violating assumption A3 of the
  derivation (`notes/berry_curvature_derivation.md`). `v_β = (h_mode − h_base)/dβ` then contained an embedding
  change, not just ∂H/∂Q.

## Fix
- **Build correctly:** build every model with `W90.model(orb_vecs=tau)` using one shared array, plus a
  post-build assertion.
- **Set τ unconditionally** in the axion path, even with external terms off, because v_β needs the shared τ.
- **Check it:** today `modules/curvature.position_terms` reads τ from the model, and `axion.dtheta` refuses a
  base/mode pair whose τ differ (`modules/axion.py:269`).

## Verification
After the fix, results are exactly τ-invariant: τ = 0, 1×, 2×, 3× and 4× the centres, and a random
non-collinear τ, all give identical totals. The int/cross/ext *split* still shifts with τ (by ±0.3281 between
τ = 0 and the centres), so always state τ when quoting a decomposition.

## Remaining
The MnBi2Te4 driver behind the paper's ~10 rad/Å (`axion_response_nk_sweep.py`) still had the bug and used no
external terms. It was deleted on 2026-09-24; MBT will be recomputed on the new pipeline.

**Sources:** [`../progress/2026-08-20_handoff_wannier_conventions.md`](../progress/2026-08-20_handoff_wannier_conventions.md)
§3a; memory `dtheta-dq-gauge-pitfall`; CLAUDE.md pitfall 3.
