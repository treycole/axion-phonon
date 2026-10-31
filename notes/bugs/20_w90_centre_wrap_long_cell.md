# 20. wannier90 centres wrap in the long MBT cell, so the WS images are wrong

**Status: FIXED** (2026-10-09) for the MBT LDA base: projection sites recentred (`w90_trial_sp_d_recentred`, job
210651). The evidence is in [`../log/2026-10-08_mbt_noncollinear_gga_lda.md`](../log/2026-10-08_mbt_noncollinear_gga_lda.md)
§11.

## Symptom

The first LDA base Wannierization (210622, 116 WFs, 4×4×4) was exact on the mesh, but its interpolated gap along
the path was **−346 meV** against DFT 198 meV.
- Bi s/p spreads were up to 56 Å². The total was 1978 Å², of which **Ω_D was 1597 Å²**.
- Some centres were ~27 Å from any equivalent atom.
- Changing the windows did nothing (four variants, 210646-210649).

## Root cause

With `num_iter = 0`, a centre is −Σ_b w_b b Im ln M_nn(k,b), and Im ln is defined only modulo 2π.
- In the rhombohedral cell (c = 81.85 Å on the axis), the shortest axial b-vector on the 4×4×4 mesh is
  (b₁+b₂+b₃)/4. Its wrap period is 4c/3 = **109.1 Å**.
- Functions at crystal x > 2/3 (z > 54.6 Å: Bi 58.32/64.46, Te 69.82/76.39) get an inconsistent centre. Example:
  Bi at 58.32 is reported at −50.4 Å (58.32 − 109.1 = −50.8).
- The **WS images are chosen from those centres**, both by wannier90's `use_ws_distance` and by our folding on the
  `.chk` centres (`wannier/wannier_io.py`). Wrong images give exact bands on the mesh and garbage between.
- **Cross-check:** `use_ws_distance = false` raised the gap to +83 meV (210650).

## Fix

In the `.win`, every atom with z > c/2 is placed at z − 81.852 (lattice vector a₁+a₂+a₃), so every projection site
has \|x\| ≤ ½. QE positions are untouched. The new `.amn` reuses the old `.mmn`/`.eig`.
- **Result (210651):**
  - total spread 382 Å²;
  - Ω_D 0.46 Å²;
  - Bi spreads 5-6 Å²;
  - every centre on its atom;
  - path gap +104 meV indirect / 185 meV smallest direct.
- **Modes:** each mode's `.win` must use the base's image choice per atom (assumption A3). The held MBT mode chain
  is updated.

## Remaining

- **Mesh resolution:** a DFT check on Γ-L (210652) shows the 4×4×4 mesh still misses the camel-back valence band near
  Γ by up to ~0.1 eV. That is a mesh-density question, not this bug.
- **Detection:** nothing in wannier90 warns. Check Ω_D and the interpolated gap for every new cell.
