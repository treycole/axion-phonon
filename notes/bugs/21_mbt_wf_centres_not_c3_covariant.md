# 21. MnBi₂Te₄: per-WF centres are not C3-covariant, so the WS fold and τ built on them break the model

**Status: FIXED** (found 2026-10-11; production option 2026-10-12): atom-pinned positions fix it (bands exact).
`load_wannier_pair(..., pin_to_atoms=True)`, set by `PIN_TO_ATOMS = True` in `run_axion.py` for MBT. The curvature "residual" first seen (8-45%) is **not** a defect: the base's
traced curvature is forced to zero by P·{T|½,½,½}, and it is zero to 1e-6 of the per-band scale for every rule.
Evidence:
`archive/diagnostics/2026-10-11_ws_positions/` (`mbt_base_full_group.py`, `mbt_eig_symmetry_fold.py` and logs).

## Symptom

LDA base Wannierization 210651 (recentred, 116 WFs, 4×4×4, projection only), interpolated with the production fold
(WS images and embedding τ on the run's own Wannier centres), checked under the 12 operations of −3m′ at 7 generic k:

| WS positions, tolerance | max \|E(gk) − E(k)\| | \|Tr Ω\| / Σₙ\|Ωₙ\| (forced zero) | relative "residual" (meaningless, see below) |
|---|---|---|---|
| WF centres, 1e-3 Å | **61.6 meV** | 3.2e-6 | 17.1% median |
| WF centres, 5e-3 Å | **64.0 meV** | 2.8e-6 | 38.4% median |
| atom-pinned, 5e-3 Å | **0.000 meV** | 1.3e-6 | 10.8% median |

The full-group table (`mbt_base_full_group.py`) normalizes each residual by max|Tr Ω|, which here is ~1e-4 Å² against
a per-band scale of ~4e2 Å² (`mbt_pt_forced_zero.py`): its percentages compare noise with noise. The base's own
curvature is therefore not a usable symmetry test for MBT; the mode (PT broken) or the eigenvalues are.

## Cause (the band part)

The Bi and Te p_x/p_y Wannier functions sit 0.06-0.15 Å off their atoms **in the plane** (one p_x orbital is not
C3-symmetric, so its centre is not pinned to the axis). C3 maps those centres to points where no centre is (max
deviation 0.144 Å; inversion and C2′ map centres onto centres to 3e-8 Å). A fold that chooses WS images from these
centres is not C3-covariant, and neither is τ. The interpolated bands then break C3 by 62 meV: the model itself is
wrong between mesh points, not only the curvature decomposition.

This is the orbital analogue of YIO's spin-partner centres (#13, 9e-4 Å), but 150× larger, so no tie tolerance can
absorb it; the MBT excess distances have no empty gap at all (`tie_window_check.py`).

## Fix

Pin each WF to its atom for the fold and τ (`POSITIONS = "pinned"`). The interpolated bands become exactly symmetric
(0.000 meV over 12 operations × 7 k).

**Checks on the pinned model:** H(R) Hermitian to 1e-15; under inversion |H(−k)| = |H(k)| with the same orbital on the
partner atom, and the phases H(−k)/H(k) are k-independent to 2e-6, i.e. H is exactly covariant with a constant
representation; an independent Kubo curvature from H alone agrees with `omega_internal`.

## Remaining

- **Production code (done 2026-10-12):** `wannier/wannier_io.py` `pin_centres_to_atoms` (nearest atom image of each base
  centre, atoms from the base `.win`; refuses a centre > 1 Å from every atom) feeds the H and A folds of both structures
  and τ; needs `shared_ws_centres=True` and a `.chk` A(R) source. The H(R) class representative is now each class's
  shortest member (the fold searches ±1 supercell around it; MBT's skewed cell needs that, as the diagnostic did).
  Check `pinned_loader_check.py`/`.log`: MBT 210651 bands 64.042 meV unpinned → **0.000 meV** pinned, same as the
  diagnostic; YIO 212615 folded H identical under the new representative rule (184/512 classes change member, |ΔH| = 0).
- **Open:** the 4×4×4 longest MBT hopping is 55.5 Å; on 8×8×8 it roughly doubles, close to the loader's
  `max_distance = 125` Å. Check it on the 8×8×8 base before dθ (a distance cut is symmetric with pinned τ, but truncates).
- **YIO:** pinned and centre-based folds agree (0.52% vs 0.52% total at 5e-3), so YIO can stay on centres.
