#!/usr/bin/env python3
"""Does use_ws_distance's closest-image selection stay consistent under base's own inversion?

_map_R_by_orbital_positions picks, per orbital pair (left,right) and per R, the
periodic image minimizing |R + tau_right - tau_left|, calling two images "tied"
only within 1e-5 A. That tolerance is far tighter than realistic Wannierization
noise (~1e-3 to 1e-2 A). So the real question isn't whether the code detects an
EXACT tie -- it almost never will -- it's whether the GAP between the 1st and
2nd closest candidate is small enough that the "closest wins" choice is fragile
to that noise, and whether the choice made for an orbital pair (s,t) is
consistent, under inversion, with the choice independently made for its exact
symmetry partner (sigma(s), sigma(t)).

Uses the 150/176 base Wannier functions whose inversion partner matched to
0.0000 A (see check_base_inversion_symmetrized.py); the other 26 are same-site-
degenerate and excluded to avoid conflating two different issues.

Run in the axion conda env: python check_ws_tie_consistency.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from wannier.wannier_io import read_wannier_checkpoint

PREFIX = "Y2Ir2O7"
BASE_DIRECTORY = (
    REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
)
MP_GRID = (8, 8, 8)
INVERSION = -np.eye(3)

# ============================================================== 1. inversion partners (as before)
checkpoint = read_wannier_checkpoint(BASE_DIRECTORY / f"{PREFIX}.chk")
centres = np.asarray(checkpoint["wannier_centres"], dtype=float)
lattice = np.asarray(checkpoint["real_lattice"], dtype=float)
n = len(centres)

shell = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
shift_cart = shell @ lattice
target = centres @ INVERSION.T
diff = target[:, None, None, :] - (centres[None, :, None, :] + shift_cart[None, None, :, :])
dist_all = np.linalg.norm(diff, axis=-1)
best_l = dist_all.argmin(axis=-1)
cost = dist_all.min(axis=-1)
row, col = linear_sum_assignment(cost)
defect = cost[row, col]
L = shift_cart[best_l[np.arange(n), col]]

clean = defect < 1e-3  # the 150 well-matched orbitals
clean_idx = np.where(clean)[0]
partner = col
print(f"Clean (unambiguous) inversion partners: {clean.sum()} / {n}")
assert np.array_equal(partner[partner], np.arange(n))  # involution, checked before

# ============================================================== 2. per-(s,t,R) WS selection, both sides
mp_grid_shell = np.array(
    [[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
)
super_shift_cart = (mp_grid_shell * np.asarray(MP_GRID)) @ lattice
TOL = 1e-5


def closest_image(s, t, R_cart):
    """Mirror _map_R_by_orbital_positions's own per-pair selection exactly."""
    displacement = centres[t] - centres[s]
    candidates = R_cart[None, :] + super_shift_cart + displacement
    distances = np.linalg.norm(candidates, axis=-1)
    dmin = distances.min()
    selected = np.where(distances <= dmin + TOL)[0]
    gap = np.sort(distances)[1] - dmin if len(distances) > 1 else np.inf
    return selected, dmin, gap


# candidates_image[i] = -(R+tau_t-tau_s) + shift_i, vs -candidates_original[j] where
# shift_j = -shift_i (worked through by hand): index i on the image side corresponds
# to the index carrying the NEGATED shift on the original side, not the same index.
neg_index = np.array([
    np.where((mp_grid_shell == -mp_grid_shell[i]).all(axis=1))[0][0]
    for i in range(len(mp_grid_shell))
])
assert (mp_grid_shell[neg_index] == -mp_grid_shell).all()


rng = np.random.default_rng(20260930)
n_r_sample = 40
n_pair_sample = 400
# sample R vectors the same way H(R) would see them: within the supercell box
R_int_sample = rng.integers(-4, 4, size=(n_r_sample, 3))
R_cart_sample = R_int_sample @ lattice

pairs = rng.choice(clean_idx, size=(n_pair_sample, 2))

rows = []
for s, t in pairs:
    sp, tp = partner[s], partner[t]
    for R_cart in R_cart_sample:
        sel_st, dmin_st, gap_st = closest_image(s, t, R_cart)
        # inversion image: R' = -R + L_t - L_s (CLAUDE.md pitfall 6), pair (sp, tp)
        R_prime_cart = -R_cart + L[t] - L[s]
        sel_sptp, dmin_sptp, gap_sptp = closest_image(sp, tp, R_prime_cart)

        # candidates_image[i] = -candidates_original[neg_index[i]] (derived by hand):
        # map the image side's selected indices through neg_index before comparing.
        mapped_sptp = {neg_index[i] for i in sel_sptp}
        consistent = set(sel_st.tolist()) == mapped_sptp
        rows.append((gap_st, gap_sptp, consistent, abs(dmin_st - dmin_sptp)))

gaps = np.array([min(r[0], r[1]) for r in rows])
consistent = np.array([r[2] for r in rows])
dmin_diff = np.array([r[3] for r in rows])
print(f"\nSanity check: |dmin_st - dmin_sptp| should be ~0 (same physical distance, both sides) "
      f"-- max={dmin_diff.max():.2e}, mean={dmin_diff.mean():.2e}")

print(f"\nTotal (s,t,R) triples tested: {len(rows)}")
print(f"Consistent under inversion: {consistent.sum()} ({consistent.mean():.2%})")
print(f"Inconsistent: {(~consistent).sum()} ({(~consistent).mean():.2%})")

print("\nGap (to 2nd-closest candidate) distribution, split by consistency:")
for label, mask in (("consistent", consistent), ("inconsistent", ~consistent)):
    if mask.sum() == 0:
        print(f"  {label}: none")
        continue
    g = gaps[mask]
    print(f"  {label:14s} n={mask.sum():5d}  gap min={g.min():.4f}  median={np.median(g):.4f}"
          f"  max={g.max():.4f} A  (fraction with gap<0.01 A: {(g<0.01).mean():.2%})")

# bucket by gap size to see if small-gap pairs are where inconsistency concentrates
print("\nInconsistency rate by gap bucket:")
buckets = [(0, 0.001), (0.001, 0.01), (0.01, 0.1), (0.1, 1.0), (1.0, np.inf)]
for lo, hi in buckets:
    mask = (gaps >= lo) & (gaps < hi)
    if mask.sum() == 0:
        continue
    rate = (~consistent[mask]).mean()
    print(f"  gap in [{lo:.4f}, {hi if hi<np.inf else 'inf':>6}) A: n={mask.sum():5d}  "
          f"inconsistent rate={rate:.2%}")
