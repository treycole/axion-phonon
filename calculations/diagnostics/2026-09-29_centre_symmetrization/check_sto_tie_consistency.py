#!/usr/bin/env python3
"""STO base analogue of check_ws_tie_consistency.py: is the WS tie-break inversion-consistent?

Cubic perovskite SrTiO3 base (trial03, 48 WFs, data/SrTiO3/base/...), inversion about
the origin (Sr site). Every orbital sits exactly on a high-symmetry Wyckoff position
(max inversion-about-origin defect 8.1e-13, vs YIO's 26/176 same-site-degenerate
orbitals with defect up to 1e-3) -- this is the true high-symmetry limit.

Mirrors check_ws_tie_consistency.py's logic exactly (partner-finding, neg_index
correspondence for the pair-dependent R shift) with BASE_DIRECTORY/PREFIX/lattice/
MP_GRID swapped to STO base.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from wannier.wannier_io import read_wannier_checkpoint

PREFIX = "SrTiO3"
BASE_DIRECTORY = (
    REPO / "data/SrTiO3/base/soc/output/181642bc889/trial_03_Sr_sp_Ti_pd_O_sp/185116bc889"
)
MP_GRID = (8, 8, 8)
INVERSION = -np.eye(3)

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

clean = defect < 1e-3
clean_idx = np.where(clean)[0]
partner = col
print(f"Clean (unambiguous) inversion partners: {clean.sum()} / {n} (max defect {defect.max():.2e})")
assert np.array_equal(partner[partner], np.arange(n))

mp_grid_shell = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
super_shift_cart = (mp_grid_shell * np.asarray(MP_GRID)) @ lattice
TOL = 1e-5


def closest_image(s, t, R_cart):
    displacement = centres[t] - centres[s]
    candidates = R_cart[None, :] + super_shift_cart + displacement
    distances = np.linalg.norm(candidates, axis=-1)
    dmin = distances.min()
    selected = np.where(distances <= dmin + TOL)[0]
    gap = np.sort(distances)[1] - dmin if len(distances) > 1 else np.inf
    return selected, dmin, gap


neg_index = np.array([
    np.where((mp_grid_shell == -mp_grid_shell[i]).all(axis=1))[0][0]
    for i in range(len(mp_grid_shell))
])
assert (mp_grid_shell[neg_index] == -mp_grid_shell).all()

rng = np.random.default_rng(20260930)
n_r_sample = 40
n_pair_sample = 400
R_int_sample = rng.integers(-4, 4, size=(n_r_sample, 3))
R_cart_sample = R_int_sample @ lattice

pairs = rng.choice(clean_idx, size=(n_pair_sample, 2))

rows = []
for s, t in pairs:
    sp, tp = partner[s], partner[t]
    for R_cart in R_cart_sample:
        sel_st, dmin_st, gap_st = closest_image(s, t, R_cart)
        R_prime_cart = -R_cart + L[t] - L[s]
        sel_sptp, dmin_sptp, gap_sptp = closest_image(sp, tp, R_prime_cart)

        mapped_sptp = {neg_index[i] for i in sel_sptp}
        consistent = set(sel_st.tolist()) == mapped_sptp
        rows.append((gap_st, gap_sptp, consistent, abs(dmin_st - dmin_sptp)))

gaps = np.array([min(r[0], r[1]) for r in rows])
consistent = np.array([r[2] for r in rows])
dmin_diff = np.array([r[3] for r in rows])
print(f"\nSanity check |dmin_st - dmin_sptp|: max={dmin_diff.max():.2e}, mean={dmin_diff.mean():.2e}")

print(f"\nTotal (s,t,R) triples tested: {len(rows)}")
print(f"Consistent under inversion: {consistent.sum()} ({consistent.mean():.2%})")
print(f"Inconsistent: {(~consistent).sum()} ({(~consistent).mean():.2%})")

print("\nInconsistency rate by gap bucket:")
buckets = [(0, 0.001), (0.001, 0.01), (0.01, 0.1), (0.1, 1.0), (1.0, np.inf)]
for lo, hi in buckets:
    mask = (gaps >= lo) & (gaps < hi)
    if mask.sum() == 0:
        continue
    rate = (~consistent[mask]).mean()
    print(f"  gap in [{lo:.4f}, {hi if hi < np.inf else 'inf':>6}) A: n={mask.sum():5d}  "
          f"inconsistent rate={rate:.2%}")
