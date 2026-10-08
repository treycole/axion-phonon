#!/usr/bin/env python3
"""Does symmetrizing c2 over mode 1's point group change dtheta?

The pointwise curvature misses covariance by up to ~15-20% (base_inversion
total residual, magnetic_symmetry_jaemo_qe76). This tests whether that noise
biases the integral dtheta = sum(c2) * d3k, by comparing it to the integral of
c2 group-averaged over mode 1's own magnetic point group (-4'3m').

No new curvature evaluation is needed: mode 1's point group maps the
reciprocal lattice to itself, so on a Gamma-centered nk^3 mesh every group
image of a mesh point is EXACTLY another mesh point (integer index
permutation, verified below) -- this reuses the already-saved c2_density grids
in nk_sweep_tif_shared_ws/*.npz directly.

Two candidate transformation laws for c2 under a magnetic point-group element
g = (R, antiunitary) are tested EMPIRICALLY against the saved data rather than
trusted from a hand derivation (the antiunitary sign either cancels in c2,
which is quadratic in the curvature, or it doesn't -- easy to get wrong on
paper, and a wrong sign would silently corrupt the result):

  (a) weight(g) = det(R)                        [sign cancels quadratically]
  (b) weight(g) = det(R) * (-1 if antiunitary)   [sign does not cancel]

Whichever fits the measured c2(g.k) vs c2(k) data is used for the projector
average c2_sym(k) = (1/|G|) sum_g weight(g) * c2(g.k).
"""

import sys
from itertools import product
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

NK_SWEEP = REPO / (
    "data/Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/192345/"
    "trial_02_Y_p_Ir_d_O_sp/no-displacement/192908bc889/nk_sweep_tif_shared_ws"
)

# Same Cartesian rotation matrices and antiunitary flags as
# validation/symmetry/check_magnetic_symmetry.ipynb's mode-1 operations.
C2Z = np.diag([-1.0, -1.0, 1.0])
C3_111 = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
MXY = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
S4Z = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
GENERATORS = [(C2Z, False), (C3_111, False), (MXY, True), (S4Z, True)]

# base's unit_cell_cart, Cartesian rows (Angstrom); same lattice for mode 1 (one cell, A3).
LATTICE = np.array([[0.0, 5.095, 5.095], [5.095, 0.0, 5.095], [5.095, 5.095, 0.0]])
RECIPROCAL = 2 * np.pi * np.linalg.inv(LATTICE).T  # rows = reciprocal vectors


# ============================================================== 1. close the group
def close_group(generators, tol=1e-9):
    elements = [(np.eye(3), False)]

    def key(R, a):
        return (tuple(np.round(R, 7).flatten()), a)

    seen = {key(*elements[0])}
    frontier = list(elements)
    while frontier:
        new_frontier = []
        for R1, a1 in frontier:
            for R2, a2 in generators:
                R, a = R1 @ R2, a1 ^ a2
                k = key(R, a)
                if k not in seen:
                    seen.add(k)
                    elements.append((R, a))
                    new_frontier.append((R, a))
        frontier = new_frontier
        if len(elements) > 96:
            raise RuntimeError("group closure did not terminate at a sane order")
    return elements


group = close_group(GENERATORS)
print(f"Closed group: {len(group)} elements")
dets = [np.linalg.det(R) for R, a in group]
n_proper = sum(1 for d in dets if d > 0)
print(f"  proper (det=+1): {n_proper}, improper (det=-1): {len(group) - n_proper}")
n_antiunitary = sum(1 for R, a in group if a)
print(f"  antiunitary: {n_antiunitary}, unitary: {len(group) - n_antiunitary}")
# cross-tab: is det always tied to antiunitary in this group?
tied = all((np.linalg.det(R) < 0) == a for R, a in group)
print(f"  det(R)<0 exactly matches antiunitary for every element: {tied}")

# ============================================================== 2. reduced-coordinate action + exact index map
def reduced_matrix(R):
    """k_reduced' = k_reduced @ G, where g acts on Cartesian k as k -> k @ R.T."""
    G = RECIPROCAL @ R.T @ np.linalg.inv(RECIPROCAL)
    return G


max_nonint = 0.0
G_reduced = []
for R, a in group:
    G = reduced_matrix(R)
    nonint = np.abs(G - np.round(G)).max()
    max_nonint = max(max_nonint, nonint)
    G_reduced.append(np.round(G).astype(int))
print(f"Max deviation from integer reduced-space matrix over the group: {max_nonint:.2e} "
      f"(should be ~0 for a genuine crystallographic operation)")
assert max_nonint < 1e-6, "reciprocal-lattice/rotation convention is wrong"


def index_map(nk, G_int, antiunitary):
    """(i1,i2,i3) -> (j1,j2,j3) such that k_reduced[j] = sign * k_reduced[i] @ G (mod 1)."""
    idx = np.stack(np.meshgrid(*[np.arange(nk)] * 3, indexing="ij"), axis=-1)  # (nk,nk,nk,3)
    sign = -1 if antiunitary else 1
    mapped = sign * np.einsum("...i,ij->...j", idx, G_int)
    return mapped % nk


# ============================================================== 3. empirically pick the weight law
def load_grid(nk):
    data = np.load(NK_SWEEP / f"nk_{nk}.npz")
    return data["c2_density"]  # (nk,nk,nk,2), complex64


print("\n=== empirical weight-law test (nk=12, base-end estimate) ===")
c2 = load_grid(12)[..., 0]  # base end
nk = c2.shape[0]
c2_flat = c2.ravel()

results = {}
for label, weight_fn in (
    ("(a) det(R) only", lambda R, a: np.linalg.det(R)),
    ("(b) det(R) * (-1)^antiunitary", lambda R, a: np.linalg.det(R) * (-1 if a else 1)),
):
    lhs, rhs = [], []
    for (R, a), G_int in zip(group, G_reduced):
        j = index_map(nk, G_int, a)
        c2_g = c2[j[..., 0], j[..., 1], j[..., 2]]
        w = weight_fn(R, a)
        lhs.append(c2_g.ravel())
        rhs.append(w * c2_flat)
    lhs = np.concatenate(lhs)
    rhs = np.concatenate(rhs)
    residual = lhs - rhs
    scale = np.sqrt(np.mean(np.abs(lhs) ** 2) + np.mean(np.abs(rhs) ** 2))
    rel_rms = np.sqrt(np.mean(np.abs(residual) ** 2)) / scale
    corr = np.abs(np.vdot(lhs, rhs)) / (np.linalg.norm(lhs) * np.linalg.norm(rhs) + 1e-300)
    results[label] = rel_rms
    print(f"  {label}: relative RMS residual={rel_rms:.4f}, |correlation|={corr:.4f}")

best = min(results, key=results.get)
print(f"\nBest-fitting law: {best}")
use_antiunitary_sign = best.startswith("(b)")


# ============================================================== 4. symmetrize and re-integrate, all nk
def weight(R, a):
    return np.linalg.det(R) * (-1 if (a and use_antiunitary_sign) else 1)


print(f"\n=== symmetrized vs. plain dtheta (weight law: {best}) ===")
print(f"{'nk':>3s}  {'endpoint':>8s}  {'plain':>12s}  {'symmetrized':>12s}  {'rel. change':>11s}")
for nk_path in sorted(NK_SWEEP.glob("nk_*.npz"), key=lambda p: int(p.stem.split("_")[1])):
    nk = int(nk_path.stem.split("_")[1])
    c2_grid = load_grid(nk)  # (nk,nk,nk,2)
    d3k = 1.0 / nk**3
    for endpoint, label in enumerate(("base", "mode")):
        c2 = c2_grid[..., endpoint]
        c2_sym = np.zeros_like(c2)
        for (R, a), G_int in zip(group, G_reduced):
            j = index_map(nk, G_int, a)
            c2_sym += weight(R, a) * c2[j[..., 0], j[..., 1], j[..., 2]]
        c2_sym /= len(group)
        plain = np.sum(c2) * d3k
        symmetrized = np.sum(c2_sym) * d3k
        rel = abs(symmetrized - plain) / abs(plain) if abs(plain) > 0 else float("nan")
        print(f"{nk:3d}  {label:>8s}  {plain.real:12.6f}  {symmetrized.real:12.6f}  {rel:10.2%}")

print(
    "\nNOTE: the symmetrized total is IDENTICAL to the plain total (0.00% at every nk, every\n"
    "endpoint) -- and that is guaranteed by construction, not a validation of anything. Every\n"
    "element of the group has weight +1 (proper&unitary: det=+1, T=0 -> +1; improper&antiunitary:\n"
    "det=-1, T=1 -> +1), so sum(weight) = |G| exactly, and group-averaging a field over its own\n"
    "orbit under a weight-preserving bijection of the mesh leaves the TOTAL SUM unchanged no\n"
    "matter how badly the field breaks covariance pointwise. This test cannot distinguish a clean\n"
    "dtheta from one biased by the known ~15% pointwise floor -- it was the wrong test for that\n"
    "question. What it DOES show: the empirical weight-law fit above (14% residual, 98%\n"
    "correlation for law (b), ~0% correlation for law (a)) is a real, independent confirmation\n"
    "that the curvature code transforms correctly under the group -- not nothing, just not an\n"
    "answer to 'does the floor bias the integral'."
)

# ============================================================== 5. a real (if rough) noise estimate
# Pointwise residual of the group relation, in ABSOLUTE units of c2 (not the ratio printed above),
# gives a per-point noise scale. Propagating it to the mesh sum as if points were independent is an
# upper bound in one sense (real noise may correlate across k for a fixed orbital pair, which would
# make the true uncertainty on dtheta SMALLER than this estimate, not larger) and optimistic in
# another (it ignores any part of the noise that is itself group-covariant, which this test cannot
# see at all, per the note above). Take it as an order-of-magnitude sanity check, not a real error bar.
print(f"\n=== rough noise-propagation sanity check (nk=12, base-end) ===")
residuals, values = [], []
for (R, a), G_int in zip(group, G_reduced):
    j = index_map(12, G_int, a)
    c2_g = c2[j[..., 0], j[..., 1], j[..., 2]]
    residuals.append((c2_g - weight(R, a) * c2).ravel())
    values.append(c2.ravel())
residuals = np.concatenate(residuals)
sigma_pointwise = np.sqrt(np.mean(np.abs(residuals) ** 2)) / np.sqrt(2)  # split between the two sides
d3k = 1.0 / 12**3
n_points = 12**3
sigma_dtheta_uncorrelated = sigma_pointwise * np.sqrt(n_points) * d3k
plain_dtheta = np.sum(c2).real * d3k
print(f"  per-point noise scale: {sigma_pointwise:.4f} (curvature units)")
print(f"  IF pointwise noise were uncorrelated across all {n_points} mesh points, propagated "
      f"1-sigma on dtheta: {sigma_dtheta_uncorrelated:.5f}")
print(f"  computed dtheta (base end, nk=12): {plain_dtheta:.5f}")
print(f"  ratio |dtheta| / this-sigma: {abs(plain_dtheta)/sigma_dtheta_uncorrelated:.2f}")
print(
    "  (a ratio >> 1 is consistent with the sign/magnitude being real under the uncorrelated-noise\n"
    "  assumption; it is NOT proof, since the noise may not be uncorrelated -- see caveats above.)"
)
