#!/usr/bin/env python3
"""Derive, from the actual DFT geometry (not literature recall), the full magnetic
point group of YIO base and the subgroup each phonon mode (1,2,3) preserves.

Procedure, for each of the 48 Cartesian operations of Oh:
  1. Crystallographic check: does g map the full 22-atom basis (4Y+4Ir+14O) onto
     itself modulo the FCC lattice (allowing a per-atom lattice shift)?  If not,
     g is not a symmetry of the lattice at this origin at all -- discard.
  2. Magnetic check (base): does g map the Ir moment pattern (pseudovector) onto
     itself WITHOUT time reversal (unprimed), or only WITH time reversal (primed),
     or neither (not a symmetry of base at all)?
  3. Mode check: does g additionally map the mode's displacement pattern (true
     vector, T-even, so independent of priming) onto itself?  If yes, g (with the
     priming found in step 2) survives in that mode's magnetic point group.

Moments used (base Y2Ir2O7.scf.out, converged AIAO, exact printed values):
  Ir1 (0,0,0):                 (-0.172951,-0.172951,-0.172951)
  Ir2 (5.095,7.6425,2.5475):   (-0.172951, 0.172951, 0.172951)
  Ir3 (7.6425,5.095,2.5475):   ( 0.172951,-0.172951, 0.172951)
  Ir4 (2.5475,2.5475,0):       ( 0.172951, 0.172951,-0.172951)
"""

import itertools
import numpy as np

np.set_printoptions(suppress=True, precision=6)

LATTICE = np.array([
    [0.000000, 5.095000, 5.095000],
    [5.095000, 0.000000, 5.095000],
    [5.095000, 5.095000, 0.000000],
])

ATOMS = {
    "Y1": (5.09500000, 5.09500000, 5.09500000),
    "Y2": (5.09500000, 2.54750000, 2.54750000),
    "Y3": (2.54750000, 5.09500000, 2.54750000),
    "Y4": (7.64250000, 7.64250000, 5.09500000),
    "Ir1": (0.00000000, 0.00000000, 0.00000000),
    "Ir2": (5.09500000, 7.64250000, 2.54750000),
    "Ir3": (7.64250000, 5.09500000, 2.54750000),
    "Ir4": (2.54750000, 2.54750000, 0.00000000),
    "O1": (3.82125000, 3.82125000, 3.82125000),
    "O2": (6.36875000, 6.36875000, 6.36875000),
    "O3": (3.44352708, 6.36875000, 6.36875000),
    "O4": (9.29397292, 6.36875000, 6.36875000),
    "O5": (6.36875000, 6.36875000, 3.44352708),
    "O6": (6.36875000, 6.36875000, 9.29397292),
    "O7": (6.36875000, 3.44352708, 6.36875000),
    "O8": (6.36875000, 9.29397292, 6.36875000),
    "O9": (3.82125000, 6.74647292, 3.82125000),
    "O10": (3.82125000, 0.89602708, 3.82125000),
    "O11": (6.74647292, 3.82125000, 3.82125000),
    "O12": (0.89602708, 3.82125000, 3.82125000),
    "O13": (3.82125000, 3.82125000, 6.74647292),
    "O14": (3.82125000, 3.82125000, 0.89602708),
}
SPECIES = {name: name[:2] if name[:2] == "Ir" else name[0] for name in ATOMS}

IR_MOMENTS = {
    "Ir1": (-0.172951, -0.172951, -0.172951),
    "Ir2": (-0.172951, 0.172951, 0.172951),
    "Ir3": (0.172951, -0.172951, 0.172951),
    "Ir4": (0.172951, 0.172951, -0.172951),
}

# mode displacements (mode2.win atoms_cart) - base.win atoms_cart; Y and O unchanged for mode2
MODE_DISP = {
    "mode1": {  # Y sublattice <111> breathing (mode1 scf.in - base.win)
        "Y1": (0.01, 0.01, 0.01), "Y2": (0.01, -0.01, -0.01),
        "Y3": (-0.01, 0.01, -0.01), "Y4": (-0.01, -0.01, 0.01),
    },
    "mode2": {  # Ir sublattice <111> breathing (mode2 scf.in - base.win)
        "Ir1": (0.01, 0.01, 0.01), "Ir2": (0.01, -0.01, -0.01),
        "Ir3": (-0.01, 0.01, -0.01), "Ir4": (-0.01, -0.01, 0.01),
    },
    "mode3": {  # O 48f sublattice <100> pairs (mode3.win atoms_cart - base.win atoms_cart)
        "O3": (-0.01, 0, 0), "O4": (0.01, 0, 0),
        "O5": (0, 0, -0.01), "O6": (0, 0, 0.01),
        "O7": (0, -0.01, 0), "O8": (0, 0.01, 0),
        "O9": (0, -0.01, 0), "O10": (0, 0.01, 0),
        "O11": (-0.01, 0, 0), "O12": (0.01, 0, 0),
        "O13": (0, 0, -0.01), "O14": (0, 0, 0.01),
    },
}

names = list(ATOMS)
positions = np.array([ATOMS[n] for n in names])
n_atoms = len(names)
lattice_inv = np.linalg.inv(LATTICE.T)  # columns are lattice vectors a1,a2,a3


def frac(cart):
    return (lattice_inv @ cart.T).T


def build_Oh():
    """48 elements of Oh as Cartesian 3x3 matrices, generated from C4z, C3_111, inversion."""
    c4z = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    c3_111 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
    inv = -np.eye(3)
    gens = [c4z, c3_111, inv]
    group = {tuple(np.round(np.eye(3), 8).flatten())}
    elements = [np.eye(3)]
    frontier = [np.eye(3)]
    while frontier:
        new_frontier = []
        for g in frontier:
            for gen in gens:
                h = np.round(g @ gen, 8)
                key = tuple(h.flatten())
                if key not in group:
                    group.add(key)
                    elements.append(h)
                    new_frontier.append(h)
        frontier = new_frontier
    assert len(elements) == 48, f"got {len(elements)} elements, expected 48"
    return elements


TAU_GRID = [n / 4.0 for n in range(4)]  # quarter-fractions: covers Fd-3m's non-primitive translations


def _try_tau(g, tau_frac, tol):
    """tau_frac: fractional-coordinate offset shared by every atom for this g.
    Returns (sigma, L) if g@r_i + tau (mod lattice) consistently maps the basis
    onto itself (species-respecting bijection), else None."""
    tau_cart = tau_frac @ LATTICE
    transformed = positions @ g.T + tau_cart
    sigma = {}
    L = {}
    for i, name in enumerate(names):
        target_frac = frac(transformed[i:i + 1])[0]
        best = None
        for j, cand in enumerate(names):
            if SPECIES[cand] != SPECIES[name]:
                continue
            cand_frac = frac(positions[j:j + 1])[0]
            shift = target_frac - cand_frac
            shift_round = np.round(shift)
            residual = np.linalg.norm((shift - shift_round) @ LATTICE)
            if residual < tol:
                if best is None or residual < best[2]:
                    best = (cand, shift_round, residual)
        if best is None:
            return None
        sigma[name] = best[0]
        L[name] = best[1] @ LATTICE
    if len(set(sigma.values())) != n_atoms:
        return None
    return sigma, L


def find_atom_mapping(g, tol=1e-3):
    """Search a quarter-fraction grid of shared non-primitive translations tau(g)
    (Fd-3m is non-symmorphic) for one that makes g a consistent symmetry of the
    22-atom basis.  Returns (sigma, L) for the first tau that works, else None."""
    for t1, t2, t3 in itertools.product(TAU_GRID, TAU_GRID, TAU_GRID):
        result = _try_tau(g, np.array([t1, t2, t3]), tol)
        if result is not None:
            return result
    return None


def test_moments(g, sigma, L, det_g):
    """Does g map the Ir moment pattern onto itself, unprimed or primed? Returns
    'unprimed', 'primed', or None."""
    max_err_unprimed = 0.0
    max_err_primed = 0.0
    for name, m in IR_MOMENTS.items():
        m = np.array(m)
        target = sigma[name]
        if target not in IR_MOMENTS:
            return None
        predicted_unprimed = det_g * (g @ m)
        predicted_primed = -det_g * (g @ m)
        actual = np.array(IR_MOMENTS[target])
        max_err_unprimed = max(max_err_unprimed, np.linalg.norm(predicted_unprimed - actual))
        max_err_primed = max(max_err_primed, np.linalg.norm(predicted_primed - actual))
    tol = 1e-3
    if max_err_unprimed < tol:
        return "unprimed"
    if max_err_primed < tol:
        return "primed"
    return None


def test_displacement(g, sigma, disp_dict, tol=1e-3):
    """Does g map the mode displacement pattern onto itself (true vector, no det factor)?"""
    for name, d in disp_dict.items():
        d = np.array(d)
        target = sigma.get(name)
        predicted = g @ d
        actual = np.array(disp_dict.get(target, (0.0, 0.0, 0.0)))
        if np.linalg.norm(predicted - actual) > tol:
            return False
        # also require atoms that DON'T move stay un-moved under sigma
    # and atoms not in disp_dict (i.e. displacement 0) must map to atoms also not in disp_dict
    for name in names:
        if name not in disp_dict:
            target = sigma[name]
            if target in disp_dict and np.linalg.norm(np.array(disp_dict[target])) > tol:
                return False
    return True


elements = build_Oh()
print(f"Built Oh: {len(elements)} elements\n")

base_ops = []  # (g, priming, sigma, L)
for g in elements:
    mapping = find_atom_mapping(g)
    if mapping is None:
        continue
    sigma, L = mapping
    det_g = np.linalg.det(g)
    priming = test_moments(g, sigma, L, det_g)
    if priming is not None:
        base_ops.append((g, priming, sigma, L))

print(f"base (m-3m' candidate): {len(base_ops)} operations survive crystallographic + "
      f"magnetic-moment check (expect 48 for the full group)\n")
n_unprimed = sum(1 for _, p, _, _ in base_ops if p == "unprimed")
n_primed = sum(1 for _, p, _, _ in base_ops if p == "primed")
print(f"  unprimed: {n_unprimed}, primed: {n_primed}\n")

mode_results = {}
for mode, disp in MODE_DISP.items():
    survivors = []
    for g, priming, sigma, L in base_ops:
        if test_displacement(g, sigma, disp):
            survivors.append((g, priming))
    mode_results[mode] = survivors
    n_u = sum(1 for _, p in survivors if p == "unprimed")
    n_p = sum(1 for _, p in survivors if p == "primed")
    print(f"{mode}: {len(survivors)} operations survive (unprimed {n_u}, primed {n_p})")

np.savez(
    "/Users/treycole/.claude/jobs/cad13fc3/tmp/point_groups.npz",
    base_ops=np.array([g for g, p, s, l in base_ops]),
    base_priming=np.array([p for g, p, s, l in base_ops]),
    **{f"{mode}_ops": np.array([g for g, p in ops]) for mode, ops in mode_results.items()},
    **{f"{mode}_priming": np.array([p for g, p in ops]) for mode, ops in mode_results.items()},
)
print("\nSaved operation matrices + priming to point_groups.npz")
