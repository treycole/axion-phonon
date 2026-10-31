r"""Step 6, first part: the tight-binding model of each structure from its Wannier90 run.

Inputs are the arrays of ``axion_phonon.io.wannier90.read_run``: V(k) = U_dis U and the overlaps M(k, b) in the
Wannier gauge from the .chk, the band energies from the .eig. Lengths in Angstrom, energies in eV, k reduced.

    hamiltonian_R        H(R) = (1/N) sum_k e^{-2 pi i k.R} V(k)^dag diag(eps_k) V(k), one block per aliasing class
    position_matrix_R    A(R) = <0i| r |Rj>, Wannier90's transl_inv_full formula from M(k, b)
    fold                 a class onto its Wigner-Seitz images: the R minimising |R + tau_j - tau_i|
    pin_to_atoms         the positions that choose those images: each base centre's nearest atom
    tb_model             H(R), A(R) and the fixed embedding tau -> a PythTB model
    pair_models          all of the above for a base and a displaced structure

A finite N1 x N2 x N3 mesh fixes H(R) and A(R) only up to the aliasing class [R] = {R + sum_i n_i N_i a_i}; how a
class is spread over its images is a choice of interpolant (CLAUDE.md pitfall 1). Here both H and A of both
structures are folded on one set of positions, the base centres' nearest atoms (A3'; bug #21), and that same set is
the embedding tau of both models (A3).
"""
from __future__ import annotations

import itertools

import numpy as np
from pythtb import Lattice, TBModel

WS_TOLERANCE = 5e-3     # A: images closer than this to the shortest are ties, weight split equally (bug #13)
ORIGIN = (0, 0, 0)


def class_representatives(mp_grid, lattice: np.ndarray) -> list[tuple[int, int, int]]:
    """One R per aliasing class: its shortest member (Cartesian). ``fold`` searches +-1 supercell around it."""
    N = np.asarray(mp_grid)
    L = np.array(list(itertools.product(range(-2, 3), repeat=3)))
    reps = []
    for n in np.ndindex(*mp_grid):
        candidates = np.asarray(n) + L * N
        reps.append(tuple(int(x) for x in candidates[np.linalg.norm(candidates @ lattice, axis=1).argmin()]))
    return reps


def hamiltonian_R(V: np.ndarray, energies: np.ndarray, kpts: np.ndarray, R_list) -> dict:
    """H(R) for each R of ``R_list``: the whole class sum (pass one R per class, then ``fold``).

    Full precision; ``_hr.dat`` prints H to 1e-6 eV, which breaks the spin-mixing symmetries (bug #22).
    """
    nk = len(kpts)
    Hk = np.einsum("kbm,kb,kbn->kmn", V.conj(), energies, V, optimize=True)
    phases = np.exp(-2j * np.pi * np.asarray(R_list, float) @ kpts.T) / nk                # (nR, nk)
    HR = (phases @ Hk.reshape(nk, -1)).reshape(len(R_list), *Hk.shape[1:])
    return {tuple(int(x) for x in R): HR[i] for i, R in enumerate(R_list)}


def fold(matrices: dict, lattice: np.ndarray, positions: np.ndarray, mp_grid, tolerance: float = WS_TOLERANCE) -> dict:
    r"""Move each (i, j) entry of each class to the image R + N L (L in {-1, 0, 1}^3) minimising
    |R + N L + tau_j - tau_i|; ties within ``tolerance`` share it equally. Matrix indices are the last two axes."""
    J = len(positions)
    input_R = list(matrices)
    R_cart = np.asarray(input_R, float) @ lattice
    images = np.array(list(itertools.product((-1, 0, 1), repeat=3)), float) * np.asarray(mp_grid, float)
    images_cart = images @ lattice
    index = {R: i for i, R in enumerate(input_R)}
    extra = []
    targets = np.empty((len(input_R), len(images)), int)
    for row, R in enumerate(input_R):
        for col, shift in enumerate(images):
            target = tuple(int(round(R[a] + shift[a])) for a in range(3))
            if target not in index:
                index[target] = len(input_R) + len(extra)
                extra.append(target)
            targets[row, col] = index[target]
    output_R = input_R + extra
    lead = matrices[input_R[0]].shape[:-2]
    out = np.zeros((len(output_R),) + lead + (J, J), complex)
    stack = np.asarray([matrices[R] for R in input_R])
    for i in range(J):
        for j in range(J):
            dist = np.linalg.norm(R_cart[:, None] + images_cart[None] + positions[j] - positions[i], axis=2)
            chosen = dist <= dist.min(axis=1, keepdims=True) + tolerance
            values = stack[..., i, j] / chosen.sum(axis=1).reshape((-1,) + (1,) * len(lead))
            rows, cols = np.nonzero(chosen)
            np.add.at(out, (targets[rows, cols], Ellipsis, i, j), values[rows])
    keep = np.abs(out).reshape(len(output_R), -1).max(axis=1) > 0.0
    return {R: out[n] for n, R in enumerate(output_R) if keep[n]}


def finite_difference_weights(b_cart: np.ndarray) -> np.ndarray:
    """w_b with sum_b w_b b_a b_c = delta_ac, one weight per shell of equal |b| (Marzari-Vanderbilt)."""
    norms = np.linalg.norm(b_cart, axis=1)
    shells, labels = [], np.empty(len(norms), int)
    for n, x in enumerate(norms):
        match = [s for s, ref in enumerate(shells) if abs(x - ref) < 1e-5]
        labels[n] = match[0] if match else len(shells)
        if not match:
            shells.append(x)
    pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    eq = np.array([[(b_cart[labels == s][:, a] * b_cart[labels == s][:, c]).sum() for s in range(len(shells))]
                   for a, c in pairs])
    target = np.array([1.0, 1, 1, 0, 0, 0])
    w, *_ = np.linalg.lstsq(eq, target, rcond=None)
    if np.abs(eq @ w - target).max() > 1e-6:
        raise ValueError("b shells are not complete")
    return w[labels]


def position_matrix_R(M, kpts, neighbours, G, reciprocal, lattice, mp_grid, centres, positions,
                      tolerance: float = WS_TOLERANCE) -> dict:
    r"""A(R) = <0i| r |Rj>, {R: (3, J, J)}: Wannier90's transl_inv_full + write_ndegen_applied formula
    (``hamiltonian_get_rmn``, Jae-Mo Lihm's branch), one b vector at a time:

    1. k-space half of the phase, M'_ij(k, b) = M_ij(k, b) e^{i b.(r_i + r_j)/2}, with this run's centres r;
    2. class sums C_b(R) = (1/N) sum_k M'(k, b) e^{-2 pi i k.R};
    3. ``fold`` on ``positions`` (shared by both structures);
    4. real-space half at the final image, i w_b b e^{-i b.R/2}, summed over b;
    5. the R = 0 diagonal set to the centres.
    """
    nk, J = len(kpts), len(centres)
    b_cart = (kpts[neighbours] + G - kpts[:, None, :]) @ reciprocal
    reference = b_cart[0]
    weights = finite_difference_weights(reference)
    order = np.empty((nk, len(reference)), int)                    # column of M holding reference b at each k
    for ik in range(nk):
        match = np.abs(b_cart[ik][None, :, :] - reference[:, None, :]).max(axis=2) < 1e-6
        if not np.all(match.sum(axis=1) == 1):
            raise ValueError(f"k point {ik}: b vectors differ from those of k point 0")
        order[ik] = match.argmax(axis=1)
    grid_index = np.rint(kpts * np.asarray(mp_grid)).astype(int) % np.asarray(mp_grid)
    if not np.allclose(kpts * np.asarray(mp_grid), np.rint(kpts * np.asarray(mp_grid)), atol=1e-7):
        raise ValueError("k points are not on the Monkhorst-Pack grid")
    total = {}
    for g, b in enumerate(reference):
        half = np.exp(0.5j * centres @ b)
        grid = np.zeros(tuple(mp_grid) + (J, J), complex)
        for ik in range(nk):
            grid[tuple(grid_index[ik])] = half[:, None] * M[ik, order[ik, g]] * half[None, :]
        grid = np.fft.fftn(grid, axes=(0, 1, 2)) / nk
        folded = fold({R: grid[R] for R in np.ndindex(*mp_grid)}, lattice, positions, mp_grid, tolerance)
        del grid
        prefactor = 1j * weights[g] * b
        for R, value in folded.items():
            term = (prefactor * np.exp(-0.5j * b @ (np.asarray(R, float) @ lattice)))[:, None, None] * value[None]
            if R in total:
                total[R] += term
            else:
                total[R] = term
    origin = total.setdefault(ORIGIN, np.zeros((3, J, J), complex))
    origin[:, np.arange(J), np.arange(J)] = centres.T
    return total


def pin_to_atoms(centres, atoms, lattice, max_offset: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Each Wannier centre replaced by its nearest atom image; returns (positions, distance of each centre to it)."""
    shifts = np.array(list(itertools.product((-1, 0, 1), repeat=3))) @ lattice
    candidates = (atoms[:, None, :] + shifts[None]).reshape(-1, 3)
    dist = np.linalg.norm(centres[:, None, :] - candidates[None], axis=2)
    nearest = dist.argmin(axis=1)
    offsets = dist[np.arange(len(nearest)), nearest]
    if offsets.max() > max_offset:
        raise ValueError(f"Wannier function {offsets.argmax()} is {offsets.max():.3f} A from every atom")
    return candidates[nearest], offsets


def tb_model(lattice: np.ndarray, tau: np.ndarray, H_R: dict, A_R: dict | None) -> TBModel:
    """PythTB model with onsite + hoppings from H(R) and the position matrix A(R); ``tau`` (J, 3) reduced.

    Every R of H must have its -R, and H's R set must lie inside A's (pitfall 1: one set of images for both).
    ``A_R=None``: H only (bands, eigenvectors); the curvature code then refuses the model.
    """
    if any(tuple(-x for x in R) not in H_R for R in H_R):
        raise ValueError("H(R) is not closed under R -> -R")
    if A_R is not None and set(H_R) - set(A_R):
        raise ValueError(f"{len(set(H_R) - set(A_R))} R vectors of H(R) are missing from A(R)")
    defect = 0.0 if A_R is None else max(np.abs(X - A_R[tuple(-x for x in R)].conj().transpose(0, 2, 1)).max()
                                         for R, X in A_R.items())
    if defect > 1e-8:
        raise ValueError(f"A(R) != A(-R)^dag by {defect:.1e} A")
    J = len(tau)
    model = TBModel(Lattice(lattice, tau, periodic_dirs=...))
    model._from_w90 = True
    model._assume_position_operator_diagonal = False
    onsite = np.diag(H_R[ORIGIN])
    if np.abs(onsite.imag).max() > 1e-9:
        raise ValueError("onsite energies are not real")
    model.set_onsite(onsite.real)
    rows, blocks, masks = [], [], []
    for R, H in H_R.items():
        if R == ORIGIN:
            keep = np.triu(np.ones((J, J), bool), k=1)
        elif R > ORIGIN:                    # one of each (R, -R) pair; PythTB adds the conjugate
            keep = np.ones((J, J), bool)
        else:
            continue
        rows.append(R)
        blocks.append(np.where(keep, H, 0))
        masks.append(keep)
    model._hops.set_rows(np.array(rows, int), np.array(blocks), np.array(masks))
    if A_R is None:
        return model
    zero = np.zeros((J, J), complex)
    model._ham_r = {R: H_R.get(R, zero) for R in A_R}
    model._ham_deg = {R: 1 for R in A_R}
    model._pos_r = {R: A_R[R] for R in sorted(A_R)}
    return model


def pair_models(base: dict, mode: dict, log=print) -> tuple[TBModel, TBModel]:
    """The base and displaced-structure models, both folded on the base centres' nearest atoms, which are also
    their one embedding tau (A3, A3'). ``base`` and ``mode`` are ``read_run`` dicts; their M is used and freed."""
    lattice, grid = base["lattice"], tuple(base["mp_grid"])
    if not np.array_equal(lattice, mode["lattice"]) or grid != tuple(mode["mp_grid"]):
        raise ValueError("base and mode differ in lattice or k mesh")
    if base["V"].shape[2] != mode["V"].shape[2]:
        raise ValueError("base and mode differ in number of Wannier functions")
    positions, offsets = pin_to_atoms(base["centres"], base["atoms"], lattice)
    log(f"pinned to atoms: centre offsets max {offsets.max():.3e} A, median {np.median(offsets):.3e} A")
    reps = class_representatives(grid, lattice)
    tau = positions @ np.linalg.inv(lattice)
    models = []
    for label, run in (("base", base), ("mode", mode)):
        H = fold(hamiltonian_R(run["V"], run["energies"], run["kpts"], reps), lattice, positions, grid)
        A = position_matrix_R(run["M"], run["kpts"], run["neighbours"], run["G"], run["reciprocal"], lattice, grid,
                              run["centres"], positions)
        run.pop("M")
        log(f"{label}: H(R) on {len(H)} R, A(R) on {len(A)} R")
        models.append(tb_model(lattice, tau, H, A))
    return models[0], models[1]


def run_model(run: dict, pin: bool = True, position_matrix: bool = False) -> TBModel:
    """The model of one Wannier run alone (band structures, parities): H(R) folded on the run's own positions, the
    centres' nearest atoms (``pin``) or the centres themselves; with A(R) if ``position_matrix`` (needs ``run["M"]``).
    For dtheta use ``pair_models``: both ends must share one set of positions."""
    lattice, grid = run["lattice"], tuple(run["mp_grid"])
    positions = pin_to_atoms(run["centres"], run["atoms"], lattice)[0] if pin else run["centres"]
    H = fold(hamiltonian_R(run["V"], run["energies"], run["kpts"], class_representatives(grid, lattice)), lattice,
             positions, grid)
    A = None
    if position_matrix:
        A = position_matrix_R(run["M"], run["kpts"], run["neighbours"], run["G"], run["reciprocal"], lattice, grid,
                              run["centres"], positions)
    return tb_model(lattice, positions @ np.linalg.inv(lattice), H, A)
