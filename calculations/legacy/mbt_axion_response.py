#!/usr/bin/env python3
"""
Compute axion response (dtheta/dQ) for all modes of a given irrep at a single nk.

This script processes all modes at once, sharing the base Hamiltonian evaluation
across modes for efficiency.

Usage example:
  python axion_response.py \
    --base-dir base/vdw_U/U_0.0/150300bc889 \
    --irrep T2- --modes mode1 mode2 mode3 mode4 \
    --amp Q_0p01 \
    --nk 20 --symmetrized
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pythtb
from pythtb import W90, TBModel
from pythtb.utils import finite_diff_coeffs, levi_civita
from tensorflow import constant as const
from tensorflow import complex64
import tensorflow as tf
import tensorflow.linalg as tfla

# Relative --base-dir / --irrep / --out-dir paths are resolved against the runs in data/.
DATA = Path(__file__).resolve().parents[2] / "data/MnBi2Te4/phonon"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute axion response for all modes of an irrep."
    )
    p.add_argument("--base-dir", type=str, required=True, help="Path to base W90 data")
    p.add_argument("--base-label", type=str, default=None, help="Label for the base")
    p.add_argument("--irrep", type=str, required=True, help="Irrep name, e.g. T2- or T1+")
    p.add_argument("--modes", type=str, nargs="+", required=True, help="List of mode names, e.g. mode1 mode2 mode3 mode4")
    p.add_argument("--amp", type=str, required=True, help="Amplitude label, e.g. Q_0p01")
    p.add_argument("--mode-subfolder", type=str, default=None, help="Optional subfolder for W90 data inside {amp}/")
    p.add_argument("--nk", type=int, default=20, help="k-mesh density (default: 20)")
    p.add_argument("--d-beta", type=float, default=0.01, help="Finite difference step (default: 0.01)")
    p.add_argument("--n-occ", type=int, default=58, help="Number of occupied states (default: 58)")
    p.add_argument("--symmetrized", action="store_true", help="Use symmetrized W90 model")
    p.add_argument("--batch-size", type=int, default=100, help="k-point batch size (default: 100)")
    p.add_argument("--out-dir", type=str, default=None, help="Override output directory")
    return p.parse_args()


def derive_base_label(base_dir: str) -> str:
    parts = Path(base_dir).parts
    if parts and parts[0] == "base":
        parts = parts[1:]
    return "-".join(parts)


def vel_fd(H_k, mu, dk_mu, order_eps, mode='central'):
    coeffs, stencil = finite_diff_coeffs(order=order_eps, mode=mode)
    print(coeffs)

    fd_sum = np.zeros_like(H_k)

    for s, c in zip(stencil, coeffs):
        fd_sum += c * np.roll(H_k, shift=-s, axis=mu)

    v = fd_sum / (dk_mu)
    return v


def berry_curvature(v_k, H_flat, occ_idxs=None):

    H_flat_tf = const(H_flat, dtype=complex64)
    v_k_tf = const(v_k, dtype=complex64)

    evals_tf, evecs_tf = tfla.eigh(H_flat_tf)

    evecs_tf = tf.transpose(evecs_tf, perm=[0, 1, 3, 2])  # (n_kpts, n_beta, n_state, n_state)
    evecs_T_tf = tf.transpose(evecs_tf, perm=[0, 1, 3, 2])  # (n_kpts, n_beta, n_state, n_state)

    evecs_conj_tf = tf.math.conj(evecs_tf)

    v_k_rot_tf = tf.matmul(
        evecs_conj_tf[None, :, :, :, :],
        tf.matmul(
            v_k_tf,
            evecs_T_tf[None, :, :, :, :]
        )
    )

    n_eigs = evals_tf.shape[-1]
    if occ_idxs is None:
        occ_idxs = np.arange(n_eigs // 2)
    elif occ_idxs == 'all':
        occ_idxs = np.arange(n_eigs)
    else:
        occ_idxs = np.array(occ_idxs)

    cond_idxs = np.setdiff1d(np.arange(n_eigs), occ_idxs)

    delta_E_tf = evals_tf[..., None, :] - evals_tf[..., :, None]
    delta_E_occ_cond_tf = tf.gather(tf.gather(delta_E_tf, occ_idxs, axis=-2), cond_idxs, axis=-1)
    delta_E_cond_occ_tf = tf.gather(tf.gather(delta_E_tf, cond_idxs, axis=-2), occ_idxs, axis=-1)
    inv_delta_E_occ_cond_tf = 1 / delta_E_occ_cond_tf
    inv_delta_E_cond_occ_tf = 1 / delta_E_cond_occ_tf

    v_occ_cond_tf = tf.gather(tf.gather(v_k_rot_tf, occ_idxs, axis=-2), cond_idxs, axis=-1)
    v_cond_occ_tf = tf.gather(tf.gather(v_k_rot_tf, cond_idxs, axis=-2), occ_idxs, axis=-1)
    v_occ_cond_tf = v_occ_cond_tf * inv_delta_E_occ_cond_tf
    v_cond_occ_tf = v_cond_occ_tf * -inv_delta_E_cond_occ_tf

    Q = tf.matmul(v_occ_cond_tf[:, None], v_cond_occ_tf[None, :])
    Q = Q.numpy()

    Omega = 1j * (Q - np.swapaxes(Q, -1, -2).conj())

    return Omega


def main():
    args = parse_args()

    base_label = args.base_label or derive_base_label(args.base_dir)
    sym_suffix = "symmetrized" if args.symmetrized else ""

    base_dir = DATA.parent / args.base_dir / sym_suffix  # base/ sits beside phonon/ in data/MnBi2Te4
    base_dir = Path(str(base_dir))

    # Build mode directories
    mode_dirs = []
    for mode_name in args.modes:
        mode_dir_parts = DATA / args.irrep / mode_name / args.amp
        if args.mode_subfolder:
            mode_dir_parts = mode_dir_parts / args.mode_subfolder
        mode_dir = mode_dir_parts / sym_suffix
        mode_dirs.append(Path(str(mode_dir)))

    if args.out_dir:
        out_dir = DATA / args.out_dir
    else:
        out_dir = DATA / args.irrep / args.amp / base_label
        if args.symmetrized:
            out_dir = out_dir / "symmetrized"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load models
    all_dirs = [base_dir] + mode_dirs
    all_labels = ["base"] + list(args.modes)

    print(f"Base label: {base_label}")
    print(f"Irrep: {args.irrep}, Modes: {args.modes}, Amp: {args.amp}")
    print()

    print("Loading Wannier90 models...")
    w90s = []
    for d, label in zip(all_dirs, all_labels):
        print(f"  {label}: {d}")
        w90 = W90(str(d), r"MnBi2Te4")
        w90s.append(w90)

    print("Building PythTB models...")
    models = []
    for w90 in w90s:
        model = w90.model(
            zero_energy=0, min_hopping_norm=1e-5,
            max_distance=125, ignorable_imaginary_part=None
        )
        models.append(model)

    nk = args.nk
    k_mesh = models[0].k_uniform_mesh([nk, nk, nk])
    d3k = 1 / (nk ** 3)

    n_occ = args.n_occ
    n_modes = len(args.modes)
    batch_size = args.batch_size

    # shape: (mode, k, displacement amp, ...)
    H = np.zeros((n_modes, k_mesh.shape[0], 2, models[0].norb, models[0].norb), dtype=np.complex64)
    # shape: (mode, cartesian dir, k, displacement amp, ...)
    v_k = np.zeros((n_modes, 3, k_mesh.shape[0], 2, models[0].norb, models[0].norb), dtype=np.complex64)

    print("Processing base (H_0)")
    for i in range(0, len(k_mesh), batch_size):
        print(f"  {i} / {len(k_mesh)} k points processed")
        k_batch = k_mesh[i:i + batch_size]

        H_batch = models[0].hamiltonian(k_batch, flatten_spin_axis=True)
        v_batch = models[0].velocity(k_batch, flatten_spin_axis=True)

        H[:, i:i + batch_size, 0, :, :] = H_batch[np.newaxis, ...]
        v_k[:, :, i:i + batch_size, 0, :, :] = v_batch[np.newaxis, ...]

    for mode_idx, model in enumerate(models[1:]):
        print(f"Processing {args.modes[mode_idx]}")

        for i in range(0, len(k_mesh), batch_size):
            print(f"  {i} / {len(k_mesh)} k points processed")
            k_batch = k_mesh[i:i + batch_size]

            H_mode_batch = model.hamiltonian(k_batch, flatten_spin_axis=True)
            v_mode_batch = model.velocity(k_batch, flatten_spin_axis=True)

            H[mode_idx, i:i + batch_size, 1, :, :] = H_mode_batch
            v_k[mode_idx, :, i:i + batch_size, 1, :, :] = v_mode_batch

    print("Computing velocity with finite difference...")
    v_beta = vel_fd(H, mu=2, dk_mu=args.d_beta, order_eps=1, mode='forward')
    v = np.concatenate((v_k, v_beta[:, np.newaxis, ...]), axis=1)

    print("Computing axion response...")
    results = {}
    for mode_idx, mode_name in enumerate(args.modes):
        print(f"  {mode_name}")
        b_curv = berry_curvature(v[mode_idx], H[mode_idx], occ_idxs=range(n_occ))
        epsilon = levi_civita(4, 4)
        chern2_density = np.einsum("ijkl, ij...mn, kl...nm->...", epsilon, b_curv, b_curv) * (1 / (16 * np.pi))
        chern2_density = chern2_density.reshape((nk, nk, nk, 2))
        dtheta = np.sum(chern2_density[:-1, :-1, :-1], axis=(0, 1, 2)) * d3k
        results[mode_name] = dtheta
        print(f"    dtheta = {dtheta}")

    out_path = out_dir / "axion_response.npz"
    print(f"\nSaving results to {out_path}")
    save_dict = {
        "H": H,
        "v_k": v_k,
        "k_mesh": k_mesh,
        "nk": np.int32(nk),
        "d3k": np.float64(d3k),
        "n_occ": np.int32(n_occ),
        "d_beta": np.float64(args.d_beta),
        "symmetrized": np.int8(1 if args.symmetrized else 0),
        "irrep": np.array([args.irrep]),
        "amp": np.array([args.amp]),
        "modes": np.array(args.modes),
        "base_dir": np.array([args.base_dir]),
        "base_label": np.array([base_label]),
    }
    for mode_name, dtheta in results.items():
        save_dict[f"dtheta_{mode_name}"] = dtheta
    np.savez_compressed(out_path, **save_dict)


if __name__ == "__main__":
    main()
