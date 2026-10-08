#!/usr/bin/env python3
"""
Compute the second Chern density on the kz=0 plane for a given base/mode pair.

Outputs:
  - <out_dir>/c2_density.npz
      contains c2_density, k_mesh, k_cart, and metadata
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import tensorflow as tf
from tensorflow import complex64, constant as const
import tensorflow.linalg as tfla
from pythtb.utils import levi_civita
from pythtb import W90

#########################
# Configuration — edit these before running

# base calculation (H_0)
base_dir_rel = "base/vdw_U/U_0.0/150300bc889"
base_label = "vdw_U-U_0.0-150300bc889"

# phonon symmetry mode
irrep = "T2-"
mode = "mode1"
amp = "Q_0p01"
mode_subfolder = None   # set if W90 data is in a subfolder inside {amp}/

# calculation parameters
nk_plane = 100          # grid density for kz=0 plane
d_beta = 0.01           # [Angstroms] finite difference step in beta
k_batch_ham = 100
k_batch_curv = 50
symmetrized = True
n_occ = 58

#########################


def berry_curvature(
        v_k: np.ndarray,
        h_flat: np.ndarray,
        occ_idxs: Iterable[int]
        ) -> np.ndarray:
    """
    Compute non-Abelian Berry curvature.

    Args:
        v_k: shape (4, nk_batch, n_state, n_state)
        h_flat: shape (nk_batch, n_state, n_state)
        occ_idxs: occupied state indices

    Returns:
        Omega: shape (4, 4, nk_batch, n_occ, n_occ)
    """
    h_flat_tf = const(h_flat, dtype=complex64)
    v_k_tf = const(v_k, dtype=complex64)

    evals_tf, evecs_tf = tfla.eigh(h_flat_tf)
    evecs_tf = tf.transpose(evecs_tf, perm=[0, 2, 1])
    evecs_t_tf = tf.transpose(evecs_tf, perm=[0, 2, 1])
    evecs_conj_tf = tf.math.conj(evecs_tf)

    v_k_rot_tf = tf.matmul(
        evecs_conj_tf[None, :, :, :],
        tf.matmul(v_k_tf, evecs_t_tf[None, :, :, :]),
    )

    n_eigs = evals_tf.shape[-1]
    occ_idxs = np.array(list(occ_idxs), dtype=np.int32)
    cond_idxs = np.setdiff1d(np.arange(n_eigs, dtype=np.int32), occ_idxs)

    E_occ = tf.gather(evals_tf, occ_idxs, axis=-1)
    E_cond = tf.gather(evals_tf, cond_idxs, axis=-1)

    delta_occ_cond = E_occ[..., :, None] - E_cond[..., None, :]
    inv_delta_occ_cond = tf.math.reciprocal(delta_occ_cond)

    rank = inv_delta_occ_cond.shape.rank
    perm = list(range(rank))
    perm[-2], perm[-1] = perm[-1], perm[-2]
    inv_delta_cond_occ = tf.transpose(inv_delta_occ_cond, perm=perm)

    v_occ_cond_tf = tf.gather(tf.gather(v_k_rot_tf, occ_idxs, axis=-2), cond_idxs, axis=-1)
    v_cond_occ_tf = tf.gather(tf.gather(v_k_rot_tf, cond_idxs, axis=-2), occ_idxs, axis=-1)
    v_occ_cond_tf *= inv_delta_occ_cond
    v_cond_occ_tf *= inv_delta_cond_occ

    q = tf.matmul(v_occ_cond_tf[:, None], v_cond_occ_tf[None, :]).numpy()
    return 1j * (q - np.swapaxes(q, -1, -2).conj())


def second_Chern_dens(
    k_mesh: np.ndarray,
    model_base,
    model_mode,
    n_occ: int,
    d_beta: float,
    k_batch_ham: int,
    k_batch_curv: int,
) -> np.ndarray:

    n_k = k_mesh.shape[0]
    n_orb = model_base.norb

    h = np.zeros((n_k, 2, n_orb, n_orb), dtype=np.complex64)
    v_cart = np.zeros((3, n_k, 2, n_orb, n_orb), dtype=np.complex64)

    for i in range(0, n_k, k_batch_ham):
        print(f"  Building Hamiltonian and velocity for k points {i} to {min(i+k_batch_ham, n_k)}...")
        j = min(i + k_batch_ham, n_k)
        kb = k_mesh[i:j]

        h0 = model_base.hamiltonian(kb, flatten_spin_axis=True).astype(np.complex64, copy=False)
        v0 = model_base.velocity(kb, flatten_spin_axis=True).astype(np.complex64, copy=False)
        h1 = model_mode.hamiltonian(kb, flatten_spin_axis=True).astype(np.complex64, copy=False)
        v1 = model_mode.velocity(kb, flatten_spin_axis=True).astype(np.complex64, copy=False)

        h[i:j, 0, :, :] = h0
        h[i:j, 1, :, :] = h1
        v_cart[:, i:j, 0, :, :] = v0
        v_cart[:, i:j, 1, :, :] = v1

    v_beta = (h[:, 1, :, :] - h[:, 0, :, :]) / d_beta

    epsilon = levi_civita(4, 4)
    chern2_flat = np.zeros((n_k, 2), dtype=np.complex64)
    occ = np.arange(n_occ, dtype=np.int32)

    for i in range(0, n_k, k_batch_curv):
        print(f"  Computing Berry curvature and Chern2 density for k points {i} to {min(i+k_batch_curv, n_k)}...")
        j = min(i + k_batch_curv, n_k)
        for beta in (0, 1):
            v4 = np.empty((4, j - i, n_orb, n_orb), dtype=np.complex64)
            v4[:3, :, :, :] = v_cart[:, i:j, beta, :, :]
            v4[3, :, :, :] = v_beta[i:j, :, :]

            b_curv = berry_curvature(v4, h[i:j, beta, :, :], occ)
            chern2_flat[i:j, beta] = (
                np.einsum("ijkl,ij...mn,kl...nm->...", epsilon, b_curv, b_curv) * (1.0 / (16.0 * np.pi))
            )

    return chern2_flat


##########################
# Resolve directories

sym_suffix = "symmetrized" if symmetrized else ""

DATA = Path(__file__).resolve().parents[2] / "data/MnBi2Te4/phonon"

base_dir = DATA.parent / base_dir_rel / sym_suffix  # base/ sits beside phonon/ in data/MnBi2Te4
base_dir = Path(str(base_dir))

mode_dir_parts = DATA / irrep / mode / amp
if mode_subfolder:
    mode_dir_parts = mode_dir_parts / mode_subfolder
mode_dir = mode_dir_parts / sym_suffix
mode_dir = Path(str(mode_dir))

out_dir = DATA / irrep / mode / amp / "c2_density" / base_label
if symmetrized:
    out_dir = out_dir / "symmetrized"
out_dir.mkdir(parents=True, exist_ok=True)

required = [
    base_dir / "MnBi2Te4_hr.dat",
    base_dir / "MnBi2Te4.win",
    mode_dir / "MnBi2Te4_hr.dat",
    mode_dir / "MnBi2Te4.win",
]
missing = [str(p) for p in required if not p.exists()]
if missing:
    raise FileNotFoundError(
        "Missing required Wannier files:\n"
        + "\n".join(missing)
        + "\nRun corresponding wannierberri.ipynb(s) first."
    )

print(f"Base label:  {base_label}")
print(f"Base dir:    {base_dir}")
print(f"Mode dir:    {mode_dir}")
print(f"Output dir:  {out_dir}")
print(f"Irrep={irrep}  mode={mode}  amp={amp}")
print()

print(f"Loading W90 models from {base_dir} and {mode_dir}")
w90_base = W90(str(base_dir), r"MnBi2Te4")
w90_mode = W90(str(mode_dir), r"MnBi2Te4")

print("Building TB models")
model_base = w90_base.model(
    zero_energy=0,
    min_hopping_norm=1e-5,
    max_distance=125,
    ignorable_imaginary_part=None,
)
model_mode = w90_mode.model(
    zero_energy=0,
    min_hopping_norm=1e-5,
    max_distance=125,
    ignorable_imaginary_part=None,
)

B = model_base.recip_lat_vecs

# gamma-centered grid in reduced coordinates
k1 = np.linspace(-0.5, 0.5, nk_plane, endpoint=False)
k2 = np.linspace(-0.5, 0.5, nk_plane, endpoint=False)

g1, g2 = np.meshgrid(k1, k2, indexing="ij")
g1_flat = g1.ravel()
g2_flat = g2.ravel()

# Solve for k3 such that kz = 0 in Cartesian coordinates
g3_flat = -(g1_flat + g2_flat)

# Keep only points in first BZ
mask = np.abs(g3_flat) <= 0.5
k_mesh = np.stack([g1_flat[mask], g2_flat[mask], g3_flat[mask]], axis=-1)
k_cart = k_mesh @ B
assert np.allclose(k_cart[:, 2], 0, atol=1e-6), "k points are not in the kz=0 plane"

c2_density = second_Chern_dens(
    k_mesh=k_mesh,
    model_base=model_base,
    model_mode=model_mode,
    n_occ=n_occ,
    d_beta=d_beta,
    k_batch_ham=k_batch_ham,
    k_batch_curv=k_batch_curv,
)

out_path = out_dir / "c2_density.npz"
print(f"Saving Chern2 density to {out_path}")
np.savez(
    out_path,
    c2_density=c2_density,
    k_mesh=k_mesh,
    k_cart=k_cart,
    nk=np.int32(nk_plane),
    symmetrized=np.int8(1 if symmetrized else 0),
    n_occ=np.int32(n_occ),
    d_beta=np.float64(d_beta),
    irrep=np.array([irrep]),
    mode=np.array([mode]),
    amp=np.array([amp]),
    base_dir=np.array([base_dir_rel]),
    base_label=np.array([base_label]),
)
