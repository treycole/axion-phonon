#!/usr/bin/env python3
"""Cross-check: does WannierBerri's own space-group symmetrizer agree with our H(R)/A(R)?

Mirrors the MBT precedent (data/MnBi2Te4/base/vdw_U/U_0.0/150300bc889/wannierberri.ipynb):
build a WannierBerri System_w90 from the same Wannier90 run (win/chk/eig -- no .mmn/.amn
needed, WannierBerri's .chk reader recovers the overlap data this project's Wannier90 fork
already writes into the checkpoint), then call systems.symmetrize(...) with YIO base's
actual magnetic point group input (AIAO Ir moments, read straight from this run's own
Y2Ir2O7.scf.out) and compare the symmetrized Ham_R/AA_R to the unsymmetrized ones.

Ir magnetic moments (atoms 5-8 in .win atom order == atoms 5-8 in scf.out), converged
values from data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/Y2Ir2O7.scf.out:
  Ir1 (-0.172951, -0.172951, -0.172951)   [-1,-1,-1]
  Ir2 (-0.172951,  0.172951,  0.172951)   [-1,+1,+1]
  Ir3 ( 0.172951, -0.172951,  0.172951)   [+1,-1,+1]
  Ir4 ( 0.172951,  0.172951, -0.172951)   [+1,+1,-1]
Y and O sites carry only tiny induced moments (~0.0013 and ~0.013 uB); treated as 0 for
the point-group search, same as the MBT notebook zeroed all non-Mn sites.

Run in the axion conda env: python yio_wannierberri_symmetrize.py
"""

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

BASE_DIRECTORY = (
    REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
)
PREFIX = "Y2Ir2O7"

import wannierberri as wberri
from wannierberri.w90files import Wannier90data
from wannierberri.w90files.win import WIN

print(f"WannierBerri {wberri.__version__}")

import os
os.chdir(BASE_DIRECTORY)

w90data = Wannier90data.from_w90_files(seedname=PREFIX, files=["win", "chk", "eig"])
# berry=False: no .mmn file exists for this run (production route is deliberately
# .mmn-free, see memory/no-mmn-production.md), and System_w90's AA_R construction
# unconditionally needs w90data.mmn when berry=True (system_w90.py:217) regardless
# of transl_inv_JM. Falls back to Ham_R-only symmetrization, matching what the MBT
# precedent (wannierberri.ipynb) actually compared -- it also only wrote/compared
# _hr.dat, never a symmetrized _r.dat.
systems = wberri.System_w90(w90data=w90data, transl_inv_JM=True, berry=False)
print(f"num_wann={systems.num_wann}  num R={len(systems.rvec.iRvec)}")

# ---- unsymmetrized Ham_R, saved before symmetrize mutates in place ----
Ham_R_before = np.array(systems.Ham_R, copy=True)
iRvec_before = np.array(systems.rvec.iRvec, dtype=int, copy=True)
centres_before = np.array(systems.wannier_centers_cart, copy=True)

# ---- atomic positions (fractional) and species labels, .win order ----
win = WIN.from_w90_file(seedname=PREFIX)
positions = np.mod(np.asarray(win["atoms_frac"], dtype=float), 1.0)
atom_name = list(win["atoms_names"])
assert atom_name == ["Y"] * 4 + ["Ir1", "Ir2", "Ir3", "Ir4"] + ["O"] * 14, atom_name

proj = ["Y:p", "Ir1:d", "Ir2:d", "Ir3:d", "Ir4:d", "O:s;p"]

magmom = (
    [[0.0, 0.0, 0.0]] * 4
    + [
        [-0.172951, -0.172951, -0.172951],
        [-0.172951, 0.172951, 0.172951],
        [0.172951, -0.172951, 0.172951],
        [0.172951, 0.172951, -0.172951],
    ]
    + [[0.0, 0.0, 0.0]] * 14
)

print("\npositions (fractional), atom_name, magmom:")
for name, pos, m in zip(atom_name, positions, magmom):
    print(f"  {name:4s} {pos} {m}")

# ---- relax WannierBerri's own point-group tolerance, as the MBT notebook did ----
import wannierberri.symmetry.point_symmetry as ps_mod
import wannierberri.system.system as sys_mod

_orig_check = ps_mod.PointGroup.check_basis_symmetry


def _check_relaxed(self, basis, tol=1e-6, rel_tol=None):
    return _orig_check(self, basis, tol=max(tol, 5e-4), rel_tol=rel_tol)


sys_mod.PointGroup.check_basis_symmetry = _check_relaxed
ps_mod.PointGroup.check_basis_symmetry = _check_relaxed

print("\nCalling systems.symmetrize(...) -- this searches for the point group first.")
systems.symmetrize(
    proj=proj,
    atom_name=atom_name,
    positions=positions.tolist(),
    soc=True,
    magmom=magmom,
    reorder_back=True,
)
print("symmetrize() finished.")

# ---- compare ----
Ham_R_after = np.array(systems.Ham_R, copy=True)
iRvec_after = np.array(systems.rvec.iRvec, dtype=int, copy=True)
centres_after = np.array(systems.wannier_centers_cart, copy=True)

print(f"\nR-vector count: before={len(iRvec_before)} after={len(iRvec_after)}")
print(f"max |centre change| = {np.abs(centres_after - centres_before).max():.3e} A")

same_R = np.array_equal(iRvec_before, iRvec_after)
print(f"same R set: {same_R}")

if same_R:
    dH = Ham_R_after - Ham_R_before
    relH = np.linalg.norm(dH) / np.linalg.norm(Ham_R_before)
    print(f"\nHam_R: relative L2 change = {relH:.4%}   max abs change = {np.abs(dH).max():.4e} eV")
else:
    print("R sets differ -- cannot do a direct elementwise comparison without a reindex.")

np.savez_compressed(
    Path(__file__).parent / "yio_wberri_symmetrize_result.npz",
    Ham_R_before=Ham_R_before, Ham_R_after=Ham_R_after,
    iRvec_before=iRvec_before, iRvec_after=iRvec_after,
    centres_before=centres_before, centres_after=centres_after,
)
print("\nSaved yio_wberri_symmetrize_result.npz")
