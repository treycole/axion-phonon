"""Step 5, displaced structures: the .amn A(k) = O(k) V_base(k) (bug #23), written as <PREFIX>.amn here.

O(k) = <psi_structure,mk| S |psi_base,nk> from step 4 (overlap/chunk_*.npz); V_base = U_dis U of the base Wannier
run with the same window (base.chk). The trial functions are then the base's Wannier functions, unshifted: in the
projection gauge this makes Lambda = 0 exactly at O(Q). Checks that the chunks cover every k once and that their k
are the .chk's, in order.
"""
from pathlib import Path

import numpy as np

from axion_phonon import paw
from axion_phonon.io import wannier90

# ---- settings (5_wannier.sh fills these in)
PREFIX = "MnBi2Te4"
# --------------
chunks = [np.load(f) for f in sorted(Path("overlap").glob("chunk_*.npz"))]
ik = np.concatenate([c["ik"] for c in chunks])
order = np.argsort(ik)
if not np.array_equal(ik[order], np.arange(1, len(ik) + 1)):
    raise SystemExit("the overlap chunks do not cover every k exactly once")
O = np.concatenate([c["O"] if "O" in c else c["O_exact"] for c in chunks])[order]   # O_exact: paw_overlap_pair.py (10-28)
xk = np.concatenate([c["xk"] for c in chunks])[order] / wannier90.BOHR       # 1/bohr -> 1/A
del chunks
chk = wannier90.read_chk("base.chk")
V = chk["v_matrix"]
if V.shape[0] != O.shape[0] or V.shape[1] != O.shape[2]:
    raise SystemExit(f"shapes: V {V.shape}, O {O.shape}")
dk = xk @ np.linalg.inv(chk["recip_lattice"]) - chk["kpt_red"]
dk = np.abs(dk - np.rint(dk)).max()
print(f"{len(ik)} k; k of the overlaps vs the base .chk: {dk:.1e}", flush=True)
if dk > 1e-6:
    raise SystemExit("the k lists differ")
A = paw.mode_amn(O, V)
s = np.linalg.svd(A, compute_uv=False)              # independent of the two NSCFs' band phases (unlike A - V)
print(f"singular values of A(k): smallest {s.min():.6f}, largest {s.max():.6f} (1 = the base Wannier functions lie "
      f"entirely in the structure's bands)", flush=True)
wannier90.write_amn(f"{PREFIX}.amn", A, "PAW overlap x base V (axion_phonon.paw.mode_amn), bug #23")
print(f"written {PREFIX}.amn", flush=True)
