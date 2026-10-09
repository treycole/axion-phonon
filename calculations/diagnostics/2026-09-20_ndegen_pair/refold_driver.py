"""DIAGNOSTIC (not production): run ``run_axion.py`` with both structures
refolded onto ONE Wigner-Seitz selection rule, the base's centres.

**H(R) only by default.**  Refolding A(R) from class sums is INVALID: under
transl_inv_full the images of an A(R) class are not copies (exp(-i b.R/2) is
applied per b at the image vector), see notes/log/2026-09-20_jaemo_ndegen_yio_mode1.md
and image_structure_check.py.  REFOLD_A=1 reproduces the retracted experiment;
its dtheta means nothing.  The H-only path is untested.

Jae-Mo's files choose the Wigner-Seitz images per run from that run's own
centres.  The aliasing-class sums of H(R) and A(R) are the mesh-determined,
centre-independent data, so:  sum each class, then redistribute it with the
in-house min|R+tau_j-tau_i| rule evaluated at the BASE centres, for base and
mode alike.  Nothing else changes: the unmodified driver runs afterwards.
"""
import os, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
import numpy as np
import pythtb
from wannier.position_matrix import _map_R_by_orbital_positions

GRID = np.array((8, 8, 8))
REFOLD_A = os.environ.get("REFOLD_A", "0") == "1"
DRIVER = str(REPO / "calculations/run_axion.py")
_original_W90 = pythtb.W90  # run_axion.py does `from pythtb import W90` after the patch below
_state = {"calls": 0, "centres": None}


def _class_sums(blocks):
    sums = {}
    for R, block in blocks.items():
        key = tuple(int(x) for x in np.mod(np.asarray(R), GRID))
        sums[key] = np.array(block, dtype=complex) if key not in sums else sums[key] + block
    return sums


def _representative(key):
    return tuple(int(((k + g // 2) % g) - g // 2) for k, g in zip(key, GRID))


def _refold(w, centres):
    H_sums = _class_sums({R: np.asarray(b["h"]) / float(b["deg"]) for R, b in w.ham_r.items()})
    A_sums = _class_sums({R: np.asarray(X) for R, X in w.pos_r.items()})  # deg == 1 for ndegen files
    H = _map_R_by_orbital_positions({_representative(k): v for k, v in H_sums.items()}, w.lat, centres, GRID)
    n = next(iter(H.values())).shape[-1]
    if not REFOLD_A:
        # H(R) images are copies, so the class-sum refold is exact for H; leave A(R) as written.
        zero_A = np.zeros((3, n, n), complex)
        for R in H:
            w.pos_r.setdefault(R, zero_A)
        w.ham_r = {R: {"h": v, "deg": 1} for R, v in H.items()}
        for R in list(w.pos_r):
            w.ham_r.setdefault(R, {"h": np.zeros((n, n), complex), "deg": 1})
        print(f"    refolded H(R) only onto base centres: {len(w.ham_r)} R vectors", flush=True)
        return
    A = _map_R_by_orbital_positions({_representative(k): v for k, v in A_sums.items()}, w.lat, centres, GRID)
    keys = set(H) | set(A)
    # Hermitian part of A, exactly as the production remap does.
    zero_A = np.zeros((3, n, n), complex)
    A_h = {}
    for R in keys:
        partner = A.get(tuple(-x for x in R), zero_A)
        A_h[R] = 0.5 * (A.get(R, zero_A) + np.conj(np.transpose(partner, (0, 2, 1))))
    zero_H = np.zeros((n, n), complex)
    w.ham_r = {R: {"h": H.get(R, zero_H), "deg": 1} for R in keys}
    w.pos_r = A_h
    print(f"    refolded onto base centres: {len(keys)} R vectors", flush=True)


def W90_refolded(directory, prefix, *args, **kwargs):
    w = _original_W90(directory, prefix, *args, **kwargs)
    _state["calls"] += 1
    if _state["calls"] == 1:  # run_axion.py builds the base first
        _state["centres"] = np.asarray(w.lattice.orb_vecs, float) @ np.asarray(w.lat, float)
    moved = np.abs(
        np.asarray(w.lattice.orb_vecs, float) @ np.asarray(w.lat, float) - _state["centres"]
    ).max()
    print(f"  refolding {Path(str(directory)).parts[-3]} (centres differ from base by {moved:.2e} A)", flush=True)
    _refold(w, _state["centres"])
    return w


pythtb.W90 = W90_refolded
# The refold applies to Jae-Mo's files, so run the driver with SOURCE = "rdat_ndegen_applied" and its own
# refolding off; its BASE/MODE settings (the production YIO mode-1 pair) are used as they are.
code = Path(DRIVER).read_text()
code, n_source = re.subn(r'^SOURCE = .*$', 'SOURCE = "rdat_ndegen_applied"', code, flags=re.M)
code, n_shared = re.subn(r'^SHARED_WS = .*$', 'SHARED_WS = False', code, flags=re.M)
assert n_source == n_shared == 1, "run_axion.py settings block changed"
exec(compile(code, DRIVER, "exec"), {"__name__": "__main__", "__file__": DRIVER})
