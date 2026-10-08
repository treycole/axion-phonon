"""Dynamic range of ws_pair_consistency on the production YIO mode-1 pair.

A. both structures mapped from their OWN centres   (what production + Jae-Mo do)
B. both structures mapped from the BASE centres    (one shared selection rule)
C. base mapped, mode left on Wannier90's min|R| set (mixed rules: the CLAUDE.md #1 failure)
D. both left on Wannier90's min|R| set             (centre-independent rule)
"""
import sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
DATA = Path(__file__).resolve().parents[3] / "data"
import numpy as np
from pythtb import W90
from wannier.position_matrix import (
    apply_orbital_dependent_R_mapping,
    ws_pair_consistency,
)

DIRS = {"base": DATA / "Y2Ir2O7/base/soc/u_3.0/output/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181538bc889", "mode": DATA / "Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/178365bc889/trial_02_Y_p_Ir_d_O_sp/no_displacement/181565bc889"}
GRID = (8, 8, 8)
t0 = time.time()


def effective(w):
    return {R: np.asarray(b["h"]) / float(b["deg"]) for R, b in w.ham_r.items()}


def centres(w):
    return np.asarray(w.lattice.orb_vecs, float) @ np.asarray(w.lat, float)


raw, own, shared = {}, {}, {}
base_centres = None
for label in ("base", "mode"):
    w = W90(str(DIRS[label]), "Y2Ir2O7")
    raw[label] = effective(w)
    c = centres(w)
    if label == "base":
        base_centres = c
    print(f"{label}: centres moved vs base = {np.abs(c - base_centres).max():.3e} A", flush=True)
    for name, store, cen in (("own", own, None), ("shared", shared, base_centres)):
        # The mapping mutates w.ham_r, so restore the raw (deg already divided out).
        w.ham_r = {R: {"h": raw[label][R], "deg": 1} for R in raw[label]}
        apply_orbital_dependent_R_mapping(w, GRID, centers_cartesian=cen, verbose=False)
        store[label] = {R: b["h"] for R, b in w.ham_r.items()}
    print(f"{label}: done at {time.time()-t0:.0f}s", flush=True)
    del w


def show(name, a, b):
    r = ws_pair_consistency(a, b, GRID)
    print(f"{name:52s} mismatched={r['entries_mismatched']:>7d}/{r['entries_compared']:<8d} "
          f"contamination={r['delta_H_contamination']:.4f}  R lists {r['n_R_base']}/{r['n_R_mode']}", flush=True)


show("A own centres (production / Jae-Mo)", own["base"], own["mode"])
show("B base centres for both (shared rule)", shared["base"], shared["mode"])
show("C base mapped, mode raw min|R| (mixed rules)", own["base"], raw["mode"])
show("D both raw min|R| (centre-independent)", raw["base"], raw["mode"])
