#!/usr/bin/env python3
"""Axion response dθ/dQ versus k-mesh for one base/mode pair: the whole calculation on one page.

Stage 1 (see ``wannier/``) already wrote ``seed_hr.dat``, ``seed_r.dat`` and ``seed_wsvec.dat`` into the two
run folders. PythTB reads them; ``modules/curvature.py`` and ``modules/axion.py`` do everything else.

Edit the settings below, then run ``python calculations/run_axion.py`` in the ``axion`` conda env (PythTB on
its ``external`` branch). The result goes into ``<MODE>/nk_sweep_<source>_<own|shared>_ws/`` (``nk_*.npz``
with dθ, the c2 density and the gap, plus ``summary.npz``/``summary.csv``), next to the data it came from.
Peak memory for the ``.mmn`` route on Y2Ir2O7 is roughly 10-12 GB.
"""

import sys
from pathlib import Path

import numpy as np
from pythtb import W90

REPO = Path("/Users/treycole/Repos/axion-phonon")
DATA = REPO / "data"
sys.path.insert(0, str(REPO))

from modules.axion import dtheta, save_dtheta_results
from modules.curvature import position_terms

# --- WS tie tolerance 1e-5 (the pre-2026-09-30 default): reset the bound defaults ---
import inspect
import wannier.position_matrix as _pm
WS_TOLERANCE = 1e-5
for _name, _f in inspect.getmembers(_pm, inspect.isfunction):
    for _arg in ("tolerance", "R_distance_tolerance"):
        if _f.__kwdefaults__ and _f.__kwdefaults__.get(_arg) == _pm._DEFAULT_WS_DISTANCE_TOLERANCE:
            _f.__kwdefaults__[_arg] = WS_TOLERANCE
        _p = [p for p in inspect.signature(_f).parameters.values() if p.default is not inspect.Parameter.empty
              and p.kind is not inspect.Parameter.KEYWORD_ONLY]
        if _f.__defaults__ and _arg in [p.name for p in _p]:
            _d = list(_f.__defaults__); _i = [p.name for p in _p].index(_arg)
            if _d[_i] == _pm._DEFAULT_WS_DISTANCE_TOLERANCE:
                _d[_i] = WS_TOLERANCE; _f.__defaults__ = tuple(_d)
print("WS tie tolerance defaults set to", WS_TOLERANCE, flush=True)

# ============================================================================== settings
PREFIX = "Y2Ir2O7"
# Wannier90 runs of the undisplaced and the displaced structure (trial_02: Y:p Ir:d O:s;p, 176 WF).
# QE 7.6 + new_ns_nc fix (2026-09); the old QE 7.5 pair was base 178382bc889/.../proj_gauge/181538bc889 and
# mode 178365bc889/.../no_displacement/181565bc889 (their text files are gzipped now).
BASE = Path("/Users/treycole/Repos/axion-phonon/data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889")
MODE = Path("/Users/treycole/Repos/axion-phonon/data/Y2Ir2O7/phonon/GM2-/mode2/Q_0.01A/output/192346/trial_02_Y_p_Ir_d_O_sp/no-displacement/193196bc889")
# .mmn A(R) caches of the two runs (used by SOURCE = "mmn"; the .mmn files themselves are deleted)
CACHE_BASE = DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_base_181538bc889.npz"
CACHE_MODE = DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_mode_181565bc889.npz"
N_OCC = 156         # occupied Wannier bands
DBETA = 0.01        # Angstrom, the phonon amplitude between the two structures
NKS = [4, 6, 8, 10, 12]

# Position matrix A(R):
#   "rdat_ndegen_applied"  Jae-Mo's write_ndegen_applied files in <run>/jaemo/
#   "mmn"                  the stock _hr.dat plus the cached .mmn A(R) (CACHE_BASE, CACHE_MODE)
#   "transl_inv_full_from_chk"  Jae-Mo's formula, rebuilt from the .chk m_matrix (production since 2026-10-07)
#   "right_centre_from_chk"     the .mmn formula, from the same .chk (used 2026-09-30 to 2026-10-06)
SOURCE = "transl_inv_full_from_chk"
# Fold both structures on the base centres (assumption A3', notes/progress/2026-09-23_ws_ties_in_beta_derivative.md).
# Needs one of the last three sources; Jae-Mo's files cannot be refolded.
SHARED_WS = True
OUT = MODE / "nk_sweep_tif_shared_ws_tol1e-5"
# ======================================================================================

OUT_NAMES = {
    ("rdat_ndegen_applied", False): "nk_sweep_new",
    ("mmn", False): "nk_sweep_mmn_own_ws",
    ("mmn", True): "nk_sweep_mmn_shared_ws",
    ("transl_inv_full_from_chk", False): "nk_sweep_tif_own_ws",
    ("transl_inv_full_from_chk", True): "nk_sweep_tif_shared_ws",
    ("right_centre_from_chk", False): "nk_sweep_chkmmn_own_ws",
    ("right_centre_from_chk", True): "nk_sweep_chkmmn_shared_ws",
}
OPTIONS = dict(zero_energy=0, min_hopping_norm=1e-5, max_distance=125, ignorable_imaginary_part=None)
if (SOURCE, SHARED_WS) not in OUT_NAMES:
    raise ValueError(f"SOURCE={SOURCE!r}, SHARED_WS={SHARED_WS}: see the options above")
OUT = Path(OUT) if OUT else MODE / OUT_NAMES[(SOURCE, SHARED_WS)]
print(f"{PREFIX}: n_occ={N_OCC}, dbeta={DBETA} A, source={SOURCE}, shared_ws={SHARED_WS}\n"
      f"  base: {BASE}\n  mode: {MODE}\n  out:  {OUT}")

if SOURCE == "rdat_ndegen_applied":
    # PythTB reads H(R) and A(R).  ONE nominal tau for both structures (assumption A3);
    # the physical motion of the Wannier centres enters through A(R) instead.
    base, mode = W90(str(BASE / "jaemo"), PREFIX), W90(str(MODE / "jaemo"), PREFIX)
    tau = base.lattice.orb_vecs.copy()
    model_base = base.model(orb_vecs=tau, **OPTIONS)
    model_mode = mode.model(orb_vecs=tau, **OPTIONS)
    del base, mode  # the models carry their own copies of H(R) and A(R)
else:
    from wannier.wannier_io import load_wannier_pair

    pair = load_wannier_pair(
        BASE,
        MODE,
        prefix=PREFIX,
        cache_directory=CACHE_BASE.parent,
        A_R_cache_files={"base": CACHE_BASE, "mode": CACHE_MODE},
        position_source=SOURCE,
        shared_ws_centres=SHARED_WS,  # tau is the base centres either way (A3)
        model_options=OPTIONS,
    )
    model_base, model_mode = pair.base_model, pair.mode_model
    del pair

# X(R) once per structure, reused for every mesh.
positions = (position_terms(model_base), position_terms(model_mode))

OUT.mkdir(parents=True, exist_ok=True)
metadata = {
    "dbeta": np.float64(DBETA),
    "n_occupied": np.int32(N_OCC),
    "base_directory": np.array([str(BASE)]),
    "mode_directory": np.array([str(MODE)]),
    "position_source": np.array([SOURCE]),
    "shared_ws_centres": np.int8(SHARED_WS),
    "frozen_gauge": np.int8(1),  # Lambda_R=None below: A_beta = 0, an explicit choice
}
results = []
for nk in NKS:
    result = dtheta(model_base, model_mode, nk, DBETA, N_OCC, positions=positions, progress=print)
    results.append(result)
    save_dtheta_results(OUT, results, metadata)
    print(f"nk={nk}:  dtheta = {result['dtheta'].real}   (base end, mode end)")
