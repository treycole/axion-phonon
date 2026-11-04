"""The Chern-Simons 3-form check, one task: partial 3-form sums of a base/structure pair on one k slice of nk^3 meshes
(axion_phonon.chern_simons; log 2026-10-14). theta per structure needs no A_beta, so [theta(structure) - theta(base)]/Q
is an estimate of dtheta/dQ independent of the frozen gauge of c2. Only for theta = 0 materials (YIO).

The cs_theta.sh job copies each Wannier run's .chk, .eig, .win and .nnkp into base/ and mode/. Both models come from
pair_models (one tau), the trials G are the SCDM columns of the base (NK_TRIAL^3 mesh), the same G for both. Writes
chunks/chunk_<TASK-1>.npz: per nk the partial sums and the smallest singular value of <u_occ|g> at both ends, and G.
"""
import os
import time
from pathlib import Path

import numpy as np

from axion_phonon import chern_simons, curvature, wannier
from axion_phonon.io import wannier90

# ---- settings (cs_theta.sh fills these in)
PREFIX = "Y2Ir2O7"
N_OCC = 156
NKS = [10, 12]           # nk^3 meshes for the zone sum
NK_TRIAL = 6             # mesh of the SCDM trial construction
N_TASKS = 2              # k slices of each mesh
BATCH = 32               # k points per call (memory)
# --------------
task = int(os.environ["TASK"]) - 1
out = Path("chunks") / f"chunk_{task:04d}.npz"
t0 = time.time()
base, mode = wannier90.read_run("base", PREFIX), wannier90.read_run("mode", PREFIX)
models = dict(zip(("base", "mode"), wannier.pair_models(base, mode)))
del base, mode
pos = {e: curvature.position_terms(m, dtype=np.complex128) for e, m in models.items()}
G = chern_simons.scdm_trials(models["base"], NK_TRIAL, N_OCC)
res = {"G": G, "nks": np.array(NKS), "n_tasks": N_TASKS, "task": task}
for e, m in models.items():
    res[f"smin_trial_{e}"] = chern_simons.smallest_sigma(m, G, chern_simons.mesh(NK_TRIAL), N_OCC)
print(f"task {task + 1}/{N_TASKS}: models and trials in {time.time() - t0:.0f} s; min sigma on the {NK_TRIAL}^3 mesh "
      f"base {res['smin_trial_base']:.3f}, structure {res['smin_trial_mode']:.3f}", flush=True)
for nk in NKS:
    k = np.array_split(chern_simons.mesh(nk), N_TASKS)[task]
    for e, m in models.items():
        res[f"sum_{e}_nk{nk}"], res[f"smin_{e}_nk{nk}"] = chern_simons.theta_sum(m, pos[e], G, k, N_OCC, BATCH)
    print(f"  nk {nk}: {len(k)} of {nk ** 3} k, partial sums base {res[f'sum_base_nk{nk}']:.6e}, "
          f"structure {res[f'sum_mode_nk{nk}']:.6e} [{time.time() - t0:.0f} s]", flush=True)
out.parent.mkdir(exist_ok=True)
np.savez(str(out) + ".part.npz", **res)
os.replace(str(out) + ".part.npz", out)
print(f"task {task + 1}: done in {time.time() - t0:.0f} s", flush=True)
