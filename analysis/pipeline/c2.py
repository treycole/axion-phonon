"""Step 6, one task: the c2 density at both ends of a base/structure pair, on the k points k[TASK-1::N_TASKS] of
klist.npz (frozen gauge A_beta = 0; production precision, complex64 in the position terms).

The 6_c2.sh job copies each Wannier run's .chk, .eig, .win and .nnkp into base/ and mode/. Writes
chunks/chunk_<TASK-1>.npz: idx (into klist.npz), c2 (n, 2) complex (base end, structure end), gap (n, 2) eV.
"""
import os
import time
from pathlib import Path

import numpy as np

from axion_phonon import axion, curvature, wannier
from axion_phonon.io import wannier90

# ---- settings (6_c2.sh fills these in)
PREFIX = "MnBi2Te4"
N_OCC = 82
Q = 0.01                 # displacement amplitude of the structure, Angstrom
N_TASKS = 10
BATCH = 50               # k points per call (memory)
# --------------
task = int(os.environ["TASK"]) - 1
out = Path("chunks") / f"chunk_{task:04d}.npz"
k_all = np.load("klist.npz")["k"]
idx = np.arange(task, len(k_all), N_TASKS)
t0 = time.time()
base, mode = wannier90.read_run("base", PREFIX), wannier90.read_run("mode", PREFIX)
mb, mm = wannier.pair_models(base, mode)
del base, mode
pb, pm = curvature.position_terms(mb), curvature.position_terms(mm)
print(f"task {task + 1}/{N_TASKS}: {len(idx)} k, models in {time.time() - t0:.0f} s", flush=True)
c2 = np.zeros((len(idx), 2), complex)
gap = np.zeros((len(idx), 2))
for s in range(0, len(idx), BATCH):
    k = k_all[idx[s:s + BATCH]]
    fb, fm = curvature.fields(mb, pb, k), curvature.fields(mm, pm, k)
    beta = axion.beta_terms(fb, fm, Q)
    for end, f in enumerate((fb, fm)):
        c = curvature.connection(f, N_OCC, beta)
        c2[s:s + BATCH, end] = axion.c2_density(curvature.omega(c))
        gap[s:s + BATCH, end] = c.E[:, N_OCC] - c.E[:, N_OCC - 1]
    print(f"  {s + len(k)}/{len(idx)} [{time.time() - t0:.0f} s]", flush=True)
out.parent.mkdir(exist_ok=True)
np.savez(str(out) + ".part.npz", idx=idx, c2=c2, gap=gap)
os.replace(str(out) + ".part.npz", out)
print(f"task {task + 1}: done in {time.time() - t0:.0f} s", flush=True)
