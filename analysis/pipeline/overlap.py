"""Step 4, one process: PAW overlaps O(k) = <psi_structure,mk| S |psi_base,nk> on the k points k[TASK-1::N_TASKS].

The 4_overlap.sh job copies both NSCF saves into base/ and mode/ and runs N_TASKS of these side by side. Writes
chunks/chunk_<TASK-1>.npz: ik (1-based, the wfc file numbers), xk (Cartesian, 1/bohr, as the wfc files store it) and
O (n_k, bands of the structure, bands of the base).
Also prints the S-orthonormality of both structures' own states at the first k (QE's davidson: ~1e-14).
"""
import os
import time
from pathlib import Path

import numpy as np

from axion_phonon import paw
from axion_phonon.io import qe

# ---- settings (4_overlap.sh fills these in)
PREFIX = "MnBi2Te4"
N_TASKS = 8
# --------------
BASE, MODE = Path("base") / f"{PREFIX}.save", Path("mode") / f"{PREFIX}.save"
task = int(os.environ["TASK"]) - 1
out = Path("chunks") / f"chunk_{task:04d}.npz"
nk = len(list(BASE.glob("wfc*.dat")))
assert nk == len(list(MODE.glob("wfc*.dat"))), "the two saves have different numbers of k points"
iks = np.arange(task, nk, N_TASKS) + 1
t0 = time.time()
systems = {}
for name, save in (("base", BASE), ("mode", MODE)):
    structure = qe.read_structure(save)
    systems[name] = paw.PAWSystem(structure, qe.read_pseudos(structure, save))
sb, sm = systems["base"], systems["mode"]
print(f"task {task + 1}/{N_TASKS}: {len(iks)} k; setup {time.time() - t0:.0f} s", flush=True)
O, xk = [], []
for i, ik in enumerate(iks):
    wb, wm = qe.read_wfc(BASE / f"wfc{ik}.dat"), qe.read_wfc(MODE / f"wfc{ik}.dat")
    if not np.allclose(wb["xk"], wm["xk"], atol=1e-8):
        raise SystemExit(f"k {ik}: the two saves have different k points")
    O.append(paw.overlap(sm, wm, sb, wb))
    xk.append(wb["xk"])
    if i == 0:
        ortho = [np.abs(paw.overlap(s, w, s, w) - np.eye(w["nbnd"])).max() for s, w in ((sb, wb), (sm, wm))]
        print(f"  S-orthonormality at k {ik}: base {ortho[0]:.1e}, structure {ortho[1]:.1e}", flush=True)
        if max(ortho) > 1e-8:
            raise SystemExit("the states are not S-orthonormal: wrong pseudopotentials or a reading error")
    print(f"  k {ik} [{time.time() - t0:.0f} s]", flush=True)
out.parent.mkdir(exist_ok=True)
np.savez(str(out) + ".part.npz", ik=iks, xk=np.array(xk), O=np.array(O))
os.replace(str(out) + ".part.npz", out)
print(f"task {task + 1}: done in {time.time() - t0:.0f} s", flush=True)
