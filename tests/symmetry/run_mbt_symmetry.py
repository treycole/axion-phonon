"""MnBi2Te4 symmetry table (paper tab:symmetry-mbt) for the base and one displaced structure, with the models of the
production dtheta: the log is the table, symmetry_base.npz / symmetry_mode.npz the numbers, and the last lines say PASS
or FAIL against symmetry_lib.THRESHOLDS (exit 1 on FAIL). What is tested: symmetry_lib.py.

Run locally (python run_mbt_symmetry.py, with the repository on PYTHONPATH and BASE_RUN / MODE_RUN set to the run
folders, BASE_SCF / MODE_SCF to the pw.x outputs of the two structures, which give the groups) or on the cluster
with run_mbt_symmetry.sh, which copies each Wannier run's .chk, .eig, .win, .nnkp and each structure's .scf.out into
base/ and mode/ of the job's scratch folder and fills in the settings below. test_mbt_symmetry.py checks the saved results.
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import symmetry_lib  # noqa: E402

# ---- settings (run_mbt_symmetry.sh fills these in)
BASE_RUN = Path("base")
MODE_RUN = Path("mode")
BASE_SCF = Path("base/MnBi2Te4.scf.out")   # pw.x output of the undistorted structure: its group
MODE_SCF = Path("mode/MnBi2Te4.scf.out")   # and of the displaced one
ANTI_TRANSLATION = True  # the base also has {1'|1/2,1/2,1/2}, which QE does not list (symmetry_lib docstring)
PREFIX = "MnBi2Te4"
N_OCC = 82
FROZEN_TOP = 10.962      # eV, dis_froz_max of the window: test 1 over the bands below it at every k (and over all)
N_GENERIC = 7            # random k
N_TUBE = 5               # k within TUBE_RADIUS of the k_z axis
TUBE_RADIUS = 0.02       # 1/A
SEED = 20260825          # as the 10-22 run: the same k
BATCH = 4                # k points per curvature call (memory)
OUT = Path(".")
# --------------
t0 = time.time()
results = symmetry_lib.run_pair(BASE_RUN, MODE_RUN, BASE_SCF, MODE_SCF, ANTI_TRANSLATION, PREFIX, N_OCC, FROZEN_TOP,
                                N_GENERIC, N_TUBE, TUBE_RADIUS, SEED, BATCH, log=lambda s: print(s, flush=True))
symmetry_lib.save(results, OUT)
bad = {label: symmetry_lib.failures(r) for label, r in results.items()}
for label, b in bad.items():
    print(f"\n{label}: " + ("PASS" if not b else "FAIL\n  " + "\n  ".join(b)))
print(f"\nall done in {time.time() - t0:.0f} s")
sys.exit(1 if any(bad.values()) else 0)
