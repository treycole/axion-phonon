"""The Chern-Simons 3-form check, the sum: theta of both structures and [theta(structure) - theta(base)]/Q per mesh,
from the cs_theta.py chunks. Writes theta_cs.npz (nks, theta_base, theta_mode, dtheta, smin)."""
import numpy as np

# ---- settings (cs_theta_sum.sh fills these in)
CHUNKS = []              # the cs_theta.py chunk files, one per task
Q = 0.01                 # displacement amplitude of the structure, Angstrom
# --------------
z = [np.load(p) for p in CHUNKS]
if sorted(int(c["task"]) for c in z) != list(range(int(z[0]["n_tasks"]))):
    raise SystemExit("not one chunk per task")
G = z[0]["G"]
if max(np.abs(c["G"] - G).max() for c in z) > 1e-10:
    raise SystemExit("the tasks built different trials G")
nks = z[0]["nks"]
rows = []
for nk in nks:
    th = {e: -sum(c[f"sum_{e}_nk{nk}"] for c in z) / nk ** 3 / (4 * np.pi) for e in ("base", "mode")}
    smin = min(float(c[f"smin_{e}_nk{nk}"]) for c in z for e in ("base", "mode"))
    rows.append((nk, th["base"], th["mode"], (th["mode"].real - th["base"].real) / Q, smin))
    print(f"nk {nk:3d}: theta base {th['base'].real:+.8f} (Im {th['base'].imag:.1e}), structure {th['mode'].real:+.8f} "
          f"(Im {th['mode'].imag:.1e}); min sigma {smin:.3f}; [theta(structure) - theta(base)] / Q = "
          f"{rows[-1][3]:+.6f} rad/A")
smin_trial = min(float(c[f"smin_trial_{e}"]) for c in z for e in ("base", "mode"))
print(f"min sigma of the trials on their mesh {smin_trial:.3f}")
print(f"dtheta/dQ (3-form) = {rows[-1][3]:.6f} rad/A (nk {nks[-1]}; Q = {Q} A)")
np.savez("theta_cs.npz", nks=nks, theta_base=np.array([r[1] for r in rows]), theta_mode=np.array([r[2] for r in rows]),
         dtheta=np.array([r[3] for r in rows]), smin=np.array([r[4] for r in rows]), Q=Q)
