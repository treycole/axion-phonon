"""Step 7: dtheta/dQ of a pair from the step-6 chunks: c2 on the whole k list, integrated (axion_phonon.tube).

Prints the convergence table (rows: in-plane refinement level m; columns: planes along the small-gap line) at both
ends, and writes c2.npz (the k-resolved integrand: k, level, plane, c2, gap) and dtheta.npz (levels, n3, dtheta).
"""
import numpy as np

from axion_phonon import tube

# ---- settings (7_dtheta.sh fills these in)
CHUNKS = []              # the step-6 chunk files, one per task
# --------------
kl = dict(np.load("klist.npz"))
c2 = np.full((len(kl["k"]), 2), np.nan, complex)
gap = np.full((len(kl["k"]), 2), np.nan)
for path in CHUNKS:
    z = np.load(path)
    if not np.isnan(c2[z["idx"], 0].real).all():
        raise SystemExit(f"{path}: k points already filled by another chunk")
    c2[z["idx"]], gap[z["idx"]] = z["c2"], z["gap"]
missing = np.isnan(c2[:, 0].real).sum()
if missing:
    raise SystemExit(f"{missing} of {len(c2)} k points have no c2")
levels, n3s, table = tube.integrate(c2, kl)
print(f"{len(c2)} k points from {len(CHUNKS)} chunks; smallest direct gap {gap.min(0)} eV (base end, structure end); "
      f"max |Im c2| / max |c2| {np.abs(c2.imag).max() / np.abs(c2).max():.1e}")
for end, label in enumerate(("base end", "structure end")):
    print(f"dtheta/dQ (rad/A), {label} (rows: in-plane level m; columns: planes along k_z)")
    print("      " + "".join(f"{n:>12d}" for n in n3s))
    for a, m in enumerate(levels):
        print(f"m={m:<3d} " + "".join(f"{v:12.5f}" for v in table[a, :, end]))
b, s = table[-1, 0]
print(f"dtheta/dQ = {b:.5f} (base end), {s:.5f} (structure end), mean {(b + s) / 2:.5f} rad/A "
      f"(m = {levels[-1]}, {n3s[0]} planes)")
np.savez("c2.npz", k=kl["k"], level=kl["level"], plane=kl["plane"], c2=c2, gap=gap)
np.savez("dtheta.npz", levels=np.array(levels), n3=np.array(n3s), dtheta=table)
