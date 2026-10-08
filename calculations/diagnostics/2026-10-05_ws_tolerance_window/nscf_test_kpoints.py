"""Reduced k-set for the NSCF symmetry test: the full m-3m' orbits (k -> g k, or -g k for primed g)
of the N_ORBITS mesh points with the largest eigenvalue asymmetry in base 192344's .eig.
Writes nscf_test_kpoints.txt (K_POINTS crystal block, weight 1/nk each)."""
import numpy as np
from pathlib import Path
HERE = Path(__file__).parent; REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
N_ORBITS = 2
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
text = (BASE / "Y2Ir2O7.nnkp").read_text().splitlines()
blk = lambda n: text[text.index(f"begin {n}") + 1:text.index(f"end {n}")]
recip = np.array([[float(x) for x in l.split()] for l in blk("recip_lattice")])
kred = np.array([[float(x) for x in l.split()[:3]] for l in blk("kpoints")[1:]])
raw = np.loadtxt(BASE / "Y2Ir2O7.eig"); eig = raw[:, 2].reshape(int(raw[:, 1].max()), int(raw[:, 0].max()))
key = lambda k: tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8); idx = {key(k): i for i, k in enumerate(kred)}
images = []
for R, pr in zip(G["base_ops"], G["base_priming"]):
    s = -1.0 if str(pr) == "primed" else 1.0
    images.append(np.array([idx[key(k)] for k in s * (kred @ recip) @ R.T @ np.linalg.inv(recip)]))
images = np.array(images)                                   # (48 ops, 512 k)
asym = np.array([np.abs(eig[images[:, i], :156] - eig[i, :156]).max() for i in range(len(kred))])
chosen, seen = [], set()
for i in np.argsort(asym)[::-1]:
    if i in seen: continue
    orbit = sorted(set(images[:, i].tolist())); seen.update(orbit); chosen.append((i, orbit))
    if len(chosen) == N_ORBITS: break
for i in np.argsort(asym)[::-1]:                            # plus the most asymmetric generic orbit
    orbit = sorted(set(images[:, i].tolist()))
    if len(orbit) == 48:
        chosen.append((i, orbit)); break
ks = sorted(set(sum((o for _, o in chosen), [])))
for i, o in chosen:
    print(f"seed k #{i + 1} {kred[i]}  max asym {1e3 * asym[i]:.2f} meV  orbit size {len(o)}")
print(f"{len(ks)} k-points total")
with open(HERE / "nscf_test_kpoints.txt", "w") as f:
    f.write("K_POINTS crystal\n%d\n" % len(ks))
    for i in ks: f.write("%15.12f %15.12f %15.12f %15.12f\n" % (*kred[i], 1.0 / len(ks)))
np.save(HERE / "nscf_test_kpoint_indices.npy", np.array(ks))
