"""Analyse the NSCF symmetry-test runs (cluster base/u_3.0/test_nscf_sym_2026-10-05/{A,B,C}).
Usage: python nscf_test_analyse.py <nscf.out> [<nscf.out> ...]
For each run: max |eps_n(g k) - eps_n(k)|, bands 1-156, per op class, over the 60-point k-set
(full orbits, so every image is in the set); Ir1-axis L point vs the other three L points; and, for
192344-save runs, the shift from the original 192344 .eig at the same k."""
import re, sys
import numpy as np
from pathlib import Path
from collections import defaultdict
HERE = Path(__file__).parent; REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
text = (BASE / "Y2Ir2O7.nnkp").read_text().splitlines()
blk = lambda n: text[text.index(f"begin {n}") + 1:text.index(f"end {n}")]
recip = np.array([[float(x) for x in l.split()] for l in blk("recip_lattice")])
kall = np.array([[float(x) for x in l.split()[:3]] for l in blk("kpoints")[1:]])
sel = np.load(HERE / "nscf_test_kpoint_indices.npy")
ksub = kall[sel]
raw = np.loadtxt(BASE / "Y2Ir2O7.eig"); eig_ref = raw[:, 2].reshape(int(raw[:, 1].max()), int(raw[:, 0].max()))[sel]
key = lambda k: tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8); idx = {key(k): i for i, k in enumerate(ksub)}
def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name
perms = []
for R, pr in zip(G["base_ops"], G["base_priming"]):
    s = -1.0 if str(pr) == "primed" else 1.0
    perms.append(np.array([idx[key(k)] for k in s * (ksub @ recip) @ R.T @ np.linalg.inv(recip)]))

def read_out(path):
    t = Path(path).read_text()
    if "JOB DONE" not in t: print(f"  WARNING: {path} has no JOB DONE")
    blocks = re.split(r"\n\s+k =", t.split("End of band structure calculation")[-1])[1:]
    ev = []
    for b in blocks:
        body = b.split("bands (ev):")[1].split("occupation numbers")[0]
        body = body.split("highest")[0].split("the Fermi")[0]
        ev.append([float(x) for x in re.findall(r"-?\d+\.\d+", body)])
    n = min(len(e) for e in ev)
    return np.array([e[:n] for e in ev])

L_axis = [i for i, k in enumerate(ksub) if np.allclose(np.mod(k, 1), 0.5)]
for path in sys.argv[1:]:
    eig = read_out(path)
    print(f"\n{path}: {eig.shape[0]} k, {eig.shape[1]} bands")
    assert eig.shape[0] == len(ksub)
    by = defaultdict(list)
    for g, (R, pr) in enumerate(zip(G["base_ops"], G["base_priming"])):
        by[(cls(R), str(pr))].append(np.abs(eig[perms[g], :156] - eig[:, :156]).max())
    for k, v in sorted(by.items(), key=lambda kv: max(kv[1])):
        print(f"  {k[0]:10s} {k[1]:9s} n={len(v):2d}  max |d eps| {1e3 * min(v):8.4f} / {1e3 * max(v):8.4f} meV")
    Ls = [i for i, k in enumerate(ksub) if np.allclose(np.abs(np.mod(k, 1) - 0.5).sum() % 1, 0) and
          sorted(np.round(np.mod(k, 1), 3)) in ([0.0, 0.0, 0.5], [0.5, 0.5, 0.5], [0.0, 0.5, 0.5])]
    if L_axis:
        ref = eig[L_axis[0], :156]
        print("  L points (Ir1 axis first):", ", ".join(f"{np.round(ksub[i], 3)} {1e3 * np.abs(eig[i, :156] - ref).max():.3f} meV" for i in Ls))
    print(f"  vs original 192344 .eig (bands 1-156): max {1e3 * np.abs(eig[:, :156] - eig_ref[:, :156]).max():.3f} meV")
