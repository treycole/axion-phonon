"""Full-mesh, full-precision DFT eigenvalue symmetry test (YIO base, m-3m', 48 ops) from an NSCF
data-file-schema.xml whose k-list is the 8^3 list of base's .nnkp (same order as the NSCF input).
Bands 1-156 (occupied). Usage: python eig_symmetry_xml_full.py <data-file-schema.xml> ..."""
import sys
import numpy as np
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict
REPO = Path("/Users/treycole/Repos/axion-phonon")
BASE = REPO / "data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889"
G = np.load(REPO / "calculations/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz", allow_pickle=True)
text = (BASE / "Y2Ir2O7.nnkp").read_text().splitlines()
blk = lambda n: text[text.index(f"begin {n}") + 1:text.index(f"end {n}")]
recip = np.array([[float(x) for x in l.split()] for l in blk("recip_lattice")])
kred = np.array([[float(x) for x in l.split()[:3]] for l in blk("kpoints")[1:]])
key = lambda k: tuple(np.round(np.mod(k, 1) * 8).astype(int) % 8); idx = {key(k): i for i, k in enumerate(kred)}
HA = 27.211386245988
def cls(R):
    d = np.linalg.det(R); P = R * d; t = round(np.trace(P))
    name = {3: "E", -1: "C2", 0: "C3", 1: "C4"}[t]
    if t == -1:
        wv, v = np.linalg.eig(P); ax = np.real(v[:, np.argmin(abs(wv - 1))])
        name = "C2<100>" if np.count_nonzero(abs(ax) > 1e-6) == 1 else "C2<110>"
    return ("" if d > 0 else "I*") + name
for path in sys.argv[1:]:
    root = ET.parse(path).getroot()
    ks = root.findall(".//band_structure/ks_energies")
    eig = np.array([[float(x) for x in k.find("eigenvalues").text.split()] for k in ks]) * HA
    assert eig.shape[0] == len(kred), (eig.shape, len(kred))
    # check k order against the nnkp list (xml k in cartesian 2pi/alat; compare via reduced coordinates)
    alat = float(root.find(".//atomic_structure").attrib["alat"])
    kc = np.array([[float(x) for x in k.find("k_point").text.split()] for k in ks])
    b = np.array([[float(x) for x in root.find(".//reciprocal_lattice").find(t).text.split()] for t in ("b1", "b2", "b3")])
    kr = kc @ np.linalg.inv(b)
    assert all(key(a) == key(c) for a, c in zip(kr, kred)), "k order differs from .nnkp"
    by = defaultdict(list)
    for R, pr in zip(G["base_ops"], G["base_priming"]):
        s = -1.0 if str(pr) == "primed" else 1.0
        perm = np.array([idx[key(k)] for k in s * (kred @ recip) @ R.T @ np.linalg.inv(recip)])
        by[(cls(R), str(pr))].append(np.abs(eig[perm, :156] - eig[:, :156]).max())
    allv = np.concatenate([v for k, v in by.items() if k[0] != "E"])
    print(f"\n{path}\n  512 k, bands 1-156; non-identity ops: max {1e3 * allv.max():.4f} meV, median {1e3 * np.median(allv):.4f} meV")
    for k, v in sorted(by.items(), key=lambda kv: max(kv[1])):
        print(f"  {k[0]:10s} {k[1]:9s} n={len(v):2d}  {1e6 * min(v):9.2f} / {1e6 * max(v):9.2f} ueV")
