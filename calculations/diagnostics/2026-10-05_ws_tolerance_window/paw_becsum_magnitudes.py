"""PAW becsum (paw.txt: becsum(ijh, na, ispin), list-directed) of Ir1-Ir4. C2x/C2y/C2z map Ir1 onto
Ir2/Ir3/Ir4 and act on real-harmonic projector pairs and on (n, mx, my, mz) by signs only, so
|becsum(:, Ir_k, comp)| must be identical for all four Ir. Usage: python ... <paw.txt> ..."""
import sys
import numpy as np
NAT, NSPIN = 22, 4
IR = [4, 5, 6, 7]
for path in sys.argv[1:]:
    v = np.array([float(x.replace("D", "E")) for x in open(path).read().split()])
    nij = v.size // (NAT * NSPIN)
    assert nij * NAT * NSPIN == v.size, v.size
    nhm = int(round((np.sqrt(8 * nij + 1) - 1) / 2))
    b = v.reshape(NSPIN, NAT, nij)          # Fortran (ijh, na, ispin) -> [ispin, na, ijh]
    print(f"\n{path}: {v.size} values -> ijh {nij} (nhm {nhm}), nat {NAT}, nspin {NSPIN}")
    for c, name in enumerate(["n", "mx", "my", "mz"]):
        ref = np.abs(b[c, IR[0]])
        diffs = [np.abs(np.abs(b[c, a]) - ref).max() for a in IR]
        print(f"  {name:2s}: max|b| {ref.max():.3e};  max | |b(Ir_k)| - |b(Ir1)| | for Ir1..Ir4: " + " ".join(f"{d:.2e}" for d in diffs))
    worst = max(np.abs(np.abs(b[c, a]) - np.abs(b[c, IR[0]])).max() for c in range(4) for a in IR)
    for a in IR[1:]:
        dd = np.abs(np.abs(b[:, a]) - np.abs(b[:, IR[0]]))
        c, j = np.unravel_index(dd.argmax(), dd.shape)
        print(f"  Ir{a - 3}: worst element comp {['n','mx','my','mz'][c]} ijh {j + 1}: Ir1 {b[c, IR[0], j]:+.6e}  Ir{a - 3} {b[c, a, j]:+.6e}")
