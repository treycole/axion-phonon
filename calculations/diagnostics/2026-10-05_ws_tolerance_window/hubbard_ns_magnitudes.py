"""Compare the final printed |n_ij| Hubbard matrices of Ir1-Ir4 (SCF output). C2<100> ops map Ir1 onto
Ir2-4 with a D(g) that is diagonal (phases only) in QE's real-harmonic x spin basis, so a covariant ns
gives four identical |n| matrices. Also prints the Hubbard moments. Usage: python ... <scf.out> ..."""
import re, sys
import numpy as np
for path in sys.argv[1:]:
    t = open(path).read()
    last = t.split("HUBBARD OCCUPATIONS")[-1]
    mats, moms = {}, {}
    for a in (5, 6, 7, 8):
        blk = last.split(f"ATOM    {a} ")[1].split("ATOM ")[0]
        occ = blk.split("occupations, | n_(i1, i2)^(sigma1, sigma2) |:")[1].split("Atomic magnetic")[0]
        mats[a] = np.array([float(x) for x in re.findall(r"-?\d+\.\d+", occ)]).reshape(10, 10)
        moms[a] = [float(x) for x in re.findall(r"-?\d+\.\d+", blk.split("mx, my, mz =")[1])[:3]]
    print(path)
    swap = np.r_[5:10, 0:5]                 # C2x, C2y rotate spin by -i sigma_x/y: up <-> down blocks
    expected = {5: mats[5], 8: mats[5], 6: mats[5][np.ix_(swap, swap)], 7: mats[5][np.ix_(swap, swap)]}
    for a, op in ((5, "E"), (6, "C2x"), (7, "C2y"), (8, "C2z")):
        print(f"  Ir{a - 4} (= {op} Ir1): Hubbard moment {moms[a]}  max | |n| - expected | = {np.abs(mats[a] - expected[a]).max():.3f}")
