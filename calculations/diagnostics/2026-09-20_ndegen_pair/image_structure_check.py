"""Are the tied Wigner-Seitz images of H(R) and A(R) equal copies of each other?

Source claim (hamiltonian_get_rmn, transl_inv_full + write_ndegen_applied):
  A(R_image) = sum_b [ndegen-weighted op_b(class)] * exp(-i b.R_image/2)
with b = +-G_i/8, so between images differing by L = 8n the per-b factor changes
by (-1)^{n_i}.  H has no such factor.  Prediction:
  * tied images of H(R): equal (to print rounding);
  * tied images of A(R): NOT equal for pairs whose n_i is odd along some axis;
  * refold-onto-own-centres error for A lives on multi-image (tie) entries only.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import numpy as np
from pythtb import W90
from wannier.position_matrix import _map_R_by_orbital_positions

GRID = np.array((8, 8, 8))
DATA = Path(__file__).resolve().parents[3] / "data"

D = DATA / "Y2Ir2O7/base/soc/u_3.0/output/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181538bc889" / "jaemo"
w = W90(str(D), "Y2Ir2O7")
centres = np.asarray(w.lattice.orb_vecs, float) @ np.asarray(w.lat, float)
H = {R: np.asarray(b["h"]) for R, b in w.ham_r.items()}
A = {R: np.asarray(X) for R, X in w.pos_r.items()}


def classes(keys):
    out = {}
    for R in keys:
        out.setdefault(tuple(int(x) for x in np.mod(np.asarray(R), GRID)), []).append(R)
    return out


def rep(key):
    return tuple(int(((k + g // 2) % g) - g // 2) for k, g in zip(key, GRID))


def refold(blocks):
    sums = {}
    for R, b in blocks.items():
        k = tuple(int(x) for x in np.mod(np.asarray(R), GRID))
        sums[k] = np.array(b, complex) if k not in sums else sums[k] + b
    return _map_R_by_orbital_positions({rep(k): v for k, v in sums.items()}, w.lat, centres, GRID)


cls = classes(H)
n = next(iter(H.values())).shape[-1]
zeroH = np.zeros((n, n), complex)
zeroA = np.zeros((3, n, n), complex)
HM, AM = refold(H), refold(A)

# --- 1. refold error split by number of images carrying weight -------------
for name, raw, mapped, zero in (("H", H, HM, zeroH), ("A", A, AM, zeroA)):
    err = {1: 0.0, 2: 0.0, 3: 0.0}
    cnt = {1: 0, 2: 0, 3: 0}
    bad = {1: 0, 2: 0, 3: 0}
    for imgs in cls.values():
        r = np.array([raw.get(R, zero) for R in imgs])
        m = np.array([mapped.get(R, zero) for R in imgs])
        axes = tuple(range(1, r.ndim - 2))  # component axis (A only)
        occupied = (np.abs(r) > 1e-6)
        if r.ndim == 4:  # A: (image, 3, n, n)
            occupied = occupied.any(axis=1)
            diff = np.abs(r - m).max(axis=1)
        else:
            diff = np.abs(r - m)
        nimg = np.minimum(occupied.sum(axis=0), 3)
        dmax = diff.max(axis=0)
        for k in (1, 2, 3):
            sel = nimg == k
            cnt[k] += int(sel.sum())
            if sel.any():
                err[k] = max(err[k], float(dmax[sel].max()))
                bad[k] += int((dmax[sel] > 5e-6).sum())
    print(f"[{name}] refold-onto-own-centres error by images carrying weight (entries with >=1 image):")
    for k in (1, 2, 3):
        print(f"      {k}{'+' if k == 3 else ' '} image(s): {cnt[k]:>9d} entries, max|diff|={err[k]:.2e}, entries >5e-6: {bad[k]}")

# --- 2. are tied images equal copies? --------------------------------------
print("\nTied pairs (exactly two images carry H weight): are the two images equal copies?")
pairs = {"even": [0, 0.0, 0.0], "odd": [0, 0.0, 0.0]}
for imgs in cls.values():
    if len(imgs) < 2:
        continue
    h = np.array([H[R] for R in imgs]); a = np.array([A[R] for R in imgs])
    occ = np.abs(h) > 1e-6
    two = occ.sum(axis=0) == 2
    if not two.any():
        continue
    i0 = np.argmax(occ, axis=0)                      # first occupied image
    for ii in range(len(imgs)):
        for jj in range(ii + 1, len(imgs)):
            sel = two & occ[ii] & occ[jj]
            if not sel.any():
                continue
            nvec = (np.asarray(imgs[ii]) - np.asarray(imgs[jj])) // 8   # L = 8 n
            parity = "odd" if (np.abs(nvec) % 2).any() else "even"
            dH = np.abs(h[ii] - h[jj])[sel]
            dA = np.abs(a[ii] - a[jj]).max(axis=0)[sel]
            scale = np.maximum(np.abs(a[ii]).max(axis=0), np.abs(a[jj]).max(axis=0))[sel]
            pairs[parity][0] += int(sel.sum())
            pairs[parity][1] = max(pairs[parity][1], float(dH.max()))
            pairs[parity][2] = max(pairs[parity][2], float((dA / np.maximum(scale, 1e-12)).max()))
for parity, (count, dh, da) in pairs.items():
    print(f"   n = (R_i - R_j)/8 has {parity} component(s): {count:>8d} tied entries;  max|H_i - H_j| = {dh:.2e} eV;  max |A_i - A_j| / |A| = {da:.2f}")
