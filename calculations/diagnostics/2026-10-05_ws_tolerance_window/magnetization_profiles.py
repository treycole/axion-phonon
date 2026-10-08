"""m(r) along real-space paths, from a QE charge-density.dat (G-space n, mx, my, mz):
f(r) = sum_G f(G) exp(i G.r), G = mill @ b (b stored in 1/bohr).  Smooth + augmentation density, not all-electron.
YIO: Ir1 -> Ir4 straight bond, and Ir1 -> bridging O -> Ir4.  MBT: Mn1 -> Mn2 along the 3-fold axis."""
import sys
import numpy as np
from scipy.io import FortranFile
BOHR = 0.529177210903
S = "/private/tmp/claude-501/-Users-treycole-Repos-axion-phonon/fd07904c-2fb6-4eec-a482-f69ca09d63a9/scratchpad/"

def load(path):
    f = FortranFile(path, "r")
    _, ngm, nspin = f.read_ints(np.int32)
    b = f.read_reals(np.float64).reshape(3, 3)
    mill = f.read_ints(np.int32).reshape(ngm, 3)
    c = np.array([f.read_reals(np.complex128) for _ in range(nspin)])
    return b, mill, c

def evaluate(b, mill, c, alat, pts_ang, chunk=20):
    G = mill @ b                                              # b is stored in 1/bohr (a_i . b_j = 2 pi)
    out = []
    for i in range(0, len(pts_ang), chunk):
        r = pts_ang[i:i + chunk] / BOHR
        ph = np.exp(1j * (r @ G.T))                            # (npt, ngm)
        out.append(np.real(ph @ c.T))                          # (npt, 4)
    return np.concatenate(out)

def path(points, n_per):
    seg = [a + (b - a) * t for a, b in zip(points[:-1], points[1:]) for t in np.linspace(0, 1, n_per, endpoint=False)]
    return np.array(seg + [points[-1]])

def angle(u, v):
    return np.degrees(np.arccos(np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-30), -1, 1)))

# ---------------- YIO ----------------
lat = np.array([[0, 5.095, 5.095], [5.095, 0, 5.095], [5.095, 5.095, 0]])
Ir1, Ir4 = np.zeros(3), np.array([2.5475, 2.5475, 0.0])
O48 = np.array([[0.53033, 0.53033, 0.53033], [0.88388, 0.88388, 0.88388]])  # placeholders, replaced below
win = open("/Users/treycole/Repos/axion-phonon/data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889/Y2Ir2O7.win").read().splitlines()
i0 = win.index("begin atoms_cart"); atoms = []
for l in win[i0 + 2:]:
    if l.strip().lower().startswith("end"): break
    p = l.split(); atoms.append((p[0], np.array(p[1:4], float)))
imgs = np.array([(i, j, k) for i in range(-1, 2) for j in range(-1, 2) for k in range(-1, 2)]) @ lat
Os = np.array([a[1] + t for a in atoms if a[0] == "O" for t in imgs])
bridge = Os[np.argmin(np.linalg.norm(Os - Ir1, axis=1) + np.linalg.norm(Os - Ir4, axis=1))]
u1 = np.array([-1, -1, -1]) / np.sqrt(3); u4 = np.array([1, 1, -1]) / np.sqrt(3)
b, mill, c = load(S + "rho_193145.dat")
alat = 13.6163
for name, pts in (("Ir1 -> Ir4 (straight, 3.60 A)", path([Ir1, Ir4], 12)),
                  (f"Ir1 -> O {np.round(bridge, 3)} -> Ir4", path([Ir1, bridge, Ir4], 6))):
    v = evaluate(b, mill, c, alat, pts)
    print(f"\nYIO {name}  (193145 density; ux = Ir1 axis (-1,-1,-1))")
    print("   s (A)     n        |m|      |m|/|m(Ir1)|   angle to Ir1 axis  angle to Ir4 axis   sign(m.ux)")
    s_along = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    m0 = np.linalg.norm(v[0, 1:])
    for sp, row in zip(s_along, v):
        m = row[1:]
        print(f"  {sp:6.2f}  {row[0]:8.4f}  {np.linalg.norm(m):9.2e}   {np.linalg.norm(m)/m0:7.3f}        {angle(m, u1):6.1f}            {angle(m, u4):6.1f}          {int(np.sign(m @ u1)):+d}")

# ---------------- MBT ----------------
b, mill, c = load(S + "rho_mbt_202833.dat")
alat = 51.7756
pts = np.array([[0, 0, z] for z in np.linspace(0, 40.926, 31)])
v = evaluate(b, mill, c, alat, pts)
print("\nMBT Mn1 (z=0) -> Mn2 (z=40.93 A) along the 3-fold axis (202833, iteration 115, not converged)")
print("  atoms on this line: Te 5.46, Te 12.03, Bi 17.39, Bi 23.54, Te 28.89, Te 35.47 (A)")
print("    z (A)     n        m_z        |m_perp|")
for p, row in zip(pts, v):
    print(f"  {p[2]:6.2f}  {row[0]:8.4f}  {row[3]:+10.2e}  {np.hypot(row[1], row[2]):9.2e}")
