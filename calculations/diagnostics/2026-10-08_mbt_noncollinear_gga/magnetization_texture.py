"""Real-space magnetization texture of a noncollinear QE run (notes/progress/2026-10-08_mbt_noncollinear_gga_lda.md §8).

Reads <prefix>.save/charge-density.dat (n, m_x, m_y, m_z in G-space), inverse-FFTs to the dense grid and reports
1. the transverse part of m: integral |m_perp|, its maximum, and where it dominates;
2. how smooth the spin density each GGA recipe would see is, measured by its total variation sum |f(i+1) - f(i)|
   over grid links: m_z (collinear), |m| (local frame, kinked at m = 0), and s|m| with s = sign(m_z) (QE's
   fixed-axis rule for lsign = .true.), plus the links where s flips at finite |m| (jumps).

Needs about 4 GB of memory for a 375^3 grid. Run on the cluster login node with python3 (numpy, scipy).
Checks: integral n = number of valence electrons, and integral |m| = QE's "absolute magnetization".
"""
import numpy as np
from scipy.io import FortranFile
from scipy import fft

# ============================================================================== settings
BASE = "/mnt/wkdk/bc889/axion/MnBi2Te4/base"
RUNS = {  # label -> charge-density.dat
    "PBE+SOC asymmetric (209025)": f"{BASE}/soc/u_4.0/209025bc889/tmp/MnBi2Te4.save/charge-density.dat",
    "PBE+SOC symmetric (test D)": f"{BASE}/soc/u_4.0/test_D_2026-10-07/209052bc889/tmp/MnBi2Te4.save/charge-density.dat",
    "LDA+SOC base (209653)": f"{BASE}/soc_lda/u_4.0/209653bc889/tmp/MnBi2Te4.save/charge-density.dat",
}
FFT_GRID = (375, 375, 375)   # dense FFT grid, from "FFT dimensions" in scf.out
OMEGA = 2994.5593            # cell volume in bohr^3, from "unit-cell volume" in scf.out
JUMP_THRESHOLD = 0.01        # a sign flip of m_z counts as a jump if |m| > this fraction of max|m| on both ends
# ======================================================================================


def read_density(path, grid):
    """Return [n, m_x, m_y, m_z] on the real-space grid (e/bohr^3, muB/bohr^3)."""
    f = FortranFile(path, "r")
    _gamma_only, ngm, nspin = f.read_ints(np.int32)
    f.read_reals(np.float64)                                  # reciprocal vectors b1, b2, b3
    mill = f.read_ints(np.int32).reshape(ngm, 3)
    fields = []
    for _ in range(nspin):
        coeff = f.read_reals(np.complex128)
        g = np.zeros(grid, complex)
        g[mill[:, 0] % grid[0], mill[:, 1] % grid[1], mill[:, 2] % grid[2]] = coeff
        fields.append((fft.ifftn(g, workers=1) * np.prod(grid)).real)
        del g
    f.close()
    return fields


def total_variation(field):
    return sum(np.abs(np.roll(field, -1, ax) - field).sum() for ax in range(3))


for label, path in RUNS.items():
    n, mx, my, mz = read_density(path, FFT_GRID)
    dv = OMEGA / n.size
    m_perp = np.hypot(mx, my)
    m_abs = np.sqrt(m_perp**2 + mz**2)
    print(f"#### {label}")
    print(f"  integral n = {n.sum() * dv:.4f} e;  integral |m| = {m_abs.sum() * dv:.4f} muB")
    print(f"  integral |m_z| = {np.abs(mz).sum() * dv:.4f};  integral |m_perp| = {m_perp.sum() * dv:.4f};"
          f"  integral (|m| - |m_z|) = {(m_abs - np.abs(mz)).sum() * dv:.4f}")
    print(f"  max |m_perp| = {m_perp.max():.4e};  max |m| = {m_abs.max():.4e} muB/bohr^3;"
          f"  grid fraction |m_perp| > |m_z| = {(m_perp > np.abs(mz)).mean():.3f}")

    s = np.sign(mz)
    fixed_axis = s * m_abs
    tv_z, tv_abs, tv_fixed = total_variation(mz), total_variation(m_abs), total_variation(fixed_axis)
    thr = JUMP_THRESHOLD * m_abs.max()
    jump_links, var_fixed_on_flips, var_z_on_flips = 0, 0.0, 0.0
    for ax in range(3):
        flip = np.roll(s, -1, ax) != s
        finite = np.minimum(m_abs, np.roll(m_abs, -1, ax)) > thr
        jump_links += int((flip & finite).sum())
        var_fixed_on_flips += np.abs(np.roll(fixed_axis, -1, ax) - fixed_axis)[flip].sum()
        var_z_on_flips += np.abs(np.roll(mz, -1, ax) - mz)[flip].sum()
    excess = tv_fixed - tv_z
    print(f"  total variation: m_z {tv_z:.1f}, |m| {tv_abs:.1f}, s|m| {tv_fixed:.1f}"
          f" (excess over m_z {excess:.1f} = {100 * excess / tv_z:.1f}%)")
    print(f"  jump links (s flips at finite |m|): {jump_links};"
          f"  share of the excess on flip links: {100 * (var_fixed_on_flips - var_z_on_flips) / max(excess, 1e-12):.0f}%")
    del n, mx, my, mz, m_perp, m_abs, s, fixed_axis
