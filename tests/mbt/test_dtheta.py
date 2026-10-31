"""MnBi2Te4 production dθ/dQ (window f10.962_w12.810_k8): the claims the paper makes about it, checked on the local
copies in data/.

- window edges: every edge at least EDGE_MIN from every DFT eigenvalue on the mesh, and the same number of bands below
  each edge at every k as the base (pipeline step 5's window_check.txt, base and every structure);
- k convergence: the dense tube list (16, 16, 5, 32) contains the production list (16, 4, 3, 8) and reproduces its
  values, and the production values are within K_CONVERGED of the dense ones (paper: "converged to 0.2%");
- the integrand is real and the gap open on every k of the production and dense lists;
- amplitude: for pattern 3, the Q -> 0 limit of a + bQ^2 through Q = 0.01 and 0.03 A is within Q_LIMIT of the
  undistorted-structure estimate at 0.01 A (paper: 0.05%).
"""
import re
from pathlib import Path

import numpy as np
import pytest

from axion_phonon import tube

# ---- settings
WINDOW = Path("/Users/treycole/Repos/axion-phonon/data/MnBi2Te4/windows/f10.962_w12.810_k8")
PATTERNS = ["T2-_p1_Q0.01", "T2-_p2_Q0.01", "T2-_p3_Q0.01", "T2-_p4_Q0.01"]
STRUCTURES = ["base"] + PATTERNS + ["T2-_p3_Q0.03"]
PRODUCTION = dict(m=4, n3=8, half=3)   # the production tube list inside the dense one
EDGE_MIN = 5.0                         # meV, as pipeline step 5
K_CONVERGED = 2e-3                     # relative, production vs dense
REPRODUCE = 1e-6                       # relative, dense list at the production setting vs dtheta.npz
IM_MAX = 1e-5                          # max |Im c2| / max |Re c2|
GAP_MIN = 0.15                         # eV, smallest direct gap (both ends) on any k
Q_LIMIT = 1e-3                         # relative
END = 0                                # undistorted-structure estimate
# --------------
pytestmark = pytest.mark.data


def need(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"not in data/: {path}")
    return path


def dense(pattern):
    D = need(WINDOW / f"dtheta_dense_{pattern}")
    kl = dict(np.load(D / "klist.npz"))
    c2 = np.full((len(kl["k"]), 2), np.nan, complex)
    gap = np.full((len(kl["k"]), 2), np.nan)
    for f in sorted((D / "c2").glob("out_*/chunks/chunk_*.npz")):
        z = np.load(f)
        c2[z["idx"]], gap[z["idx"]] = z["c2"], z["gap"]
    assert not np.isnan(c2.real).any(), f"{D}: missing chunks"
    return kl, c2, gap


def production(pattern):
    return np.load(need(WINDOW / f"dtheta_{pattern}" / "dtheta.npz"))["dtheta"][2, 0]   # [m 4, 8 planes], both ends


@pytest.mark.parametrize("structure", STRUCTURES)
def test_window_edges(structure):
    lines = [l for l in need(WINDOW / structure / "window_check.txt").read_text().splitlines() if "dis_" in l]
    assert len(lines) == 4, lines
    for l in lines:
        assert l.rstrip().endswith("OK"), l
        assert float(re.search(r"closest eigenvalue\s+([\d.]+) meV", l).group(1)) >= EDGE_MIN, l
        count = re.search(r"different count than the base: (\d+)", l)
        assert count is None or int(count.group(1)) == 0, l


@pytest.mark.parametrize("pattern", PATTERNS)
def test_k_converged(pattern):
    kl, c2, gap = dense(pattern)
    prod = production(pattern)
    at_prod = tube.integrate_at(c2, kl, **PRODUCTION)
    assert np.all(np.abs(at_prod / prod - 1) < REPRODUCE), (at_prod, prod)
    top = tube.integrate_at(c2, kl, m=int(kl["M"]), n3=int(kl["n3"]), half=int(kl["half"]))
    assert np.all(np.abs(prod / top - 1) < K_CONVERGED), f"production {prod}, dense {top}"


@pytest.mark.parametrize("pattern", PATTERNS)
def test_integrand_real_and_gap_open(pattern):
    z = np.load(need(WINDOW / f"dtheta_{pattern}" / "c2.npz"))
    _, c2_dense, gap_dense = dense(pattern)
    for name, c2, gap in (("production", z["c2"], z["gap"]), ("dense", c2_dense, gap_dense)):
        assert np.abs(c2.imag).max() / np.abs(c2.real).max() < IM_MAX, name
        assert gap.min() > GAP_MIN, (name, gap.min())


def test_amplitude_limit():
    a, b = production("T2-_p3_Q0.01")[END], production("T2-_p3_Q0.03")[END]
    limit = (9 * a - b) / 8                                   # a + bQ^2 through Q = 0.01 and 0.03
    assert abs(a / limit - 1) < Q_LIMIT, (a, b, limit)
