"""SrTiO3 Ti Born effective charge from the mixed (k, u) curvature agrees with DFPT (the paper's check of the mixed
planes). Reads what analysis/validation/sto_born_charge.ipynb saved; run the notebook to refresh it.

Each trial set's Z*_{Ti,zz} on the densest mesh must agree with DFPT within AGREE, have converged in the mesh within
CONVERGED, and the three sets (different splits between ionic and electronic parts) must agree within SPREAD.
"""
from pathlib import Path

import numpy as np
import pytest

# ---- settings
RESULTS = Path("/Users/treycole/Repos/axion-phonon/data/SrTiO3/born_charge/born_charge.csv")
AGREE = 0.03                  # relative to DFPT (2026-10-31: +1.6%, +1.1%, +0.6%)
CONVERGED = 1e-3              # relative change between the two densest meshes
SPREAD = 0.02                 # relative, max - min over the trial sets
# --------------
pytestmark = pytest.mark.data


def test_born_charge():
    pd = pytest.importorskip("pandas")
    if not RESULTS.exists():
        pytest.skip(f"not in data/: {RESULTS} (run analysis/validation/sto_born_charge.ipynb)")
    r = pd.read_csv(RESULTS)
    finals = []
    for label, g in r.groupby("set"):
        g = g.sort_values("nk")
        z, z_prev = g["Z* (Wannier)"].iloc[-1], g["Z* (Wannier)"].iloc[-2]
        assert abs(g["vs DFPT"].iloc[-1]) < AGREE, f"{label}: {g['vs DFPT'].iloc[-1]:+.2%} from DFPT"
        assert abs(z / z_prev - 1) < CONVERGED, f"{label}: not converged in the mesh ({z_prev:.4f} -> {z:.4f})"
        assert g["|Im|"].max() < 1e-8, f"{label}: imaginary part {g['|Im|'].max():.1e}"
        finals.append(z)
    assert (max(finals) - min(finals)) / np.mean(finals) < SPREAD, f"trial sets disagree: {finals}"
