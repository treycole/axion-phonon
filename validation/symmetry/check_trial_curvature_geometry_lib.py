"""Source-node and symmetry checks for the two Ti_Q_2A interpolants.

The raw-MMN check establishes equality of the discrete occupied subspaces and
Wilson loops.  This asks two separate downstream questions:

* how different are the continuous curvatures evaluated at every source-mesh
  node; and
* does each interpolant separately obey the surviving crystal symmetries away
  from those nodes?

The traced Cartesian curvature is represented as the axial vector
``(Omega_yz, Omega_zx, Omega_xy)``.  For an orthogonal spatial operation ``O``
it transforms as ``det(O) O Omega``; time reversal sends ``Omega(k)`` to
``-Omega(-k)``.

Pure library: no argparse, no ``__main__``. Used directly by
``check_trial_curvature_geometry.ipynb`` and imported by
``check_rdat_ndegen_applied_symmetry.py`` (``symmetry_checks``).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import RESULTS  # noqa: E402
from validation.symmetry.compare_trial_curvature_lib import (  # noqa: E402
    PREFIX,
    TRIALS,
    build_trial,
    comparison_metrics,
    evaluate_trial,
)

OUTPUT = RESULTS / "trial_geometry"


def relative_metrics(residual: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    """Return absolute and scale-aware residual norms."""

    denominator = float(np.linalg.norm(reference))
    return {
        "max_abs_A2": float(np.max(np.abs(residual))),
        "rms_A2": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
        "relative_l2": float(
            np.linalg.norm(residual) / denominator if denominator else np.nan
        ),
    }


def symmetry_checks(model, pos, points: np.ndarray) -> dict:
    """Test T, C4z, and Mx on generic paired off-grid points."""

    operations = {
        "time_reversal": {
            "matrix": -np.eye(3),
            "axial_matrix": -np.eye(3),
        },
        "C4z": {
            "matrix": np.array(
                [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
            ),
        },
        "Mx": {
            "matrix": np.diag([-1.0, 1.0, 1.0]),
        },
    }
    original = evaluate_trial(model, pos, points)["trace"]
    result = {}
    for name, operation in operations.items():
        spatial = operation["matrix"]
        transformed_points = (points @ spatial.T) % 1.0
        transformed = evaluate_trial(model, pos, transformed_points)["trace"]
        axial = operation.get(
            "axial_matrix", np.linalg.det(spatial) * spatial
        )
        expected = original @ axial.T
        result[name] = relative_metrics(transformed - expected, expected)
    return result
