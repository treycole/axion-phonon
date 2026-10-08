"""Functions for comparing the occupied traced Berry curvature of two Ti_Q_2A trials.

Both Wannierizations describe the same 38-band occupied manifold (the lowest
Ti-semicore Kramers pair is outside the target manifold).  The comparison uses
the center-aware position matrix reconstructed from each trial's MMN overlaps,
with the same orbital-dependent real-space mapping applied to H(R) and A(R).

Pure library: no argparse, no ``__main__``. Used directly by
``compare_trial_curvature.ipynb`` and imported by
``external_curvature/check_vs_wilson_flux.py``,
``symmetry/check_trial_curvature_geometry.py``,
``symmetry/check_rdat_ndegen_applied_symmetry.py``, and
``symmetry/compare_position_matrix_sources.py``.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pythtb import W90

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import (
    RESULTS,
    STO_AA_CACHE,
    STO_TI_Q_2A_RUN,
)
from validation._paths import DATA
from modules import curvature as cv  # noqa: E402
from wannier.position_matrix import build_model_with_position_matrix  # noqa: E402
from wannier.wannier_io import check_mmn_compatibility, read_wannier_checkpoint
from wannier.position_matrix import (
    apply_orbital_dependent_R_mapping,
    build_position_matrix_from_mmn,
)


PREFIX = "SrTiO3"
RUN = STO_TI_Q_2A_RUN
MMN_DONOR = RUN / "trial_03_Sr_sp_Ti_pd_O_sp/proj_gauge"
CACHE_DIR = STO_AA_CACHE
OUTPUT_DIR = RESULTS / "trial_curvature"
N_OCC = 38

TRIALS = {
    "trial03": {
        "directory": MMN_DONOR,
        "label": "trial 03: Sr sp, Ti pd, O sp (48 WFs)",
        "color": "C0",
        "cache": CACHE_DIR / "AA_Ti_Q_2A_trial03_48wf_mmn_pair_ws1em05.npz",
    },
    "trial04": {
        "directory": DATA / "SrTiO3/Ti_Q_2A/output/188914bc889/trial_04_Sr_sp_Ti_p_O_sp/proj_gauge/188927bc889",
        "label": "trial 04: Sr sp, Ti p, O sp (38 WFs)",
        "color": "C3",
        "cache": CACHE_DIR / "AA_Ti_Q_2A_trial04_38wf_mmn_pair_ws1em05.npz",
    },
}

COMPONENTS = ((1, 2, "yz"), (2, 0, "zx"), (0, 1, "xy"))


def to_cartesian(model, curvature_reduced: np.ndarray) -> np.ndarray:
    """Convert both reciprocal indices of a Berry-curvature two-form."""

    reciprocal_inverse = np.linalg.inv(
        np.asarray(model.recip_lat_vecs, dtype=float)
    )
    return np.einsum(
        "ia,jb,ab...->ij...",
        reciprocal_inverse,
        reciprocal_inverse,
        curvature_reduced,
    )


def build_trial(key: str, trials: dict = TRIALS):
    """Build one consistently mapped H(R)/A(R) model.

    ``trials[key]["source"]`` picks the position-matrix source (default
    ``"mmn"``, the production route); any other key of
    ``wannier.position_matrix.POSITION_SOURCES`` -- e.g.
    ``"rdat_ndegen_applied"`` -- goes through the shared entry point instead,
    with no ``.mmn`` read at all.
    """

    setup = trials[key]
    directory = Path(setup["directory"])
    source = setup.get("source", "mmn")
    checkpoint = read_wannier_checkpoint(directory / f"{PREFIX}.chk")
    if checkpoint["num_wann"] < N_OCC:
        raise ValueError(
            f"{key} has only {checkpoint['num_wann']} WFs for {N_OCC} occupied bands"
        )

    w90 = W90(str(directory), PREFIX)
    if source == "mmn":
        mmn_directory = Path(setup.get("mmn_directory", MMN_DONOR))
        check_mmn_compatibility(directory, mmn_directory, PREFIX)
        apply_orbital_dependent_R_mapping(w90, checkpoint["mp_grid"])
        model = w90.model(zero_energy=0.0, min_hopping_norm=0.0)
        diagnostics = build_position_matrix_from_mmn(
            model,
            directory,
            PREFIX,
            mmn_directory=mmn_directory,
            cache=setup["cache"],
            translational_invariant_JM=False,
            translational_invariant_MV=False,
            R_distance_tolerance=1e-5,
        )
        print(
            f"{key}: {model.norb} WFs, {N_OCC} occupied; source=mmn; "
            f"A(R) Hermiticity={diagnostics['mmn_pair_hermiticity_relative']:.3e}, "
            f"outside-H(R) weight={diagnostics['weight_outside_ham_R']:.3e}"
        )
    else:
        model, result = build_model_with_position_matrix(
            w90,
            directory=directory,
            prefix=PREFIX,
            source=source,
            model_options=dict(zero_energy=0.0, min_hopping_norm=0.0),
        )
        print(
            f"{key}: {model.norb} WFs, {N_OCC} occupied; source={source}; "
            f"reads {result.reads}; diagnostics={result.diagnostics}"
        )
    pos = cv.position_terms(model, dtype=np.complex128)
    return model, pos


def evaluate_trial(model, pos, kpoints: np.ndarray) -> dict[str, np.ndarray]:
    """Return energies and the three Cartesian components of Tr Omega."""

    f = cv.fields(model, pos, kpoints)
    curvature = to_cartesian(model, cv.omega(cv.connection(f, N_OCC)))
    trace = np.column_stack(
        [
            np.trace(curvature[i, j], axis1=-2, axis2=-1).real
            for i, j, _ in COMPONENTS
        ]
    )
    energies = np.linalg.eigvalsh(f.H)[..., :N_OCC]
    return {"trace": trace, "energies": energies, "curvature": curvature}


def comparison_metrics(results: dict[str, dict[str, np.ndarray]]) -> dict:
    """Compute pointwise and norm-scaled differences between the trials."""

    left = results["trial03"]["trace"]
    right = results["trial04"]["trace"]
    difference = right - left
    metrics: dict[str, object] = {
        "n_occupied": N_OCC,
        "max_abs_energy_difference_eV": float(
            np.max(
                np.abs(
                    results["trial04"]["energies"]
                    - results["trial03"]["energies"]
                )
            )
        ),
        "components": {},
    }
    for index, (_, _, name) in enumerate(COMPONENTS):
        scale = max(float(np.max(np.abs(left[:, index]))), 1e-30)
        reference_l2 = float(np.linalg.norm(left[:, index]))
        metrics["components"][name] = {
            "trial03_max_abs_A2": float(np.max(np.abs(left[:, index]))),
            "trial04_max_abs_A2": float(np.max(np.abs(right[:, index]))),
            "max_abs_difference_A2": float(np.max(np.abs(difference[:, index]))),
            "rms_difference_A2": float(
                np.sqrt(np.mean(np.abs(difference[:, index]) ** 2))
            ),
            "max_difference_over_trial03_peak": float(
                np.max(np.abs(difference[:, index])) / scale
            ),
            "relative_l2_difference": float(
                np.linalg.norm(difference[:, index]) / reference_l2
                if reference_l2 > 0
                else np.nan
            ),
            "pearson_correlation": float(
                np.corrcoef(left[:, index], right[:, index])[0, 1]
                if np.std(left[:, index]) > 0 and np.std(right[:, index]) > 0
                else np.nan
            ),
        }
    vector_scale = max(float(np.linalg.norm(left)), 1e-30)
    metrics["all_components_relative_l2_difference"] = float(
        np.linalg.norm(difference) / vector_scale
    )
    metrics["all_components_max_abs_difference_A2"] = float(
        np.max(np.abs(difference))
    )
    return metrics


def ab_initio_grid_subspace_metrics(trials: dict = TRIALS) -> dict[str, float]:
    """Verify that the checkpoints span the same frozen occupied projector."""

    trial03_dir = Path(trials["trial03"]["directory"])
    trial04_dir = Path(trials["trial04"]["directory"])
    checkpoint03 = read_wannier_checkpoint(trial03_dir / f"{PREFIX}.chk")
    checkpoint04 = read_wannier_checkpoint(trial04_dir / f"{PREFIX}.chk")

    raw_eigenvalues = np.loadtxt(trial03_dir / f"{PREFIX}.eig")
    eigenvalues = np.empty(
        (checkpoint03["num_kpts"], checkpoint03["num_bands"]), dtype=float
    )
    eigenvalues[
        raw_eigenvalues[:, 1].astype(int) - 1,
        raw_eigenvalues[:, 0].astype(int) - 1,
    ] = raw_eigenvalues[:, 2]

    minimum_singular_value = 1.0
    maximum_projector_difference = 0.0
    for ik in range(checkpoint03["num_kpts"]):
        V03 = checkpoint03["v_matrix"][ik]
        V04 = checkpoint04["v_matrix"][ik]
        H03 = V03.conj().T @ (eigenvalues[ik, :, None] * V03)
        _, eigenvectors03 = np.linalg.eigh(H03)
        occupied03 = V03 @ eigenvectors03[:, :N_OCC]
        singular_values = np.linalg.svd(
            occupied03.conj().T @ V04, compute_uv=False
        )
        minimum_singular_value = min(
            minimum_singular_value, float(singular_values.min())
        )
        projector_difference = np.linalg.norm(
            occupied03 @ occupied03.conj().T - V04 @ V04.conj().T
        )
        maximum_projector_difference = max(
            maximum_projector_difference, float(projector_difference)
        )

    return {
        "minimum_principal_overlap_singular_value": minimum_singular_value,
        "maximum_projector_frobenius_difference": maximum_projector_difference,
    }


def save_results(
    t: np.ndarray,
    kpoints: np.ndarray,
    results: dict[str, dict[str, np.ndarray]],
    metrics: dict,
    output_dir: Path = OUTPUT_DIR,
) -> None:
    """Write machine-readable data, summary metrics, and an overlay plot."""

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "generic_line_trace_curvature.csv"
    npz_path = output_dir / "generic_line_trace_curvature.npz"
    json_path = output_dir / "generic_line_trace_curvature_summary.json"
    plot_path = output_dir / "generic_line_trace_curvature.png"

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        header = ["t", "kx_red", "ky_red", "kz_red"]
        for _, _, component in COMPONENTS:
            header.extend(
                [
                    f"trial03_trace_omega_{component}_A2",
                    f"trial04_trace_omega_{component}_A2",
                    f"trial04_minus_trial03_{component}_A2",
                ]
            )
        writer.writerow(header)
        for point_index in range(len(t)):
            row = [t[point_index], *kpoints[point_index]]
            for component_index in range(len(COMPONENTS)):
                left = results["trial03"]["trace"][point_index, component_index]
                right = results["trial04"]["trace"][point_index, component_index]
                row.extend([left, right, right - left])
            writer.writerow(row)

    np.savez_compressed(
        npz_path,
        t=t,
        kpoints_reduced=kpoints,
        components=np.array([name for _, _, name in COMPONENTS]),
        trial03_trace=results["trial03"]["trace"],
        trial04_trace=results["trial04"]["trace"],
        trial03_occupied_energies=results["trial03"]["energies"],
        trial04_occupied_energies=results["trial04"]["energies"],
    )
    json_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    print(f"saved {csv_path}")
    print(f"saved {npz_path}")
    print(f"saved {json_path}")


def plot_comparison(t: np.ndarray, results: dict[str, dict[str, np.ndarray]], trials: dict = TRIALS):
    """Build (but do not save) the overlay figure; returns (fig, axes)."""

    fig, axes = plt.subplots(3, 1, figsize=(8.2, 8.4), sharex=True, dpi=180)
    for component_index, (i, j, component) in enumerate(COMPONENTS):
        axis = axes[component_index]
        for key, setup in trials.items():
            axis.plot(
                t,
                results[key]["trace"][:, component_index],
                color=setup["color"],
                linewidth=1.6,
                label=setup["label"],
            )
        axis.axhline(0.0, color="0.35", linewidth=0.7)
        axis.set_ylabel(rf"$\mathrm{{Tr}}\,\Omega_{{{component}}}$ [$\AA^2$]")
        axis.grid(alpha=0.2)
        if component_index == 0:
            axis.legend(frameon=False, fontsize=8)
    axes[-1].set_xlabel(r"$t$ in $\mathbf{k}(t)=(t,0.30,0.15)$")
    fig.suptitle(r"Ti-displaced SrTiO$_3$: occupied traced Berry curvature")
    fig.tight_layout()
    return fig, axes


def load_results(output_dir: Path = OUTPUT_DIR) -> tuple[dict, dict]:
    """Reload a previous run's arrays and metrics without recomputing."""

    npz_path = output_dir / "generic_line_trace_curvature.npz"
    json_path = output_dir / "generic_line_trace_curvature_summary.json"
    data = np.load(npz_path)
    results = {
        "trial03": {"trace": data["trial03_trace"], "energies": data["trial03_occupied_energies"]},
        "trial04": {"trace": data["trial04_trace"], "energies": data["trial04_occupied_energies"]},
    }
    metrics = json.loads(json_path.read_text())
    return {"t": data["t"], "kpoints": data["kpoints_reduced"], "results": results}, metrics
