"""Write dtheta results (one nk_<nk>.npz per mesh, summary.npz, summary.csv). No physics here."""
from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np


def save_dtheta_results(
    output_directory: str | Path,
    results: Sequence[Mapping[str, object]],
    metadata: Mapping[str, object],
) -> tuple[Path, Path]:
    """Save completed mesh results and rebuild the NPZ/CSV sweep summaries."""

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    ordered = sorted(results, key=lambda result: int(result["nk"]))
    if not ordered:
        raise ValueError("No dtheta results to save.")

    for result in ordered:
        nk = int(result["nk"])
        np.savez_compressed(
            output_directory / f"nk_{nk}.npz",
            nk=np.int32(nk),
            d3k=np.float64(result["d3k"]),
            dtheta=result["dtheta"],
            dtheta_simpson=result["dtheta_simpson"],
            minimum_gap=result["minimum_gap"],
            minimum_gap_kpoint=result["minimum_gap_kpoint"],
            c2_density=result["c2_density"],
            # File-level aliases keep existing analysis notebooks readable.
            min_gap=result["minimum_gap"],
            min_gap_k=result["minimum_gap_kpoint"],
            chern2_density=result["c2_density"],
            **metadata,
        )

    nks = np.array([int(result["nk"]) for result in ordered], dtype=np.int32)
    dtheta = np.array([result["dtheta"] for result in ordered], dtype=np.complex64)
    dtheta_simpson = np.array(
        [result["dtheta_simpson"] for result in ordered], dtype=np.complex64
    )
    minimum_gap = np.array(
        [result["minimum_gap"] for result in ordered], dtype=np.float64
    )
    minimum_gap_kpoint = np.array(
        [result["minimum_gap_kpoint"] for result in ordered], dtype=np.float64
    )
    d3k = np.array([result["d3k"] for result in ordered], dtype=np.float64)

    summary_npz = output_directory / "summary.npz"
    np.savez_compressed(
        summary_npz,
        nks=nks,
        dtheta=dtheta,
        dtheta_simpson=dtheta_simpson,
        minimum_gap=minimum_gap,
        minimum_gap_kpoint=minimum_gap_kpoint,
        d3k=d3k,
        min_gap=minimum_gap,
        min_gap_k=minimum_gap_kpoint,
        **metadata,
    )

    summary_csv = output_directory / "summary.csv"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "nk",
                "dtheta_base_real",
                "dtheta_base_imag",
                "dtheta_mode_real",
                "dtheta_mode_imag",
                "dtheta_simpson_base_real",
                "dtheta_simpson_mode_real",
                "minimum_gap_base",
                "minimum_gap_mode",
                "minimum_gap_kpoint_base_x",
                "minimum_gap_kpoint_base_y",
                "minimum_gap_kpoint_base_z",
                "minimum_gap_kpoint_mode_x",
                "minimum_gap_kpoint_mode_y",
                "minimum_gap_kpoint_mode_z",
                "d3k",
            ]
        )
        for result in ordered:
            dtheta_endpoint = np.asarray(result["dtheta"])
            dtheta_simp = np.asarray(result["dtheta_simpson"])
            gaps = np.asarray(result["minimum_gap"])
            gap_k = np.asarray(result["minimum_gap_kpoint"])
            writer.writerow(
                [
                    int(result["nk"]),
                    dtheta_endpoint[0].real,
                    dtheta_endpoint[0].imag,
                    dtheta_endpoint[1].real,
                    dtheta_endpoint[1].imag,
                    dtheta_simp[0].real,
                    dtheta_simp[1].real,
                    gaps[0],
                    gaps[1],
                    *gap_k[0].tolist(),
                    *gap_k[1].tolist(),
                    result["d3k"],
                ]
            )

    return summary_npz, summary_csv
