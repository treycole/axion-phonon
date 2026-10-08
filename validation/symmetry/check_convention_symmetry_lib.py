"""Compare four internally consistent SrTiO3 interpolation conventions.

The pipelines are deliberately kept separate:

``mmn_pair``
    Build the centre-aware position matrix from ``.mmn``/``.chk`` and reassign
    both it and ``H(R)`` by the pair-dependent physical distance
    ``|R + tau_t - tau_s|``.

``rdat_direct``
    Read ``A(R)`` directly from ``_r.dat`` and ``H(R)`` directly from
    ``_hr.dat``.  No Wigner--Seitz reassignment is applied to either matrix.

``rdat_pair``
    Read ``A(R)`` from the ordinary ``_r.dat``, strip Wannier90's scalar
    Wigner--Seitz weights, and re-distribute both it and ``H(R)`` by
    ``|R + tau_t - tau_s|`` -- i.e. apply the ``use_ws_distance`` rule in house.
    Needs no ``.mmn`` and no extra files, so this is the cheap route.

``rfull_pair``
    Read ``A(R)`` from postw90's ``_r_full.dat`` (``transl_inv_full=T`` with
    ``use_ws_distance=T``, so the blocks already carry pair-dependent effective
    weights) and reassign ``H(R)`` by the same physical distance used above.
    This is the native Wannier90 ordering: it keeps every neighbour link ``b``
    separate through the Fourier transform and maps to Wigner--Seitz afterwards,
    whereas ``mmn_pair`` sums the links first.  The two do not commute at finite
    mesh density, so this pipeline is a genuine third convention, not a
    reimplementation of the first.

For every model, the Hamiltonian and position matrices are asserted to have
exactly the same R-vector set. Tr Omega_xy is then evaluated on:

* two paths in the undisplaced P+T-symmetric structure (zero everywhere),
* a path lying in each Ti-z or Sr-z displaced structure's kx=1/2 mirror plane
  (zero everywhere), and
* a generic displaced path k=(t, 0.30, 0.15), whose intersections with mirror
  planes at t=0, 0.3, 0.5, 0.7 and 1 are symmetry-enforced zeros.

Pure library: no argparse, no ``__main__``. Used directly by
``check_convention_symmetry.ipynb`` and imported by
``check_pr702_rdat_symmetry.py`` (as ``ccs``).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pythtb import W90

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import (
    RESULTS,
    STO_AA_CACHE,
    STO_BASE_RUN,
    STO_SR_Q_01A_RUN,
    STO_TI_Q_01A_RUN,
)
from validation._paths import DATA
from modules import curvature as cv
from wannier.position_matrix import build_model_with_position_matrix

PREFIX = "SrTiO3"
BASE = STO_BASE_RUN
MODE = STO_TI_Q_01A_RUN
SR_MODE = STO_SR_Q_01A_RUN
AA_CACHE = STO_AA_CACHE
OUTPUT_DIR = RESULTS / "convention_symmetry"

TRIALS = {
    "a": {
        "base": DATA / "SrTiO3/base/soc/output/181642bc889/trial_01_Sr_spd_Ti_spd_O_sp/182494bc889",
        "mode": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_01_Sr_spd_Ti_spd_O_sp/no_displacement/182499bc889",
        "base_mmn": BASE,   # the NSCF root now holds the base .mmn
        "mode_mmn": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889",
        "n_occ": 30,
        "label": "set (a): 60 WFs, 30 occupied",
    },
    "b": {
        "base": DATA / "SrTiO3/base/soc/output/181642bc889/trial_02_Sr_p_O_sp/183813bc889",
        "mode": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_02_Sr_p_O_sp/no_displacement/183965bc889",
        "base_mmn": BASE,   # the NSCF root now holds the base .mmn
        "mode_mmn": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889",
        "n_occ": 30,
        "label": "set (b): 30 WFs, occupied-only",
    },
    "c": {
        "base": DATA / "SrTiO3/base/soc/output/181642bc889/trial_03_Sr_sp_Ti_pd_O_sp/185116bc889",
        "mode": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889",
        "base_mmn": BASE,   # the NSCF root now holds the base .mmn
        "mode_mmn": DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889",
        "n_occ": 38,
        "label": "set (c): 48 WFs, 38 occupied",
    },
    "sr_b": {
        "base": DATA / "SrTiO3/base/soc/output/181642bc889/trial_02_Sr_p_O_sp/183813bc889",
        "mode": DATA / "SrTiO3/Sr_Q_01A/output/184367bc889/trial_02_Sr_p_O_sp/no_displacement/188438bc889",
        "base_mmn": BASE,
        "mode_mmn": SR_MODE,
        "n_occ": 30,
        "label": "Sr-z set (b): 30 WFs, occupied-only",
    },
    "sr_c": {
        "base": DATA / "SrTiO3/base/soc/output/181642bc889/trial_03_Sr_sp_Ti_pd_O_sp/185116bc889",
        "mode": DATA / "SrTiO3/Sr_Q_01A/output/184367bc889/trial_03_Sr_sp_Ti_pd_O_sp/no_displacement/188444bc889",
        "base_mmn": BASE,
        "mode_mmn": SR_MODE,
        "n_occ": 38,
        "label": "Sr-z set (c): 48 WFs, 38 occupied",
    },
}

CONVENTIONS = {
    "mmn_pair": {
        "label": ".mmn, centre-aware; pair-dependent H/A",
        "color": "C0",
        "linestyle": "-",
    },
    "rdat_direct": {
        "label": "direct _r.dat + _hr.dat",
        "color": "C1",
        "linestyle": "--",
    },
    "rfull_pair": {
        "label": "postw90 _r_full.dat (transl_inv_full); pair-dependent H/A",
        "color": "C2",
        "linestyle": ":",
    },
    "rdat_pair": {
        "label": "_r.dat + use_ws_distance remap; tau = Wannier centres",
        "color": "C3",
        "linestyle": "-.",
    },
    "rdat_pair_tau0": {
        "label": "_r.dat + use_ws_distance remap; tau = 0",
        "color": "C4",
        "linestyle": (0, (3, 1, 1, 1, 1, 1)),
    },
}

FORCED_ZEROS = np.array([0.0, 0.3, 0.5, 0.7, 1.0])


def to_cartesian(model, curvature_reduced: np.ndarray) -> np.ndarray:
    """Convert both reduced reciprocal indices of a curvature two-form."""
    reciprocal_inverse = np.linalg.inv(np.asarray(model.recip_lat_vecs, dtype=float))
    return np.einsum(
        "ia,jb,ab...->ij...",
        reciprocal_inverse,
        reciprocal_inverse,
        curvature_reduced,
    )


def trace_omega_xy(model, pos: cv.PositionTerms, kpoints, n_occ: int) -> np.ndarray:
    """Occupied-band trace of the total Cartesian Omega_xy."""
    kpoints = np.asarray(kpoints, dtype=float)
    f = cv.fields(model, pos, kpoints)
    curvature_total = to_cartesian(model, cv.omega(cv.connection(f, n_occ)))
    return np.trace(curvature_total[0, 1], axis1=-2, axis2=-1).real


def _cache_path(trial: str, structure: str, directory: Path, aa_cache: Path = AA_CACHE) -> Path:
    if trial in {"sr_b", "sr_c"}:
        trial_dir = {
            "sr_b": "trial_02_Sr_p_O_sp",
            "sr_c": "trial_03_Sr_sp_Ti_pd_O_sp",
        }[trial]
        return aa_cache / (
            f"AA_{trial_dir}_{structure}_{directory.name}_wbnone_ws1em05.npz"
        )
    return aa_cache / (
        f"AA_symmetry_{trial}_{structure}_{directory.name}_wbnone_ws1em05.npz"
    )


def convention_available(trial: str, structure: str, convention: str, trials: dict = TRIALS) -> bool:
    """Whether the inputs this convention needs are actually on disk.

    Only ``rfull_pair`` can be missing: ``_r_full.dat`` has to be produced by a
    separate postw90 run per Wannierization, and not every trial set has one.
    ``rdat_pair`` needs nothing beyond the ``_r.dat`` every run already writes.
    """
    if convention in ("rdat_pair", "rdat_pair_tau0"):
        return True
    if convention == "rfull_pair":
        return (
            Path(trials[trial][structure]) / f"{PREFIX}_r_full.dat"
        ).is_file()
    return True


def build_model(
    trial: str,
    structure: str,
    convention: str,
    trials: dict = TRIALS,
) -> tuple[object, cv.PositionTerms, dict]:
    """Build one model under one convention, via the shared entry point.

    Every route lives in ``wannier/position_matrix.py``; this function only
    supplies the SrTiO3-specific paths and the two local variations that are
    not properties of the source itself:

    * ``mmn_pair`` borrows its ``.mmn`` from a donor trial (same nscf,
      different projections) and uses a named cache;
    * ``rdat_pair_tau0`` pins the embedding to tau = 0 rather than the Wannier
      centres.  That is a ``model_options`` choice, not a different A(R).
    """
    setup = trials[trial]
    directory = Path(setup[structure])
    w90 = W90(str(directory), PREFIX)

    # Convention name -> position-matrix source.  The two `rdat_pair*` names
    # differ only in the embedding, so they share one source.
    source = {
        "mmn_pair": "mmn",
        "rdat_direct": "rdat_direct",
        "rfull_pair": "rfull",
        "rdat_pair": "rdat_pair",
        "rdat_pair_tau0": "rdat_pair",
    }.get(convention)
    if source is None:
        raise ValueError(f"Unknown convention {convention!r}")

    options = dict(zero_energy=0.0, min_hopping_norm=0.0)
    if convention == "rdat_pair_tau0":
        # Set tau on the model, never on w90.lattice: TBModel.lattice returns
        # a copy, so mutating it is a silent no-op that would leave H and v at
        # the Wannier centres while only the external terms moved.
        options["orb_vecs"] = np.zeros_like(
            np.asarray(w90.lattice.orb_vecs, dtype=float)
        )

    donor = Path(setup[f"{structure}_mmn"]) if convention == "mmn_pair" else None
    model, result = build_model_with_position_matrix(
        w90,
        directory=directory,
        prefix=PREFIX,
        source=source,
        model_options=options,
        mmn_directory=donor,
        cache=_cache_path(trial, structure, directory) if donor else None,
        R_distance_tolerance=1e-5,
    )

    if convention == "rdat_pair_tau0" and np.abs(model.orb_vecs).max() > 1e-14:
        raise RuntimeError("tau = 0 embedding did not take effect")

    info = result.diagnostics
    pos = cv.position_terms(model, dtype=np.complex128)
    provenance = {
        "position_source": result.reads,
        "ws_rule": result.ws_rule,
        "n_R_H": result.n_R,
        "n_R_A": result.n_R,
        "same_R_set": True,
        "zero_H_blocks_added": int(info.get("zero_H_blocks_added", 0)),
        "mmn_pair_hermiticity_relative": info.get(
            "mmn_pair_hermiticity_relative", np.nan
        ),
        "rdat_pair_hermiticity_relative": info.get(
            "rdat_pair_hermiticity_relative", np.nan
        ),
    }
    return model, pos, provenance


def calculate(n_points: int, trials=None, conventions=None, trial_set: dict = TRIALS,
              convention_set: dict = CONVENTIONS):
    trials = tuple(trials or trial_set)
    conventions = tuple(conventions or convention_set)
    t = np.linspace(0.0, 1.0, n_points)
    paths = {
        "generic": np.column_stack((t, np.full_like(t, 0.30), np.full_like(t, 0.15))),
        "mirror_kx_half": np.column_stack(
            (np.full_like(t, 0.50), t, np.full_like(t, 0.15))
        ),
    }
    curves = {}
    provenance = {}
    for trial in trials:
        setup = trial_set[trial]
        for convention in conventions:
            missing = [
                structure
                for structure in ("base", "mode")
                if not convention_available(trial, structure, convention, trial_set)
            ]
            if missing:
                print(
                    f"{trial}/{convention}: skipped, no {PREFIX}_r_full.dat for "
                    f"{', '.join(missing)}",
                    flush=True,
                )
                continue
            for structure in ("base", "mode"):
                print(f"{trial}/{structure}/{convention}", flush=True)
                model, pos, prov = build_model(trial, structure, convention, trial_set)
                provenance[(trial, structure, convention)] = prov
                for path_name, kpoints in paths.items():
                    curves[(trial, structure, convention, path_name)] = trace_omega_xy(
                        model, pos, kpoints, setup["n_occ"]
                    )
    return t, paths, curves, provenance


def write_outputs(output_dir: Path, t, paths, curves, provenance,
                   trial_set: dict = TRIALS, convention_set: dict = CONVENTIONS) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    # Only the (trial, convention) pairs calculate() could actually build.
    done = sorted({(trial, conv) for trial, _, conv in provenance})
    trials_done = [trial for trial in trial_set if any(a == trial for a, _ in done)]
    conventions_done = [
        conv for conv in convention_set if any(b == conv for _, b in done)
    ]

    arrays = {"t": t, "forced_zeros": FORCED_ZEROS}
    for name, kpoints in paths.items():
        arrays[f"k_{name}"] = kpoints
    for key, values in curves.items():
        arrays["__".join(key)] = values
    np.savez_compressed(output_dir / "symmetry_curves.npz", **arrays)

    metrics_path = output_dir / "symmetry_metrics.csv"
    with metrics_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "trial", "structure", "convention", "diagnostic", "t",
                "value_A2", "abs_value_A2", "n_R_H", "n_R_A",
                "same_R_set", "zero_H_blocks_added", "position_source", "ws_rule",
            ]
        )
        for trial, convention in done:
            for structure in ("base", "mode"):
                prov = provenance[(trial, structure, convention)]
                for path_name in paths:
                    values = curves[(trial, structure, convention, path_name)]
                    idx = int(np.argmax(np.abs(values)))
                    writer.writerow(
                        [trial, structure, convention, f"max_{path_name}",
                         f"{t[idx]:.8f}", f"{values[idx]:.16e}",
                         f"{abs(values[idx]):.16e}", prov["n_R_H"],
                         prov["n_R_A"], int(prov["same_R_set"]),
                         prov["zero_H_blocks_added"],
                         prov["position_source"], prov["ws_rule"]]
                    )
                if structure == "mode":
                    generic = curves[(trial, structure, convention, "generic")]
                    for zero in FORCED_ZEROS:
                        idx = int(round(zero * (len(t) - 1)))
                        if not np.isclose(t[idx], zero, atol=1e-12):
                            raise ValueError(
                                f"n_points={len(t)} does not sample t={zero}"
                            )
                        writer.writerow(
                            [trial, structure, convention, "generic_forced_zero",
                             f"{zero:.8f}", f"{generic[idx]:.16e}",
                             f"{abs(generic[idx]):.16e}", prov["n_R_H"],
                             prov["n_R_A"], int(prov["same_R_set"]),
                             prov["zero_H_blocks_added"],
                             prov["position_source"], prov["ws_rule"]]
                        )

    fig = plot_conventions(t, curves, trials_done, conventions_done, trial_set, convention_set)
    fig.savefig(output_dir / "berry_curvature_convention_comparison.png", bbox_inches="tight")
    plt.close(fig)

    print(f"saved {output_dir / 'symmetry_curves.npz'}")
    print(f"saved {metrics_path}")
    print(f"saved {output_dir / 'berry_curvature_convention_comparison.png'}")


def plot_conventions(t, curves, trials_done, conventions_done,
                      trial_set: dict = TRIALS, convention_set: dict = CONVENTIONS):
    """Build (but do not save) the per-trial, per-path comparison figure."""
    n_rows = len(trials_done)
    fig, axes = plt.subplots(
        n_rows, 3, figsize=(13.5, 3.25 * n_rows), dpi=180, sharex="col",
        squeeze=False,
    )
    columns = (
        ("base", "generic", "Undisplaced: generic path (P+T zero)"),
        ("mode", "mirror_kx_half", "Displaced: kx=1/2 mirror plane"),
        ("mode", "generic", "Displaced: generic path"),
    )
    for row, trial in enumerate(trials_done):
        setup = trial_set[trial]
        for col, (structure, path_name, title) in enumerate(columns):
            ax = axes[row, col]
            ax.axhline(0.0, color="0.35", linewidth=0.7)
            for convention in conventions_done:
                style = convention_set[convention]
                key = (trial, structure, convention, path_name)
                if key not in curves:
                    continue
                values = curves[key]
                ax.plot(
                    t,
                    values,
                    color=style["color"],
                    linestyle=style["linestyle"],
                    linewidth=1.35,
                    label=style["label"],
                )
            if col < 2:
                ax.set_yscale("symlog", linthresh=1e-9, linscale=0.8)
            else:
                for zero in FORCED_ZEROS:
                    ax.axvline(zero, color="0.8", linestyle=":", linewidth=0.7)
            if row == 0:
                ax.set_title(title, fontsize=10)
            if col == 0:
                ax.set_ylabel(setup["label"] + "\nTr $\\Omega_{xy}$ [$\\AA^2$]")
            if row == n_rows - 1:
                ax.set_xlabel("path parameter $t$")
            ax.grid(alpha=0.15)
    axes[0, 2].legend(frameon=False, fontsize=8)
    fig.suptitle("SrTiO$_3$: consistent interpolation-convention comparison", y=0.995)
    fig.tight_layout()
    return fig


def load_results(output_dir: Path = OUTPUT_DIR):
    """Reload a previous run's curves; provenance/per-point metrics stay in the CSV only."""
    data = np.load(output_dir / "symmetry_curves.npz")
    t = data["t"]
    paths = {name[2:]: data[name] for name in data.files if name.startswith("k_")}
    curves = {
        tuple(name.split("__")): data[name]
        for name in data.files
        if name not in ("t", "forced_zeros") and not name.startswith("k_")
    }
    return t, paths, curves
