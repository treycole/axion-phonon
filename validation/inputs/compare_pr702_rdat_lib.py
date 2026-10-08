"""Check that Wannier90 PR #702's ``_r.dat`` is postw90's ``A(R)``, rounded.

Two routines now claim to build the translation-equivariant position matrix.

``postw90 _r_full.dat``
    ``get_AA_R`` with ``transl_inv_full=T``, exported by the local ``write_aa_r``
    patch.  Each neighbour link is multiplied by ``exp(i b.(r_n + r_m)/2)``
    before the Fourier transform and by ``exp(-i b.R/2)`` after it, and the
    links are kept separate throughout.  Written at ``ES24.16E3``.

``wannier90.x _r.dat``
    ``plot_write_rmn`` with ``transl_inv_full=T``, added by upstream PR #702.
    Both phases are combined into one ``exp(i b.(r_n + r_m - R)/2)`` applied
    inside the k-sum, using each k-point's own ``b`` rather than a canonical
    neighbour ordering.  Written at ``6F12.6``, unchanged from before.

If the two are the same matrix, then every component of ``_r.dat`` must be
*exactly* the six-decimal rendering of the corresponding ``_r_full.dat`` value.
That is the test: not "agrees to tolerance" but bit-exact after rounding, which
distinguishes an algebraic identity from two formulas that merely happen to be
close. Also reports what the ``6F12.6`` format costs, since PR #702 does not
change it -- see notes/progress/2026-08-27_rfull_vs_mmn_position_matrix.md section 0.

Pure library: no argparse, no ``__main__``. Used directly by
``compare_pr702_rdat.ipynb`` and imported by ``inputs/compare_jaemo_ndegen_rdat.ipynb``
(``read_matrix``, ``fortran_round``).
"""

from __future__ import annotations

import csv
import sys
from itertools import islice
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import RESULTS, STO_TI_Q_2A_RUN  # noqa: E402

#: 6F12.6 resolves 1e-6; Fortran rounds half away from zero.
PRINT_STEP = 1e-6

#: Colours follow the entity across every PR #702 figure (dataviz reference
#: palette, slots 1-5, adjacent pairs validated); see
#: validation/symmetry/check_pr702_rdat_symmetry.py.
SERIES = {
    "stock": ("#1baf7a", "stock _r.dat"),
    "PR702": ("#eda100", "PR #702 _r.dat, 6F12.6"),
    "PR702hp": ("#e87ba4", "PR #702 _r.dat, ES24.16E3"),
}
INK, INK_SECONDARY, MUTED = "#0b0b0b", "#52514e", "#898781"
SURFACE, GRID, AXIS = "#fcfcfb", "#e1e0d9", "#c3c2b7"

#: Differences below this are drawn at the floor, exact zeros included.
FLOOR = 1e-17

#: A full-precision _r.dat must match postw90 to summation-order roundoff.
FULL_PRECISION_TOL = 1e-12

#: Default candidate/reference pair used when no explicit case list is given.
DEFAULT_RUN_DIR = STO_TI_Q_2A_RUN / "trial_03_Sr_sp_Ti_pd_O_sp/proj_gauge"
DEFAULT_CANDIDATE = DEFAULT_RUN_DIR / "SrTiO3_r.dat"
DEFAULT_REFERENCE = DEFAULT_RUN_DIR / "SrTiO3_r_full.dat"


def read_matrix(path: Path, chunk: int = 200_000) -> dict:
    """Stream a ``_r.dat``/``_r_full.dat`` into index and value arrays."""
    idx_parts: list[np.ndarray] = []
    val_parts: list[np.ndarray] = []
    with open(path) as handle:
        header = handle.readline().rstrip()
        num_wann = int(handle.readline())
        nrpts = int(handle.readline())
        while True:
            lines = list(islice(handle, chunk))
            if not lines:
                break
            block = np.array("".join(lines).split(), dtype=np.float64)
            block = block.reshape(len(lines), 11)
            idx_parts.append(block[:, :5].astype(np.int64))
            val_parts.append(block[:, 5:])
    idx = np.concatenate(idx_parts)
    val = np.concatenate(val_parts)
    expected = nrpts * num_wann * num_wann
    if len(idx) != expected:
        raise ValueError(f"{path}: {len(idx)} rows, header implies {expected}")
    return dict(path=path, header=header, num_wann=num_wann, nrpts=nrpts,
                idx=idx, val=val)


def fortran_round(values: np.ndarray, step: float = PRINT_STEP) -> np.ndarray:
    """Render as Fortran F12.6 would: round half *away from zero*."""
    return np.sign(values) * np.floor(np.abs(values) / step + 0.5) * step


def print_format(data: dict) -> str:
    """``6F12.6`` if every value sits on the 1e-6 grid, else full precision."""
    on_grid = np.abs(data["val"] - fortran_round(data["val"])) < 1e-12
    return "6F12.6" if on_grid.all() else "ES24.16E3"


def hermiticity(data: dict) -> float:
    """max |A_nm(R) - conj(A_mn(-R))| / max|A|, the relative residual."""
    idx, val = data["idx"], data["val"]
    lookup = {tuple(row): i for i, row in enumerate(idx)}
    partner = np.empty(len(idx), dtype=np.int64)
    for i, row in enumerate(idx):
        key = (-row[0], -row[1], -row[2], row[4], row[3])
        partner[i] = lookup.get(key, -1)
    if (partner < 0).any():
        missing = int((partner < 0).sum())
        raise ValueError(f"{data['path']}: {missing} rows have no -R partner")
    other = val[partner]
    residual = np.abs(val[:, 0::2] - other[:, 0::2]).max()
    residual = max(residual, np.abs(val[:, 1::2] + other[:, 1::2]).max())
    return residual / np.abs(val).max()


def compare(candidate: dict, reference: dict) -> dict:
    if not np.array_equal(candidate["idx"], reference["idx"]):
        raise ValueError(
            "row (R, n, m) layout differs -- postw90 reorders its R list when "
            "use_ws_distance=T, so build the reference with use_ws_distance=F"
        )
    cand, ref = candidate["val"], reference["val"]
    rounded = fortran_round(ref)
    exact = np.abs(cand - rounded) < 1e-12
    diff = np.abs(cand - ref)
    big = np.abs(ref) > 1e-5
    killed = (ref != 0) & (rounded == 0)
    return dict(
        candidate_format=print_format(candidate),
        components=int(cand.size),
        bit_exact=int(exact.sum()),
        bit_exact_fraction=float(exact.mean()),
        max_residual_after_rounding=float(np.abs(cand - rounded).max()),
        max_abs_diff=float(diff.max()),
        components_above_1e5=int(big.sum()),
        mean_relative_diff_above_1e5=float((diff[big] / np.abs(ref[big])).mean()) if big.any() else 0.0,
        max_relative_diff_above_1e5=float((diff[big] / np.abs(ref[big])).max()) if big.any() else 0.0,
        rounded_to_zero_fraction=float(killed.mean()),
        hermiticity_candidate=hermiticity(candidate),
        hermiticity_reference=hermiticity(reference),
    )


def passes(result: dict) -> bool:
    """Bit-exact after rounding for a 6F12.6 file; machine epsilon otherwise."""
    if result["candidate_format"] == "6F12.6":
        return result["bit_exact"] == result["components"]
    return result["max_abs_diff"] < FULL_PRECISION_TOL


def report(result: dict, candidate: dict, reference: dict) -> bool:
    passed = passes(result)
    print(f"candidate  {candidate['path']}")
    print(f"           {candidate['header']}")
    print(f"reference  {reference['path']}")
    print(f"           {reference['header']}")
    print(f"\n  num_wann={reference['num_wann']}  nrpts={reference['nrpts']}  "
          f"components={result['components']}")
    print(f"  bit-exact after 6-decimal rounding  {result['bit_exact']} / "
          f"{result['components']}  ({result['bit_exact_fraction']:.6%})")
    print(f"  max residual after rounding         {result['max_residual_after_rounding']:.3e}")
    print(f"  max |candidate - reference|         {result['max_abs_diff']:.3e}")
    print(f"  hermiticity, candidate / reference  {result['hermiticity_candidate']:.3e}"
          f" / {result['hermiticity_reference']:.3e}")
    print("\n  what the 6F12.6 format costs:")
    print(f"    reference values rounded to zero  {result['rounded_to_zero_fraction']:.2%}")
    print(f"    components with |ref| > 1e-5      {result['components_above_1e5']}")
    print(f"    max relative error there          {result['max_relative_diff_above_1e5']:.3e}")
    if passed and result["candidate_format"] != "6F12.6":
        print("\n  -> IDENTICAL: PR #702 and postw90 agree to machine precision; "
              "they are the same computation.")
    elif passed:
        print("\n  -> IDENTICAL: PR #702 writes postw90's A(R) to the precision "
              "6F12.6 can carry.")
    else:
        print("\n  -> DIFFERENT: the two constructions are not the same matrix.")
    return passed


# --------------------------------------------------------------------------- #
# Several systems, one table and one figure
# --------------------------------------------------------------------------- #


def per_R_max(data: dict, diff: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """|R| in lattice units, and max |diff| over the (n, m, x) block at that R."""
    unique, inverse = np.unique(data["idx"][:, :3], axis=0, return_inverse=True)
    worst = np.zeros(len(unique))
    np.maximum.at(worst, inverse.ravel(), diff.max(axis=1))
    return np.linalg.norm(unique, axis=1), worst


def run_case(label, candidate_path, reference_path, control_path=None) -> dict:
    reference = read_matrix(reference_path)
    case = dict(label=label, reference=reference, rows=[])
    entries = [("candidate", candidate_path)]
    if control_path is not None:
        entries.append(("control", control_path))
    for role, path in entries:
        data = read_matrix(path)
        series = "stock" if role == "control" else (
            "PR702" if print_format(data) == "6F12.6" else "PR702hp")
        print(f"\n=== {label}: {SERIES[series][1]}\n")
        result = compare(data, reference)
        passed = report(result, data, reference)
        diff = np.abs(data["val"] - reference["val"])
        case["rows"].append(dict(
            role=role, series=series, passed=passed, result=result,
            diff=diff, per_R=per_R_max(data, diff),
            num_wann=reference["num_wann"], nrpts=reference["nrpts"],
        ))
    return case


def _rounded_only(fmt):
    """Bit-exactness after rounding only means something for a 6F12.6 file."""
    return lambda r: fmt(r) if r["candidate_format"] == "6F12.6" else "n/a (full precision)"


TABLE_COLUMNS = [
    ("bit-exact vs 6F12.6(ref)", _rounded_only(lambda r: f"{r['bit_exact_fraction']:.6%}")),
    ("max abs diff", lambda r: f"{r['max_abs_diff']:.2e}"),
    ("max resid. after rounding",
     _rounded_only(lambda r: f"{r['max_residual_after_rounding']:.1e}")),
    ("mean rel. diff, |A|>1e-5", lambda r: f"{r['mean_relative_diff_above_1e5']:.2e}"),
    ("max rel. diff, |A|>1e-5", lambda r: f"{r['max_relative_diff_above_1e5']:.2e}"),
    ("Hermiticity (file)", lambda r: f"{r['hermiticity_candidate']:.2e}"),
]


def write_tables(output_dir: Path, cases: list[dict]) -> list[Path]:
    csv_path = output_dir / "file_level_summary.csv"
    keys = list(cases[0]["rows"][0]["result"])
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["system", "file", "num_wann", "nrpts", *keys])
        for case in cases:
            for row in case["rows"]:
                writer.writerow([case["label"], SERIES[row["series"]][1],
                                 row["num_wann"], row["nrpts"],
                                 *(row["result"][k] for k in keys)])

    md_path = output_dir / "file_level_summary.md"
    lines = [
        "`_r.dat` vs postw90 `write_aa_r` `_r_full.dat` "
        "(`transl_inv_full=T`, `use_ws_distance=F`), matched row-for-row",
        "",
        "| system | file | WF | nrpts | components | "
        + " | ".join(name for name, _ in TABLE_COLUMNS) + " |",
        "|---|---|---|---|---|" + "---|" * len(TABLE_COLUMNS),
    ]
    for case in cases:
        for row in case["rows"]:
            r = row["result"]
            lines.append(
                f"| {case['label']} | {SERIES[row['series']][1]} | {row['num_wann']} "
                f"| {row['nrpts']} | {r['components']:,} | "
                + " | ".join(fmt(r) for _, fmt in TABLE_COLUMNS) + " |")
    reference = cases[0]["rows"][0]["result"]
    lines += ["", f"Reference Hermiticity: {reference['hermiticity_reference']:.1e}.  "
              "Rounded-to-zero fraction of the reference under 6F12.6: "
              + ", ".join(f"{c['label']} {c['rows'][0]['result']['rounded_to_zero_fraction']:.1%}"
                          for c in cases) + "."]
    md_path.write_text("\n".join(lines) + "\n")
    return [csv_path, md_path]


def style_axis(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID, linewidth=0.6, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=8)


def plot(cases: list[dict]):
    """Build (but do not save) the residual-distribution figure."""
    plt.rcParams.update({"font.family": "sans-serif", "text.color": INK,
                         "axes.labelcolor": INK_SECONDARY, "axes.titlecolor": INK})
    grid = np.logspace(np.log10(FLOOR), -1, 400)
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 2.7 * len(cases)), dpi=180,
                             squeeze=False, facecolor=SURFACE,
                             gridspec_kw=dict(width_ratios=(1.25, 1)))
    for row_index, case in enumerate(cases):
        tail_ax, r_ax = axes[row_index]
        for row in case["rows"]:
            color, label = SERIES[row["series"]]
            ordered = np.sort(row["diff"].ravel())
            above = 1.0 - np.searchsorted(ordered, grid, side="right") / ordered.size
            tail_ax.plot(grid, np.where(above > 0, above, np.nan), color=color,
                         linewidth=2, label=label, solid_capstyle="round")
            norm, worst = row["per_R"]
            r_ax.plot(norm, np.clip(worst, FLOOR, None), "o", markersize=4.5,
                      markerfacecolor=color, markeredgecolor=SURFACE,
                      markeredgewidth=1.0, label=label, alpha=0.9)
        for ax in (tail_ax, r_ax):
            style_axis(ax)
            ax.set_yscale("log")
        tail_ax.set_xscale("log")
        tail_ax.axvline(PRINT_STEP / 2, color=MUTED, linewidth=0.8)
        if row_index == 0:
            tail_ax.annotate("half the 6F12.6 step", (PRINT_STEP / 2, 0.04),
                             xycoords=("data", "axes fraction"), xytext=(4, 0),
                             textcoords="offset points", ha="left", va="bottom",
                             fontsize=7.5, color=INK_SECONDARY)
        tail_ax.set_xlim(FLOOR, 1e-1)
        tail_ax.set_ylim(1e-7, 2)
        r_ax.set_ylim(FLOOR * 0.5, 1e-1)
        tail_ax.set_ylabel(r"fraction with $|\Delta A| > x$")
        r_ax.set_ylabel(r"max $|\Delta A|$ in block [$\AA$]")
        first = case["rows"][0]["result"]
        tail_ax.set_title(
            f"{case['label']}  ({case['rows'][0]['num_wann']} WF, "
            f"{case['rows'][0]['nrpts']} R, {first['components']:,} components)",
            loc="left", fontsize=9.5)
        r_ax.set_title("worst component at each R", loc="left", fontsize=9.5)
        if row_index == len(cases) - 1:
            tail_ax.set_xlabel(r"$x$,  $|A_{\rm file} - A_{\rm postw90}|$ per component [$\AA$]")
            r_ax.set_xlabel("|R| [lattice units]")
        tail_ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.suptitle("Wannier90 PR #702 _r.dat vs our postw90 write_aa_r export "
                 "(same checkpoint, transl_inv_full=T, use_ws_distance=F)",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    return fig


def save_outputs(output_dir: Path, cases: list[dict], make_plot: bool = True) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    saved = write_tables(output_dir, cases)
    if make_plot:
        fig = plot(cases)
        png_path = output_dir / "pr702_file_level.png"
        fig.savefig(png_path, bbox_inches="tight", facecolor=SURFACE)
        saved.append(png_path)
    return saved
