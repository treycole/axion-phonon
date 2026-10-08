#!/usr/bin/env python3
"""Does our Python ``transl_inv_full`` A(R) reproduce Jae-Mo's ``_r.dat`` 1:1?

``wannier.position_matrix.position_matrix_transl_inv_full`` rebuilds
``<0i|r|Rj>`` from the checkpoint's ``m_matrix`` with Wannier90's
``transl_inv_full`` + ``write_ndegen_applied`` formula, on the run's own
centres.  Jae-Mo's ``wannier90.x`` (``restart = plot`` on the same ``.chk``,
all three flags) wrote ``_r.dat`` from the same data.  The file is printed at
``F12.6``, so agreement is judged after rounding ours to six decimals.

    python compare_chk_transl_inv_full.py <run_dir> <prefix> [<r.dat dir>] [--out results.json]

``<r.dat dir>`` defaults to ``<run_dir>`` (SrTiO3 trial04 keeps both together);
for Y2Ir2O7 it is ``<run_dir>/jaemo``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from wannier.position_matrix import position_matrix_transl_inv_full  # noqa: E402
from wannier.wannier_io import read_wannier_checkpoint  # noqa: E402


def compare(run_dir: Path, prefix: str, rdat_dir: Path) -> dict:
    start = time.time()
    checkpoint = read_wannier_checkpoint(run_dir / f"{prefix}.chk", include_m_matrix=True)
    for key in ("v_matrix",):
        checkpoint.pop(key, None)  # not needed; frees memory
    ours = position_matrix_transl_inv_full(checkpoint, run_dir / f"{prefix}.nnkp")
    del checkpoint
    print(f"  built ours: {len(ours)} R vectors in {time.time() - start:.0f} s")

    keys = sorted(ours)
    index = {R: i for i, R in enumerate(keys)}
    stack = np.stack([ours[R] for R in keys])  # (nR, 3, J, J)
    del ours
    hermiticity = max(
        float(np.abs(stack[index[R]] - np.conj(stack[index[tuple(-v for v in R)]]).swapaxes(1, 2)).max())
        for R in keys
        if tuple(-v for v in R) in index
    )

    path = rdat_dir / f"{prefix}_r.dat"
    with open(path) as handle:
        handle.readline()
        num_wann = int(handle.readline())
        nrpts = int(handle.readline())
    stats = dict(rows=0, rows_R_missing_in_ours=0, max_abs_diff=0.0, exact_after_rounding=0,
                 components=0, diff_sq=0.0, norm_sq=0.0, nonzero_file_rows=0)
    file_R = set()
    reader = pd.read_csv(path, sep=r"\s+", header=None, skiprows=3, chunksize=4_000_000,
                         dtype=np.float64, engine="c")
    for chunk in reader:
        values = chunk.to_numpy()
        R = values[:, :3].astype(int)
        m = values[:, 3].astype(int) - 1
        n = values[:, 4].astype(int) - 1
        theirs = values[:, 5::2] + 1j * values[:, 6::2]  # (rows, 3)
        file_R.update(map(tuple, np.unique(R, axis=0)))
        rows = np.array([index.get(tuple(r), -1) for r in map(tuple, R)])
        present = rows >= 0
        stats["rows_R_missing_in_ours"] += int((~present).sum() and np.abs(theirs[~present]).max() > 0)
        mine = np.zeros_like(theirs)
        mine[present] = stack[rows[present], :, m[present], n[present]]
        rounded = np.round(mine.real, 6) + 1j * np.round(mine.imag, 6)
        diff = np.abs(rounded - theirs)
        stats["rows"] += len(values)
        stats["components"] += diff.size
        stats["exact_after_rounding"] += int((diff < 1.5e-6).sum())  # one print step of slack
        stats["max_abs_diff"] = max(stats["max_abs_diff"], float(np.abs(mine - theirs).max()))
        stats["diff_sq"] += float((np.abs(mine - theirs) ** 2).sum())
        stats["norm_sq"] += float((np.abs(theirs) ** 2).sum())
        stats["nonzero_file_rows"] += int((np.abs(theirs).max(axis=1) > 0).sum())
    ours_nonzero = {R for R in keys if np.abs(stack[index[R]]).max() > 5e-7}
    return {
        "num_wann": num_wann,
        "nrpts_file": nrpts,
        "n_R_ours": len(keys),
        "n_R_ours_above_print": len(ours_nonzero),
        "R_in_ours_above_print_but_not_in_file": len(ours_nonzero - file_R),
        "rows": stats["rows"],
        "exact_fraction": stats["exact_after_rounding"] / stats["components"],
        "max_abs_diff": stats["max_abs_diff"],
        "relative_l2": float(np.sqrt(stats["diff_sq"] / stats["norm_sq"])),
        "hermiticity_max_ours": hermiticity,
        "seconds": time.time() - start,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("prefix")
    parser.add_argument("rdat_dir", type=Path, nargs="?")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = compare(args.run_dir, args.prefix, args.rdat_dir or args.run_dir)
    print(json.dumps(result, indent=2))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
