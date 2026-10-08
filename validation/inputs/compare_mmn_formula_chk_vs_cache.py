#!/usr/bin/env python3
"""Does the ``.mmn`` A(R) rebuilt from the ``.chk`` ``m_matrix`` reproduce the cache?

``wannier.position_matrix.position_matrix_mmn_formula`` applies the production
``.mmn`` formula to the checkpoint's rotated overlaps instead of the ``.mmn``
file.  The cache in ``_AA_cache/`` was built from the ``.mmn`` itself (since
deleted).  Agreement to the ``.mmn``'s 12-digit text precision means the two
are the same input, so the two position-matrix formulas can be compared on
one ``.chk`` with nothing else different.

    python compare_mmn_formula_chk_vs_cache.py <run_dir> <prefix> <cache.npz> [--out results.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from wannier.position_matrix import position_matrix_mmn_formula  # noqa: E402
from wannier.wannier_io import read_wannier_checkpoint  # noqa: E402


def compare(run_dir: Path, prefix: str, cache: Path) -> dict:
    checkpoint = read_wannier_checkpoint(run_dir / f"{prefix}.chk", include_m_matrix=True)
    checkpoint.pop("v_matrix", None)
    ours = position_matrix_mmn_formula(checkpoint, run_dir / f"{prefix}.nnkp")
    del checkpoint
    with np.load(cache, allow_pickle=False) as archive:
        reference_R = [tuple(int(v) for v in R) for R in archive["iRvec"]]
        reference = archive["AA"]
        reference_centres = archive["wcc_cart"]
    index = {tuple(int(v) for v in R): i for i, R in enumerate(ours["iRvec"])}
    diff_sq = norm_sq = worst = 0.0
    for i, R in enumerate(reference_R):
        mine = ours["AA"][index[R]] if R in index else 0.0
        difference = np.abs(mine - reference[i])
        worst = max(worst, float(difference.max()))
        diff_sq += float((difference ** 2).sum())
        norm_sq += float((np.abs(reference[i]) ** 2).sum())
    return {
        "n_R_ours": len(index),
        "n_R_cache": len(reference_R),
        "R_in_cache_not_ours": sum(R not in index for R in reference_R),
        "max_abs_diff_A": worst,
        "relative_l2": (diff_sq / norm_sq) ** 0.5,
        "centres_max_diff_A": float(np.abs(ours["wcc_cart"] - reference_centres).max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("prefix")
    parser.add_argument("cache", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = compare(args.run_dir, args.prefix, args.cache)
    print(json.dumps(result, indent=2))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
