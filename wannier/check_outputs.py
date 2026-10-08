#!/usr/bin/env python3
"""Stage 1 output checks: did Wannier90 write what stage 2 expects?

    python wannier/check_outputs.py RUN_DIR [MODE_DIR] [--prefix SEED] [--quick]

One directory checks that structure.  Two directories (base, mode) also check that
they are compatible and report the base/mode Wigner-Seitz consistency.  ``--quick``
reads only file headers (seconds); without it the files are loaded through PythTB
(a minute or two, several GB of memory for 176 Wannier functions).

Exit status is non-zero if any check fails.  The pair consistency number is
information, not a pass/fail check (see notes/progress/2026-09-20_jaemo_ndegen_yio_mode1.md).
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FAILURES: list[str] = []


def report(ok: bool, message: str) -> bool:
    print(f"  {'PASS' if ok else 'FAIL'}  {message}")
    if not ok:
        FAILURES.append(message)
    return ok


def find_prefix(directory: Path, prefix: str | None) -> str:
    if prefix:
        return prefix
    wins = sorted(directory.glob("*.win"))
    if len(wins) != 1:
        raise SystemExit(f"{directory}: need exactly one *.win (or pass --prefix)")
    return wins[0].stem


def plain_or_gz(path: Path) -> Path:
    return path if path.is_file() else Path(f"{path}.gz")


def read_header_line(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        return handle.readline()


def read_hr_degeneracies(path: Path) -> tuple[int, int, list[int]]:
    """num_wann, nrpts and the ndegen block of seed_hr.dat (header only)."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        handle.readline()
        num_wann, nrpts = int(handle.readline()), int(handle.readline())
        values: list[int] = []
        for _ in range(-(-nrpts // 15)):
            values += [int(x) for x in handle.readline().split()]
    return num_wann, nrpts, values


def win_keyword(text: str, key: str) -> str | None:
    match = re.search(rf"^\s*{key}\s*[=:]\s*(\S+)", text, re.I | re.M)
    return match.group(1).lower().strip(".") if match else None


def check_files(directory: Path, prefix: str, quick: bool) -> dict:
    print(f"\n{directory}  (seed {prefix})")
    info: dict = {}
    for suffix in ("_hr.dat", "_r.dat", "_wsvec.dat", "_centres.xyz"):
        path = plain_or_gz(directory / f"{prefix}{suffix}")
        report(path.is_file(), f"{prefix}{suffix} present")

    header = read_header_line(plain_or_gz(directory / f"{prefix}_wsvec.dat"))
    report("use_ws_distance=.true." in header, "_wsvec.dat header: use_ws_distance=.true.")
    report("write_ndegen_applied=.true." in header,
           "_wsvec.dat header: write_ndegen_applied=.true.")

    win = (directory / f"{prefix}.win").read_text(errors="replace")
    report(win_keyword(win, "transl_inv_full") in ("true", "t"), ".win: transl_inv_full = true")
    report(win_keyword(win, "use_ws_distance") in ("true", "t"), ".win: use_ws_distance = true")
    match = re.search(r"^\s*mp_grid\s*[=:]?\s*(\d+)\s+(\d+)\s+(\d+)", win, re.I | re.M)
    info["mp_grid"] = tuple(int(x) for x in match.groups()) if match else None
    report(info["mp_grid"] is not None, f".win: mp_grid = {info['mp_grid']}")

    num_wann, nrpts, degeneracies = read_hr_degeneracies(plain_or_gz(directory / f"{prefix}_hr.dat"))
    info["num_wann"], info["nrpts"] = num_wann, nrpts
    report(set(degeneracies) == {1}, f"_hr.dat degeneracy block all ones ({nrpts} R vectors)")
    return info


def check_loaded(directory: Path, prefix: str):
    """Load through PythTB; check R=0 diagonal of A(R) vs the centres, and hermiticity."""
    from pythtb import W90

    w90 = W90(str(directory), prefix)
    centres = np.asarray(w90.lattice.orb_vecs, float) @ np.asarray(w90.lat, float)
    diagonal = np.real(np.einsum("ass->sa", np.asarray(w90.pos_r[(0, 0, 0)])))
    error = float(np.abs(diagonal - centres).max())
    report(error < 1e-5, f"A(R=0) diagonal vs Wannier centres: {error:.1e} A (limit 1e-5; print floor 5e-7)")

    defect = norm = 0.0
    for R, block in w90.pos_r.items():
        partner = w90.pos_r.get(tuple(-x for x in R))
        if partner is None:
            defect = np.inf
            break
        difference = block - np.conj(np.transpose(partner, (0, 2, 1)))
        defect += float(np.vdot(difference, difference).real)
        norm += float(np.vdot(block, block).real)
    relative = float(np.sqrt(defect / norm)) if norm else np.nan
    report(relative < 1e-8, f"A(R) = A(-R)^dagger, relative defect {relative:.1e}")
    return w90


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directories", nargs="+", type=Path, help="one directory, or base and mode")
    parser.add_argument("--prefix", help="seedname (default: the only .win)")
    parser.add_argument("--quick", action="store_true", help="headers only; do not load through PythTB")
    args = parser.parse_args()
    if len(args.directories) > 2:
        parser.error("give one directory, or a base and a mode directory")

    infos, w90s = [], []
    for directory in args.directories:
        prefix = find_prefix(directory, args.prefix)
        infos.append(check_files(directory, prefix, args.quick))
        if not args.quick:
            w90s.append(check_loaded(directory, prefix))

    if len(args.directories) == 2:
        print("\nbase / mode pair")
        base, mode = infos
        report(base["mp_grid"] == mode["mp_grid"], f"same mp_grid ({base['mp_grid']} vs {mode['mp_grid']})")
        report(base["num_wann"] == mode["num_wann"], f"same num_wann ({base['num_wann']} vs {mode['num_wann']})")
        if w90s:
            from wannier.position_matrix import ws_pair_consistency

            H = [{R: b["h"] for R, b in w.ham_r.items()} for w in w90s]  # deg == 1 (checked above)
            pair = ws_pair_consistency(H[0], H[1], base["mp_grid"])
            print(
                f"  info  base/mode Wigner-Seitz images: {pair['entries_mismatched']}/"
                f"{pair['entries_compared']} entries changed support; the change of "
                f"interpolant is {pair['delta_H_contamination']:.2%} of |H_mode - H_base|"
            )

    print("\n" + ("ALL CHECKS PASSED" if not FAILURES else f"{len(FAILURES)} CHECK(S) FAILED"))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    raise SystemExit(main())
