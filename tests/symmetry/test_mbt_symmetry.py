"""MnBi2Te4 production models obey their magnetic point groups (paper tab:symmetry-mbt): symmetry_lib.failures() is
empty for the base (-3m' with the anti-translation partners: 24 operations) and each displaced structure (-3'm',
12), at the production window. Each group is read from the structure's pw.x output.

For each structure, in this order:
1. the saved results of a cluster run (run_mbt_symmetry.sh), copied to SAVED/<structure>/ (symmetry_base.npz,
   symmetry_mode.npz);
2. otherwise computed here from the local Wannier runs, if both .chk are in data/ (base and pattern 3 are: ~2.5 min,
   4 GB);
3. otherwise skipped, with the reason.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import symmetry_lib  # noqa: E402

# ---- settings
WINDOW = Path("/Users/treycole/Repos/axion-phonon/data/MnBi2Te4/windows/f10.962_w12.810_k8")
STRUCTURES = ["T2-_p1_Q0.01", "T2-_p2_Q0.01", "T2-_p3_Q0.01", "T2-_p4_Q0.01"]
SAVED = WINDOW            # cluster results: SAVED / symmetry_<structure> / out / symmetry_{base,mode}.npz
BASE_RUN = WINDOW / "base/out"
DFT = Path("/Users/treycole/Repos/axion-phonon/data/MnBi2Te4/dft")   # DFT / <structure> / scf/out/<PREFIX>.scf.out
ANTI_TRANSLATION = True   # the base also has {1'|1/2,1/2,1/2}, which QE does not list (symmetry_lib docstring)
PREFIX = "MnBi2Te4"
N_OCC = 82
FROZEN_TOP = 10.962       # eV, dis_froz_max of the window
N_GENERIC = 7
N_TUBE = 5
TUBE_RADIUS = 0.02        # 1/A
SEED = 20260825
BATCH = 4
# --------------


def has_chk(run: Path) -> bool:
    return (run / f"{PREFIX}.chk").exists() or (run / f"{PREFIX}.chk.gz").exists()


@pytest.fixture(scope="module", params=STRUCTURES)
def results(request):
    s = request.param
    saved = SAVED / f"symmetry_{s}" / "out"
    if (saved / "symmetry_mode.npz").exists():
        return s, symmetry_lib.load(saved)
    mode_run = WINDOW / s / "out"
    if not (has_chk(BASE_RUN) and has_chk(mode_run)):
        pytest.skip(f"{s}: no saved results in {saved} and no local .chk in {mode_run} (run run_mbt_symmetry.sh)")
    scf = (DFT / "base/scf/out" / f"{PREFIX}.scf.out", DFT / s / "scf/out" / f"{PREFIX}.scf.out", ANTI_TRANSLATION)
    settings = (N_OCC, FROZEN_TOP, N_GENERIC, N_TUBE, TUBE_RADIUS, SEED, BATCH)
    pair = request.getfixturevalue("mbt_pair")
    if pair["base_run"] == BASE_RUN and pair["mode_run"] == mode_run:          # the session's pair (conftest.py)
        return s, symmetry_lib.pair_results(pair["base"], pair["mode"], pair["base_model"], pair["mode_model"],
                                            BASE_RUN, mode_run, *scf, *settings, log=lambda line: None)
    return s, symmetry_lib.run_pair(BASE_RUN, mode_run, *scf, PREFIX, *settings, log=lambda line: None)


@pytest.mark.data
@pytest.mark.parametrize("label", ["base", "mode"])
def test_symmetry(results, label):
    s, r = results
    bad = symmetry_lib.failures(r[label])
    assert not bad, f"{s}, {label}:\n  " + "\n  ".join(bad)
