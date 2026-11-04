"""Puts the repository root on sys.path (for axion_phonon), registers the `data` marker (tests that read the runs in
data/, git-ignored; they skip, naming the missing path, where a run is not there; `pytest tests -m "not data"` runs
only the unit tests), and builds the production models the data tests share, once per session.

Fixtures (each a dict: the read_run dict of each run without its overlaps M, and the models):
  mbt_pair   MnBi2Te4 base and pattern 3 (Q = 0.01 A), production window: pair_models (~2 min, 4 GB)
  yio_base   Y2Ir2O7 base, production window (wannier90 220455): run_model with A(R), pinned (the base model of pair_models)
"""
import sys
from pathlib import Path

import pytest

# ---- settings
DATA = Path("/Users/treycole/Repos/axion-phonon/data")
MBT_WINDOW = DATA / "MnBi2Te4/windows/f10.962_w12.810_k8"
MBT_BASE_RUN = MBT_WINDOW / "base/out"
MBT_MODE_RUN = MBT_WINDOW / "T2-_p3_Q0.01/out"
MBT_N_OCC = 82
YIO_WINDOW = DATA / "Y2Ir2O7/windows/f12.422_w22.644_k8"
YIO_BASE_RUN = YIO_WINDOW / "base/out"
YIO_N_OCC = 156
# --------------

REPOSITORY = Path(__file__).resolve().parent.parent
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))


def pytest_configure(config):
    config.addinivalue_line("markers", "data: reads the runs in data/; skips where a run is not there")


def has_run(folder: Path, prefix: str) -> bool:
    return any((folder / f"{prefix}.{e}").exists() for e in ("chk", "chk.gz")) and (folder / f"{prefix}.eig").exists()


@pytest.fixture(scope="session")
def mbt_pair():
    from axion_phonon import wannier
    from axion_phonon.io import wannier90
    for run in (MBT_BASE_RUN, MBT_MODE_RUN):
        if not has_run(run, "MnBi2Te4"):
            pytest.skip(f"not in data/: {run}")
    base, mode = wannier90.read_run(MBT_BASE_RUN, "MnBi2Te4"), wannier90.read_run(MBT_MODE_RUN, "MnBi2Te4")
    mb, mm = wannier.pair_models(base, mode, log=lambda s: None)     # pops M from both runs
    return dict(base=base, mode=mode, base_model=mb, mode_model=mm, n_occ=MBT_N_OCC, dbeta=0.01,
                base_run=MBT_BASE_RUN, mode_run=MBT_MODE_RUN)


@pytest.fixture(scope="session")
def yio_base():
    from axion_phonon import wannier
    from axion_phonon.io import wannier90
    if not has_run(YIO_BASE_RUN, "Y2Ir2O7"):
        pytest.skip(f"not in data/: {YIO_BASE_RUN}")
    run = wannier90.read_run(YIO_BASE_RUN, "Y2Ir2O7")
    model = wannier.run_model(run, pin=True, position_matrix=True)
    return dict(run=run, model=model, n_occ=YIO_N_OCC, run_dir=YIO_BASE_RUN)
