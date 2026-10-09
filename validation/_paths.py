"""Every data location the validation scripts read, named once.

The first-principles runs live under ``data/`` (``DATA`` below) in folders named
by cluster job ID.  Their full paths are written out here, once, so a check
imports a constant instead of spelling the path.  To repoint a check at a new
run, edit the path here (or override it in the check's config cell).

Importing this module also puts the repository root on ``sys.path``, which is
what makes ``from modules.curvature import ...`` work when a script is
run directly from any directory.

Layout
------
``REPOSITORY``  repository root
``RESULTS``     ``validation/results`` -- every check writes its summary here
``STO_*``       SrTiO3: convention, position-matrix and Born-charge work
``YIO_*``       Y2Ir2O7: the production axion-response pair

Nothing here is guaranteed to exist: the large runs are not in version control.
Call :func:`require` to fail early with a readable message instead of a
``FileNotFoundError`` from somewhere deep in a reader.
"""

from __future__ import annotations

import sys
from pathlib import Path

VALIDATION = Path(__file__).resolve().parent
REPOSITORY = VALIDATION.parent

if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))


RESULTS = VALIDATION / "results"
DATA = REPOSITORY / "data"

# --------------------------------------------------------------------------- #
# SrTiO3
# --------------------------------------------------------------------------- #

STO = DATA / "SrTiO3"
STO_PREFIX = "SrTiO3"

#: Shared ``A(R)`` archive.  Replaces the ``.mmn`` as stored input; see
#: notes/log/2026-08-27_rfull_vs_mmn_position_matrix.md.
STO_AA_CACHE = STO / "_AA_cache"

#: Undisplaced cubic reference (also holds SrTiO3.ph.out, the DFPT Z* source).
STO_BASE_RUN = DATA / "SrTiO3/base/soc/output/181642bc889"

#: Ti displaced by +0.01 A along z -- the Born-effective-charge finite difference.
STO_TI_Q_01A_RUN = DATA / "SrTiO3/Ti_Q_01A/output/181839bc889"

#: Sr displaced by +0.01 A along z -- the second Born-charge finite difference.
STO_SR_Q_01A_RUN = DATA / "SrTiO3/Sr_Q_01A/output/184367bc889"

#: Ti displaced by 2 A -- the large-displacement trial-comparison study.
STO_TI_Q_2A = STO / "Ti_Q_2A"
STO_TI_Q_2A_RUN = DATA / "SrTiO3/Ti_Q_2A/output/188914bc889"

# --------------------------------------------------------------------------- #
# Y2Ir2O7 -- the production axion-response pair
# --------------------------------------------------------------------------- #

YIO = DATA / "Y2Ir2O7/phonon"
YIO_PREFIX = "Y2Ir2O7"
YIO_AA_CACHE = YIO / "_AA_cache"

#: trial_02 (Y:p, Ir:d, O:s;p), 176 Wannier functions, n_occ = 156, U = 3.0 eV.
#: This is the set the production sweep uses.  The legacy 196-orbital pair
#: (base 173260bc889, mode 173195bc889, n_occ = 104) is a different Wannier
#: space and does not test the models in use.
YIO_TRIAL = "trial_02_Y_p_Ir_d_O_sp"
YIO_N_OCCUPIED = 156

YIO_BASE = DATA / "Y2Ir2O7/base/soc/u_3.0/output/178382bc889/trial_02_Y_p_Ir_d_O_sp/proj_gauge/181538bc889"
YIO_MODE1 = DATA / "Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/178365bc889/trial_02_Y_p_Ir_d_O_sp/no_displacement/181565bc889"

#: Frozen-phonon amplitude separating the two Wannierizations, in Angstrom.
YIO_D_BETA = 0.01

YIO_AA_CACHE_FILES = {
    "base": DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_base_181538bc889.npz",
    "mode": DATA / "Y2Ir2O7/phonon/_AA_cache/AA_trial_02_Y_p_Ir_d_O_sp_mode_181565bc889.npz",
}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def require(path: Path, what: str) -> Path:
    """Return ``path``, or raise with the name of the missing input."""
    if not path.exists():
        raise FileNotFoundError(
            f"{what} not found at {path}\n"
            "The first-principles runs are not in version control.  Edit "
            f"{Path(__file__).relative_to(REPOSITORY)} to point at your copy."
        )
    return path


def result_dir(name: str) -> Path:
    """Return (and create) ``validation/results/<name>``."""
    directory = RESULTS / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory
