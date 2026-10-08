r"""Ti Born effective charge in cubic SrTiO3: DFPT versus Wannier mixed curvature.

Compares the row
``(Z*_{Ti,zx}, Z*_{Ti,zy}, Z*_{Ti,zz})`` from Quantum ESPRESSO DFPT against a
one-sided finite difference between the undisplaced structure and ``Ti_Q_01A``,
where Ti alone moves by ``+0.01 A`` along z.

Method
------
With ``u = u_{Ti,z}`` the electronic response is the BZ integral of the trace of
the mixed ``(k, u)`` curvature (note "Wannier mixed-curvature result"),

    Z^el_{Ti,z beta} = (Omega / (2 pi)^3) Int d^3k_cart Tr B_{k_beta, u}.

``curvature.fields`` and ``TBModel.velocity`` both work in **reduced** reciprocal
axes, so the Cartesian prefactor must not be applied a second time.  Writing
``B = model.recip_lat_vecs`` (rows b_1,b_2,b_3, carrying the 2 pi) the mixed
curvature carries one lower k index, hence

    B^cart_{beta,u} = sum_a (B^-1)_{beta,a} B^red_{a,u},
    Int d^3k_cart   = |det B| Int_[0,1)^3 d^3k_red,    det B = (2 pi)^3 / Omega,

and the two volume factors cancel exactly:

    Z^el_{Ti,z beta} = sum_a (B^-1)_{beta,a} <Tr B^red_{a,u}>_BZ,

with ``<.>_BZ`` the plain mean over a Gamma-centred uniform mesh of [0,1)^3.
The curvature is split as ``B = B_int + B_cross + B_ext`` and each piece is
integrated separately; only the sum is physical (note "Wannier mixed-curvature
result").  Adding the ionic term to the zz component alone gives the total to
compare with ``ph.x``.

The ionic term is **not** simply ``z_valence``.  ``dis_win_min`` leaves the
Ti 3s and Ti 3p semicore bands outside the Wannier space.  They translate
rigidly with the displaced atom, so their response is folded into
``Z^ion,eff = z_valence - N_EXCLUDED_ON_DISPLACED``.  See that constant.

Sign convention: the code stores ``B_{k_beta, u} = -B_{u, k_beta}``, i.e.
direction index 3 is ``u`` and ``B[a, 3]`` is ``B_{k_a, u}``, which is the
convention the note's ``+`` prefactor assumes.  The physical check is that
``Z^el_{Ti,zz}`` comes out near
``Z*_QE - Z^ion,eff = 7.29 - 4 = +3.29``.

Pure library: no argparse, no ``__main__``. Every argparse option of the old
CLI has a same-named, same-default config variable in ``born_effective_charge.ipynb``;
this module holds only the parsers, the curvature math, and the driver helpers
they call. Also imported by ``validation/observable/audit_mixed_born_response.ipynb``
(``centre_cross_check``, ``parse_ph_born_charges``).

Requires the displaced Wannier90 run to be present at ``MODE_W90`` below (see
``preflight`` output for the exact file list).
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

from pythtb import W90

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation._paths import (
    STO,
    STO_AA_CACHE,
    STO_BASE_RUN,
    STO_SR_Q_01A_RUN,
    STO_TI_Q_01A_RUN,
)
from validation._paths import DATA
from modules import curvature as cv
from modules.axion import beta_terms, parametric_connection, prepare_lambda
from wannier.wannier_io import check_mmn_compatibility, load_Lambda_R
from wannier.position_matrix import (
    apply_orbital_dependent_R_mapping,
    build_position_matrix_from_mmn,
    install_postw90_position_matrix,
)


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

PREFIX = "SrTiO3"

# Undisplaced reference.  The W90 paths may point either at the directory
# holding PREFIX.win or at a parent of it: _resolve_w90_dir() descends into the
# newest subdirectory that has one, and prints what it picked.
# TRIAL selects the projection set.  The base and displaced runs must use the
# SAME trial set, otherwise the two Wannier spaces are not comparable and the
# finite difference measures the change of projection, not the displacement.
#
#   trial_01_Sr_spd_Ti_spd_O_sp : 60 WFs, spans occupied + empty (n_occ = 30)
#   trial_02_Sr_p_O_sp          : 30 WFs, spans the occupied manifold alone
#   trial_03_Sr_sp_Ti_pd_O_sp   : 48 WFs, wider window, only Ti 3s excluded
#
# With an occupied-only set there are no conduction states inside the model, so
# B_int and B_cross vanish identically and the entire response is the external
# (position-matrix) term.  That is not a special case in the code -- the empty
# conduction block makes both terms zero on their own -- but it does mean the
# answer stands or falls on prefix_tb.dat / prefix_r.dat, and it enables an
# exact cross-check against the Wannier centres (see centre_cross_check).
TRIAL = "trial_03_Sr_sp_Ti_pd_O_sp"

BASE_RUN = STO_BASE_RUN
BASE_W90 = BASE_RUN / TRIAL
PH_OUT = BASE_RUN / f"{PREFIX}.ph.out"          # DFPT reference, undisplaced
DISPLACED_UPF = BASE_RUN / "Ti-sp_r.upf"         # source of Z^ion, per atom

# Displaced structure.  MODE_W90 is the directory holding the .amn produced by
# generate_rigid_shift_amn.sh together with the Wannier90 output built from it.
# "no_displacement" keeps the *unshifted* base Wannier functions as trial
# orbitals, which is what makes Lambda(R) = 0 (the frozen Wannier gauge) the
# self-consistent choice below.  Point this at "rigid_shift" instead only if
# you also supply LAMBDA_FILE, since rigidly shifted trial orbitals move the
# Wannier functions with u and Lambda(R) is then not zero.
MODE_RUN = STO_TI_Q_01A_RUN
MODE_W90 = MODE_RUN / TRIAL / "no_displacement"   # must match TRIAL

# --position-source mmn rebuilds <0s|r|Rt> from the .mmn overlaps instead of
# reading prefix_r.dat.  The two routes are different finite-mesh
# discretizations; wannier/wannier_io.py explains the centre-aware link
# phasing now implemented in-house.  The .mmn depends only on the nscf, not on the projections, so
# a trial set whose own run kept no .mmn can borrow one from any other trial
# built on the same nscf; the borrow is checked, not assumed.  Keyed by trial,
# value is the directory to take the .mmn from (None = that run's own).
# The base .mmn now sits at the NSCF root (BASE_RUN) rather than inside any one
# trial, which is where it belongs: it is a property of the NSCF, not of the
# projections.  check_mmn_compatibility accepts a donor that is a *parent* of
# the borrowing directory -- containment establishes the shared NSCF, which the
# root directory cannot demonstrate with its own .nnkp because it has none.
MMN_DONOR = {}   # set from DISPLACEMENT[atom]["mmn_donor"] in the notebook config

# Cache for the in-house A(R) assembly; the .mmn parse is ~30 s and ~1 GB.
AA_CACHE_DIR = STO_AA_CACHE

# Displacement: Ti z, forward difference only (note "One-sided finite difference").
D_BETA = 0.01                       # [Angstrom]
DISPLACED_ATOM = "Ti"
DISPLACED_AXIS = 2                  # 0=x, 1=y, 2=z

# Everything that depends on *which* atom was displaced.  Selected with
# --displaced-atom (now the DISPLACED_ATOM config variable); the Ti entry is
# the default.
#
# ``excluded_on_atom`` is the part of the excluded occupied manifold that sits
# on the displaced atom and therefore moves into the ionic term (see the
# Z^ion,eff derivation below).  It is a property of the *pair* (window, atom),
# not of the window alone, which is why it lives here and not in TRIAL_SETUP.
#
# ``semicore_n_occ`` is the --semicore-test null.  That test integrates the
# lowest states and expects ~0 because they carry no character of the displaced
# atom.  States 0-11 are Sr 4p + O 2s: 2.5% Ti, so the null is calibrated for
# Ti -- but 49% Sr, so for an Sr displacement the very same group *must*
# respond and None disables the test rather than reporting a meaningless FAIL.
DISPLACEMENT = {
    "Ti": dict(
        mode_run=STO_TI_Q_01A_RUN,
        upf="Ti-sp_r.upf",
        default_trial="trial_03_Sr_sp_Ti_pd_O_sp",
        # Ti 3s (2) + Ti 3p (6) excluded by dis_win_min = -10 eV; trial_03's
        # wider window leaves only Ti 3s out.
        excluded_on_atom={"trial_01_Sr_spd_Ti_spd_O_sp": 8,
                          "trial_02_Sr_p_O_sp": 8,
                          "trial_03_Sr_sp_Ti_pd_O_sp": 2},
        semicore_n_occ={"trial_01_Sr_spd_Ti_spd_O_sp": 12,
                        "trial_02_Sr_p_O_sp": 12,
                        "trial_03_Sr_sp_Ti_pd_O_sp": None},
        # Only trial_03's displaced run kept an .mmn; the others borrow it.
        mmn_donor={
            "trial_01_Sr_spd_Ti_spd_O_sp": dict(
                base=BASE_RUN,
                mode=DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889"),
            "trial_02_Sr_p_O_sp": dict(
                base=BASE_RUN,
                mode=DATA / "SrTiO3/Ti_Q_01A/output/181839bc889/trial_03_Sr_sp_Ti_pd_Osp/no_displacement/185119bc889"),
            "trial_03_Sr_sp_Ti_pd_O_sp": dict(base=BASE_RUN),
        },
    ),
    "Sr": dict(
        mode_run=STO_SR_Q_01A_RUN,
        upf="Sr-sp_r.upf",
        default_trial="trial_02_Sr_p_O_sp",
        # In trial_02, of the ten bands below dis_win_min, only Sr 4s (2) sits
        # on Sr; Ti 3s and Ti 3p are on the undisplaced atom.  Trial_03's wider
        # window includes Sr 4s and leaves only Ti 3s out, so no excluded
        # electron follows the displaced Sr ion.
        excluded_on_atom={"trial_02_Sr_p_O_sp": 2,
                          "trial_03_Sr_sp_Ti_pd_O_sp": 0},
        # States 0-11 are 49% Sr: they respond, so there is no null here.
        semicore_n_occ={"trial_02_Sr_p_O_sp": None,
                        "trial_03_Sr_sp_Ti_pd_O_sp": None},
        # The Sr displaced .mmn is a property of the shared NSCF run and lives
        # at its DFT output root, one level above the trial directories.
        mmn_donor={"trial_02_Sr_p_O_sp": dict(base=BASE_RUN,
                                               mode=STO_SR_Q_01A_RUN),
                   "trial_03_Sr_sp_Ti_pd_O_sp": dict(
                       base=BASE_RUN,
                       mode=STO_SR_Q_01A_RUN)},
    ),
}

# Occupied states *of the Wannier model*.  Left as None this is counted as the
# bands lying below the SCF highest-occupied level and cross-checked against
# N_OCC_EXPECTED; set an int to force it.
#
# This is NOT the electron count.  The pseudopotentials put 40 electrons in the
# cell (Sr-sp 10 + Ti-sp 12 + 3 x O 6), but ``dis_win_min = -10 eV`` throws the
# deepest ten occupied spinor bands out of the Wannier window.  From the QE
# eigenvalues at the first nscf k-point:
#
#     -46.34 x2   Ti 3s          |  excluded, below dis_win_min
#     -23.05 x2   Sr 4s          |
#     -22.65 x2 } Ti 3p (SOC)    |
#     -22.40 x4 }                |
#     ------------------------------------------------------------
#      -7.31 x2 } Sr 4p (SOC)    |  inside the window
#      -6.33 x4 }                |
#      -5.42 x2 } O 2s           |
#      -4.29 x4 }                |
#       7.18 .. 10.06  x18  O 2p |
#
# so the Wannier occupied manifold is 6 + 6 + 18 = 30 states, ending at the SCF
# VBM of 10.056 eV.  The in-window assignment is not guesswork: projecting the
# base eigenvectors onto the Sr / Ti / O blocks of the Wannier basis gives
# states 0-11 as 49% Sr, 2.5% Ti, 48% O -- Sr 4p + O 2s, with no Ti 3p.  Both
# Ti semicore shells are therefore outside the window.
N_OCC = None

# Per-trial band accounting. n_occ_expected is a cross-check on the detected
# filling; n_excluded_total is verified at run time against (sum z_valence) -
# n_occ, and n_excluded_on_displaced is the part of it sitting on the displaced
# atom, which moves into the ionic term.
TRIAL_SETUP = {
    # dis_win_min = -10 eV: Ti 3s, Ti 3p and Sr 4s all excluded.
    "trial_01_Sr_spd_Ti_spd_O_sp": dict(n_occ=30, excluded_total=10, excluded_on_atom=8,
                                        semicore_n_occ=12),
    # dis_win = dis_froz = [-20, 11]: same exclusions, no disentanglement.
    "trial_02_Sr_p_O_sp":          dict(n_occ=30, excluded_total=10, excluded_on_atom=8,
                                        semicore_n_occ=12),
    # dis_win_min = -40 eV reaches below Sr 4s and Ti 3p, so only Ti 3s is left out.
    # That puts Ti 3p *inside* the window, so the lowest states are no longer a
    # Ti-free group and --semicore-test has no calibrated null here: states 0-11
    # sit partly on the displaced atom and must respond.  semicore_n_occ = None
    # disables the test rather than reporting a meaningless FAIL.
    "trial_03_Sr_sp_Ti_pd_O_sp":   dict(n_occ=38, excluded_total=2,  excluded_on_atom=2,
                                        semicore_n_occ=None),
}

# Ionic charge to add to the zz component.
#
# The note says to use z_valence = +12, which is right only when the Wannier
# space carries *every* valence electron.  It does not: the ten bands above are
# missing, and eight of them -- Ti 3s (2) and Ti 3p (6) -- sit on the displaced
# atom.  Both shells are core-like and translate rigidly with the Ti nucleus,
# contributing exactly -8 to dP/du_Ti, which the Wannier integral cannot see.
# Folding that into the ionic term,
#
#     Z^ion,eff_Ti = z_valence - N_EXCLUDED_ON_DISPLACED = 12 - 8 = +4,
#
# which is also the chemically sensible Ti core charge (nucleus screened by the
# 3s3p shells, i.e. nominal Ti 4+ with d0).  The remaining electronic response
# is then *positive*: Z*_Ti = 7.29 exceeds the nominal +4 precisely because of
# the anomalous Ti 3d - O 2p charge transfer.  N_EXCLUDED_TOTAL is checked at
# run time against (sum of z_valence) - n_occ, so widening or narrowing
# dis_win_min cannot change the accounting unnoticed.


# Overall sign relating Tr B_{k_beta,u} to dP_beta/du.  +1 is the note's
# convention, and it is correct.  Derivation: r_bar_n = (Omega/(2pi)^3) Int
# d^3k A_nn (Blount) with A = i<u|d u>, so P^el = -(e/Omega) Sum_n r_bar_n =
# -(e/(2pi)^3) Int d^3k Tr A.  Then
#     Z^el = (Omega/e) dP^el/du = -(Omega/(2pi)^3) Int d^3k d_u Tr A_{k_beta},
# and since Tr[A_u, A_k] = 0 and Int d^3k d_{k_beta} Tr A_u = 0 on the periodic
# BZ, Int Tr B_{u,k_beta} = Int d_u Tr A_{k_beta}, giving
#     Z^el_{Ti,z beta} = -(Omega/(2pi)^3) Int Tr B_{u,k_beta}
#                      = +(Omega/(2pi)^3) Int Tr B_{k_beta,u}.
# The code stores B[a,3] = B_{k_a,u}, so the prefactor is +1.
CURVATURE_SIGN = +1

# --semicore-test: the Sr 4p + O 2s group (states 0-11, isolated by a 9.3 eV
# gap) is deep, atomic-like, and carries essentially no Ti character (2.5%).
# None of it is bound to the displaced atom, so it must give Z^el close to
# ZERO -- a null test of the whole pipeline that needs no DFPT reference and
# no knowledge of the answer.  A large value here would mean the machinery is
# manufacturing a response where there is none.
SEMICORE_N_OCC = 12
SEMICORE_EXPECTED = 0.0
SEMICORE_TOL = 0.5

NKS = [4, 6, 8, 10, 12, 14, 16]     # BZ mesh sweep, nk^3 points each
K_BATCH = 512
MIN_HOPPING_NORM = 1e-5
DTYPE = np.complex128               # 60 orbitals is cheap; do not use complex64

_AXIS_NAMES = "xyz"

# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #


def parse_ph_born_charges(path: Path, *, asr: bool) -> dict[int, tuple[str, np.ndarray]]:
    r"""Born effective-charge tensors from a ``ph.x`` output.

    ``ph.x`` prints rows labelled by the **field** direction, ``Ex (...)`` being
    the forces along (x,y,z) for a field along x.  The printed matrix is
    therefore ``M[beta][alpha]``, and this returns ``Z[alpha][beta] = M.T`` so
    that ``Z[a]`` is the row of the note, ``Z*_{kappa, alpha beta}``.  For ideal
    cubic SrTiO3 every tensor here is diagonal, so the transpose is immaterial;
    it matters for lower-symmetry structures.

    Parameters
    ----------
    asr : bool
        Select the block with or without the acoustic sum rule applied.  The
        note calls for ``asr=False``: displacing Ti alone cannot test a sum rule
        over all atoms.
    """
    text = path.read_text(errors="ignore")
    marker = (
        "Effective charges (d Force / dE) in cartesian axis with asr applied:"
        if asr
        else "Effective charges (d Force / dE) in cartesian axis without acoustic sum rule applied (asr)"
    )
    start = text.find(marker)
    if start < 0:
        raise ValueError(f"No '{'with' if asr else 'without'} asr' Z* block in {path}")
    block = text[start + len(marker):]
    # Stop before the next "Effective charges" header, otherwise the non-ASR
    # block runs on into the ASR one and the later atoms overwrite the earlier.
    nxt = block.find("Effective charges")
    if nxt >= 0:
        block = block[:nxt]

    atom_re = re.compile(r"atom\s+(\d+)\s+(\S+)\s+Mean Z\*:")
    row_re = re.compile(r"E\*?[xyz]\s*\(([^)]*)\)")

    out: dict[int, tuple[str, np.ndarray]] = {}
    entries = list(atom_re.finditer(block))
    for n, m in enumerate(entries):
        stop = entries[n + 1].start() if n + 1 < len(entries) else len(block)
        rows = row_re.findall(block[m.end():stop])
        if len(rows) != 3:
            break  # ran past the end of this block into the next section
        M = np.array([[float(x) for x in r.split()] for r in rows], dtype=float)
        out[int(m.group(1))] = (m.group(2), M.T)
    if not out:
        raise ValueError(f"Could not parse any Z* tensors from {path}")
    return out


def parse_upf_zvalence(path: Path) -> float:
    """``z_valence`` from a UPF pseudopotential (the ionic charge)."""
    head = path.read_text(errors="ignore")[:20000]
    m = re.search(r'z_valence\s*=\s*"?\s*([0-9.eEdD+-]+)', head, re.I)
    if m is None:
        m = re.search(r"([0-9.]+)\s+Z valence", head, re.I)
    if m is None:
        raise ValueError(f"No z_valence found in {path}")
    return float(m.group(1).replace("D", "E").replace("d", "e"))


def parse_win_atoms(win_path: Path) -> tuple[list[str], np.ndarray]:
    """Species and Cartesian positions (Angstrom) from a ``.win atoms_cart`` block."""
    lines = win_path.read_text(errors="ignore").splitlines()
    try:
        i0 = next(i for i, s in enumerate(lines)
                  if s.strip().lower().startswith("begin atoms_cart"))
    except StopIteration:
        raise ValueError(
            f"No atoms_cart block in {win_path} (atoms_frac is not handled here)."
        ) from None
    i1 = next(i for i in range(i0 + 1, len(lines))
              if lines[i].strip().lower().startswith("end atoms_cart"))
    body = [s.strip() for s in lines[i0 + 1:i1] if s.strip()]
    scale = 1.0
    if body and body[0].lower() in ("ang", "angstrom", "bohr"):
        if body[0].lower() == "bohr":
            scale = 0.529177210903
        body = body[1:]
    species, pos = [], []
    for s in body:
        t = s.split()
        species.append(t[0])
        pos.append([float(x) * scale for x in t[1:4]])
    return species, np.array(pos, dtype=float)


def parse_mp_grid(win_path: Path) -> tuple[int, int, int]:
    """``mp_grid`` from a ``.win``; the k-mesh the Wannierization was built on."""
    m = re.search(r"^\s*mp_grid\s*[=:]?\s*(\d+)\s+(\d+)\s+(\d+)",
                  win_path.read_text(errors="ignore"), re.I | re.M)
    if m is None:
        raise ValueError(f"No mp_grid in {win_path}")
    return tuple(int(g) for g in m.groups())


def parse_win_setup(win_path: Path) -> dict:
    """num_wann / num_bands and the projections block from a ``.win``."""
    text = win_path.read_text(errors="ignore")
    out = {}
    for key in ("num_wann", "num_bands"):
        m = re.search(rf"^\s*{key}\s*[=:]?\s*(\d+)", text, re.I | re.M)
        if m:
            out[key] = int(m.group(1))
    block = re.search(r"begin projections(.*?)end projections", text, re.I | re.S)
    out["projections"] = (
        tuple(ln.strip().lower() for ln in block.group(1).splitlines() if ln.strip())
        if block else ()
    )
    return out


def parse_scf_species(scf_path: Path) -> list[str]:
    """Species symbols in ``ATOMIC_SPECIES`` order, with their UPF file names."""
    lines = scf_path.read_text(errors="ignore").splitlines()
    i0 = next(i for i, s in enumerate(lines) if s.strip().upper().startswith("ATOMIC_SPECIES"))
    out = []
    for s in lines[i0 + 1:]:
        t = s.split()
        if len(t) < 3 or not re.match(r"^[A-Za-z]", t[0]):
            break
        out.append((t[0], t[2]))
    return out


def parse_scf_atoms(scf_path: Path) -> tuple[list[str], np.ndarray]:
    """Species and positions from a pw.x ``ATOMIC_POSITIONS angstrom`` block."""
    lines = scf_path.read_text(errors="ignore").splitlines()
    i0 = next(i for i, s in enumerate(lines) if s.strip().upper().startswith("ATOMIC_POSITIONS"))
    if "angstrom" not in lines[i0].lower():
        raise ValueError(f"{scf_path}: expected ATOMIC_POSITIONS in angstrom")
    species, pos = [], []
    for s in lines[i0 + 1:]:
        t = s.split()
        if len(t) < 4 or not re.match(r"^[A-Za-z]", t[0]):
            break
        species.append(t[0])
        pos.append([float(x) for x in t[1:4]])
    return species, np.array(pos, dtype=float)


# --------------------------------------------------------------------------- #
# Curvature
# --------------------------------------------------------------------------- #


def mixed_curvature_traces(
    model_base,
    model_mode,
    k_red: np.ndarray,
    n_occ: int,
    d_beta: float,
    pos_base: cv.PositionTerms,
    pos_mode: cv.PositionTerms,
    lambda_prepared=None,
) -> dict[str, np.ndarray]:
    r"""Per-k traces ``Tr B^red_{k_a, u}`` for the int / cross / ext pieces.

    Everything is evaluated at ``u = 0`` (the undisplaced structure): the base
    eigenstates set the gauge, and the displaced model enters only through the
    forward differences ``dH/du`` and ``dA_i/du`` (``axion.beta_terms``).

    ``lambda_prepared`` is ``axion.prepare_lambda(Lambda_R, dtype)``; ``None`` is the
    frozen Wannier gauge ``A_u = 0``.  ``n_occ`` bands ``0 .. n_occ-1`` are occupied.

    Returns arrays of shape ``(3, Nk)``, complex; the imaginary parts are
    round-off and are reported by the caller as a diagnostic.
    """
    f_base = cv.fields(model_base, pos_base, k_red)
    f_mode = cv.fields(model_mode, pos_mode, k_red)
    if lambda_prepared is None:
        beta = beta_terms(f_base, f_mode, d_beta)
    else:
        A_beta, dA_beta = parametric_connection(pos_base, lambda_prepared, k_red)
        beta = beta_terms(f_base, f_mode, d_beta, A_beta, dA_beta)
    c = cv.connection(f_base, n_occ, beta)

    trace = lambda matrix: np.trace(matrix[:3, 3], axis1=-2, axis2=-1)
    return {
        "int": trace(cv.omega_internal(c)),
        "cross": trace(cv.omega_cross(c)),
        "ext": trace(cv.omega_external(c)),
    }


def z_electronic(traces: dict[str, np.ndarray], n_k_total: int, B_inv: np.ndarray) -> dict[str, np.ndarray]:
    r"""Reduced-axis trace sums -> Cartesian ``Z^el_{Ti,z beta}``.

    ``Z^el_beta = sum_a (B^-1)_{beta,a} * (1/N_k) sum_k Tr B^red_{a,u}``.  The
    caller accumulates ``sum_k`` over k-batches and passes the total ``N_k``.
    """
    return {
        key: CURVATURE_SIGN * (B_inv @ (val / n_k_total))
        for key, val in traces.items()
    }


# --------------------------------------------------------------------------- #
# Driver helpers
# --------------------------------------------------------------------------- #


def resolve_w90_dir(d: Path, label: str, prefix: str = PREFIX) -> Path:
    """``d`` itself if it holds ``prefix.win``, else its newest such subdirectory.

    The run directories here gain a job-id level (``w90/182494bc889``) whenever a
    Wannierization is repeated, so resolving one level keeps the configuration
    above from going stale on every rerun. The choice is printed, never guessed
    silently.
    """
    if (d / f"{prefix}.win").exists():
        return d
    if not d.is_dir():
        return d
    subs = sorted(
        (s for s in d.iterdir() if s.is_dir() and (s / f"{prefix}.win").exists()),
        key=lambda s: s.stat().st_mtime,
    )
    if not subs:
        return d
    if len(subs) > 1:
        print(f"  {label}: {len(subs)} candidate runs under {d}, "
              f"using the newest: {subs[-1].name}")
    else:
        print(f"  {label}: resolved to {d.name}/{subs[-1].name}")
    return subs[-1]


def position_source_available(d: Path, prefix: str = PREFIX) -> tuple[bool, str]:
    """Whether ``d`` holds a position matrix pythtb can actually read.

    pythtb reads ``prefix_tb.dat`` (``write_tb``) or the legacy ``prefix_r.dat``,
    and does **not** decompress: a directory holding only the ``.gz`` that
    ``send_job.sh`` leaves behind loads without a position matrix, silently
    disabling every external term.
    """
    for name in (f"{prefix}_tb.dat", f"{prefix}_r.dat"):
        if (d / name).exists():
            return True, name
    gz = [f"{n}.gz" for n in (f"{prefix}_tb.dat", f"{prefix}_r.dat")
          if (d / f"{n}.gz").exists()]
    if gz:
        return False, f"only compressed {', '.join(gz)} -- run: gunzip {d}/{gz[0]}"
    return False, f"no {prefix}_tb.dat or {prefix}_r.dat"


def preflight(self_test: bool, base_w90: Path, mode_w90: Path,
              position_source: str = "rdat", donor: dict | None = None,
              *, prefix: str = PREFIX, ph_out: Path = PH_OUT,
              displaced_upf: Path = DISPLACED_UPF) -> list[str]:
    """Files that must exist before anything is loaded."""
    missing = []
    if position_source == "mmn":
        donor = donor or {}
        pairs = [("base", base_w90, donor.get("base") or base_w90)]
        if not self_test:
            pairs.append(("mode", mode_w90, donor.get("mode") or mode_w90))
        for label, d_chk, d_mmn in pairs:
            if not (d_chk / f"{prefix}.chk").exists():
                missing.append(
                    f"{d_chk / f'{prefix}.chk'}  ({label} .chk -- required by "
                    "position_source='mmn'; it holds the U matrices and cannot "
                    "be borrowed from another trial)")
            if not (d_mmn / f"{prefix}.mmn").exists():
                missing.append(
                    f"{d_mmn / f'{prefix}.mmn'}  ({label} .mmn -- may be borrowed "
                    "from any trial sharing this nscf; set MMN_DONOR)")
    for p in (base_w90 / f"{prefix}.win", ph_out, displaced_upf):
        if not p.exists():
            missing.append(str(p))
    dirs = [("base", base_w90)] + ([] if self_test else [("mode", mode_w90)])
    for label, d in dirs:
        if not d.is_dir():
            missing.append(f"{d}  ({label} Wannier90 directory)")
            continue
        if label == "mode" and not (d / f"{prefix}.win").exists():
            missing.append(str(d / f"{prefix}.win"))
        if position_source == "rfull":
            rfull = d / f"{prefix}_r_full.dat"
            if not rfull.is_file():
                missing.append(
                    f"{rfull}  ({label} pair-weighted postw90 position matrix)"
                )
        else:
            ok, why = position_source_available(d, prefix)
            if not ok:
                missing.append(f"{d}  (position matrix: {why})")
    return missing


def parse_scf_fermi(scf_out: Path) -> float:
    """Highest occupied level (or Fermi energy) in eV from a pw.x output."""
    text = scf_out.read_text(errors="ignore")
    m = None
    for pat in (r"highest occupied level \(ev\):\s*([-\d.]+)",
                r"highest occupied, lowest unoccupied level \(ev\):\s*([-\d.]+)",
                r"the Fermi energy is\s*([-\d.]+)"):
        found = re.findall(pat, text)
        if found:
            m = float(found[-1])
            break
    if m is None:
        raise ValueError(f"No Fermi/highest-occupied level in {scf_out}")
    return m


def parse_centres(path: Path) -> np.ndarray:
    """Wannier centres (Angstrom) from ``prefix_centres.xyz``, in file order."""
    lines = path.read_text(errors="ignore").splitlines()
    n = int(lines[0].split()[0])
    out = [[float(x) for x in t[1:4]]
           for t in (ln.split() for ln in lines[2:2 + n]) if t[0] == "X"]
    return np.array(out, dtype=float)


def centre_cross_check(base_w90: Path, mode_w90: Path, d_beta: float,
                        prefix: str = PREFIX) -> np.ndarray | None:
    r"""``Z^el`` from the Wannier centres alone, exact when n_occ == num_wann.

    With the Wannier space equal to the occupied manifold and Lambda = 0, the
    trace is invariant under the band rotation and the whole curvature is the
    external term, so

        Z^el_beta = -d_u Int d^3k_red Tr A^cart_beta = -d_u Sum_s r_bar_{s,beta}

    because the nominal positions tau are held fixed across the two models.
    That is the textbook Wannier-centre formula, and it uses nothing but
    ``prefix_centres.xyz`` -- no curvature, no k-mesh.  Returns ``None`` if
    either file is absent.
    """
    fb, fm = (d / f"{prefix}_centres.xyz" for d in (base_w90, mode_w90))
    if not (fb.exists() and fm.exists()):
        return None
    cb, cm = parse_centres(fb), parse_centres(fm)
    if cb.shape != cm.shape:
        raise ValueError(f"centre files disagree in size: {cb.shape} vs {cm.shape}")
    return -(cm - cb).sum(axis=0) / d_beta


def detect_n_occ(model, e_ref: float, tol: float = 0.1, nk: int = 4) -> int:
    """Occupied states of the Wannier model: those lying entirely below ``e_ref``.

    The Wannier space spans occupied *and* empty bands, and its occupied count
    is not the electron count -- bands below ``dis_win_min`` are simply not
    represented.  Anchoring to the SCF highest-occupied level is unambiguous;
    picking the largest gap is not, since in SrTiO3 the O 2s-2p ionic gap
    (9.3 eV) is much wider than the fundamental one (1.9 eV).
    """
    k = model.k_uniform_mesh([nk, nk, nk], include_endpoints=False)
    e = np.linalg.eigvalsh(model.hamiltonian(k, flatten_spin_axis=True))
    n = int(np.sum(e.max(axis=0) <= e_ref + tol))
    if n == 0:
        raise ValueError(
            f"No band group lies below the reference energy {e_ref:.3f} eV; "
            "is the model's zero_energy consistent with the SCF output?"
        )
    vbm = float(e[:, n - 1].max())
    if n == e.shape[1]:
        # Every Wannier function is occupied: the Wannierization spans the
        # valence manifold alone.  There is no CBM inside the model.
        print(f"  all {n} Wannier states lie below E_ref = {e_ref:.3f} eV "
              f"(VBM {vbm:.3f}) -- occupied-only Wannierization")
        return n
    cbm = float(e[:, n].min())
    print(f"  {n} states lie below E_ref = {e_ref:.3f} eV  "
          f"(VBM {vbm:.3f}, CBM {cbm:.3f}, indirect gap {cbm - vbm:.4f} eV)")
    if cbm <= vbm:
        print("  !! WARNING: the Wannier model is metallic at this filling; the "
              "Kubo denominators below assume a gap.")
    return n


def check_structures(self_test: bool, base_w90: Path, mode_w90: Path,
                      *, prefix: str = PREFIX, base_run: Path = BASE_RUN,
                      mode_run: Path = MODE_RUN, displaced_atom: str = DISPLACED_ATOM,
                      displaced_axis: int = DISPLACED_AXIS, d_beta: float = D_BETA) -> None:
    """Verify the two runs differ only by the intended displacement."""
    print("Consistency checks")
    b_spec, b_scf = parse_scf_atoms(base_run / f"{prefix}.scf.in")
    b_wspec, b_win = parse_win_atoms(base_w90 / f"{prefix}.win")

    # The .win atom block sets the projection centres.  If it disagrees with the
    # SCF geometry the Wannier functions are built around the wrong sites.
    if b_wspec != b_spec or not np.allclose(b_win, b_scf, atol=1e-6):
        d = np.abs(b_win - b_scf).max()
        print(f"  !! WARNING: base {prefix}.win positions differ from the base SCF "
              f"by up to {d:.4f} A.")
        for i, (s, w, c) in enumerate(zip(b_spec, b_win, b_scf)):
            if not np.allclose(w, c, atol=1e-6):
                print(f"     atom {i + 1} {s}: win {w}  scf {c}")
        print("     The base projections are not centred on the base atoms; the "
              "reference gauge is not the cubic one and the zx/zy components "
              "need not vanish.")
    else:
        print("  base .win projection centres match the base SCF geometry.")

    if self_test:
        return

    m_spec, m_scf = parse_scf_atoms(mode_run / f"{prefix}.scf.in")
    if m_spec != b_spec:
        raise ValueError("Base and displaced runs list different species/order.")
    delta = m_scf - b_scf
    moved = np.where(np.abs(delta).max(axis=1) > 1e-8)[0]
    print(f"  displaced atoms: {[f'{b_spec[i]}#{i + 1}' for i in moved] or 'none'}")
    if len(moved) != 1 or b_spec[moved[0]] != displaced_atom:
        raise ValueError(
            f"Expected exactly one displaced {displaced_atom}; got "
            f"{[(b_spec[i], delta[i]) for i in moved]}"
        )
    d_vec = delta[moved[0]]
    got = d_vec[displaced_axis]
    if not np.allclose(np.delete(d_vec, displaced_axis), 0.0, atol=1e-8):
        raise ValueError(f"{displaced_atom} displacement is not along "
                         f"{_AXIS_NAMES[displaced_axis]}: {d_vec}")
    print(f"  d{displaced_atom},{_AXIS_NAMES[displaced_axis]} = {got:+.6f} A "
          f"(configured d_beta = {d_beta:+.6f} A)")
    if not np.isclose(got, d_beta, rtol=0, atol=1e-8):
        raise ValueError(
            f"SCF displacement {got:+.6f} A does not match d_beta = {d_beta:+.6f} A. "
            "The finite difference would be divided by the wrong step."
        )

    # Catch a trial-set mismatch here, before the expensive model load.
    b_set = parse_win_setup(base_w90 / f"{prefix}.win")
    m_set = parse_win_setup(mode_w90 / f"{prefix}.win")
    if b_set.get("num_wann") != m_set.get("num_wann") or \
            b_set["projections"] != m_set["projections"]:
        raise ValueError(
            "Base and displaced runs use different projection sets:\n"
            f"  base  num_wann={b_set.get('num_wann')}  "
            f"projections={list(b_set['projections'])}\n"
            f"  mode  num_wann={m_set.get('num_wann')}  "
            f"projections={list(m_set['projections'])}\n"
            "The two Wannier spaces must be the same, or the finite difference "
            "measures\na change of projection rather than the displacement. "
            "Point the mode directory at the\ndisplaced run built with the same "
            "trial."
        )
    print(f"  both runs use num_wann={b_set['num_wann']}, "
          f"projections={list(b_set['projections'])}")

    # The displaced .win must place its projections on the displaced atoms.
    m_wspec, m_win = parse_win_atoms(mode_w90 / f"{prefix}.win")
    if m_wspec != m_spec or not np.allclose(m_win, m_scf, atol=1e-6):
        print(f"  !! WARNING: displaced {prefix}.win positions differ from the "
              f"displaced SCF by up to {np.abs(m_win - m_scf).max():.4f} A.")
    else:
        print("  displaced .win projection centres match the displaced SCF geometry.")
