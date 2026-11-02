"""Read and write Wannier90 files: .chk, .eig, .win, .nnkp, .amn. No physics here.

``read_run`` gathers everything step 6 needs from one Wannier90 run folder into one dict of arrays; the physics in
``axion_phonon/wannier.py`` takes that dict. Lengths in Angstrom, energies in eV, k in reduced coordinates. A file
gzipped by the cluster epilog (``name.gz``) is read in place of ``name``.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

BOHR = 0.529177210903       # Angstrom


def _existing(path: str | Path) -> Path:
    """``path``, or ``path.gz`` when only that exists."""
    path = Path(path)
    if path.is_file():
        return path
    gz = Path(f"{path}.gz")
    if gz.is_file():
        return gz
    raise FileNotFoundError(path)


def _text(path: str | Path) -> str:
    path = _existing(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="ignore") as f:
        return f.read()


class _GzipFortranFile:
    """Sequential little-endian Fortran records from a gzipped file (records over 2 GiB come as subrecords)."""

    def __init__(self, path: Path):
        self._f = gzip.open(path, "rb")

    def read_record(self, dtype):
        chunks = []
        while True:
            marker = int(np.frombuffer(self._f.read(4), "<i4")[0])
            chunks.append(self._f.read(abs(marker)))
            self._f.read(4)
            if marker >= 0:
                break
        return np.frombuffer(b"".join(chunks), np.dtype(dtype)).copy()

    def close(self):
        self._f.close()


def read_chk(path: str | Path, m_matrix: bool = False) -> dict:
    """Wannier90 checkpoint. ``v_matrix`` (nk, num_bands, num_wann) = U_dis U, the Wannier functions in the Bloch basis.

    ``m_matrix=True`` also returns ``m_matrix`` (nk, nntot, num_wann, num_wann) = <u_mk|u_n,k+b> in the Wannier gauge,
    b in the .nnkp neighbour order (1.2 GB for MnBi2Te4 at 8x8x8).
    """
    path = _existing(path)
    f = _GzipFortranFile(path) if path.suffix == ".gz" else FortranFile(str(path), "r")
    ints, floats = (lambda: f.read_record("<i4")), (lambda: f.read_record("<f8"))

    def complex_():
        x = floats()
        return x[::2] + 1j * x[1::2]

    f.read_record("c")                                     # header
    nb = int(ints()[0])
    ints()                                                 # num_exclude_bands
    ints()                                                 # exclude_bands
    real_lattice = floats().reshape(3, 3, order="F")
    recip_lattice = floats().reshape(3, 3, order="F")
    nk = int(ints()[0])
    mp_grid = tuple(int(x) for x in ints())
    kpts = floats().reshape(nk, 3)
    nntot = int(ints()[0])
    nw = int(ints()[0])
    f.read_record("c")                                     # checkpoint label
    disentangled = bool(ints()[0])
    if disentangled:
        floats()                                           # omega_invariant
        lwindow = ints().reshape(nk, nb).astype(bool)
        ndimwin = ints()
        U_dis = complex_().reshape(nk, nw, nb).swapaxes(1, 2)
    U = complex_().reshape(nk, nw, nw).swapaxes(1, 2)
    m = f.read_record("<f8")
    M = m.view(np.complex128).reshape(nk, nntot, nw, nw).swapaxes(2, 3) if m_matrix else None
    del m
    centres = floats().reshape(nw, 3)
    spreads = floats()
    f.close()
    if disentangled:
        V = np.zeros((nk, nb, nw), complex)
        for ik in range(nk):
            V[ik, lwindow[ik]] = U_dis[ik][:ndimwin[ik]] @ U[ik]
    else:
        V = U
    out = dict(num_bands=nb, num_wann=nw, num_kpts=nk, nntot=nntot, mp_grid=mp_grid, kpt_red=kpts,
               real_lattice=real_lattice, recip_lattice=recip_lattice, v_matrix=V, wannier_centres=centres,
               wannier_spreads=spreads)
    if m_matrix:
        out["m_matrix"] = M
    return out


def read_eig(path: str | Path, num_kpts: int) -> np.ndarray:
    """``.eig`` as (nk, num_bands), eV."""
    e = np.loadtxt(_existing(path))
    nb = int(e[:, 0].max())
    if e.shape[0] != num_kpts * nb or not np.array_equal(e[:, 1].astype(int), np.repeat(np.arange(1, num_kpts + 1), nb)):
        raise ValueError(f"{path}: not {num_kpts} k x {nb} bands in order")
    return e[:, 2].reshape(num_kpts, nb)


def _block(text: str, name: str) -> list[list[str]]:
    m = re.search(rf"^\s*begin\s+{name}\s*$(.*?)^\s*end\s+{name}\s*$", text, re.I | re.M | re.S)
    if m is None:
        raise ValueError(f"no {name} block")
    lines = [line.split("!")[0].split("#")[0].split() for line in m.group(1).splitlines()]
    return [line for line in lines if line]


def read_win(path: str | Path) -> dict:
    """``.win``: ``lattice`` (3, 3) rows a_i (A), ``atoms`` (nat, 3) Cartesian (A), ``symbols``, ``mp_grid``,
    ``windows`` {dis_win_min, dis_win_max, dis_froz_min, dis_froz_max: eV} (those the file sets)."""
    text = _text(path)
    rows = _block(text, "unit_cell_cart")
    scale = 1.0
    if len(rows[0]) == 1:
        scale = BOHR if rows[0][0].lower().startswith("bohr") else 1.0
        rows = rows[1:]
    lattice = np.array(rows, float) * scale
    for name in ("atoms_cart", "atoms_frac"):
        try:
            rows = _block(text, name)
        except ValueError:
            continue
        unit = 1.0
        if len(rows[0]) == 1:
            unit = BOHR if rows[0][0].lower().startswith("bohr") else 1.0
            rows = rows[1:]
        symbols = [r[0] for r in rows]
        xyz = np.array([r[1:4] for r in rows], float)
        atoms = xyz * unit if name == "atoms_cart" else xyz @ lattice
        break
    else:
        raise ValueError(f"{path}: no atoms_cart or atoms_frac block")
    m = re.search(r"^\s*mp_grid\s*[=:]?\s*(\d+)\s+(\d+)\s+(\d+)", text, re.I | re.M)
    windows = {key: float(v) for key, v in re.findall(r"^\s*(dis_(?:win|froz)_(?:min|max))\s*[=:]?\s*(\S+)", text,
                                                         re.I | re.M)}
    return dict(lattice=lattice, atoms=atoms, symbols=symbols, mp_grid=tuple(int(x) for x in m.groups()),
                windows=windows)


def read_nnkp_neighbours(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """``.nnkp`` nnkpts block: neighbour index (nk, nntot), 0-based, and its reciprocal-lattice shift G (nk, nntot, 3)."""
    rows = _block(_text(path), "nnkpts")
    nntot = int(rows[0][0])
    table = np.array(rows[1:], int).reshape(-1, nntot, 5)
    if not np.array_equal(table[:, :, 0], np.repeat(np.arange(1, len(table) + 1)[:, None], nntot, axis=1)):
        raise ValueError(f"{path}: nnkpts rows are not ordered by k point")
    return table[:, :, 1] - 1, table[:, :, 2:]


def read_run(folder: str | Path, prefix: str, overlaps: bool = True) -> dict:
    """Everything step 6 needs from one Wannier90 run: the .chk (with m_matrix), .eig, .nnkp neighbours and .win.

    Keys: ``lattice`` (rows a_i, A), ``reciprocal`` (rows b_i, 1/A), ``mp_grid``, ``kpts`` (nk, 3) reduced,
    ``V`` (nk, num_bands, J), ``energies`` (nk, num_bands), ``M`` (nk, nntot, J, J), ``neighbours`` (nk, nntot),
    ``G`` (nk, nntot, 3), ``centres`` (J, 3) Cartesian, ``atoms`` (nat, 3) Cartesian, ``symbols``.
    ``overlaps=False`` skips M, neighbours and G (and the .nnkp): enough for H(R), e.g. band structures.
    """
    folder = Path(folder)
    chk = read_chk(folder / f"{prefix}.chk", m_matrix=overlaps)
    win = read_win(folder / f"{prefix}.win")
    if tuple(win["mp_grid"]) != tuple(chk["mp_grid"]):
        raise ValueError(f"{folder}: .win and .chk k meshes differ")
    if np.abs(win["lattice"] - chk["real_lattice"]).max() > 1e-8:
        raise ValueError(f"{folder}: .win and .chk lattices differ")
    run = dict(lattice=chk["real_lattice"], reciprocal=chk["recip_lattice"], mp_grid=chk["mp_grid"],
               kpts=chk["kpt_red"], V=chk["v_matrix"], energies=read_eig(folder / f"{prefix}.eig", chk["num_kpts"]),
               centres=chk["wannier_centres"], atoms=win["atoms"], symbols=win["symbols"])
    if overlaps:
        neighbours, G = read_nnkp_neighbours(folder / f"{prefix}.nnkp")
        if len(neighbours) != chk["num_kpts"]:
            raise ValueError(f"{folder}: .nnkp does not match the .chk")
        run.update(M=chk["m_matrix"], neighbours=neighbours, G=G)
    return run


def write_amn(path: str | Path, A: np.ndarray, header: str) -> None:
    """``.amn`` from A (nk, num_bands, num_wann)."""
    nk, nb, nw = A.shape
    k, n, m = np.meshgrid(np.arange(nk), np.arange(nw), np.arange(nb), indexing="ij")
    rows = np.column_stack([m.ravel() + 1, n.ravel() + 1, k.ravel() + 1, A[k, m, n].real.ravel(),
                            A[k, m, n].imag.ravel()])
    with open(path, "w") as f:
        f.write(header + "\n" + f"{nb:12d}{nk:12d}{nw:12d}\n")
        np.savetxt(f, rows, fmt="%5d%5d%5d%18.12f%18.12f")
