import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
import re
from pathlib import Path

def kpath_distance(
        k_frac: np.ndarray, 
        b1: np.ndarray, 
        b2: np.ndarray, 
        b3: np.ndarray) -> np.ndarray:
    """
    Build 1D cumulative k-path distance (in 1/Å) from fractional k-points.

    Parameters
    ----------
    k_frac : (nks, 3)
        Fractional k-points (crystal coords).
    b1,b2,b3 : (3,) in 1/Å
        Reciprocal lattice basis vectors in Cartesian coords.

    Returns
    -------
    x : (nks,)
        Cumulative distance along the path.
    """
    B = np.vstack([b1, b2, b3]).T     # 3x3, columns are basis vectors
    k_cart = (k_frac @ B.T)           # (nks,3) Cartesian k
    dk = np.linalg.norm(np.diff(k_cart, axis=0), axis=1)
    x = np.zeros(len(k_cart), dtype=float)
    x[1:] = np.cumsum(dk)
    return x

def read_bands_dat_filband(path: str):
    """
    Parse QE bands.x filband output (e.g. '*_bands.dat').

    Returns
    -------
    E : (nks, nbnd) float
        Energies in eV (row = k-point, col = band).
    k_frac : (nks, 3) float
        Fractional (crystal) k-points for each row of E.
    meta : dict
        Parsed metadata (e.g., {'nbnd': 240, 'nks': 251}) if present.
    """
    # Try to grab nbnd/nks from header
    meta = {}
    m = re.search(r"nbnd\s*=\s*(\d+).+nks\s*=\s*(\d+)", open(path).read(5000), re.I|re.S)
    if m:
        meta["nbnd"] = int(m.group(1))
        meta["nks"]  = int(m.group(2))

    def is_k_marker(s: str) -> bool:
        # k-marker line has exactly three floats
        try:
            vals = [float(x) for x in s.split()]
            return len(vals) == 3
        except ValueError:
            return False

    klist = []          # list of [kx, ky, kz] (fractional)
    energies_rows = []  # list of list-of-energies per k
    ebuf = []

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if is_k_marker(s):
                # starting a new k-point: flush previous energies if any
                if ebuf:
                    energies_rows.append(ebuf)
                    ebuf = []
                kx, ky, kz = (float(x) for x in s.split())
                klist.append([kx, ky, kz])
            else:
                # energy values (possibly many per line)
                try:
                    vals = [float(x) for x in s.split()]
                except ValueError:
                    continue
                ebuf.extend(vals)

    # flush last k
    if ebuf:
        energies_rows.append(ebuf)

    # Convert
    E_raw = np.array(energies_rows, dtype=float)   # ragged → rectangular next
    k_frac = np.array(klist, dtype=float)

    # Infer dimensions if header absent / inconsistent
    nks = meta.get("nks", E_raw.shape[0])
    nbnd = meta.get("nbnd", int(max(len(row) for row in E_raw)))

    # Normalize to (nks, nbnd)
    E = np.full((nks, nbnd), np.nan, dtype=float)
    for i in range(min(nks, len(E_raw))):
        row = E_raw[i]
        E[i, :min(nbnd, len(row))] = row[:nbnd]

    return E, k_frac, meta


def read_w90_band_dat(path: str, return_labels: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """
    Read Wannier90 'bands' output '*_band.dat' (blocks per band, each row: k  E).
    Returns:
        k : (nks,) cumulative k-distance along the path (already metric-adjusted)
        E : (nband, nks) energies in eV
    """
    blocks, cur = [], []
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                if cur:
                    blocks.append(cur)
                    cur = []
                continue
            cur.append([float(x) for x in s.split()])
    if cur: 
        blocks.append(cur)
    if not blocks: 
        raise ValueError(f"No data blocks found in {path}")

    nks = len(blocks[0])
    k = np.array([row[0] for row in blocks[0]], float)
    E = np.array([[row[1] for row in blk] for blk in blocks], float)  # (nband, nks)

    if return_labels:
        path_labels = read_labelinfo(path.replace("_band.dat", "_band.labelinfo.dat"))
        return k, E, path_labels
    
    return k, E

def read_labelinfo(path: str) -> List[Tuple[str, float]]:
    """Read (label, x_position) from '*_band.labelinfo.dat'."""
    out = []
    with open(path, "r") as f:
        for ln in f:
            if not ln.strip(): 
                continue
            parts = ln.split()
            out.append((parts[0], float(parts[2])))  # label, cumulative distance
    return out

def plot_bands(
        k: np.ndarray, 
        E: np.ndarray, 
        label_ticks: List[Tuple[str, float]],
        ef: float | None = None, 
        title: str = "Band structure"
        ) -> Tuple[plt.Figure, plt.Axes]:
    """Plot bands vs 1D k-path; returns output filepath."""
    fig, ax = plt.subplots()

    # plot bands
    Y = E if ef is None else (E - ef)
    for band in Y:
        ax.plot(k, band, color='b', lw=0.7)

    # labels and ticks
    ax.set_xticks([pos for lab, pos in label_ticks])
    ax.set_xticklabels([lab for lab, pos in label_ticks])
    for pos in [pos for lab, pos in label_ticks]:
        ax.axvline(pos, color='gray', lw=0.5, ls='--')

    # Fermi level
    if ef is not None:
        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.set_ylabel(r"Energy − $E_F$ (eV)")
    else:
        ax.set_ylabel("Energy (eV)")   
        
    ax.set_xlim(k[0], k[-1])
    ax.set_xlabel(r"$k$-path distance")
    ax.set_title(title)

    return fig, ax

def read_dos(filename):
    data = np.loadtxt(filename)
    energies = data[:, 0]
    dos = data[:, 1]
    return energies, dos

def read_bands_gnu(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    blocks = []
    block = []

    for line in lines:
        if line.strip() == "":
            if block:
                blocks.append(np.array(block, dtype=float))
                block = []
        else:
            block.append([float(x) for x in line.strip().split()])

    if block:
        blocks.append(np.array(block, dtype=float))

    # Each block corresponds to one band, all with same k-path x
    blocks = np.array(blocks)  # shape: (nbands, nkpts, 2)

    k_path = blocks[0][:, 0]    # use k-distance from first block
    bands = np.array([b[:, 1] for b in blocks])  # shape: (nbands, nkpts)

    return k_path, bands

def get_symmetry_kpath_labels(distances, segments, labels):
    """
    distances: list of k-point cumulative distances (x-axis)
    segments: list of number of k-points between symmetry points (e.g., [20, 20, 20])
    labels: list of high-symmetry point labels (e.g., ['Γ', 'X', 'M', 'Γ'])
    """
    x_ticks = [0]
    x_labels = [labels[0]]
    idx = 0
    for seg, label in zip(segments, labels[1:]):
        idx += seg
        x_ticks.append(distances[idx])
        x_labels.append(label)
    return x_ticks, x_labels

def read_amn(path):
    """
    Robust parser for Wannier90 .amn files.

    Returns
    -------
    A : np.ndarray, complex128, shape = (nk, nbnd, nproj)
        A[k, m, n] = <psi_{m,k} | g_n>.
    kpoints : list[tuple[float,float,float]] or None
        k-points in crystal coords (only for block format that lists them).
    meta : dict
        {'header': (nbnd_h, nproj_h, nk_h), 'inferred': (nbnd, nproj, nk), 'layout': 'row'|'block'}
    """
    def decomment(s: str) -> str:
        return s.split('!')[0].split('#')[0].strip()

    def is_int(s: str) -> bool:
        try: int(s); return True
        except: return False

    def is_float(s: str) -> bool:
        try: float(s); return True
        except: return False

    lines = [decomment(L) for L in Path(path).read_text().splitlines()]
    lines = [L for L in lines if L]

    # Find the first line with three positive ints (wannier header), but don't fully trust it
    nbnd_h = nproj_h = nk_h = None
    header_idx = None
    for i, L in enumerate(lines):
        parts = L.split()
        if len(parts) >= 3:
            try:
                a, b, c = map(int, parts[:3])
                if a > 0 and b > 0 and c > 0:
                    nbnd_h, nproj_h, nk_h = a, b, c
                    header_idx = i
                    break
            except ValueError:
                pass
    if header_idx is None:
        raise ValueError("Could not locate the 'nbnd nproj nk' header in .amn")

    # Detect layout using a few lines after the header
    sample = []
    for L in lines[header_idx+1 : header_idx+1+20]:
        parts = L.split()
        if parts:
            sample.append(parts)

    has_row_layout = any(
        len(p) >= 5 and is_int(p[0]) and is_int(p[1]) and is_int(p[2]) and is_float(p[3]) and is_float(p[4])
        for p in sample
    )

    if has_row_layout:
        # ---- ROW LAYOUT: each data line = m n ik Re Im ----
        rows = []
        max_m = max_n = max_ik = 0
        for L in lines[header_idx+1:]:
            parts = L.split()
            if len(parts) >= 5 and is_int(parts[0]) and is_int(parts[1]) and is_int(parts[2]) and is_float(parts[3]) and is_float(parts[4]):
                m = int(parts[0]); n = int(parts[1]); ik = int(parts[2])
                re = float(parts[3]); im = float(parts[4])
                rows.append((ik, m, n, re, im))
                if m   > max_m:  max_m  = m
                if n   > max_n:  max_n  = n
                if ik  > max_ik: max_ik = ik

        # Infer true dims from the data (this resolves swapped nk/nproj)
        nbnd, nproj, nk = max_m, max_n, max_ik

        A = np.zeros((nk, nbnd, nproj), dtype=np.complex128)
        for ik, m, n, re, im in rows:
            A[ik-1, m-1, n-1] = re + 1j*im

        return A, None, {"header": (nbnd_h, nproj_h, nk_h), "inferred": (nbnd, nproj, nk), "layout": "row"}

    # ---- BLOCK LAYOUT: optional 'kx ky kz' line, then nbnd*nproj lines of m n Re Im ----
    nbnd, nproj, nk = nbnd_h, nproj_h, nk_h  # block format usually has correct header
    A = np.zeros((nk, nbnd, nproj), dtype=np.complex128)
    kpoints = []
    i = header_idx + 1
    k = 0
    count = 0

    while i < len(lines) and k < nk:
        parts = lines[i].split()

        # optional k-vector line (3 floats, optional weight)
        if count == 0 and len(parts) in (3, 4) and all(is_float(x) for x in parts[:3]):
            kpoints.append(tuple(float(x) for x in parts[:3]))
            i += 1
            continue

        # projection line: m n Re Im
        if len(parts) >= 4 and is_int(parts[0]) and is_int(parts[1]) and is_float(parts[2]) and is_float(parts[3]):
            m = int(parts[0]); n = int(parts[1]); re = float(parts[2]); im = float(parts[3])
            if not (1 <= m <= nbnd and 1 <= n <= nproj):
                raise ValueError(f"Index out of range at k={k+1}: m={m}, n={n} (nbnd={nbnd}, nproj={nproj})")
            A[k, m-1, n-1] = re + 1j*im
            count += 1
            if count == nbnd * nproj:
                k += 1
                count = 0
            i += 1
            continue

        i += 1

    if k != nk:
        raise ValueError(f"Incomplete .amn: parsed {k} of {nk} k-point blocks.")

    if len(kpoints) != nk:
        kpoints = None

    return A, kpoints, {"header": (nbnd_h, nproj_h, nk_h), "inferred": (nbnd, nproj, nk), "layout": "block"}

def read_bands_out(filepath):
    """k-points (nk, 3) and eigenvalues (nbands, nk) in eV from a pw.x ``calculation='bands'`` output."""
    import re

    # Matches floats, including cases like "0.6835-0.0000" (missing space)
    float_re = re.compile(r'[+-]?\d+\.\d+')
    with open(filepath, 'r') as f:
        lines = f.readlines()
    kpts, all_eigs = [], []
    i = 0
    while i < len(lines):
        line = lines[i]
        if 'bands (ev)' in line:
            k_str = line.split('k =')[1].split('(')[0]
            kpts.append([float(x) for x in float_re.findall(k_str)])
            i += 2  # blank line after the header
            eigs = []
            while i < len(lines) and lines[i].strip():
                eigs.extend(float(x) for x in float_re.findall(lines[i]))
                i += 1
            all_eigs.append(eigs)
        else:
            i += 1
    return np.array(kpts), np.array(all_eigs).T


def fermi_energy(scf_out):
    """Fermi energy (metals) or highest occupied level (insulators) in eV from a pw.x output."""
    import re

    value = None
    for line in open(scf_out):
        if "the Fermi energy is" in line or "highest occupied" in line:
            # "the Fermi energy is  11.2531 ev", "highest occupied level (ev):  10.0560",
            # "highest occupied, lowest unoccupied level (ev):  10.0491  11.9568"
            value = float(re.findall(r"-?\d+\.\d+", line.split(":")[-1] if ":" in line else line)[0])
    if value is None:
        raise ValueError(f"no Fermi energy or highest occupied level in {scf_out}")
    return value


def read_kpath(bands_in):
    """Reduced k-points (nk, 3) and [(index, label)] of a ``K_POINTS crystal_b`` block (labels after ``!``).

    The points are the ones pw.x lists in the matching ``bands.out``, so a Wannier model evaluated at
    them lines up with the DFT bands point for point."""
    lines = open(bands_in).read().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().upper().startswith("K_POINTS"))
    nodes, counts, ticks, index = [], [], [], 0
    for line in lines[start + 2: start + 2 + int(lines[start + 1].split()[0])]:
        values, _, label = line.partition("!")
        values = values.split()
        nodes.append([float(x) for x in values[:3]])
        counts.append(int(values[3]))
        label = label.strip()
        ticks.append((index, r"$\Gamma$" if label.upper() in ("G", "GAMMA") else label))
        index += int(values[3])
    nodes = np.array(nodes)
    points = [nodes[i] + (nodes[i + 1] - nodes[i]) * j / counts[i]
              for i in range(len(nodes) - 1) for j in range(counts[i])]
    return np.array(points + [nodes[-1]]), ticks


def read_pdos_by_element(directory, prefix):
    """{(element, orbital): (energy, dos)} summed over atoms, from ``projwfc.x`` ``pdos_atm#`` files.

    Element is the species label with trailing digits removed (``Ir1`` -> ``Ir``)."""
    import glob
    import os
    import re

    summed = {}
    for path in sorted(glob.glob(os.path.join(directory, f"{prefix}.pdos_atm#*"))):
        match = re.search(r'pdos_atm#\d+\((\w+?)\d*\)_wfc#\d+\((\w+)\)', os.path.basename(path))
        if match:
            energy, dos = read_dos(path)
            key = (match.group(1), match.group(2))
            summed[key] = (energy, summed[key][1] + dos if key in summed else dos)
    return summed
