"""Read Quantum ESPRESSO files: wfcN.dat, the save's XML, PAW UPFs, pw.x and projwfc.x outputs. No physics here.

Units as QE stores them: bohr, bohr^-1 (the wfc files store k and the reciprocal vectors in bohr^-1, Cartesian).
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


def read_wfc(path: str | Path) -> dict:
    """``wfcN.dat`` (io_base.write_wfc). ``evc`` (nbnd, npol * igwx): row n is [spin-up G..., spin-down G...];
    ``mill`` (igwx, 3) the Miller indices of the G vectors, ``b`` rows b1, b2, b3, ``xk`` the k point."""
    f = FortranFile(str(path), "r")
    rec = f.read_record(np.dtype([("ik", "<i4"), ("xk", "<f8", 3), ("ispin", "<i4"), ("gamma", "<i4"),
                                  ("scalef", "<f8")]))[0]
    ngw, igwx, npol, nbnd = f.read_ints(np.int32)
    b = f.read_reals(np.float64).reshape(3, 3)
    mill = f.read_ints(np.int32).reshape(igwx, 3)
    evc = np.empty((nbnd, npol * igwx), complex)
    for n in range(nbnd):
        evc[n] = f.read_reals(np.float64).view(complex)
    f.close()
    if rec["gamma"]:
        raise NotImplementedError("gamma_only wavefunctions")
    return dict(ik=int(rec["ik"]), xk=np.array(rec["xk"]), npol=int(npol), nbnd=int(nbnd), igwx=int(igwx),
                b=b, mill=mill, evc=evc)


def read_structure(save: str | Path) -> dict:
    """Cell (rows, bohr), atoms ``tau`` (bohr), their species ``names``, and ``species`` [(name, UPF file)] of a save."""
    root = ET.parse(Path(save) / "data-file-schema.xml").getroot()
    out = root.find("output")
    species = [(s.get("name"), s.find("pseudo_file").text.strip()) for s in out.find("atomic_species")]
    struct = out.find("atomic_structure")
    atoms = struct.find("atomic_positions")
    names = [a.get("name") for a in atoms]
    tau = np.array([[float(x) for x in a.text.split()] for a in atoms])
    cell = np.array([[float(x) for x in struct.find("cell").find(f"a{i}").text.split()] for i in (1, 2, 3)])
    return dict(species=species, names=names, tau=tau, cell=cell, omega=abs(np.linalg.det(cell)),
                alat=float(struct.get("alat")))


def _floats(node) -> np.ndarray:
    return np.array(node.text.split(), float)


@dataclass
class Pseudo:
    """The parts of a PAW UPF (v2) that build S and the partial-wave terms. Radial functions are r * f(r).
    A norm-conserving UPF (no PP_AUGMENTATION, e.g. ONCV for SrTiO3) reads with ``augmented = False`` and zero Q:
    S = 1 and the overlap is the plane-wave part alone."""
    r: np.ndarray
    rab: np.ndarray
    beta: list                                      # projectors r * beta(r)
    lll: np.ndarray
    jjj: np.ndarray | None                          # None without spin-orbit
    kkbeta: int
    qqq: np.ndarray                                 # PP_Q, nbeta x nbeta
    qfunc0: np.ndarray                              # L = 0 augmentation functions, (nbeta, nbeta, mesh)
    aewfc: list = field(default_factory=list)       # all-electron partial waves (large component)
    pswfc: list = field(default_factory=list)       # pseudo partial waves
    iraug: int = 0
    augmented: bool = True


def read_upf(path: str | Path) -> Pseudo:
    text = Path(path).read_text()
    text = re.sub(r"&(?!amp;|lt;|gt;|quot;|apos;)", "&amp;", text)      # bare '&' in PP_INFO
    root = ET.fromstring(text[text.index("<UPF"):])
    mesh = root.find("PP_MESH")
    r, rab = _floats(mesh.find("PP_R")), _floats(mesh.find("PP_RAB"))
    nl = root.find("PP_NONLOCAL")
    nbeta = int(root.find("PP_HEADER").get("number_of_proj"))
    beta, lll, kbeta = [], [], []
    for i in range(1, nbeta + 1):
        node = nl.find(f"PP_BETA.{i}")
        beta.append(_floats(node))
        lll.append(int(node.get("angular_momentum")))
        kbeta.append(int(node.get("cutoff_radius_index")))
    aug = nl.find("PP_AUGMENTATION")
    qfunc0 = np.zeros((nbeta, nbeta, len(r)))
    if aug is None:                                                     # norm-conserving: no augmentation
        if root.find("PP_HEADER").get("pseudo_type", "").strip().upper() != "NC":
            raise ValueError(f"{path}: no PP_AUGMENTATION, but pseudo_type is not NC")
        iraug, qqq = 0, np.zeros((nbeta, nbeta))
    else:
        iraug = int(aug.get("cutoff_r_index"))
        qqq = _floats(aug.find("PP_Q")).reshape(nbeta, nbeta)
        if aug.get("q_with_l", "F").strip().lower() not in ("t", "true", ".true."):
            raise NotImplementedError("augmentation without q_with_l")
        for i in range(1, nbeta + 1):
            for j in range(i, nbeta + 1):
                node = aug.find(f"PP_QIJL.{i}.{j}.0")
                if node is not None:
                    qfunc0[i - 1, j - 1] = qfunc0[j - 1, i - 1] = _floats(node)
        qfunc0[:, :, iraug:] = 0.0                                      # read_upf_new: zero beyond iraug
    so = root.find("PP_SPIN_ORB")
    jjj = None
    if so is not None:
        jjj = np.array([float(so.find(f"PP_RELBETA.{i}").get("jjj")) for i in range(1, nbeta + 1)])
    full = root.find("PP_FULL_WFC")
    aewfc = [_floats(full.find(f"PP_AEWFC.{i}")) for i in range(1, nbeta + 1)] if full is not None else []
    pswfc = [_floats(full.find(f"PP_PSWFC.{i}")) for i in range(1, nbeta + 1)] if full is not None else []
    return Pseudo(r=r, rab=rab, beta=beta, lll=np.array(lll), jjj=jjj,
                  kkbeta=max(max(kbeta), iraug), qqq=qqq, qfunc0=qfunc0, aewfc=aewfc, pswfc=pswfc, iraug=iraug,
                  augmented=aug is not None)


def read_pseudos(structure: dict, pseudo_dir: str | Path) -> dict[str, Pseudo]:
    """The UPF of every species of ``read_structure``'s output, by species name."""
    return {name: read_upf(Path(pseudo_dir) / upf) for name, upf in structure["species"]}


def read_pw_eigenvalues(path: str | Path) -> np.ndarray:
    """Band energies (nk, nbnd), eV, from a pw.x NSCF output (one ``bands (ev)`` block per k)."""
    eig, cur, reading = [], [], False
    for line in open(path):
        if "bands (ev)" in line:
            if cur:
                eig.append(cur)
            cur, reading = [], True
            continue
        if reading:
            if "occupation numbers" in line or "highest occupied" in line or "Fermi energy" in line \
                    or "Writing" in line:
                reading = False
                continue
            cur += [float(x) for x in re.findall(r"-?\d+\.\d{4}", line)]
    if cur:
        eig.append(cur)
    if len({len(e) for e in eig}) != 1:
        raise ValueError(f"{path}: inconsistent band counts per k")
    return np.array(eig)


def read_bands_out(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """k points (nk, 3; Cartesian, 2 pi / alat) and eigenvalues (nbands, nk), eV, of a pw.x ``calculation='bands'``
    output."""
    number = re.compile(r"[+-]?\d+\.\d+")                       # also splits "0.6835-0.0000"
    lines = open(path).readlines()
    kpts, eigs, i = [], [], 0
    while i < len(lines):
        if "bands (ev)" in lines[i]:
            kpts.append([float(x) for x in number.findall(lines[i].split("k =")[1].split("(")[0])])
            i += 2
            e = []
            while i < len(lines) and lines[i].strip():
                e.extend(float(x) for x in number.findall(lines[i]))
                i += 1
            eigs.append(e)
        else:
            i += 1
    return np.array(kpts), np.array(eigs).T


def fermi_energy(path: str | Path) -> float:
    """Fermi energy (metals) or highest occupied level (insulators), eV, from a pw.x output (the last one printed)."""
    value = None
    for line in open(path):
        if "the Fermi energy is" in line or "highest occupied" in line:
            value = float(re.findall(r"-?\d+\.\d+", line.split(":")[-1] if ":" in line else line)[0])
    if value is None:
        raise ValueError(f"no Fermi energy or highest occupied level in {path}")
    return value


def read_symmetry(path: str | Path) -> dict:
    """The symmetry operations a pw.x output lists (verbosity 'high'): ``names`` (QE's description), ``cart``
    (n, 3, 3) the Cartesian rotation matrices, ``time_reversal`` (n,) bool (QE's 'Time Reversal 1': the operation is
    combined with time reversal), ``fractional`` (n,) bool (it carries a fractional translation), and
    ``translations_disabled`` (QE's 'fractional translations are disabled', e.g. for a supercell)."""
    text = Path(path).read_text()
    if "Sym. Ops." not in text:
        raise ValueError(f"{path}: no symmetry listing (run pw.x with verbosity = 'high')")
    block = text.split("Sym. Ops.")[1].split("Cartesian axes")[0]
    names, cart, trev, frac = [], [], [], []
    for entry in re.split(r"\n\s+isym =", block)[1:]:
        names.append(entry.split("\n")[0].split(None, 1)[1].strip())
        found = re.search(r"Time Reversal\s+(\d)", entry)
        trev.append(bool(int(found.group(1))) if found else False)
        rows = re.findall(r"\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)", re.split(r"\n\s*cart\.\s+s\(", entry)[1])[:3]
        cart.append(np.array(rows, float))
        frac.append("f =(" in entry)
    return dict(names=names, cart=np.array(cart), time_reversal=np.array(trev), fractional=np.array(frac),
                translations_disabled="fractional translations are disabled" in text)


def read_kpath(path: str | Path) -> tuple[np.ndarray, list]:
    """Reduced k points (nk, 3) and [(index, label)] of a pw.x input's ``K_POINTS crystal_b`` path (labels after
    ``!``), expanded exactly as pw.x lists them in the bands output."""
    lines = open(path).read().splitlines()
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
    points = [nodes[i] + (nodes[i + 1] - nodes[i]) * j / counts[i] for i in range(len(nodes) - 1) for j in range(counts[i])]
    return np.array(points + [nodes[-1]]), ticks


def read_pdos_by_element(directory: str | Path, prefix: str) -> dict:
    """{(element, orbital): (energy, dos)} summed over atoms (and over j with spin-orbit), from projwfc.x
    ``<prefix>.pdos_atm#`` files. The element is the species label without trailing digits (``Mn1`` -> ``Mn``); the
    orbital is s, p, d or f, with the wavefunction number added when an element has two of the same l (``s#2``: Mn 4s)."""
    files = sorted(Path(directory).glob(f"{prefix}.pdos_atm#*"))
    pat = re.compile(r"pdos_atm#\d+\(([A-Za-z]+?)\d*\)_wfc#(\d+)\(([spdf])")
    shells = {}                                                     # (element, l) -> wfc numbers present
    for path in files:
        m = pat.search(path.name)
        if m:
            shells.setdefault((m.group(1), m.group(3)), set()).add(int(m.group(2)))
    summed = {}
    for path in files:
        m = pat.search(path.name)
        if not m:
            continue
        el, n, l = m.group(1), int(m.group(2)), m.group(3)
        ls = sorted(shells[(el, l)])
        # spin-orbit splits one shell into j = l -/+ 1/2 with consecutive wfc numbers; a second shell of the same l
        # (Mn 3s and 4s) has its own pair: label it by its rank
        label = l
        if l == "s" and len(ls) > 1:
            label = f"s#{ls.index(n) + 1}"
        elif l != "s" and len(ls) > 2:
            label = f"{l}#{ls.index(n) // 2 + 1}"
        data = np.loadtxt(path)
        key = (el, label)
        summed[key] = (data[:, 0], summed[key][1] + data[:, 1] if key in summed else data[:, 1])
    return summed
