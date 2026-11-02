"""Input files of a displaced structure from the base's, and the window check. Called by the step scripts:

    python inputs.py nscf   <base nscf.in> <structure scf.in>                       -> nscf.in on stdout
    python inputs.py win    <base .win> <base scf.in> <structure scf.in>            -> .win on stdout
    python inputs.py window <nscf.out> <.win> <tolerance eV> [<base nscf.out>]      -> report; exit 1 if any edge fails

nscf: the base's NSCF input with the structure's ATOMIC_POSITIONS. The cell must be the base's (a displaced structure
differs from the base only in its atomic positions).

win: the base's .win with each atom moved by its displacement (structure scf.in - base scf.in). The base .win may put
an atom at another lattice image than its scf.in (MBT: sites at |x| <= 1/2, bug W90 centre wrap), so the .win keeps
the base's image and only adds the displacement.

window: every edge (dis_win_min/max, dis_froz_min/max) at least <tolerance> from every eigenvalue at every k and, with
a base NSCF, the same number of bands below each edge at each k as the base (an edge near a band splits a k-star; a
band crossing an edge between base and structure builds the two models differently).
"""
import re
import sys

import numpy as np

from axion_phonon.io.qe import read_pw_eigenvalues


def _block(text, name, n):
    """(start, end) character span of the n lines after the line starting with ``name``, and those lines."""
    m = re.search(rf"^{name}\b.*\n", text, re.M)
    if m is None:
        sys.exit(f"no {name} card")
    lines = text[m.end():].splitlines(keepends=True)[:n]
    return m.start(), m.end() + sum(len(x) for x in lines), m.group(0) + "".join(lines)


def _nat(text):
    return int(re.search(r"^\s*nat\s*=\s*(\d+)", text, re.M | re.I).group(1))


def _cell(text):
    """The cell as written: the CELL_PARAMETERS card, or (ibrav != 0, e.g. SrTiO3) the ibrav / celldm / A..cosC lines."""
    if re.search(r"^CELL_PARAMETERS\b", text, re.M):
        return _block(text, "CELL_PARAMETERS", 3)[2].split()
    return sorted(re.findall(r"^\s*((?:ibrav|celldm\(\d\)|a|b|c|cosab|cosac|cosbc)\s*=\s*\S+)", text, re.M | re.I))


def _atoms(text):
    _, _, blk = _block(text, "ATOMIC_POSITIONS", _nat(text))
    head, *rows = blk.splitlines()
    return head, [(r.split()[0], np.array(r.split()[1:4], float)) for r in rows]


def nscf(base_nscf, scf):
    t, s = open(base_nscf).read(), open(scf).read()
    if _nat(t) != _nat(s):
        sys.exit("nat differs")
    if _cell(t) != _cell(s):
        sys.exit("the cell differs from the base")
    (h1, a1), (h2, a2) = _atoms(t), _atoms(s)
    if h1.split() != h2.split() or [x for x, _ in a1] != [x for x, _ in a2]:
        sys.exit("ATOMIC_POSITIONS units or species differ from the base")
    i, j, _ = _block(t, "ATOMIC_POSITIONS", _nat(t))
    sys.stdout.write(t[:i] + _block(s, "ATOMIC_POSITIONS", _nat(s))[2] + t[j:])


def win(base_win, base_scf, scf):
    w = open(base_win).read()
    m = re.search(r"^begin atoms_cart\s*\n(\s*ang\w*\s*\n)?(.*?)^end atoms_cart", w, re.M | re.S | re.I)
    if m is None or (m.group(1) or "ang").strip().lower()[:3] != "ang":
        sys.exit("base .win needs atoms_cart in Angstrom")
    (hb, ab), (hs, as_) = _atoms(open(base_scf).read()), _atoms(open(scf).read())
    if "angstrom" not in hb.lower() or "angstrom" not in hs.lower():
        sys.exit("scf.in ATOMIC_POSITIONS must be in angstrom")
    rows = [r.split() for r in m.group(2).strip().splitlines()]
    if [r[0] for r in rows] != [x for x, _ in ab] or [x for x, _ in ab] != [x for x, _ in as_]:
        sys.exit("atom order differs between the base .win, the base scf.in and the structure scf.in")
    out = []
    for r, (_, xb), (name, xs) in zip(rows, ab, as_):
        x = np.array(r[1:4], float) + (xs - xb)
        out.append(f"{name:4s} {x[0]:16.12f} {x[1]:16.12f} {x[2]:16.12f}\n")
    sys.stdout.write(w[:m.start(2)] + "".join(out) + w[m.end(2):])


def window(nscf_out, win_file, tol, base_out=None):
    edges = {}
    for line in open(win_file):
        m = re.match(r"\s*(dis_(?:win|froz)_(?:min|max))\s*=\s*(\S+)", line)
        if m:
            edges[m.group(1)] = float(m.group(2))
    E = read_pw_eigenvalues(nscf_out)
    base = read_pw_eigenvalues(base_out) if base_out else None
    print(f"{nscf_out}: {E.shape[0]} k, {E.shape[1]} bands" + (f"; base {base_out}" if base_out else ""))
    ok = True
    for name, x in sorted(edges.items()):
        d = np.abs(E - x)
        k, b = np.unravel_index(d.argmin(), d.shape)
        below = (E < x).sum(1)
        good = d.min() >= tol
        text = (f"  {name:13s} {x:9.4f} eV: closest eigenvalue {d.min() * 1e3:8.1f} meV away (k {k + 1}, band {b + 1});"
                f" bands below: {below.min()}-{below.max()}")
        if base is not None:
            diff = np.flatnonzero(below != (base < x).sum(1))
            good &= len(diff) == 0
            text += f"; k with a different count than the base: {len(diff)}"
        ok &= good
        print(text + ("  OK" if good else "  FAIL"))
    if edges.get("dis_win_max") is not None and (E < edges["dis_win_max"]).sum(1).max() >= E.shape[1]:
        print("  the top NSCF band is below dis_win_max  FAIL")
        ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    what, args = sys.argv[1], sys.argv[2:]
    if what == "nscf":
        nscf(*args)
    elif what == "win":
        win(*args)
    elif what == "window":
        window(args[0], args[1], float(args[2]), args[3] if len(args) > 3 else None)
    else:
        sys.exit(f"unknown: {what}")
