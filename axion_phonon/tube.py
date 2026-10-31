r"""Step 6, last part: the k list for the BZ integral of c2 in MnBi2Te4, and the integral itself.

c2 is large (~1e3) wherever the direct gap is ~0.2 eV, which in MnBi2Te4 is the whole k_z line through Gamma (the
bands hardly disperse along k_z, c = 81.85 A). In the model's reduced coordinates (b1, b2, b3; a1 + a2 + a3 = c)
that line is k = t (1, 1, 1). The integral over the unit cube is unchanged by the unimodular change of basis to
(g1, g2, g3) = (b1, b2, b1 + b2 + b3), in which the tube is u1 = u2 = 0:

    k = (u1 + u3, u2 + u3, u3) mod 1.

Sampling: n3 planes u3 = i / n3; in each, the nk x nk grid in (u1, u2), with every cell within ``half`` cells of
the tube replaced by an m x m subgrid. The subgrid is nested (levels 1, 2, 4, .., m), so one list gives the
integral at every level and at n3, n3/2, n3/4, n3/8 planes: the convergence table.
"""
from __future__ import annotations

import numpy as np


def klist(nk: int = 16, m: int = 8, half: int = 3, n3: int = 32) -> dict:
    """``k`` (N, 3) reduced; ``level`` (N,) in-plane level of each point (0 = coarse, outside the box); ``plane`` (N,);
    ``cell`` (N,) Chebyshev distance in cells of the point's coarse cell from the tube."""
    h = 1.0 / nk
    levels = [2 ** i for i in range(int(np.log2(m)) + 1)]
    j = np.arange(-m // 2, m // 2)
    sub = np.stack(np.meshgrid(j, j, indexing="ij"), -1).reshape(-1, 2)
    level_of_sub = np.array([min(lv for lv in levels if np.all(s % (m // lv) == 0)) for s in sub])
    uv, lev, cell = [], [], []
    for a in range(nk):
        for b in range(nk):
            p = np.array([a, b]) * h
            c = int(np.abs(np.rint(((p + 0.5) % 1.0 - 0.5) / h)).max())
            if c <= half:
                uv.append(p[None] + sub * h / m)
                lev.append(level_of_sub)
                cell.append(np.full(len(sub), c))
            else:
                uv.append(p[None])
                lev.append(np.array([0]))
                cell.append(np.array([c]))
    uv, lev, cell = np.concatenate(uv), np.concatenate(lev), np.concatenate(cell)
    for lv in levels:                                            # each level tiles the plane exactly once
        area = (lev == 0).sum() * h ** 2 + ((lev >= 1) & (lev <= lv)).sum() * (h / lv) ** 2
        assert abs(area - 1.0) < 1e-12, (lv, area)
    u3 = np.arange(n3) / n3
    k = np.concatenate([np.column_stack([uv[:, 0] + t, uv[:, 1] + t, np.full(len(uv), t)]) % 1.0 for t in u3])
    return dict(k=k, level=np.tile(lev, n3), plane=np.repeat(np.arange(n3), len(uv)), cell=np.tile(cell, n3),
                nk=nk, M=m, half=half, n3=n3)


def integrate(c2: np.ndarray, kl: dict) -> tuple[list, list, np.ndarray]:
    r"""dtheta(m, n3) = (1/n3) sum_{planes} [sum_{coarse} c2 h^2 + sum_{box, level <= m} c2 (h/m)^2].

    ``c2`` (N, ...) on the points of ``kl`` (the output of ``klist``). Returns (levels, plane counts, table with
    shape (levels, plane counts) + c2.shape[1:]); entries whose points are missing (NaN) stay NaN.
    """
    level, plane = kl["level"], kl["plane"]
    nk, m, n3_all = int(kl["nk"]), int(kl["M"]), int(kl["n3"])
    h = 1.0 / nk
    levels = [2 ** i for i in range(int(np.log2(m)) + 1)]
    n3s = [n3_all // 2 ** i for i in range(4) if n3_all // 2 ** i >= 2]
    table = np.full((len(levels), len(n3s)) + c2.shape[1:], np.nan)
    for a, lv in enumerate(levels):
        for b, n3 in enumerate(n3s):
            on = plane % (n3_all // n3) == 0
            out, box = on & (level == 0), on & (level >= 1) & (level <= lv)
            if np.isnan(c2[out | box].real).any():
                continue
            table[a, b] = ((c2[out] * h ** 2).sum(0) + (c2[box] * (h / lv) ** 2).sum(0)).real / n3
    return levels, n3s, table
