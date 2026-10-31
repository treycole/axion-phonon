"""Figures for the 2026-10-27 meeting note (notes/meetings/2026-10-27.md).

New plots of results that so far existed only as tables: MBT tube k convergence and the Q dependence of mode 3,
the MBT window study, and the curvature symmetry with atom-pinned positions (YIO base + modes, MBT modes).
Also copies the figures made since 10-21 by their notebooks. Writes PNGs next to this file.
Run with the `axion` env: python make_figures.py
"""
import re
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- settings ----
REPO = Path("/Users/treycole/Repos/axion-phonon")
OUT = REPO / "notes/meetings/2026-10-27_figures"
T2 = REPO / "data/MnBi2Te4/archive/phonon/T2-"
TUBES = {  # label: production tube folder (holds box_H5/ and convergence_M8_N16/)
    "1 (Mn)": T2 / "mode1/Q_0.01A/lda/dtheta_tube_conv1e-10_214630",
    "2 (Bi)": T2 / "mode2/Q_0.01A/lda/dtheta_tube_conv1e-10_214633",
    "3 (Te1)": T2 / "mode3/Q_0.01A/lda/dtheta_tube_conv1e-10_214627",
    "4 (Te2)": T2 / "mode4/Q_0.01A/lda/dtheta_tube_conv1e-10_214636",
}
TUBE_M3_Q03 = T2 / "mode3/Q_0.03A/lda/dtheta_tube_conv1e-10_214980"
WINDOWS = T2 / "mode3/Q_0.01A/lda/window_study_2026-10-22"     # <set>/tube_dtheta_tube.log
WINDOW_PROD = (10.278, 15.169)                                   # production frozen / outer top (eV)
VBM = 7.268                                                      # eV, base NSCF 214290 (setup_mbt_window_study.sh)
YSYM = REPO / "archive/diagnostics/2026-10-26_yio_symmetry_conv1e-10"
YIO_MODELS = {"base 214249": "base_{}_214249.npz", "mode 1 214147": "mode1_{}_214147.npz",
              "mode 2 214151": "mode2_{}_214151.npz", "mode 3 214155": "mode3_{}_214155.npz"}
MSYM = REPO / "archive/diagnostics/2026-10-22_mbt_mode_symmetry"   # mode{n}_symmetry.npz
COPIES = {  # meeting copy: source (made by the notebook named in README.md)
    "c2_distribution.png": REPO / "paper/figures/c2_distribution.png",
    "mbt_c2_patterns.png": T2 / "c2_distribution_modes.png",
    "yio_c2_patterns.png": REPO / "data/Y2Ir2O7/phonon/GM2-/c2_distribution_modes.png",
    "pbe_vs_lda_bands.png": REPO / "data/pbe_vs_lda_bands.png",
}
FLOOR = 1e-5   # %; residuals below this are drawn at the floor (exact to machine precision)
# ------------------

# reference categorical palette (dataviz skill), light mode; text in ink tokens, never series colors
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cd1"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
PARTS = ("internal", "cross", "external", "total")
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.dpi": 160,
                     "legend.frameon": False, "lines.linewidth": 1.6, "lines.markersize": 6})


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight")
    print("wrote", OUT / name)


def converged(tube):
    """Box half-width 5 (m 4, 8 planes) + [(m 8, 16 planes) - (m 4, 8 planes)] at half-width 3; (base, mode) end."""
    box = np.load(tube / "box_H5/dtheta_convergence.npz")
    conv = np.load(tube / "convergence_M8_N16/dtheta_convergence.npz")
    b = box["dtheta"][list(box["halves"]).index(5), list(box["levels"]).index(4), list(box["n3"]).index(8)]
    lv, n3 = list(conv["levels"]), list(conv["n3"])
    d = conv["dtheta"]
    return b + d[lv.index(8), n3.index(16)] - d[lv.index(4), n3.index(8)]


# ---------------------------------------------------------------- 1. MBT tube k convergence + mode 3 vs Q
fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
for i, (label, tube) in enumerate(TUBES.items()):
    ref = converged(tube)[0]
    conv = np.load(tube / "convergence_M8_N16/dtheta_convergence.npz")
    box = np.load(tube / "box_H5/dtheta_convergence.npz")
    lv, n3 = list(conv["levels"]), list(conv["n3"])
    rel = lambda v: 100 * (v / ref - 1)
    ax[0, 0].plot(conv["levels"], rel(conv["dtheta"][:, n3.index(16), 0]), "o-", color=C[i], label=f"pattern {label}")
    ax[0, 1].plot(conv["n3"], rel(conv["dtheta"][lv.index(8), :, 0]), "o-", color=C[i])
    h = box["halves"]
    ax[1, 0].plot(h, rel(box["dtheta"][:, list(box["levels"]).index(4), list(box["n3"]).index(8), 0]), "o-",
                  color=C[i])
for a, xl, t in ((ax[0, 0], "in-plane refinement level m (16 k_z planes, box half-width 3)", "(a) in-plane level"),
                 (ax[0, 1], "k_z planes (m = 8, box half-width 3)", "(b) k_z planes"),
                 (ax[1, 0], "refined-box half-width (cells; m = 4, 8 planes)", "(c) refined box")):
    a.set_yscale("symlog", linthresh=0.1)
    a.axhspan(-0.2, 0.2, color=GRID, lw=0)
    a.axhline(0, color=INK2, lw=0.8)
    a.set_xlabel(xl)
    a.set_ylabel("∂θ/∂Q vs converged (%)")
    a.set_title(t, loc="left", color=INK)
ax[0, 0].set_xscale("log", base=2); ax[0, 1].set_xscale("log", base=2)
ax[0, 0].legend(fontsize=8, loc="lower right")
ax[1, 0].text(2.05, 0.13, "±0.2%", color=INK2, fontsize=8)
ax[1, 0].text(4.6, 0.3, "production:\nhalf-width 3, m 4, 8 planes", color=INK2, fontsize=8, ha="right")

# (d) mode 3: both estimates at Q 0.01 and 0.03 with a + b Q^2 fits
a = ax[1, 1]
qs = np.array([0.01, 0.03])
ends = np.array([converged(TUBES["3 (Te1)"]), converged(TUBE_M3_Q03)])          # (Q, end)
qq = np.linspace(0, 0.035, 50)
for e, name, col in ((0, "undistorted-structure estimate", C[0]), (1, "displaced-structure estimate", C[1])):
    b, a0 = np.polyfit(qs ** 2, ends[:, e], 1)
    a.plot(qq, a0 + b * qq ** 2, "--", color=col, lw=1)
    a.plot(qs, ends[:, e], "o", color=col, label=f"{name}: Q → 0 {a0:.3f}")
    print(f"mode 3 {name}: {ends[:, e]} -> Q^2 fit intercept {a0:.4f}")
a.set_xlabel("Q (Å)"); a.set_ylabel("∂θ/∂Q (rad/Å)")
a.set_title("(d) pattern 3 (Te1): the two estimates vs Q", loc="left", color=INK)
a.legend(fontsize=8, loc="upper left")
fig.suptitle("MnBi₂Te₄ T₂⁻, ∂θ/∂Q on the k_z tube (conv_thr 1e-10, pinned, Q = 0.01 Å unless stated)", color=INK)
fig.tight_layout()
save(fig, "mbt_dtheta_convergence.png")


# ---------------------------------------------------------------- 2. MBT window study
def tube_m4(log):
    """(base, mode) end at m 4, 8 planes from a dtheta_tube.log."""
    rows = re.findall(r"^m=4\s+(\S+)", log.read_text(), re.M)
    return float(rows[0]), float(rows[1])


sets = {}
for d in sorted(WINDOWS.iterdir()):
    froz, outer = (float(x[1:].replace("p", ".")) for x in d.name.split("_"))
    sets[froz, outer] = tube_m4(d / "tube_dtheta_tube.log")
sets[WINDOW_PROD] = tube_m4(TUBES["3 (Te1)"] / "dtheta_tube.log")
prod = sets[WINDOW_PROD][0]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for a, fixed, idx, xl in ((ax[0], WINDOW_PROD[1], 0, "frozen-window top − VBM (eV); outer top VBM+7.90"),
                          (ax[1], WINDOW_PROD[0], 1, "outer-window top − VBM (eV); frozen top VBM+3.01")):
    pts = sorted((k[idx], v[0]) for k, v in sets.items() if k[1 - idx] == fixed)
    x, y = np.array(pts).T
    a.axhspan(-1, 1, color=GRID, lw=0)
    a.axhline(0, color=INK2, lw=0.8)
    a.plot(x - VBM, 100 * (y / prod - 1), "o-", color=C[0])
    a.plot(WINDOW_PROD[idx] - VBM, 0, "o", color=C[1], label="production")
    a.set_xlabel(xl)
ax[0].set_ylabel("∂θ/∂Q vs production (%)")
ax[0].text(1.6, -2.4, "±1%", color=INK2, fontsize=8)
ax[0].legend(fontsize=8, loc="upper right")
ax[0].set_title("(a) frozen top", loc="left", color=INK); ax[1].set_title("(b) outer top", loc="left", color=INK)
fig.suptitle("MnBi₂Te₄ pattern 3, Q = 0.01 Å: window study (undistorted-structure estimate, m 4, 8 planes)", color=INK)
fig.tight_layout()
save(fig, "mbt_window_study.png")

# ---------------------------------------------------------------- 3. symmetry with atom-pinned positions
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
a = ax[0]
rng = np.random.default_rng(0)
for j, part in enumerate(PARTS):
    for s, (kind, mk) in enumerate((("full_group", "o"), ("pinned", "o"))):
        vals = np.concatenate([np.load(YSYM / f.format(kind))["table"][1:, j] * 100 for f in YIO_MODELS.values()])
        x = j * 3 + s + rng.uniform(-0.25, 0.25, vals.size)
        a.scatter(x, np.maximum(vals, FLOOR), s=10, color=C[s], alpha=0.6, edgecolor="none",
                  label=("Wannier centres" if s == 0 else "atom-pinned") if j == 0 else None)
# x carries no value: each part has two clusters (Wannier-centre / pinned positions), jittered so the dots don't overlap
a.set_xticks([j * 3 + s for j in range(4) for s in (0, 1)], ["centres", "pinned"] * 4, fontsize=8)
for j, part in enumerate(PARTS):
    a.text(j * 3 + 0.5, -0.12, part, transform=a.get_xaxis_transform(), ha="center", va="top", color=INK)
a.grid(False, axis="x")
a.set_yscale("log"); a.set_ylim(FLOOR / 2, 2)
a.set_ylabel("relative violation (%)")
a.set_title("(a) Y₂Ir₂O₇: base (47 ops) + modes 1–3 (23 ops each)", loc="left", color=INK)
a.legend(fontsize=8, loc="upper right", markerscale=2)
a.text(-0.4, FLOOR * 1.4, "floor = exact (< 1e-5 %)", color=INK2, fontsize=8)

a = ax[1]
for m in range(1, 5):
    z = np.load(MSYM / f"mode{m}_symmetry.npz", allow_pickle=True)
    g, names = z["generic"] * np.r_[1, 100, 100, 100, 100], list(z["names"])
    ops = [i for i, n in enumerate(names) if n not in ("E", "-E") and not n.startswith("control")]
    for j in range(4):
        x = j * 2 + (m - 2.5) * 0.15
        a.scatter(np.full(len(ops), x), np.maximum(g[ops, j + 1], FLOOR), s=12, color=C[m - 1], edgecolor="none",
                  label=f"pattern {m}" if j == 0 else None)
        a.scatter(x + 0.9, g[names.index("control C4z"), j + 1], marker="x", s=18, color=C[m - 1])
# x carries no value: per part, the four patterns side by side (−3′m′ operations), then the C4z control (×)
a.set_xticks([j * 2 + d for j in range(4) for d in (0, 0.9)], ["−3′m′", "C4z"] * 4, fontsize=8)
for j, part in enumerate(PARTS):
    a.text(j * 2 + 0.45, -0.12, part, transform=a.get_xaxis_transform(), ha="center", va="top", color=INK)
a.grid(False, axis="x")
a.set_yscale("log"); a.set_ylim(FLOOR / 2, 100)
a.set_title("(b) MnBi₂Te₄ modes: −3′m′ (dots), C4z control (×)", loc="left", color=INK)
a.legend(fontsize=8, loc="center left", bbox_to_anchor=(0, 0.3), ncol=2)
fig.suptitle("Curvature symmetry per operation (conv_thr 1e-10 runs, H from .chk, WS tol 5e-3 Å)", color=INK)
fig.tight_layout()
save(fig, "symmetry_pinned.png")

for name, src in COPIES.items():
    shutil.copy2(src, OUT / name)
    print("copied", src, "->", OUT / name)
