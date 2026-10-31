"""Figures for the 2026-10-12 meeting note (notes/meetings/2026-10-12.md).

Reads the saved symmetry tables, the nk-sweep summaries and the Gamma2- force data; writes PNGs next to this file.
Run with the `axion` env: python make_figures.py
"""
import runpy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# ---- settings ----
OUT = Path("/Users/treycole/Repos/axion-phonon/notes/meetings/2026-10-12_figures")
SYM = Path("/Users/treycole/Repos/axion-phonon/archive/diagnostics/2026-10-05_ws_tolerance_window")
MODESYM = Path("/Users/treycole/Repos/axion-phonon/archive/diagnostics/2026-10-12_yio_lda_modes")
GROUPS = Path("/Users/treycole/Repos/axion-phonon/archive/diagnostics/2026-09-29_centre_symmetrization/logs/point_groups.npz")
DYNMAT = Path("/Users/treycole/Repos/axion-phonon/archive/diagnostics/2026-10-12_yio_lda_modes/dynmat/gm2_dynmat.py")
GM2 = Path("/Users/treycole/Repos/axion-phonon/data/Y2Ir2O7/phonon/GM2-")
SWEEPS = {  # (mode, functional): nk_sweep summary.csv
    (1, "LDA"): GM2 / "mode1/Q_0.01A/lda/w90_trial02_212982/nk_sweep_tif_shared_ws/summary.csv",
    (2, "LDA"): GM2 / "mode2/Q_0.01A/lda/w90_trial02_212994/nk_sweep_tif_shared_ws/summary.csv",
    (3, "LDA"): GM2 / "mode3/Q_0.01A/lda/w90_trial02_212995/nk_sweep_tif_shared_ws/summary.csv",
    (1, "PBE"): GM2 / "mode1/Q_0.01A/pbe/w90_trial02_212996/nk_sweep_tif_shared_ws/summary.csv",
    (3, "PBE"): GM2 / "mode3/Q_0.01A/pbe/w90_trial02_213002/nk_sweep_tif_shared_ws/summary.csv",
}
# YIO base curvature symmetry, one saved table per step of the note's §2 table: (label, file, key)
STEPS = [
    ("last week\n(PBE,\nunpatched)", "control_192344_2026-10-10.npz", "tol_1e-03"),
    ("QE #16\npatched;\nedge on band", "base_full_group_YP_210600.npz", "tol_1e-03"),
    ("window\nedges\nin gaps", "base_full_group_YP_frozfix_212616.npz", "tol_1e-03"),
    ("LDA", "base_full_group_YL2_212615.npz", "tol_1e-03"),
    ("tie tol.\n5e-3 Å", "base_full_group_YL2_212615_tol5e-3.npz", "tol_5e-03"),
    ("full-prec.\nH, no cut", "base_full_group_YL2_212615_hchk_nocut_tol5e-3.npz", "table"),
]
MODE_TABLES = {"base (48 ops)": SYM / "base_full_group_YL2_212615_hchk_nocut_tol5e-3.npz",
               "mode 1 (24 ops)": MODESYM / "mode1_full_group_212982.npz",
               "mode 2 (24 ops)": MODESYM / "mode2_full_group_212994.npz",
               "mode 3 (24 ops)": MODESYM / "mode3_full_group_212995.npz"}
# DFT eigenvalue symmetry, worst |Δε| in μeV (log 10-07 §5, log 10-08 §12.4, log 10-12 §12)
EIG = [("PBE, unpatched QE", 7590.0), ("PBE, patched (#16)", 137.0), ("LDA base", 1.18), ("LDA modes 1–3", 1.0)]
# ΔE / (−½ F·Δr): YIO LDA (log 10-12 §2), MBT LDA after #19 (log 10-11 §5)
EVF = [("YIO 1 (Y)", 1.6133, 1.6142), ("YIO 2 (Ir)", 1.8093, 1.8096), ("YIO 3 (O)", 0.5528, 0.5486),
       ("MBT 1 (Mn)", 0.0245, 0.0197), ("MBT 2", 0.1138, 0.1099), ("MBT 3", 0.0757, 0.0749), ("MBT 4", 0.0470, 0.0467)]
# ------------------

# reference categorical palette (dataviz skill), light mode; text in ink tokens, never series colors
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
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


# ---- 1. YIO symmetry: DFT eigenvalues and the curvature path ----
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"width_ratios": [1, 1.9]})
labels, vals = zip(*EIG)
y = np.arange(len(vals))[::-1]
a.barh(y, vals, color=C[0], height=0.55)
a.set_xscale("log"); a.set_yticks(y, labels); a.set_xlim(0.3, 3e4)
a.set_xlabel("worst |ε(gk) − ε(k)|  (μeV)"); a.grid(axis="y", visible=False)
for yi, v in zip(y, vals):
    a.text(v * 1.3, yi, ("≤ 1 (print)" if v == 1.0 else f"{v:g}"), va="center", color=INK, fontsize=9)
a.set_title("(a) DFT eigenvalues", loc="left", color=INK)

med = np.array([np.median(np.load(SYM / f, allow_pickle=True)[k][1:], axis=0) * 100 for _, f, k in STEPS])
x = np.arange(len(STEPS))
markers = ["o", "s", "^", "D"]
for j, part in enumerate(PARTS):
    b.plot(x, med[:, j], marker=markers[j], color=C[j], label=part)
    nudge = {"internal": 0.6, "external": 1.6}.get(part, 1.0)  # separate the two end labels near 1e-2
    b.text(x[-1] + 0.15, med[-1, j] * nudge, f"{part} {med[-1, j]:.2g}%", va="center", color=INK, fontsize=9)
b.set_yscale("log"); b.set_xticks(x, [s for s, _, _ in STEPS], fontsize=8.5); b.set_xlim(-0.3, len(x) + 0.6)
b.set_ylabel("curvature residual, median over 47 ops (%)"); b.grid(axis="x", visible=False)
b.legend(loc="lower left", ncols=4, fontsize=9)
b.set_title("(b) Berry curvature, YIO base", loc="left", color=INK)
fig.tight_layout(); save(fig, "yio_symmetry.png")

# ---- 2. Production residual per operation: which operations carry the cross part ----
groups = np.load(GROUPS, allow_pickle=True)
ops = {"base (48 ops)": groups["base_ops"], "mode 1 (24 ops)": groups["mode1_ops"],
       "mode 2 (24 ops)": groups["mode2_ops"], "mode 3 (24 ops)": groups["mode3_ops"]}
fig, axs = plt.subplots(1, 4, figsize=(11, 3.3), sharey=True)
for ax, (name, path) in zip(axs, MODE_TABLES.items()):
    t = np.load(path, allow_pickle=True)["table"][1:] * 100
    keeps = np.array([abs(abs(r[2, 2]) - 1) < 1e-8 for r in ops[name]])[1:]
    rng = np.random.default_rng(0)
    for j in range(4):
        for sel, col, lab in ((keeps, C[0], "keeps z axis"), (~keeps, C[1], "moves z axis")):
            ax.scatter(j + rng.uniform(-0.18, 0.18, sel.sum()), np.maximum(t[sel, j], 1e-5), s=16, color=col,
                       edgecolor="#fcfcfb", linewidth=0.6, label=lab if j == 0 else None)
    ax.set_yscale("log"); ax.set_xticks(range(4), PARTS, rotation=30, fontsize=8.5); ax.grid(axis="x", visible=False)
    ax.set_title(name, loc="left", color=INK, fontsize=10)
axs[0].set_ylabel("residual per operation (%)")
axs[0].legend(loc="upper left", fontsize=8.5)
fig.suptitle("YIO LDA, production settings: the parts' residuals sit on operations that move the z axis; the total does not",
             x=0.01, ha="left", color=INK, fontsize=10.5)
fig.tight_layout(); save(fig, "yio_symmetry_per_operation.png")

# ---- 3. dθ/dQ vs nk ----
fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
for ax, mode, atoms in zip(axs, (1, 2, 3), ("Y", "Ir", "O 48f")):
    for func, col in (("LDA", C[0]), ("PBE", C[1])):
        if (mode, func) not in SWEEPS:
            continue
        d = np.genfromtxt(SWEEPS[(mode, func)], delimiter=",", names=True)
        ax.plot(d["nk"], d["dtheta_base_real"], marker="o", color=col, label=f"{func}, base end")
        ax.plot(d["nk"], d["dtheta_mode_real"], marker="o", ls="--", mfc="#fcfcfb", color=col, label=f"{func}, mode end")
    ax.set_xticks([4, 6, 8, 10, 12]); ax.set_xlabel("nk (nk³ mesh)")
    ax.set_title(f"mode {mode} ({atoms})", loc="left", color=INK)
axs[0].set_ylabel("dθ/dQ (rad/Å)")
axs[0].legend(fontsize=8.5, loc="lower right")
fig.tight_layout(); save(fig, "dtheta_vs_nk.png")

# ---- 4. Energy vs force check ----
fig, ax = plt.subplots(figsize=(7.5, 3.0))
names, ratios = [n for n, *_ in EVF], [e / f for _, e, f in EVF]
cols = [C[0] if n.startswith("YIO") else C[1] for n in names]
ax.axhline(1, color=INK2, lw=1)
ax.scatter(range(len(EVF)), ratios, color=cols, s=50, zorder=3, edgecolor="#fcfcfb", linewidth=1)
for i, r in enumerate(ratios):
    ax.text(i, r + 0.015, f"{r:.3f}", ha="center", va="bottom", color=INK, fontsize=8.5)
ax.set_xticks(range(len(EVF)), names, fontsize=8.5); ax.set_ylim(0.9, 1.33); ax.grid(axis="x", visible=False)
ax.set_ylabel("ΔE / (−½ F·Δr)")
ax.text(len(EVF) - 0.6, 1.17, "MBT before QE #19: ΔE < 0 (ratio ≈ −130 at Q = 0.002 Å)", ha="right", va="top", color=INK2, fontsize=8.5)
ax.set_title("Energy vs force, Q = 0.01 Å (LDA+U)", loc="left", color=INK)
fig.tight_layout(); save(fig, "energy_vs_force.png")

# ---- 5. Physical Γ2- modes ----
g = runpy.run_path(str(DYNMAT))
Cm, V, w2, dth = g["C"], g["V"], g["w2"], np.array(g["DTHETA_DQ"])
freq = np.sqrt(w2) * g["conv"] / (2 * np.pi * 2.99792458e10)
share = np.array([[np.sum((Cm[i, n] * V[i]) ** 2) for i in range(3)] for n in range(3)])
share /= share.sum(1, keepdims=True)
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 2.8), gridspec_kw={"width_ratios": [1.4, 1]})
yy = np.arange(3)[::-1]
left = np.zeros(3)
for i, (atom, col) in enumerate(zip(("Y", "Ir", "O 48f"), C)):
    a.barh(yy, share[:, i], left=left, color=col, height=0.55, edgecolor="#fcfcfb", linewidth=2, label=atom)
    for yi, l, s in zip(yy, left, share[:, i]):
        if s > 0.12:
            a.text(l + s / 2, yi, f"{atom} {s:.2f}", ha="center", va="center", color="#fcfcfb", fontsize=8.5)
    left += share[:, i]
a.set_yticks(yy, [f"{f:.0f} cm⁻¹" for f in freq]); a.set_xlim(0, 1); a.grid(visible=False)
a.set_xlabel("share of |Δr|²"); a.legend(ncols=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
dq = Cm.T @ dth
b.barh(yy, dq, color=C[0], height=0.55)
b.axvline(0, color=INK2, lw=1)
for yi, v in zip(yy, dq):
    b.text(v + (0.0002 if v > 0 else -0.0002), yi, f"{v:+.5f}", va="center", ha="left" if v > 0 else "right",
           color=INK, fontsize=8.5)
b.set_yticks(yy, [f"{f:.0f} cm⁻¹" for f in freq]); b.set_xlim(-0.0060, 0.0030); b.grid(axis="y", visible=False)
b.set_xlabel("dθ/dq  (rad Å⁻¹ amu⁻½; sign of each mode arbitrary)")
fig.tight_layout(); save(fig, "physical_modes.png")
