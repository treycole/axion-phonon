"""Y2Ir2O7 symmetry table (paper tab:symmetry-yio) and Fig. 6(b) for the base and one displaced structure, with the
models of the production dtheta (axion_phonon.wannier.pair_models: H and A(R) from the .chk, pinned to the base
centres' nearest atoms, one tau for both). The log is the table; symmetry_base.npz / symmetry_mode.npz hold every
operation's residuals; the last lines say PASS or FAIL against THRESHOLDS (exit 1 on FAIL).

The YIO counterpart of run_mbt_symmetry.py. AIAO breaks time reversal, so the band-traced Berry curvature is not zero,
and the test is that pseudovector itself (as in test_yio_symmetry.py and the paper):
1. DFT eigenvalues of the lowest EIG_BANDS bands on the Wannier mesh (.eig): e_n(g k) = e_n(k), k -> -g k if primed.
2. Interpolated bands at N_GENERIC random k.
3. omega_a = Re Tr F_a (F_a = eps_abc Omega_bc / 2, Cartesian) of the N_OCC occupied bands, for the internal, cross
   and external parts and the total: omega(g k) = det(R) R omega(k), and -det(R) R omega(k) for primed g (k -> -g k).
The group of each structure is read from its pw.x output (QE lists the whole group for YIO). Controls, which must
fail: the fourfold C4z and S4z with whichever priming is not in the group. For a displaced structure, a control that is
in the undistorted group (C4z T: the operation the Gamma2- displacement breaks, at O(Q)) only has to exceed
CONTROL_OVER_RESIDUAL times the largest residual of the group ("broken"); the others are non-symmetries (~200%).

Run on the cluster with run_yio_symmetry.sh (it fills the settings below), or locally with BASE_RUN / MODE_RUN set to
the Wannier run folders (.chk, .eig, .win, .nnkp) and BASE_SCF / MODE_SCF to the pw.x outputs.
"""
import sys
import time
from pathlib import Path

import numpy as np

from axion_phonon import curvature, wannier
from axion_phonon.io import qe, wannier90

# ---- settings (run_yio_symmetry.sh fills these in)
BASE_RUN = Path("base")
MODE_RUN = Path("mode")
BASE_SCF = Path("base/Y2Ir2O7.scf.out")
MODE_SCF = Path("mode/Y2Ir2O7.scf.out")
PREFIX = "Y2Ir2O7"
N_OCC = 156
EIG_BANDS = 260          # lowest bands in test 1 (paper)
N_GENERIC = 7            # random k: Cartesian components uniform in +-0.23 1/A
SEED = 20260825          # as the 10-26 table: the same k
BATCH = 4                # k points per curvature call (memory)
OUT = Path(".")
# --------------
PARTS = ("internal", "cross", "external", "total")
THRESHOLDS = {"dft_ueV": 2.0, "bands_meV": 1e-3, "internal": 1e-5, "cross": 1e-3, "external": 1e-3, "total": 1e-4,
              "control_min": 0.5,   # maxima over the operations other than the identity; each non-symmetry control at least
              "control_over_residual": 10}   # each broken-symmetry control: at least this many times the largest residual
C4Z = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
S4Z = np.array([[0.0, 1, 0], [-1, 0, 0], [0, 0, -1]])
CONTROLS = (("C4z", C4Z, False), ("C4z T", C4Z, True), ("S4z", S4Z, False), ("S4z T", S4Z, True))


def group(scf_out):
    s = qe.read_symmetry(scf_out)
    ops, primed = s["cart"], np.asarray(s["time_reversal"], bool)
    assert np.allclose(ops[0], np.eye(3)) and not primed[0], f"{scf_out}: first operation is not the identity"
    names = [n + (" T" if p else "") for n, p in zip(s["names"], primed)]
    ctrl = [c for c in CONTROLS if not any(np.allclose(R, c[1]) and p == c[2] for R, p in zip(ops, primed))]
    return ops, primed, names, ctrl


def to_reduced(K, reciprocal, op, anti):
    return ((-1.0 if anti else 1.0) * K @ op.T) @ np.linalg.inv(reciprocal)


def dft_test(run, ops, primed):
    kred, eig, grid = run["kpts"], run["energies"][:, :EIG_BANDS], np.array(run["mp_grid"])
    index = {tuple(np.rint(k * grid).astype(int) % grid): i for i, k in enumerate(kred)}
    out = []
    for op, anti in zip(ops, primed):
        n = to_reduced(kred @ run["reciprocal"], run["reciprocal"], op, anti) * grid
        assert np.allclose(n, np.rint(n), atol=1e-6), "operation does not map the mesh onto itself"
        out.append(np.abs(eig[[index[tuple(v % grid)] for v in np.rint(n).astype(int)]] - eig).max() * 1e6)
    return np.array(out)


def trace_vector(model, pos, k, reciprocal):
    rinv = np.linalg.inv(reciprocal)
    out = {p: [] for p in PARTS}
    for s in range(0, len(k), BATCH):
        c = curvature.connection(curvature.fields(model, pos, k[s:s + BATCH]), N_OCC)
        om = {"internal": curvature.omega_internal(c), "cross": curvature.omega_cross(c),
              "external": curvature.omega_external(c)}
        om["total"] = om["internal"] + om["cross"] + om["external"]
        for p in PARTS:
            tr = np.trace(np.einsum("ia,jb,ab...->ij...", rinv, rinv, om[p], optimize=True), axis1=-2, axis2=-1).real
            out[p].append(np.stack((tr[1, 2], tr[2, 0], tr[0, 1]), axis=-1))
    return {p: np.concatenate(v) for p, v in out.items()}


def model_tests(label, model, reciprocal, ops, primed, names, ctrl, log):
    """(n_ops + n_controls, 5) rows: max band change (meV), then internal, cross, external, total residuals."""
    pos = curvature.position_terms(model, dtype=np.complex128)
    K0 = np.random.default_rng(SEED).uniform(-0.23, 0.23, size=(N_GENERIC, 3))
    k0 = K0 @ np.linalg.inv(reciprocal)
    e0 = np.linalg.eigvalsh(model.hamiltonian(k0, flatten_spin_axis=True))
    w0 = trace_vector(model, pos, k0, reciprocal)
    log(f"\n   {label}, {N_GENERIC} random k: |omega| scale "
        + ", ".join(f"{p} {np.abs(v).max():.3e}" for p, v in w0.items()))
    cases = list(zip(names, ops, primed)) + [("control " + c[0], c[1], c[2]) for c in ctrl]
    width = max(len(c[0]) for c in cases) + 1
    log(f"   {'operation':{width}s} {'W de (meV)':>11s} {'internal':>9s} {'cross':>9s} {'external':>9s} {'total':>9s}")
    rows = []
    for name, op, anti in cases:
        kt = to_reduced(K0, reciprocal, op, anti)
        dE = np.abs(np.linalg.eigvalsh(model.hamiltonian(kt, flatten_spin_axis=True)) - e0).max() * 1e3
        w = trace_vector(model, pos, kt, reciprocal)
        factor = (-1.0 if anti else 1.0) * np.linalg.det(op)
        res = []
        for p in PARTS:
            pred = factor * w0[p] @ op.T
            res.append(np.abs(w[p] - pred).max() / max(np.abs(np.concatenate((w[p], pred))).max(), 1e-30))
        log(f"   {name:{width}s} {dE:11.6f} " + " ".join(f"{v:8.4%}" for v in res))
        rows.append([dE] + res)
    rows = np.array(rows)
    sym = rows[1:len(ops)]
    log(f"   {label}, excl. identity ({len(ops) - 1} operations): max W de {sym[:, 0].max():.2e} meV; median (max) "
        + ", ".join(f"{p} {np.median(sym[:, i + 1]):.4%} ({sym[:, i + 1].max():.4%})" for i, p in enumerate(PARTS))
        + "; controls total " + ", ".join(f"{c[0]} {v:.1%}" for c, v in zip(ctrl, rows[len(ops):, 4])))
    return rows


def failures(r):
    n = r["n_ops"]
    sym, ctl = r["rows"][1:n], r["rows"][n:]
    bad = [f"DFT {r['dft_ueV'][1:].max():.2f} ueV"] if r["dft_ueV"][1:].max() > THRESHOLDS["dft_ueV"] else []
    if sym[:, 0].max() > THRESHOLDS["bands_meV"]:
        bad.append(f"bands {sym[:, 0].max():.2e} meV")
    for i, p in enumerate(PARTS):
        if sym[:, i + 1].max() > THRESHOLDS[p]:
            bad.append(f"{p} {sym[:, i + 1].max():.2e}")
    floor = THRESHOLDS["control_over_residual"] * sym[:, 4].max()
    for nm, v, broken in zip(r["names"][n:], ctl[:, 4], r["control_broken"]):
        if v < (floor if broken else THRESHOLDS["control_min"]):
            bad.append(f"{nm} ({'broken symmetry' if broken else 'non-symmetry'}) total {v:.2e}")
    return bad


if __name__ == "__main__":
    log = lambda s: print(s, flush=True)   # noqa: E731
    t0 = time.time()
    base, mode = wannier90.read_run(BASE_RUN, PREFIX), wannier90.read_run(MODE_RUN, PREFIX)
    mb, mm = wannier.pair_models(base, mode, log=log)
    results, base_keys = {}, set()
    for label, run, model, scf in (("base", base, mb, BASE_SCF), ("mode", mode, mm, MODE_SCF)):
        ops, primed, names, ctrl = group(scf)
        if label == "base":
            base_keys = {(tuple(np.round(R, 6).ravel()), bool(p)) for R, p in zip(ops, primed)}
        broken = np.array([(tuple(np.round(c[1], 6).ravel()), bool(c[2])) in base_keys for c in ctrl], bool)
        dft = dft_test(run, ops, primed)
        log(f"\n{label}: group from {scf}: {len(ops)} operations ({int(primed.sum())} primed); 1. DFT .eig, lowest "
            f"{EIG_BANDS} bands on the {tuple(run['mp_grid'])} mesh: max |de| {dft[1:].max():.2f} ueV")
        rows = model_tests(label, model, base["reciprocal"], ops, primed, names, ctrl, log)
        results[label] = {"names": np.array(names + ["control " + c[0] for c in ctrl]),
                          "primed": np.r_[primed, [c[2] for c in ctrl]], "n_ops": len(ops), "dft_ueV": dft,
                          "rows": rows, "control_broken": broken}
        np.savez(OUT / f"symmetry_{label}.npz", **results[label])
    bad = {label: failures(r) for label, r in results.items()}
    for label, b in bad.items():
        print(f"\n{label}: " + ("PASS" if not b else "FAIL\n  " + "\n  ".join(b)))
    print(f"\nall done in {time.time() - t0:.0f} s")
    sys.exit(1 if any(bad.values()) else 0)
