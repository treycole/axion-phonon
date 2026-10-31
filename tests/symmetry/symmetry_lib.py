"""Symmetry of the MnBi2Te4 models of the paper dtheta: the functions behind run_mbt_symmetry.py (the paper's table,
tab:symmetry-mbt) and test_mbt_symmetry.py (the pass/fail thresholds).

The group of each structure is read from its pw.x output (axion_phonon.io.qe.read_symmetry: QE's Cartesian matrices
and 'Time Reversal' flags), not typed in. The undistorted structure also has the anti-translation {1'|1/2,1/2,1/2},
which exchanges the septuple layers; QE does not list it, because it disables fractional translations for the
magnetic supercell. With anti_translation, each listed operation is added again with the opposite priming, so the
base is tested against all 24 operations of its magnetic point group in k space. Test 1 checks that the added ones
are symmetries (if they were not, the DFT eigenvalues would differ by meV).

For one base/pattern pair, built exactly as in the production c2 jobs (axion_phonon.wannier.pair_models: H and A(R)
from the .chk, folded on the base centres' nearest atoms, one tau for both), and for each of the two models:
1. DFT eigenvalues on the Wannier mesh (.eig): e_n(g k) = e_n(k), k -> -g k for primed operations; the bands below the
   frozen-window top at every k, and all bands.
2. Interpolated bands at random k and at k within a small radius of the k_z axis (the tube where dtheta/dQ sits).
3. Curvature: Tr F is PT-forced to zero, so the gauge-invariant tensor T_ab = Re Tr[F_a F_b] of the occupied bands
   (F_a = eps_abc Omega_bc / 2, Cartesian), T(g k) = R T(k) R^T, for the internal, cross, external parts and the total.
Controls, which must fail: C4z (not a symmetry of the lattice) and unprimed inversion k -> -k (the operation the
displacement breaks), each where it is not in the group.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from axion_phonon import curvature, wannier
from axion_phonon.io import qe, wannier90

PARTS = ("internal", "cross", "external", "total")
# Pass/fail (test_mbt_symmetry.py, and the PASS/FAIL lines of run_mbt_symmetry.py). Residuals are maxima over the
# operations other than the identity; curvature ones are relative to the largest |T_ab| of the same part (10-31 local
# test, base and pattern 3: total 1.9e-4, cross 4.4e-4, external 1.5e-3, internal 0, DFT 0.06 ueV, bands < 1e-4 meV).
THRESHOLDS = {
    "dft_ueV": 1.0,                 # DFT eigenvalues below the frozen-window top
    "bands_meV": 1e-3,              # interpolated bands
    "internal": 1e-6,
    "cross": 5e-3,
    "external": 5e-3,
    "total": 1e-3,
    "trace": 1e-4,                  # max |Tr F| / |F| scale: PT forces Tr F = 0
    "control_min": 5e-3,            # each control's total residual at random k, at least
    "control_over_residual": 10,    # and, at every k set, at least this many times the largest residual of the group
}
CONTROLS = (("C4z", np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]]), False),
            ("-E unprimed", -np.eye(3), False))


def group(scf_out: Path, anti_translation: bool):
    """(operations (n, 3, 3) Cartesian, primed (n,), names) of a structure, from its pw.x output; with
    anti_translation, every operation again with the opposite priming (its product with {1'|t})."""
    s = qe.read_symmetry(scf_out)
    ops, primed = s["cart"], s["time_reversal"]
    names = [n + (" T" if p else "") for n, p in zip(s["names"], primed)]
    if not (np.allclose(ops[0], np.eye(3)) and not primed[0]):
        raise ValueError(f"{scf_out}: the first operation is not the identity")
    if anti_translation:
        ops, primed = np.concatenate([ops, ops]), np.concatenate([primed, ~primed])
        names = names + [n + " {1'|t}" for n in names]
    keys = [(tuple(np.round(R, 6).ravel()), bool(p)) for R, p in zip(ops, primed)]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{scf_out}: repeated operations (anti_translation={anti_translation})")
    return ops, primed, names


def controls(ops, primed):
    """(name, operation, primed) of CONTROLS that are not in the group."""
    return [c for c in CONTROLS if not any(np.allclose(R, c[1]) and p == c[2] for R, p in zip(ops, primed))]


def to_reduced(K, reciprocal, op, anti):
    """Reduced k of g K (or -g K for a primed operation); K Cartesian rows, reciprocal rows b_i."""
    return ((-1.0 if anti else 1.0) * K @ op.T) @ np.linalg.inv(reciprocal)


def dft_test(run: dict, ops, primed, frozen_top: float):
    """(bands below frozen_top at every k, (n_ops, 2) max |e_n(g k) - e_n(k)| in micro-eV: those bands, all bands)."""
    kred, eig, grid = run["kpts"], run["energies"], np.array(run["mp_grid"])
    index = {tuple(np.rint(k * grid).astype(int) % grid): i for i, k in enumerate(kred)}
    above = eig.max(0) > frozen_top
    if not above.any():
        raise ValueError(f"no band rises above the frozen-window top {frozen_top} eV")
    n_low = int(np.argmax(above))
    out = []
    for op, anti in zip(ops, primed):
        n = to_reduced(kred @ run["reciprocal"], run["reciprocal"], op, anti) * grid
        assert np.allclose(n, np.rint(n), atol=1e-6), "operation does not map the mesh onto itself"
        j = [index[tuple(v % grid)] for v in np.rint(n).astype(int)]
        d = np.abs(eig[j] - eig)
        out.append((d[:, :n_low].max() * 1e6, d.max() * 1e6))
    return n_low, np.array(out)


def curvature_tensor(model, pos, k, reciprocal, n_occ: int, batch: int):
    """T_ab = Re Tr[F_a F_b] (Cartesian, (nk, 3, 3)) per part, and max |Tr F_a| over k."""
    rinv = np.linalg.inv(reciprocal)
    out = {name: [] for name in PARTS}
    trace = 0.0
    for s in range(0, len(k), batch):
        c = curvature.connection(curvature.fields(model, pos, k[s:s + batch]), n_occ)
        om = {"internal": curvature.omega_internal(c), "cross": curvature.omega_cross(c),
              "external": curvature.omega_external(c)}
        om["total"] = om["internal"] + om["cross"] + om["external"]
        for name in out:
            cart = np.einsum("ia,jb,ab...->ij...", rinv, rinv, om[name], optimize=True)   # reduced -> Cartesian
            F = np.stack((cart[1, 2], cart[2, 0], cart[0, 1]))                           # (3, nk, v, v)
            out[name].append(np.einsum("akmn,bknm->kab", F, F, optimize=True).real)
            if name == "total":
                trace = max(trace, np.abs(np.einsum("akmm->ak", F)).max())
    return {name: np.concatenate(v) for name, v in out.items()}, trace


def k_sets(reciprocal, n_generic: int, n_tube: int, tube_radius: float, seed: int):
    """Reduced k: n_generic random, and n_tube within tube_radius (1/A) of the k_z axis between Gamma and Z."""
    rinv = np.linalg.inv(reciprocal)
    rng = np.random.default_rng(seed)
    generic = rng.uniform(-0.23, 0.23, size=(n_generic, 3)) @ rinv
    axis = reciprocal.sum(0)                                        # b1 + b2 + b3, along c; Z at its half
    phi = rng.uniform(0, 2 * np.pi, n_tube)
    perp = np.stack([np.cos(phi), np.sin(phi), np.zeros(n_tube)], 1) * rng.uniform(0.3, 1, (n_tube, 1)) * tube_radius
    tube = (perp + rng.uniform(0, 1, (n_tube, 1)) * axis / 2) @ rinv
    return {"generic": generic, "tube": tube}


def model_tests(label, model, reciprocal, ops, primed, names, ctrl, ksets, n_occ, batch, log=print):
    """Bands and curvature tensor under every operation and control: {k set: (n_ops + n_controls, 5) rows of
    (max band change meV, internal, cross, external, total)}, and {k set: max |Tr F| / |F| scale}."""
    pos = curvature.position_terms(model, dtype=np.complex128)
    cases = list(zip(names, ops, primed)) + [("control " + c[0], c[1], c[2]) for c in ctrl]
    width = max(len(c[0]) for c in cases) + 1
    rows, trace = {}, {}
    for kset, k0 in ksets.items():
        K0 = k0 @ reciprocal
        e0 = np.linalg.eigvalsh(model.hamiltonian(k0, flatten_spin_axis=True))
        T0, tr0 = curvature_tensor(model, pos, k0, reciprocal, n_occ, batch)
        trace[kset] = tr0 / np.sqrt(np.abs(T0["total"]).max())
        log(f"\n   {label}, {kset} k ({len(k0)}): |T| scale "
            + ", ".join(f"{n} {np.abs(v).max():.3e}" for n, v in T0.items())
            + f"; max |Tr F| / |F| scale {trace[kset]:.1e} (PT-forced zero)")
        log(f"   {'operation':{width}s} {'W de (meV)':>11s} {'internal':>9s} {'cross':>9s} {'external':>9s} {'total':>9s}")
        r = []
        for name, op, anti in cases:
            kt = to_reduced(K0, reciprocal, op, anti)
            dE = np.abs(np.linalg.eigvalsh(model.hamiltonian(kt, flatten_spin_axis=True)) - e0).max() * 1e3
            T, _ = curvature_tensor(model, pos, kt, reciprocal, n_occ, batch)
            res = []
            for part in PARTS:
                pred = np.einsum("ab,kbc,dc->kad", op, T0[part], op)
                res.append(np.abs(T[part] - pred).max() / max(np.abs(pred).max(), np.abs(T[part]).max()))
            log(f"   {name:{width}s} {dE:11.4f} " + " ".join(f"{v:8.4%}" for v in res))
            r.append([dE] + res)
        rows[kset] = np.array(r)
        sym = rows[kset][1:len(ops)]
        log(f"   {label}, excl. identity ({len(ops) - 1} operations), {kset}: max W de {sym[:, 0].max():.4f} meV; "
            "median (max) "
            + ", ".join(f"{p} {np.median(sym[:, i + 1]):.4%} ({sym[:, i + 1].max():.4%})" for i, p in enumerate(PARTS))
            + "; controls total " + ", ".join(f"{c[0]} {v:.2%}" for c, v in zip(ctrl, rows[kset][len(ops):, 4])))
    return rows, trace


def run_pair(base_run: Path, mode_run: Path, base_scf: Path, mode_scf: Path, anti_translation: bool, prefix: str,
             n_occ: int, frozen_top: float, n_generic: int, n_tube: int, tube_radius: float, seed: int, batch: int,
             log=print) -> dict:
    """Every test on the base and the displaced model of one pair: {"base": results, "mode": results}, each a dict
    that save() writes and failures() checks. base_scf, mode_scf: the pw.x outputs the groups are read from;
    anti_translation: the base has {1'|t} (module docstring)."""
    base, mode = wannier90.read_run(base_run, prefix), wannier90.read_run(mode_run, prefix)
    mb, mm = wannier.pair_models(base, mode, log=log)
    return pair_results(base, mode, mb, mm, base_run, mode_run, base_scf, mode_scf, anti_translation, n_occ,
                        frozen_top, n_generic, n_tube, tube_radius, seed, batch, log)


def pair_results(base, mode, mb, mm, base_run, mode_run, base_scf, mode_scf, anti_translation, n_occ, frozen_top,
                 n_generic, n_tube, tube_radius, seed, batch, log=print) -> dict:
    """run_pair on runs already read (their kpts, energies) and models already built by pair_models."""
    reciprocal = base["reciprocal"]
    results, groups = {}, {}
    for label, run, folder, scf, anti in (("base", base, base_run, base_scf, anti_translation),
                                          ("mode", mode, mode_run, mode_scf, False)):
        ops, primed, names = group(scf, anti)
        ctrl = controls(ops, primed)
        groups[label] = (ops, primed, names, ctrl)
        n_low, dft = dft_test(run, ops, primed, frozen_top)
        log(f"{label}: {folder}; group from {scf}: {len(ops)} operations"
            f"{' (with the anti-translation partners)' if anti else ''}; 1. DFT .eig on the {tuple(run['mp_grid'])} "
            f"mesh: {n_low} bands below {frozen_top} eV at every k; max |de| over the {len(ops) - 1} operations: "
            f"{dft[1:, 0].max():.2f} ueV (those bands), {dft[1:, 1].max():.2f} ueV (all {run['energies'].shape[1]})")
        results[label] = {"names": np.array(names + ["control " + c[0] for c in ctrl]),
                          "primed": np.r_[primed, [c[2] for c in ctrl]], "n_ops": len(ops), "n_low": n_low,
                          "dft_ueV": dft}
    ksets = k_sets(reciprocal, n_generic, n_tube, tube_radius, seed)
    for label, model in (("base", mb), ("mode", mm)):
        ops, primed, names, ctrl = groups[label]
        rows, trace = model_tests(label, model, reciprocal, ops, primed, names, ctrl, ksets, n_occ, batch, log)
        results[label].update({kset: v for kset, v in rows.items()},
                              **{f"trace_{kset}": v for kset, v in trace.items()})
    return results


def save(results: dict, out: Path) -> None:
    for label, r in results.items():
        np.savez(out / f"symmetry_{label}.npz", **r)


def load(folder: Path) -> dict:
    return {label: dict(np.load(folder / f"symmetry_{label}.npz")) for label in ("base", "mode")}


def failures(r: dict) -> list[str]:
    """Every threshold of THRESHOLDS that one model's results (an entry of run_pair's output, or load()'s) violate."""
    t, n_ops, bad = THRESHOLDS, int(r["n_ops"]), []
    if r["dft_ueV"][1:, 0].max() > t["dft_ueV"]:
        bad.append(f"DFT eigenvalues: {r['dft_ueV'][1:, 0].max():.2f} ueV > {t['dft_ueV']}")
    for kset in ("generic", "tube"):
        rows = r[kset]
        sym, ctrl = rows[1:n_ops], rows[n_ops:]
        if sym[:, 0].max() > t["bands_meV"]:
            bad.append(f"{kset}: interpolated bands {sym[:, 0].max():.1e} meV > {t['bands_meV']}")
        for i, part in enumerate(PARTS):
            if sym[:, i + 1].max() > t[part]:
                bad.append(f"{kset}: {part} {sym[:, i + 1].max():.2e} > {t[part]}")
        if float(r[f"trace_{kset}"]) > t["trace"]:
            bad.append(f"{kset}: |Tr F| / |F| {float(r[f'trace_{kset}']):.1e} > {t['trace']}")
        for name, v in zip(r["names"][n_ops:], ctrl[:, 4]):
            if v < t["control_over_residual"] * sym[:, 4].max() or (kset == "generic" and v < t["control_min"]):
                bad.append(f"{kset}: {name} total {v:.2e} does not fail clearly (group max {sym[:, 4].max():.2e})")
    return bad
