"""Our curvature against PythTB's, on the production models (port of
analysis/validation/external_curvature/check_vs_pythtb.ipynb, the old suite's primary gate).

- k-k block: omega_cross + omega_external (axion_phonon.curvature) equals PythTB's own per-k implementation,
  berry_curvature(include_external=True) - (include_external=False), itself validated against postw90's J0 + J1;
  the internal part equals PythTB's include_external=False. MnBi2Te4 and Y2Ir2O7 bases, non-Abelian (full
  n_occ x n_occ matrices: the MnBi2Te4 trace is zero by PT, the matrices are not).
- mixed k-beta block (no published implementation to compare with): on the MnBi2Te4 base/pattern-3 pair with a
  nonzero A_beta = i Lambda, the 4x4 assembly keeps the k-k block unchanged and its (l, beta) entries equal a
  hand assembly of all four cross terms, the mixed curl and -i[A_l, A_beta]; every block is Hermitian.
"""
import numpy as np
import pytest

from axion_phonon import axion, curvature

# ---- settings
K_TEST = np.vstack([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.25, 0.25, 0.25],
                    np.random.default_rng(7).uniform(0, 1, size=(3, 3))])
VS_PYTHTB = 1e-8              # max |ours - PythTB| / max |PythTB|
HAND = 1e-10                  # max |assembly - hand| / max |assembly|
# --------------
pytestmark = pytest.mark.data


@pytest.fixture(params=["mbt", "yio"])
def base(request):
    if request.param == "mbt":
        p = request.getfixturevalue("mbt_pair")
        return p["base_model"], p["n_occ"]
    p = request.getfixturevalue("yio_base")
    return p["model"], p["n_occ"]


def test_k_k_block_vs_pythtb(base):
    model, n_occ = base
    occ = np.arange(n_occ)
    ref_tot = model.berry_curvature(K_TEST, occ_idxs=occ, non_abelian=True, include_external=True)
    ref_int = model.berry_curvature(K_TEST, occ_idxs=occ, non_abelian=True, include_external=False)
    c = curvature.connection(curvature.fields(model, curvature.position_terms(model, dtype=np.complex128), K_TEST), n_occ)
    internal, ext = curvature.omega_internal(c), curvature.omega_cross(c) + curvature.omega_external(c)
    scale = np.abs(ref_tot).max()
    assert np.abs(internal - ref_int).max() < VS_PYTHTB * scale, "internal part != PythTB include_external=False"
    err = np.abs(ext - (ref_tot - ref_int)).max()
    assert err < VS_PYTHTB * scale, f"cross + external vs PythTB: {err:.2e} (scale {scale:.2e})"


def test_mixed_block_by_hand(mbt_pair):
    mb, mm, n_occ, dbeta = mbt_pair["base_model"], mbt_pair["mode_model"], mbt_pair["n_occ"], mbt_pair["dbeta"]
    pb, pm = curvature.position_terms(mb, dtype=np.complex128), curvature.position_terms(mm, dtype=np.complex128)
    f, f_mode = curvature.fields(mb, pb, K_TEST), curvature.fields(mm, pm, K_TEST)
    lam = 1j * np.diag(np.linspace(0.01, 0.02, mb.norb))           # anti-Hermitian Lambda(0): A_beta = i Lambda
    lam[0, 1], lam[1, 0] = 0.03, -0.03
    A_beta, dA_beta = axion.parametric_connection(pb, axion.prepare_lambda({(0, 0, 0): lam}, np.complex128), K_TEST)
    beta = axion.beta_terms(f, f_mode, dbeta, A_beta, dA_beta)
    assert np.abs(beta.A).max() > 0

    c3, c4 = curvature.connection(f, n_occ), curvature.connection(f, n_occ, beta)
    ext3 = curvature.omega_cross(c3) + curvature.omega_external(c3)
    ext4 = curvature.omega_cross(c4) + curvature.omega_external(c4)
    assert ext4.shape == (4, 4, len(K_TEST), n_occ, n_occ)
    assert np.allclose(ext4[:3, :3], ext3), "k-k block changed by the mixed-plane path"

    # hand assembly in the band gauge: d = internal connection (31), a / A_vv = rotated A (28), mixed curl
    E, U = np.linalg.eigh(f.H)
    adj = lambda X: np.conj(np.swapaxes(X, -1, -2))
    rot = lambda X: adj(U) @ X @ U
    occ, cond = np.arange(n_occ), np.arange(n_occ, mb.norb)
    gap = E[:, occ][:, None, :] - E[:, cond][:, :, None]
    dH = np.concatenate([f.dH, beta.dH[None]])
    d = np.stack([rot(X)[:, cond[:, None], occ] / gap for X in dH])                # (4, nk, c, v)
    a = np.stack([rot(X) for X in f.A] + [rot(beta.A)])
    a_cv, a_vv = a[:, :, cond[:, None], occ], a[:, :, occ[:, None], occ]
    mixed_curl = np.stack([rot(X)[:, occ[:, None], occ] for X in beta.curl])
    scale = np.abs(ext4).max()
    for i in range(3):
        x_i, x_b = -1j * a_cv[i], -1j * a_cv[3]
        cross = 1j * (adj(d[i]) @ x_b + adj(x_i) @ d[3] - adj(d[3]) @ x_i - adj(x_b) @ d[i])
        pure = mixed_curl[i] - 1j * (a_vv[i] @ a_vv[3] - a_vv[3] @ a_vv[i])
        assert np.abs(ext4[i, 3] - (cross + pure)).max() < HAND * scale, f"(k_{i}, beta) entry != hand assembly"
    assert np.abs(ext4 - np.conj(np.swapaxes(ext4, -1, -2))).max() < HAND * scale, "curvature blocks not Hermitian"
