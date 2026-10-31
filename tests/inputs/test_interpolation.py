"""The interpolation inputs of the production models, before any curvature (port of
analysis/validation/inputs/check_interpolation_inputs.ipynb, which checked the .mmn route on SrTiO3; production builds
H(R) and A(R) from the .chk, axion_phonon/wannier.py). On the MnBi2Te4 and Y2Ir2O7 production bases:

- H(R) reproduces the Wannier-gauge Hamiltonian of the .chk, V^dag E V, at every mesh k (folding chooses images
  within each aliasing class, so mesh values are unchanged);
- off the mesh, H(k), A(k) and curl A(k) are Hermitian, and covariant under k -> k + G with the convention-I phase
  diag(exp(2 pi i tau.G));
- analytic dH/dk and curl A (axion_phonon.curvature.fields) match central differences of the same interpolants.
"""
import numpy as np
import pytest

from axion_phonon import curvature

# ---- settings
H_MESH_EV = 1e-8              # max |eigenvalue difference| on the mesh, eV
HERMITIAN = 1e-12             # relative Frobenius
PERIODIC = 2e-11              # relative Frobenius
DERIVATIVE = 2e-7             # relative Frobenius, best step
STEPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)   # reduced
PLANES = ((1, 2), (2, 0), (0, 1))                     # curl[l] is the plane (l+1, l+2)
G = np.array([1, -1, 2])
# --------------
pytestmark = pytest.mark.data


@pytest.fixture(params=["mbt", "yio"])
def base(request):
    if request.param == "mbt":
        p = request.getfixturevalue("mbt_pair")
        return p["base"], p["base_model"]
    p = request.getfixturevalue("yio_base")
    return p["run"], p["model"]


def rel(a, b):
    return np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300)


def test_hamiltonian_on_mesh(base):
    run, model = base
    HW = np.einsum("kbi,kb,kbj->kij", run["V"].conj(), run["energies"], run["V"], optimize=True)
    H = model.hamiltonian(run["kpts"], flatten_spin_axis=True)
    err = np.abs(np.linalg.eigvalsh(H) - np.linalg.eigvalsh(HW)).max()
    assert err < H_MESH_EV, f"H(R) vs .chk on the mesh: {err:.2e} eV"


def test_hermitian_and_periodic(base):
    _, model = base
    k = np.random.default_rng(910247).uniform(-0.4, 1.4, size=(5, 3))
    pos = curvature.position_terms(model, dtype=np.complex128)
    f, fG = curvature.fields(model, pos, k), curvature.fields(model, pos, k + G)
    phase = np.exp(2j * np.pi * (np.asarray(model.orb_vecs) @ G))
    shift = lambda X: np.conj(phase)[:, None] * X * phase[None, :]          # D^dag X D
    for name, X, XG in (("H", f.H, fG.H), ("A", f.A, fG.A), ("curl A", f.curl, fG.curl)):
        assert rel(X, np.conj(np.swapaxes(X, -1, -2))) < HERMITIAN, f"{name} not Hermitian"
        assert rel(XG, shift(X)) < PERIODIC, f"{name}(k + G) != D^dag {name}(k) D: {rel(XG, shift(X)):.1e}"


def test_derivatives(base):
    _, model = base
    k = np.random.default_rng(240827).uniform(-0.35, 1.35, size=(4, 3))
    pos = curvature.position_terms(model, dtype=np.complex128)
    f = curvature.fields(model, pos, k)
    best_dH, best_curl = np.inf, np.inf
    for h in STEPS:
        dH = np.empty_like(f.dH)
        dA = np.empty((3,) + f.A.shape, complex)                            # dA[mu][nu] = d_mu A_nu
        for mu in range(3):
            e = np.zeros(3)
            e[mu] = h
            fp, fm = curvature.fields(model, pos, k + e), curvature.fields(model, pos, k - e)
            dH[mu] = (fp.H - fm.H) / (2 * h)
            dA[mu] = (fp.A - fm.A) / (2 * h)
        curl = np.stack([dA[mu, nu] - dA[nu, mu] for mu, nu in PLANES])
        best_dH, best_curl = min(best_dH, rel(dH, f.dH)), min(best_curl, rel(curl, f.curl))
    assert best_dH < DERIVATIVE, f"dH/dk vs central difference: {best_dH:.1e}"
    assert best_curl < DERIVATIVE, f"curl A vs central difference: {best_curl:.1e}"
