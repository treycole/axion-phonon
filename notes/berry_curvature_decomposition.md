# Implementation notes: external Berry-curvature terms

**The derivation lives in [`berry_curvature_derivation.md`](berry_curvature_derivation.md).**
Equation numbers below refer to that note. This file holds only what is specific
to the code: the symbol map, the two things that are easy to get wrong, and the
validation checklist.

---

## 1. Symbol map

| Derivation | [VS] | `modules/curvature.py` / `modules/axion.py` |
| --- | --- | --- |
| $\boldsymbol\tau_s$ (fixed across $\beta$) | footnote 7 | `PositionTerms.tau` (= `model.orb_vecs`, set by `W90.model(orb_vecs=tau)`) |
| $X_{\mu,st}(\mathbf R)$, Eq. (8) | Eq. (11) | `model._pos_r` -> `PositionTerms.X` |
| $\mathbf b_{st}(\mathbf R)$, Eq. (7) | Eq. (7) | `bond` in `position_terms` (from `R` and `tau`) |
| $\mathbb{A}_\mu(\mathbf k)$, Eq. (8) | Eq. (18) | `Fields.A` (`curvature.fields`) |
| $\mathbb{\Omega}_{\mu\nu}(\mathbf k)$, Eq. (9) | Eq. (21) | `Fields.curl`, weights `PositionTerms.curl_weight` |
| $\Lambda_{st}(\mathbf R)$, Eq. (33) | Eq. (71) | `Lambda_R` argument of `axion.dtheta` (`axion.prepare_lambda`) |
| $\mathbb{A}_\beta(\mathbf k)$, Eq. (33) | Eq. (70) | `BetaTerms.A` (`axion.parametric_connection`) |
| $\mathbb{\Omega}_{i\beta}(\mathbf k)$, Eq. (35) | Eq. (72) | `BetaTerms.curl` (`axion.beta_terms`) |
| $\partial_\beta\mathbb{A}_i$ | — | `dA` inside `axion.beta_terms` |
| $d_\mu=-iA^{{\rm int},cv}_\mu$, Eq. (28), (31) | Eq. (53) | `Connection.d` |
| $a_\mu=A^{{\rm ext},cv}_\mu$, Eq. (15) | Eq. (50) | `Connection.a` |
| $A^{{\rm ext},vv}_\mu$, Eq. (16) | Eq. (33) | `Connection.A_vv` |
| $\Omega^{{\rm ext},vv}_{\mu\nu}$, Eq. (16) | Eq. (62) | `Connection.Omega_vv` |
| $\mathcal{B}^{\rm cross}$, Eq. (30) | part of Eq. (61) | `curvature.omega_cross`, pythtb `BX` |
| $\mathcal{B}^{\rm ext}$, Eq. (27) | part of Eq. (61) | `curvature.omega_external`, pythtb `BE` |

`omega_cross` and `omega_external` are the two position-matrix pieces (Eqs. 30, 27);
`omega_internal` is the Kubo term $\mathcal{B}^{\rm int}$, Eq. (29); `omega` sums
the three (Eq. 24).

---

## 2. Two things that are easy to get wrong

**$\Lambda$ is optional; $\mathbb{A}_\beta=0$ is never silent.**
`connection(f, n_occ, beta)` takes a `BetaTerms` for four directions — it carries
the mixed curl $\mathbb{\Omega}_{i\beta}$, built from $\partial_\beta\mathbb{A}_i$
(the finite difference of the two structures' connections, `axion.beta_terms`).
Whether it *also* carries a nonzero $\mathbb{A}_\beta$ is optional: pass
`Lambda_R=...` to `axion.dtheta` for the full treatment (Eq. 33–36), or omit it
for the frozen gauge (§8.5), where `BetaTerms.A is None` (`Connection.frozen_gauge`
is `True`) and every term containing it is skipped rather than multiplied by zeros.

**Fixed orbital positions.** Assumption (A3) requires $\boldsymbol\tau_s$ held fixed
across the two structures. Build *both* PythTB models with
`W90.model(orb_vecs=tau)` using that same array: `position_terms` reads it from the
model, and `axion.dtheta` refuses a mismatch.
Both then form $X=\langle\mathbf 0s|\hat r|\mathbf Rt\rangle-\tau_s\delta_{\mathbf R\mathbf 0}\delta_{st}$,
retaining the Wannier-center offset. Leaving each structure at its own Wannier
centers zeroes that diagonal and silently discards the center motion that
$\partial_\beta X$ is meant to carry.

---

## 3. Checks

1. $\mathbb{A}^\dagger_\mu=\mathbb{A}_\mu$ for all four $\mu$ including $\beta$; $\mathbb{\Omega}^\dagger_{\mu\nu}=\mathbb{\Omega}_{\mu\nu}$.
2. $\Lambda^{\rm II}(\mathbf k)$ antihermitian at every $\mathbf k$ **before** the Fourier transform, and $\Lambda(\mathbf R)=-\Lambda^\dagger(-\mathbf R)$ after. If only the second holds, the $\mathbf R$-space symmetrization is laundering a $k$-space error.
3. $\Lambda$ **zero test**: run Eq. (36) at $\beta=0$; $\tfrac12(M-M^\dagger)$ must vanish to machine precision.
4. $\mathcal{B}^\dagger_{\mu\nu}=\mathcal{B}_{\mu\nu}$ and $\mathcal{B}_{\nu\mu}=-\mathcal{B}_{\mu\nu}$, separately for int, cross, ext, and including the mixed planes.
5. **Gauge covariance**: rotate the occupied states by a random unitary $G$; each piece must map to $G^\dagger\mathcal{B}G$, and ${\rm tr}(\mathcal{B}_{\mu\nu}\mathcal{B}_{\rho\sigma})$ must be unchanged. Note $\mathcal{B}$ itself is *not* invariant — comparing matrix elements will look like a failure when nothing is wrong.
6. $X_\mu(\mathbf R)=0\Rightarrow\mathcal{B}^{\rm cross}=\mathcal{B}^{\rm ext}=0$ exactly.
7. Feeding $\Lambda=0$ must reproduce the frozen-gauge expressions of §8.5.
8. **Cancellation** (§5.3): compute $(A^{{\rm ext},cv}_\mu)^\dagger A^{{\rm ext},cv}_\nu$ explicitly and the $q$-slice of $-iU^\dagger[\mathbb{A}_\mu,\mathbb{A}_\nu]U$; they must be equal and opposite. This is the one step where a sign slip survives silently.
9. **Convention**: $\mathbb{A}_\mu$ and $\mathbb{H}$ in the same convention. Recompute in convention II and confirm $\mathcal{B}$ is unchanged.
10. **$\tau$ independence**: $\boldsymbol\tau_s$ is a pure gauge choice, so every invariant of $\mathcal{B}$ must be independent of it. Verified to ~1e-12 for hBN and Fe.
11. $\mathcal{B}_{ij}$ against PythTB `berry_curvature(include_external=True, non_abelian=True)`, itself validated against `postw90` $J_0+J_1$.
12. Convergence with $J$ of the *total* $\mathcal{B}$ and of $\partial_\beta\theta$, while the three pieces individually shift.
13. **Endpoint spread**: the two $\partial_\beta\theta$ estimates (spatial fields at $\beta$ and at $\beta+\Delta\beta$, one shared secant) agree only if $c_2$ is linear in $\beta$, so their spread measures curvature of the response, **not** validity; there is no right endpoint. Their mean is the $O(\Delta\beta^2)$ estimate. Only a sign flip between them is diagnostic. See the `axion.py` docstring.

Checks 11 and 1–10 are the only ones with an external reference. Per §10.4 of the
derivation note, the $k$–$\beta$ blocks have **no** published implementation to
compare against, so checks 3, 7 and 8, and the analytic moving-frame tests in
`validation/unit/test_curvature.py`, carry the weight there.

### Where each check runs

| Check | Script |
| --- | --- |
| 1, 4, 5, 6 | [`validation/unit/test_curvature.py`](../validation/unit/test_curvature.py) |
| 9, 10 | [`validation/symmetry/check_convention_symmetry.ipynb`](../validation/symmetry/check_convention_symmetry.ipynb) |
| 11 (k–k vs pythtb/postw90) | [`validation/external_curvature/check_vs_pythtb.ipynb`](../validation/external_curvature/check_vs_pythtb.ipynb) |
| 11 (vs postw90 directly) | [`check_vs_postw90_line.ipynb`](../validation/external_curvature/check_vs_postw90_line.ipynb), [`check_vs_postw90_full_mesh.ipynb`](../validation/external_curvature/check_vs_postw90_full_mesh.ipynb) |
| 12 (mesh convergence) | [`validation/observable/born_effective_charge.ipynb`](../validation/observable/born_effective_charge.ipynb) |
| 13 (endpoint spread, a linearity diagnostic) | [`validation/observable/audit_mixed_born_response.ipynb`](../validation/observable/audit_mixed_born_response.ipynb) |
| 2, 3, 7, 8 | **not automated** — see [`validation/README.md`](../validation/README.md) |

The full index, including the input-layer and symmetry checks that this list does
not cover, is [`validation/README.md`](../validation/README.md).
