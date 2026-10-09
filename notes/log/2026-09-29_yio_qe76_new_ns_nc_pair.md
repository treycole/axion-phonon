# YIO base/mode-1 pair on QE 7.6 + the `new_ns_nc` fix, and a centre-symmetrization test (2026-09-29)

Follows [`2026-09-21_qe_new_ns_nc_bug.md`](2026-09-21_qe_new_ns_nc_bug.md) (the bug and
patch) and [`2026-09-21_yio_symmetry.md`](2026-09-21_yio_symmetry.md) (the 4-14% floor
this rerun was meant to test) and [`2026-09-23_ws_ties_in_beta_derivative.md`](2026-09-23_ws_ties_in_beta_derivative.md)
(the shared-centre machinery this rerun uses). First real base/mode pair computed with
the patched DFT.

## 0. The pair

`pw.x` = `~/q-e-nsfix` on the cluster: QE 7.6 + the `new_ns_nc` two-line patch, **not**
the PAW `v_rad`/`g_rad` fix (`qe_paw_nc_uninit_vrad`, found later the same week). **This
turns out not to matter for YIO**: a direct comparison run (192344 vs 193145, same input,
patched vs unpatched PAW) gave iteration-1 energy agreeing to 2e-5 Ry and only a 1.5% Ir
moment shift that can't even be confidently pinned on the PAW fix over the 2.5 months of
other upstream changes between the two builds. The PAW bug's large effect was specific to
MBT (Bi/Te sit at exactly zero moment, where the bug bites hardest); YIO's Ir moments are
never near zero. Not a reason to hold off on this pair's numbers.

- Base: `data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/trial_02_Y_p_Ir_d_O_sp/192473bc889`
- Mode 1: `data/Y2Ir2O7/phonon/GM2-/mode1/Q_0.01A/output/192345/trial_02_Y_p_Ir_d_O_sp/no-displacement/192908bc889`
- Both Wannierized with `use_ws_distance`, `transl_inv_full`, `write_ndegen_applied` all
  true, `num_iter = 0` (no further localization beyond the initial projection — this
  matters, see §3). `_hr.dat`/`_r.dat`/`_wsvec.dat` came from a separate `restart = plot`
  pass (the in-job plot stage silently died, cause not found — see `wannier90-dev` build
  notes; a dedicated compute-node job worked).
- Minimum gap: 1.04 eV, up from the unpatched pair's 0.72 eV — consistent with pitfall 7's
  "Hubbard spin channel too weak" diagnosis of the bug's effect.

## 1. Symmetry (`check_magnetic_symmetry`, `jaemo` route)

Results: `validation/results/magnetic_symmetry_jaemo_qe76/` (old: `magnetic_symmetry_jaemo/`).

The part of the curvature built from H(R) alone, and the band energies, both improved
2-4x:

| test | H(R)-only rel., old -> new | band residual (meV), old -> new |
|---|---|---|
| base inversion | 0.045 -> 0.013 | 3.4 -> 1.9 |
| mode1 C2z | 0.052 -> 0.023 | 16.1 -> 4.8 |
| mode1 C3[111] | 0.040 -> 0.032 | 8.7 -> 1.8 |
| mode1 MxyT | 0.031 -> 0.011 | 3.5 -> 1.8 |
| mode1 S4zT | 0.038 -> 0.021 | 18.7 -> 5.5 |

So the DFT bug was real and the fix helps the Hamiltonian. But the **total** curvature
got worse on the test that matters most (base inversion: 0.097 -> 0.152), because the
position-matrix (A(R)) parts now dominate (external relative residual 1.6, up from 0.79).
S4z controls stay at ~2, so the test itself is still discriminating real symmetry from
noise correctly. **Corrects `2026-09-21_yio_symmetry.md`'s "4-14% floor, common to both
position-matrix routes": the floor moved and is no longer flat across parts — it now
sits almost entirely in A(R), not H(R).**

## 2. `dtheta/dQ` for mode 1

`transl_inv_full_from_chk` + `SHARED_WS=True` (rebuilds Jae-Mo's formula from the `.chk`, no
`.mmn`; base and mode share base's Wigner-Seitz images, per `ws_ties_in_beta_derivative`).
`calculations/run_axion.py` now points at this pair with this source.

| nk | base end | mode end |
|---|---|---|
| 4 | -0.01334 | -0.01316 |
| 6 | -0.01835 | -0.01835 |
| 8 | -0.01803 | -0.01801 |
| 10 | -0.01751 | -0.01760 |
| 12 | **-0.01752** | **-0.01761** |

Converged (nk 10 -> 12: 0.04%), endpoints agree to 0.5%, 0 base/mode WS-image
mismatches. About 5x the unpatched pair's `-0.0028/-0.0037` on the same route (which was
itself unconverged), and ~3x the unpatched `.mmn`-formula value (`-0.0053`). Sign
unchanged. `right_centre_from_chk` (the other `.chk` formula) OOM-killed on an 18 GB laptop before
nk=4 (77 GB peak footprint) — no second-route cross-check yet.

**Not yet known: how the §1 symmetry residual propagates into this integral.** The
BZ sum could average pointwise symmetry-breaking down, or could carry a comparable bias.
Untested.

## 3. Does symmetrizing the Wannier centres fix the A(R) residual? No — and why not

Hypothesis going in: `transl_inv_full` puts a per-image phase `exp(-i b.R/2)` on each
periodic copy, so tied images differ by `(-1)^n`; if the raw `.chk` centres used to break
ties are only *approximately* symmetric, A(R) would fail to be exactly covariant even
where H(R) is fine. Testable in isolation on base alone (no pairing): refold both H(R)
(`apply_orbital_dependent_R_mapping`) and A(R) (`position_matrix_transl_inv_full`) on a
centre array forced to satisfy `INVERSION @ e_s - L_s = e_{partner(s)}` exactly, instead
of the raw centres, and rerun `base_inversion`'s covariance test.

Script: `calculations/diagnostics/2026-09-29_centre_symmetrization/check_base_inversion_symmetrized.py`.
Partner found by bipartite (Hungarian) nearest-centre matching mod lattice shifts;
symmetrized centre = `c_s + r_s/2` where `r_s` is the pairwise defect (derivation and an
exactness self-check are in the script).

**Result: every part got worse** (internal 0.013 -> 0.079, cross 0.47 -> 1.50, external
1.59 -> 1.62, total 0.152 -> 0.192) — the opposite of the hypothesis.

**Why, on inspection:** 150 of 176 WFs pair to *exactly* 0.0000 A under naive nearest-centre
inversion matching (double-checked against the `.win` atoms, which pair exactly at this
same origin — not an origin-choice bug). But 26 WFs (13 pairs) all show the *identical*
defect, 6.2401 A, sitting at permutations of the Ir sublattice coordinates
`{0, 2.5475, 5.095, 7.6425}` — a same-site orbital-degeneracy collision, not noise: several
WFs crowd near the same point, and pure geometric nearest-neighbor matching (no orbital-
character information) picks the wrong partner for that subset, moving some centres by up
to 2.55 A and corrupting their WS-image assignment for real. This contaminates the test:
the 85% that were already exact had no room to improve, and the 15% that were mismatched
were actively broken further.

**Conclusion:** this test is negative evidence against the tie-phase hypothesis, not
confirmation of it — the raw centres are *already* essentially exact for the large
majority of WFs, so a centre-based fix isn't the likely explanation for the residual.

**`num_iter = 0` is not the explanation either, and `num_iter > 0` is not a valid test
of it.** Projection-only WFs (`num_iter = 0`) are the *best*-case setup for symmetry, not
a source of the floor: trial orbitals are pinned at atom centres by construction, so their
centres should already track the point group closely (matching the 85%-exact finding
above). Maximal localization (`num_iter > 0`) would not be a clean test of "does more
localization reduce the floor" even if it were tried: (a) it applies a k-dependent unitary
rotation among the WFs that the external-curvature formula's derivation does not
accommodate, so its output can't be fed through the existing pipeline at all, and (b) it's
a numerical spread minimization with no symmetry constraint of its own, so it isn't even
guaranteed to preserve symmetry better than projection — plausibly worse. **Genuinely
open**: checked whether the trial projections lack per-site local axes (`.win`'s
`begin projections` block has plain `Ir1: d` / `Ir2: d` / ... with no `zaxis=`/`xaxis=`,
so Wannier90 uses the global Cartesian frame for every site) — this is likely *not* a bug
given the lattice's own Cartesian axes already coincide with the cubic crystallographic
axes the space group operations are simple in, but it wasn't verified computationally.
Also open: a proper orbital-character-aware pairing (using the `.chk`'s disentanglement
matrices, not geometry) for the 15% mismatched subset. Neither settled.

## Symmetrizing the BZ integral (2026-09-29, later)

Tested whether §1's residual biases §2's `dtheta`, by group-averaging the already-saved
`c2_density` grids over mode 1's magnetic point group (script:
`calculations/diagnostics/2026-09-29_centre_symmetrization/symmetrize_c2_test.py`). Since
the group maps the reciprocal lattice to itself, every mesh point's group image is exactly
another mesh point (integer index permutation, no interpolation) — cheap, reuses
`nk_sweep_tif_shared_ws/*.npz` directly, no new curvature evaluation.

Closing the group from the notebook's own generators gives exactly 24 elements
(`-4'3m'`), with **every improper element paired with time reversal** (never appears
unitary alone) — forced group structure, not a fit. This matters: `c2` is quadratic in
the curvature, so the naive guess "the antiunitary sign cancels, weight = det(R) alone"
is tempting but wrong — checked empirically against the saved data rather than trusted
from a hand derivation, since a group with an equal proper/improper split would force
`weight=det(R)` to make the *total* integral vanish identically (it doesn't, so that
law is wrong). `weight(g) = det(R) * (-1)^antiunitary` fits the data at 14% residual
(matching the known floor) with 98% correlation; `det(R)` alone has ~0% correlation.
Real, independent confirmation the curvature code transforms correctly under the group.

**But the actual test — symmetrized vs. plain `dtheta` — is a trivial null result, not
a validation.** Every element of this group has `weight = +1` (proper+unitary and
improper+antiunitary both give `+1`), so `sum(weight) = |G|` exactly, and group-averaging
over a measure-preserving bijection of the mesh with that property leaves the *total sum*
identical no matter how badly the field breaks covariance pointwise (confirmed: 0.00%
change at every nk, both endpoints). This test cannot distinguish a clean `dtheta` from
one biased by the floor — wrong test for the question, not a wrong answer to it. A rough
substitute (propagate the measured pointwise noise scale to the mesh sum assuming
point-to-point independence) gives `|dtheta|/sigma ~ 226` at nk=12, consistent with the
sign/magnitude being real, but the independence assumption is unverified and could go
either way if the noise is correlated across k (plausible if it traces to specific
orbitals rather than to k itself).

## Next

1. `right_centre_from_chk` cross-check needs a machine with more memory than this laptop (~77 GB peak,
   traced to `_map_R_by_orbital_positions` allocating a dense array over a 27x-expanded
   R-set, called 4x for the shared-centre `right_centre_from_chk` route vs. `transl_inv_full_from_chk`'s 8
   small per-b-vector calls), or a cheaper implementation of that function.
2. The A(R) floor's real-space origin is still unexplained; no cheap analysis-of-existing-
   data test has found it (§3, and the symmetrized-integral test above). The next genuine
   test needs new data or code, not another trick on this pair's outputs.
