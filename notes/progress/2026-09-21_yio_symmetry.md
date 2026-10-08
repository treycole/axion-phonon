# YIO curvature symmetry at base and displaced, two position-matrix routes (2026-09-21)

`validation/symmetry/check_magnetic_symmetry.ipynb` on Y2Ir2O7 (176 WF, 156 occupied), base
`m-3m'` and mode 1 `-4'3m'`, once with the `.mmn` position matrix and once with Jae-Mo's
`write_ndegen_applied` files (`POSITION_SOURCE = "jaemo"`, read through PythTB exactly as written).
Results: `validation/results/magnetic_symmetry/` (mmn) and `magnetic_symmetry_jaemo/`.

## Short answer

The curvature obeys the magnetic point group to a floor of **about 4 to 14% of the typical curvature**
(|omega| ~ 0.7 A^2), in both structures and in both routes. A non-symmetry is violated by about 190%,
so the test resolves the difference. The floor comes from the shared H(R), not from the position matrix.

## What changed in the notebook

* The mode-1 z axis: all three components vanish. C2z gives `omega_x = omega_y = 0`, and the antiunitary
  `S4z T`, which also fixes the axis, leaves no allowed direction (checked numerically: intersection of
  the two constraints is empty). The old test expected `omega_z` allowed and normalised by it, which
  produced the 78 to 193% rows. It now mirrors the base test: all components against the generic curvature.
* `mode1_S4zT` covariance added (the operation the corrected statement depends on).
* Calibration: `base_S4z_control` and `mode1_S4z_control` apply S4z **without** T. Both structures are nearly
  inversion symmetric, so the true relation is `omega(S4z k) = +S4z omega(k)` and the unitary test
  predicts `-S4z omega(k)`: the residual should be about 2. My first control, inversion on mode 1, was a
  mistake: the 0.01 A displacement breaks inversion far below the floor, so it could not discriminate.
* `POSITION_SOURCE = "jaemo"` loads `<run>/jaemo/` with no remapping.

## Results (total curvature, |violation| / generic scale)

| test | .mmn | Jae-Mo |
|---|---|---|
| base inversion | 0.079 | 0.097 |
| base z axis, mirror intersection (omega = 0) | 0.049 | 0.076 |
| mode 1 C2z | 0.133 | 0.138 |
| mode 1 C3 [111] | 0.131 | 0.126 |
| mode 1 M_xy T | 0.049 | 0.066 |
| mode 1 S4z T | 0.063 | 0.081 |
| mode 1 z axis, and Gamma (omega = 0) | 0.065, 0.051 | 0.081, 0.042 |
| mode 1 [111] line, [1,-1,0] line | 0.085, 0.038 | 0.093, 0.054 |
| **control, base S4z without T** | **1.93** | **1.91** |
| **control, mode 1 S4z without T** | **1.96** | **1.92** |
| mode 1 inversion (broken at O(Q)) | 0.102 | 0.214 |

Control over the worst true symmetry: 15x (mmn), 14x (Jae-Mo). The predicted value from the saved arrays
was 1.975. Along the z axes |omega| is at most 0.047 (mmn) and 0.059 (Jae-Mo) in every component, for
both structures.

Band energies `max |E(Rk) - E(k)|` agree between the routes to 0.004 meV (same H(k)): 3.4 meV for
inversion, 16 meV for C2z, 18.7 meV for S4z T. The S4z controls give 14.4 and 19.1 meV, so the band
energies alone cannot separate a broken symmetry from a kept one here; the curvature can.

## Where the floor comes from

* The residual at each k is reproduced by both routes: correlation +0.83 to +0.96 across the six covariance
  tests, so roughly 83 to 96% of the variance is common. The routes share H(R) and differ only in how A(R) is
  built, so the floor is not a property of the interpolant. Together with the raw QE eigenvalue
  violations (4 to 17 meV) and the unequal Hubbard `ns` on the four Ir sites (see the memory notes of
  2026-08-25) it points at the DFT input; the QE symmetrization bug itself is still unconfirmed.
* The Jae-Mo route is a little less symmetric than `.mmn` (rms residual 17 to 27% larger; the totals of the
  true symmetries differ by at most 2.7 percentage points), and both are far better than the stock `_r.dat` route saved earlier
  (123% for inversion).
* The stock `_r.dat` comparison in `results/magnetic_symmetry_rdat` is from an older run and is not a
  fair like-for-like with this table.

## Two things that matter for dtheta/dQ

1. **The floor is not common-mode between the two structures.** Inversion residuals of base and mode 1 at
   the same seven k: band energies correlation +0.01 (rms difference 0.80 meV = sqrt(0.45^2 + 0.67^2),
   the independent-noise value), curvature +0.26 (mmn) and +0.29 (Jae-Mo). So the symmetry-violating
   error does not cancel in `H_mode - H_base` or in a curvature difference. Whether the BZ integral
   suppresses it enough is **not measured**; a symmetric mesh removes the non-invariant part of the error
   to first order but not a symmetric part.
2. **At Q = 0.01 A the inversion breaking is at or below the floor.** On the `.mmn` route `mode1_inversion`
   (0.102) lies inside the range of the true symmetries (0.049 to 0.133). On the Jae-Mo route it is
   0.214, above that range, but the two routes disagree about it much more (rms difference 0.034 A^2)
   than about any true symmetry (0.006 to 0.018). That is suggestive of the route dependence already
   seen in dtheta (about 0.5 to 0.6x), but it rests on seven k-points and one displacement, and is not
   established.

## Not established

* How the floor propagates into dtheta/dQ (no test symmetrizes the inputs or repeats with an independent
  DFT realisation).
* The physical size of the inversion-odd curvature at 0.01 A. A second, larger amplitude (the parked
  second-amplitude test) would raise it above the floor and also test linearity.
* The internal, cross and external pieces are individually far less covariant (external 33 to 84% of its
  own scale on the true symmetries, `.mmn` route) than the total (4 to 14%). Only the total needs to be covariant when
  the Wannier set is not symmetry adapted.

## DFT-side diagnosis of the floor (added 2026-09-21, from the base SCF outputs)

Inputs: `base/soc/u_3.0/178382bc889/Y2Ir2O7.scf.out` (`HUBBARD (ortho-atomic)`) and the first steps of a
`HUBBARD (atomic)` rerun pasted by the user. The four Ir are four QE species (`Ir1`-`Ir4`).

* **The input is symmetric.** The first printed ns is the analytic AIAO guess in both runs (eigenvalues exactly
  0.4 x5 and 1.0 x5, moments exactly +-sqrt(3) along the four <111>, |m| = 3.00000, tilt 0.00 deg).
* **The asymmetry is produced by the first computed ns and then stays.** Ortho-atomic, Hubbard-projected |m| on
  Ir1..Ir4: 0.36372 / 0.35709 / 0.35766 / 0.34792 at the first computed ns (spread 4.4%, tilts 1.1 / 0.4 / 1.9 /
  1.1 deg), 9.0% and 1.4 / 1.4 / 5.0 / 2.5 deg at convergence. Not an instability growing out of noise.
* **The projector is not the cause.** `atomic` gives the same pattern: |Tr up - Tr down| = 0.0311, 0.0314
  (Ir1, Ir2) against 0.0356, 0.0360 (Ir3, Ir4), ratio 1.16 (ortho-atomic converged: 0.0306, 0.0312 against
  0.0336, 0.0342, ratio 1.12). Spin-summed Tr[ns] is identical on the four sites to 5 digits in both.
* **The site magnetisation printed by `report = 1` is exactly symmetric.** Sphere-integrated |m| = 0.1225080 on
  all four Ir to 7 digits, tilt 0.00 deg, charge 4.308703 on all four, at every iteration. QE symmetrises the
  density itself, so this shows the density symmetriser works and the ns symmetriser does not; it does not
  show the Hamiltonian was symmetric before symmetrisation.
* **Order of magnitude closes.** The Hubbard potential built from the asymmetric ns, U x (delta ns) = 3 eV x about
  5e-3, is about 15 meV, the size of the 4 to 17 meV band-energy violations.

**Cause found and verified afterwards:** [`2026-09-21_qe_new_ns_nc_bug.md`](2026-09-21_qe_new_ns_nc_bug.md).
`new_ns_nc` applies the rotation in the wrong orientation relative to `nb = irt(isym, na)`, which is wrong for
every operation of order 3, 4 or 6 that permutes Hubbard atoms; `symmetry_with_labels` is what made such
operations appear. The suspects listed above other than this one (time-reversal spin matrix, `find_u`,
ortho-atomic projector) are cleared. A patched `pw.x` gives Ir moments identical to the printed precision.

