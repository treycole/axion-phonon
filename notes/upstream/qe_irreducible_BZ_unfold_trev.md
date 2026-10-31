# QE MR draft: `irreducible_BZ` drops −k stars when lattice rotations appear only with time reversal

**Target:** gitlab.com/QEF/q-e, merge request into `develop` + issue. **Patch:** `patches/qe_irreducible_BZ_unfold_trev.patch`
(branch `fix/irreducible_BZ-unfold-trev`). **cc:** Jae-Mo Lihm (author of 1ce38add3).

---

**Title:** irreducible_BZ: k-point list misses whole stars for magnetic groups with nsym == nrot and time reversal
(7.5CN, 7.6, develop)

**Summary.** Since 1ce38add3 ("Simplify irreducible_BZ", 2025-09-17; in qe-7.5CN and qe-7.6), the unfold step of
`PW/src/irreducible_BZ.f90` maps k → −S k for every lattice operation `irot` with `t_rev(irot) == 1`:

```fortran
xks(:) = invs(:,1,irot) * xkg(1) + invs(:,2,irot) * xkg(2) + invs(:,3,irot) * xkg(3)
IF (magnetic_sym .AND. (t_rev(irot) == 1)) xks = -xks        ! line 117
```

The input k points come from `kpoint_grid`, which reduces with the *plain* rotations of the Bravais lattice
(`t_rev = 0`), so the unfold has to use plain rotations too. When the magnetic group contains every lattice rotation
(`nsym == nrot`) but some only combined with time reversal, the set {S : unitary} ∪ {−S : with T} can be a proper
subgroup of the lattice group, and part of each star is never generated. Before 1ce38add3, `irrek_nc` applied `t_rev`
only in the reduction step. The fix removes the one line; the reduction keeps its `t_rev` handling.

**Example.** Rhombohedral MnBi2Te4 with a zone-boundary (T2⁻) distortion: noncollinear + SOC, 4×4×4 MP grid, magnetic
group with E, C3±, C2′ unitary and i, S6±, σd only with T (12 = nrot operations). The unfold actions are E, C3±,
C2′, −i = E, −S6 = C3, −σd = C2′: only D3, so −k is never generated.
- pw.x prints **13 k points instead of 16**. Three −k orbits are missing and their partners carry double weight
  (12/64 instead of 6/64); the orbits cover 52 of the 64 grid points.
- The SCF converges to a non-variational state: total energy 0.31 mRy **below** the undistorted structure (a
  displacement from a force-free structure must cost +½kQ²), |m|abs 10.27 instead of 9.99 μB. Forces are nearly
  right (1e-4 Ry/bohr).
- With `nosym = .true.`, or with the patch, the same input converges to E_base + 0.0027 mRy and |m|abs 9.99; the
  patched run agrees with `nosym` to 1e-7 Ry and prints 16 k points.

**Who is affected.** Noncollinear magnetic runs (`magnetic_sym`) whose magnetic point group contains every rotation
of the lattice point group, at least one only with time reversal. Runs with nsym < nrot, or with all lattice
rotations unitary, get correct lists (checked: pyrochlore Y2Ir2O7, m-3m′ and three 24-operation subgroups; MnBi2Te4
undistorted, −3m′); for those the patched and unpatched k lists are byte-identical.

**How to check a k list.** Compare each weight with the size of its orbit under the magnetic group (k → S k, or −S k
when the operation carries T) and check that the orbits cover the grid. A Python port of `kpoint_grid` +
`irreducible_BZ` reproduces the 13-point list exactly before the patch and gives 16 after it.

**Test.** Proposed test-suite case: a small noncollinear magnet with a group of that type (e.g. a rhombohedral or
hexagonal antiferromagnet with moments along c and a distortion that makes inversion antiunitary), checking the
number of k points and weights.

**Attachments to prepare:** the input of the MnBi2Te4 case (`data/MnBi2Te4/archive/qe_patches/irreducible_BZ/`), the outputs
210645 (unpatched), 210656 (`nosym`), 212618 (patched), and `kpoint_weights_check.py`
(`archive/diagnostics/2026-10-10_qe_irreducible_bz/`).
