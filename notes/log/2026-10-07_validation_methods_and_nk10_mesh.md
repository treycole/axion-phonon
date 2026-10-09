# How the YIO symmetry checks work, and the 10×10×10 ab-initio mesh test (2026-10-07)

Follows [`2026-10-06_segni_fix_validated_and_reruns.md`](2026-10-06_segni_fix_validated_and_reruns.md). This note
covers four things:
- the 10×10×10 ab-initio mesh test on the base (§1);
- a fix to the gzipped-checkpoint reader that the 10×10×10 `.chk` needed (§2);
- why the SCF density can be exactly symmetric while the NSCF eigenvalues are not (§3);
- how each symmetry comparison and validation is actually done (§4), for bugs #7 and #16 and the curvature tables.

**Bottom line.** A denser ab-initio mesh does **not** reduce the YIO symmetry floor. On the pre-segni-fix DFT the
total curvature residual got worse (`.chk` route 1.57% → 3.60% median), and the band-energy residual stayed at
~5 meV on every mesh and route. That ~5 meV is bug #16, which a k-mesh can't remove. The test was run on 192344's
DFT, which predates the #16 fix, so it says nothing yet about the patched reruns.

## 1. The 10×10×10 mesh test

**Setup.**
- **DFT:** base 192344 (QE 7.6 + `new_ns_nc` fix, **without** the segni fix). New NSCF on a 10×10×10 mesh
  (1000 k, `nosym`) from 192344's SCF density, then pw2wannier90 and wannier90 (job 202714; cluster folder
  `phonon/base/u_3.0/192344bc889/mesh_10x10x10/`). The `.win` is identical to the 8×8×8 trial_02 run apart from
  `mp_grid`: 176 WF, `num_iter = 0`, `dis_num_iter = 0`, `use_ws_distance`, `transl_inv_full`,
  `write_ndegen_applied`.
- **Wannier90 rerun:** 202714's final `wannier90.x` (64 MPI ranks) was killed by SGE `h_vmem` (maxvmem 233 GB vs
  64 × 3.6 GB) right after disentanglement. Rerun as **202972** in `mesh_10x10x10/wannier_np1/` with
  `mpirun -np 1` on a fully reserved 64-slot node, reusing 202714's `.amn`/`.mmn`/`.eig`. It finished in 67 min with
  a peak of 29.6 GB. (wannier90.x copies its arrays on every rank, and with zero iterations there is almost
  nothing to parallelize: always run it on one rank.)
- **Spreads:** Ω_total 196.8 Å² against 188.6 Å² at 8×8×8. Almost all of the change is in Ω_I (162.5 vs 155.0),
  which Wannier90 evaluates by finite differences on the mesh and so depends on it. Ω_D and Ω_OD barely move.
- **Local copy:** `data/Y2Ir2O7/base/soc/u_3.0/output/192344bc889/mesh_10x10x10/202972bc889/` (`_tb.dat` not
  copied, since it duplicates `_hr`/`_r`). H(R) has 1891 R vectors against 1097 at 8×8×8.

**Check.** The full m-3m′ table (48 operations; method in §4.4) on two position-matrix routes:
- **`chk`:** H(R) from `_hr.dat` folded on the `.chk` centres, and A(R) from `position_matrix_mmn_formula` on the
  `.chk` `m_matrix` (= `right_centre_from_chk`), at the default WS tolerance (1e-3 Å).
- **`jaemo`:** `_hr.dat`/`_r.dat` from the `write_ndegen_applied` run, read as written
  (`build_model_with_position_matrix(source="rdat_ndegen_applied")`).

The 8×8×8 `chk` row is the existing base table in
`2026-09-29_centre_symmetrization/logs/full_symmetry_tables.log`, computed with the same metric and k-points.
Script, log and arrays: `calculations/diagnostics/2026-10-05_yio_nk10_symmetry/base_symmetry_vs_mesh.{py,log,npz}`.

**Result (median over the 47 non-identity operations):**

| | `chk` 8³ | `chk` 10³ | `jaemo` 8³ | `jaemo` 10³ |
|---|---|---|---|---|
| internal (H only) | 1.06% | 1.68% | 2.25% | 4.55% |
| cross | 16.0% | 12.3% | 57.6% | 26.2% |
| external (A(R)) | 25.9% | **51.3%** | 103% | 101% |
| **total** | **1.57%** | **3.60%** | **8.10%** | **11.9%** |
| total, worst op | 2.42% | 4.16% | 17.2% | 19.2% |
| inversion, total | 0.24% | 0.22% | 15.2% | 13.5% |
| band energies (max \|Δε\|) | — | 4.55 meV | 4.91 meV | 4.64 meV |

**Reading.**
- **The mesh doesn't help.** Total goes up on both routes, and only the cross part improves. So the floor isn't
  Fourier-interpolation error from too coarse a mesh.
- **The band residual is mesh-independent, at ~5 meV.** That is the NSCF Hamiltonian's own asymmetry (bug #16:
  4.83 meV on the same DFT). It is in the DFT, so no Wannier-side mesh change can remove it.
- **The `chk` route is 2-5× more symmetric than the Jae-Mo files** at both meshes. For inversion it is 0.2% against
  13-15%. The internal column differs between the routes too, although both read the same `_hr.dat`. The
  difference must be the WS fold: the in-house fold (with the 1e-3 Å tie tolerance, #13) versus Wannier90's
  `use_ws_distance` choice as written. Probably Wannier90's tie test behaves like our old tight tolerance; not
  verified.
- **Not yet known:** whether the denser mesh helps once #16 is fixed. That needs the test redone on the
  segni-fixed reruns. There is no mode-1 10×10×10 run, so dθ can't use this mesh.

## 2. Checkpoint reader: Fortran subrecords

`wannier/wannier_io.py::_SequentialFortranReader` (used for gzipped `.chk`) assumed one Fortran record per write.
The 10×10×10 `m_matrix` record is 8 × 1000 × 176² × 16 B = 3.96 GB. That exceeds 2 GiB, so ifort splits it into
subrecords, with a negative leading length marker meaning "continued". The reader raised "inconsistent Fortran
record markers". Every 8×8×8 checkpoint is below the limit.

**Fix (2026-10-05):**
- **Subrecords:** `read_record` now follows the subrecords (negative leading marker = more to come; trailing
  markers compared by magnitude) and joins them into one preallocated buffer.
- **No copy:** `m_matrix` is a `complex128` view of the record. Fortran stores it as (re, im) pairs, which is
  already the complex128 layout. The old `re + 1j*im` made about two extra full-size copies (~8 GB at 10³).

**Checks:**
- **8×8×8:** every array bit-identical to the old reader. That checkpoint is uncompressed, so this tests the
  view change only.
- **10×10×10:** reads in 19 s. The centres match `_centres.xyz` to 5e-9 Å, and the overlap blocks at the first
  and the last (k, b) are near-unitary (singular values 0.991-0.9999). So the subrecords were joined in order.
- **Backup:** the pre-edit file was kept in the session scratchpad (`wannier/` is not in git).

## 3. Why the SCF density is symmetric but the NSCF eigenvalues are not

The density is symmetric because QE enforces it at every SCF step, not because the Hamiltonian is symmetric.
1. **The step from occupations to potential doesn't respect symmetry.** In each PAW sphere, the GGA input is
   n↑,↓ = ½(n ± s|m|), where s is taken from a fixed global axis `ux` (Ir1's starting moment), and `ux` doesn't
   rotate with the crystal. Exactly symmetric sphere occupations (`becsum`) therefore still produce a PAW
   potential term D_ij on Ir1 that is not the image of Ir2-4's. Symmetric input, a function that doesn't respect
   the symmetry, asymmetric output.
2. **The SCF hides it.** Every iteration QE symmetrizes ρ, m, ns and `becsum`, and also averages D_ij over the 48
   operations (`PAW_symmetrize_ddd`). The SCF Hamiltonian and everything stored are exactly symmetric. But they
   contain the group average of the spurious term, so the converged state is shifted, and the SCF must be rerun
   too.
3. **The NSCF exposes it.** The NSCF builds D_ij once from the stored `becsum` and never averages it (with
   `nosym` only the identity exists anyway). Its Hamiltonian keeps Ir1's own term, so ε(k) ≠ ε(gk), Ir1's on-site
   levels sit 12 meV off Ir2-4, and Wannier90, the curvature and dθ inherit it.

The "energy" that breaks symmetry is therefore the band eigenvalues compared at k and gk. A total energy is a
single number and can't be tested this way.

## 4. How each check is done

### 4.1 `new_ns_nc` (bug #7)

1. **Python port of the symmetrizer** (`calculations/diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py`).
   - **Port:** QE's own `irt`, `sr`, `t_rev`, `find_u` and `d_matrix`, fed with the operation data that run
     178382 printed (`verbosity = 'high'`).
   - **Test:** start from QE's analytic AIAO ns, check it is invariant under all 48 forward actions (6e-15), then
     apply the symmetrizer.
   - **As written:** it does not return the input (|m| 3.0 → 1.0), so it isn't a valid symmetrizer.
   - **With the pre-image `irt(invs(isym), na)`:** it returns symmetric input unchanged and is idempotent, to 1e-15.
2. **Local pw.x, patched vs unpatched** (`calculations/qe_patches/verification/compare_hubbard_moments.py`).
   - **Setup:** the YIO input, the same Davidson settings and seeds, so both binaries receive identical raw
     occupations and only the symmetrizer differs.
   - **Compared:** the four Ir Hubbard moments per iteration, as the spread in |m| and each moment's tilt from its
     own ⟨111⟩ axis.
   - **Result:** 2.5% / 0.7° unpatched, 0.00% / 0.00° patched.
   - **Caveat:** only the comparison counts; that local build's SCF was numerically broken.
3. **Cluster, production input** (192271 patched vs 178382 unpatched).
   - **Symmetry:** 48 operations found (36 with fractional translations), so the operations that expose the bug
     were active.
   - **Result:** patched, the four Ir are equal and on their own axes at every iteration; unpatched, the Hubbard
     moment collapses (0.058 converged, against 0.48 patched at iteration 20).
   - **Caveat:** stopped at 20 iterations, unconverged.
4. **Minimal Fe4 case** (QE develop vs `fix/new_ns_nc-preimage`, 2 iterations, 1 rank, 9 s).
   - **Result:** the four Fe moments spread 7.9-9.6% unpatched and 0.00% patched.
   - **Pattern:** the unpatched split pairs exactly the atoms that the image/pre-image mix-up predicts.
5. **Stored quantities of the converged patched SCFs** (192344, 193145; scripts in
   `calculations/diagnostics/2026-10-05_ws_tolerance_window/`).
   - **`density_symmetry.py`:** ρ and m in G-space under each {R|t}: ρ(RG) = e^{−iRG·t} ρ(G) and
     m(RG) = e^{−iRG·t} s·det(R)·R m(G), with s = −1 for primed operations. t is recovered from the atom basis.
     Exact to 1e-15.
   - **`hubbard_ns_covariance.py`:** ns(g·a) = W ns(a) W† with W = D_d(R) ⊗ U(R), and time reversal on the spin
     part for primed operations. The conventions (D vs Dᵀ, U vs U*, spin-block order) are not assumed: they are
     fitted on Ir1's 12 site operations, then applied to all 48. Exact to 3e-10, phases included, once QE's
     Condon-Shortley signs on xz/yz are used.
   - **`paw_becsum_magnitudes.py`:** magnitudes only, under the three C2 operations that map Ir1 onto Ir2-4 (they
     act on real-harmonic projector pairs by signs only). Equal to 1e-14. This is the weakest of the checks.
   - **Scope:** these show the stored SCF state is symmetric. By §3 they could not detect #16.

**PAW `v_rad` comparison** (bug #12; 192344 vs 193145, same input). Compared the total energy, the iteration count
and the Ir moments: 2.3e-5 Ry, 411 vs 252 iterations, |m| +1.5%. All of these are symmetrized SCF quantities. That
is why "the PAW fix doesn't matter for YIO" was right about #12 but blind to #16.

### 4.2 Bug #16: finding it

1. **DFT eigenvalues** (`eig_symmetry.py`; full precision from `data-file-schema.xml` in `nscf_test_analyse_xml.py`).
   - **Data:** the NSCF `.eig` on the full 8³ mesh with `nosym`, bands 1-156.
   - **For each operation g:** take k to Cartesian (k_red·B, with B the reciprocal lattice), apply R (−R for primed
     operations), go back to reduced coordinates, wrap into [0, 1), and look the point up in the mesh.
   - **Reported:** max |ε_n(gk) − ε_n(k)| per operation.
   - **Fractional translations** don't move k (they only rephase the wavefunctions), so they don't enter.
2. **Why every gk is on the fixed grid.**
   - **The grid is Γ-centred:** every k is exactly n/8 with integer n (max deviation 0), and Γ is included.
   - **Every m-3m′ operation is a symmetry of the fcc lattice:** in reduced reciprocal coordinates each is an
     integer matrix (all 48 integer to 1e-17). Example: C3[111], the Cartesian permutation (x, y, z) → (z, x, y),
     acts on reduced coordinates as a cyclic permutation.
   - **Consequence:** an integer matrix maps n/8 to m/8, a grid point after wrapping, which is allowed since
     ε(k+G) = ε(k). Time reversal maps n/8 to −n/8, also on a Γ-centred grid. Each operation therefore permutes
     the 512 k-points, and the script's dictionary lookup would fail if any gk fell off the grid.
   - **When it would fail:** for a shifted Monkhorst-Pack grid (a ½-offset grid is not closed under all
     rotations), or for an operation that isn't a lattice symmetry.
3. **On-site levels** (`onsite_spectra.py`).
   - **Assignment:** each Wannier function goes to the atom nearest its centre. Each Ir has 10 WFs (d, × 2 spinor).
   - **Levels:** the 10×10 on-site block H_ij(R=0) of one Ir holds the trigonal crystal field, the SOC and the
     exchange splitting. Its eigenvalues are what the notes call that atom's "on-site levels".
   - **Test:** g maps Ir1's WF set onto Ir2's, so H₀(Ir2) = U H₀(Ir1) U† and the spectra must coincide whatever
     the gauge within each set.
   - **Result:** Ir2-4 agree to ~0.3 meV and Ir1 is 11.6-11.9 meV off. Controls: Y ×4 ≤ 0.1 meV, O(8b) 0.04 meV.
     The 12 O(48f) split into two groups of 6 (bonded to Ir1 or not), 1.6 meV apart.
4. **Isolating the cause** (cluster NSCF tests, `base/u_3.0/test_nscf_sym_2026-10-05`; table in the 10-05 note).
   - **k-set:** `nscf_test_kpoints.py` builds complete symmetry orbits (all gk, −gk for primed operations) of the 2
     most asymmetric mesh points plus one generic 48-point orbit, 60 k in all. Each image is in the set by
     construction.
   - **Changes, one at a time:** binary, process count, no Hubbard U, atomic projectors, LDA, and `ux` moved by a
     tiny Y seed moment. Only LDA (0.0 meV) and moving `ux` (0.1-0.2 meV) removed the breaking.

### 4.3 Bug #16: validating the fix (run C2)

- **Setup:** C2 (job 205711) is test C with only the binary changed: the same 193145 SCF save, the same 60 k, a
  byte-identical `nscf.in`. Eigenvalues are compared at full precision from each run's `data-file-schema.xml`.
- **Result:** worst operation 4.83 meV → 46 μeV, and Ir1's site group no longer stands out.
- **Limit:** the SCF underneath is still unpatched, so this validates the NSCF Hamiltonian, not a fully patched
  chain. The production reruns (207546-207549, plus 209009 for mode 2) do the full chain.

### 4.4 Curvature covariance tables (all YIO tables, including §1)

- **k-points:** 7 generic points from a fixed seed (`default_rng(20260825)`), uniform in ±0.23 Å⁻¹ Cartesian.
  They are off-grid; the Wannier-interpolated model can be evaluated at any k, and so can their images.
- **Quantity:** ω(k) = (Tr Ω_yz, Tr Ω_zx, Tr Ω_xy) over the 156 occupied bands, split into internal, cross and
  external by `modules/curvature.py`.
- **For each operation:** evaluate ω at k′ = R k (k′ = −R k for primed operations) and compare with the
  prediction ±det(R)·R·ω(k), with the minus sign for primed operations.
- **Reported:** max |ω(k′) − prediction| divided by max |ω| over the evaluated and predicted values, per part.
  Also max |ε_n(k′) − ε_n(k)| over all bands.
- **Calibration:** the S4z control (not a symmetry) gives ~199%, so the test separates kept from broken
  symmetry. A relative number is only comparable within one part, since external is ~25× smaller in absolute size.

## 5. Open
1. (#16) Redo §1 on the segni-fixed base once its Wannier stage is done. Run a mode-1 10×10×10 only if the patched
   base improves with the mesh.
2. (#16) Delete the pre-segni 10×10×10 data once agreed: ~36 GB on the cluster (`mesh_10x10x10/202714bc889/`
   nscf `tmp/`, `.mmn.gz`, `.amn.gz`; `wannier_np1/` holds hard links to the same data) and ~13 GB locally.
