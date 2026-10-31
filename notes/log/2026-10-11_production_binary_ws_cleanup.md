# Production pw.x `db4b99b4`, WS tolerance 5e-3, MBT per-WF centres, YIO mode 2 restart, cluster cleanup (2026-10-11)

## Bottom line

- **Production pw.x is now `db4b99b4`** (user, 2026-10-11): `~/q-e-kfix`, branch `production/axion-phonon` =
  new_ns_nc (#7) + segni (#16) + irreducible_BZ (#19) on QE develop (the v_rad fix #12 is upstream). Its k lists are
  byte-identical to the previous `e270b05b` for the YIO base, YIO modes 1-3 and the MBT base, so those runs stay
  valid; only the MBT T₂⁻ modes change (13 → 16 k). The four T₂⁻ SCFs were resubmitted (212708-212711).
- **WS tie tolerance default 5e-3 Å** (user OK). For YIO it sits in an empty gap (near-ties ≤ 1.7e-3 Å, next
  excess 2.1e-2 Å). Unit tests pass.
- **MnBi2Te4 cannot use per-WF centres for the fold** (new bug #21): the Bi/Te p_x, p_y centres sit 0.06-0.15 Å
  off-axis, so the fold breaks C3 and the interpolated bands by 62 meV. Atom-pinned positions make them exact. The
  MBT base's traced curvature is forced to zero (P·{T|½,½,½}), so relative curvature-symmetry tables there are
  noise.
- **YIO curvature symmetry, all variants in one table** (§3): the clean runs at 5e-3 give 0.5-0.7% median total;
  what is left is internal, on the 32 spin-mixing operations. Recommendation: accept.
- **YIO mode 2 at Q = 0.01 restarted from the converged Q = 0.005 density** (212707), one amplitude for all modes.
- **Cluster: 6.5 TB used.** A cleanup script (`scripts/cluster_cleanup_2026-10-11.py`) frees ~5-6 TB; running it
  from here was blocked by the permission check, so it waits for the user.
- **Upstream drafts** for every external bug: `notes/upstream/`.

## 1. Production binary

**Source.** Cluster `~/q-e-kfix` was `~/q-e-segnifix` (branch `fix/paw-segni-lsign`, 4b480ac68) with an uncommitted
edit of `PW/src/irreducible_BZ.f90`. Committed there as fa9b19c8f on a new branch `production/axion-phonon`, and
fetched into `~/Repos/q-e` with `fix/paw-segni-lsign`. The binary is unchanged (`db4b99b4`). Ancestry checked: the
v_rad fix 87ac195cf and the new_ns_nc commit 74bc68eb3 are in it.

**GitLab fork.** Pushed 2026-10-11 (user): `production/axion-phonon`, plus the earlier `fix/new_ns_nc-preimage`. Single-commit branches on upstream develop (c5a3cf95b) for
the MRs: `fix/paw-segni-lsign-develop` (5af5489f0), `fix/irreducible_BZ-unfold-trev` (rebased, 308a2f8b2). Upstream
develop still has all three bugs (checked the lines).

**Does it change YIO?** Dry runs of both binaries at low cutoff (ecutwfc 15, one electron step, killed after the
k list) on every current input, `~/kcheck_2026-10-11`:

| input | ops | k points e270b05b → db4b99b4 | k list and weights |
|---|---|---|---|
| YIO base YP | 48 | 8 → 8 | identical |
| YIO modes 1, 2, 3 | 24 each | 8 → 8 | identical |
| MBT base | 12 | 13 → 13 | identical |
| MBT T₂⁻ modes 1-4 | 12 each | **13 → 16** | different |

The patch removes a line that acts only when a lattice rotation is in the group with time reversal and the unfold
then misses part of a star; for every YIO group the unfold already reached the whole star.

**MBT reruns.** `phonon_lda/T2-/mode_{1..4}/Q_0p01_kfix_2026-10-11/`: same input and start (base state,
`startingpot = 'file'`), pw.x `db4b99b4`, SCF only (the NSCF waits for the mesh decision, §4). All four print 16 k.
The held Wannier chains on the invalid runs (MBT 210633-210640) and on the unconverged YIO mode-2 state (210624-210627)
were removed; their folders stay as templates.

## 2. WS tie tolerance 5e-3, and why 1e-3 failed

**What the fold does.** For each pair of Wannier functions (m, n) and each aliasing class R (mod the 8×8×8 mesh),
the fold picks the image R + N·T with the shortest bond |R + N·T + c_n − c_m| and puts H(R), A(R) there; images
within `tolerance` of the shortest share it equally. Symmetry-related bonds must get the same treatment.

**Why 1e-3 failed for YIO.** The spinor projections are quantized along z. The 16 operations of 4/mmm about z keep
each function's spin character, so they map functions to functions and centres to centres. The other 32 (the
threefold axes and the x/y fourfold axes) mix a function with its same-site spin partner, whose centre differs by up
to 8.6e-4 Å. Bonds that symmetry makes exactly equal therefore differ by up to ~1.7e-3 Å. At 1e-3, some of those
near-ties were recognized and their symmetry partners were not, so one bond's H(R) was split over two images while
its partner's sat on one: a non-covariant interpolant, 18% in the external term.

**The window** (`archive/diagnostics/2026-10-05_ws_tolerance_window/tie_window_check.py`, all pairs, all classes,
27 images; log next to it):

| run | largest near-tie excess | smallest genuine excess | entries 1e-3 < x ≤ 5e-3 |
|---|---|---|---|
| YIO LDA 212615 | 1.675e-3 Å | 2.233e-2 Å | 5560 |
| YIO PBE 212616 | 1.721e-3 Å | 2.134e-2 Å | 6992 |
| MBT LDA 210651 | none: continuous through 5e-3 (4.91e-3, 5.09e-3) | | 3448 |

Any tolerance in (1.7e-3, 2.1e-2) gives the same YIO fold; 5e-3 is in the middle. Applied:
`_DEFAULT_WS_DISTANCE_TOLERANCE = 5e-3` (`wannier/position_matrix.py:501`); cache provenance records the tolerance,
so 1e-3 caches are not reused. `validation/unit`: 55 passed.

**Pinning.** "Pinned" replaces each function's position by its atom's (nearest atom image) in the fold and in τ. It
is a choice, not a correction: the image rule only selects the interpolant between mesh points (CLAUDE.md pitfall 1)
and τ is a free embedding shared by base and mode (A3); H(R) and A(R) still use the same rule. Atom positions are
covariant by construction, so ties are exact. For YIO it gives the same totals as 5e-3 (0.52% vs 0.52%), so YIO stays
on centres + 5e-3. For MBT it is needed (§4).

## 3. YIO curvature symmetry: every saved table

`archive/diagnostics/2026-10-05_ws_tolerance_window/symmetry_comparison_table.md` (generated by
`symmetry_comparison_table.py`; per class and per part). Median (max) over the 47 non-identity operations of m-3m′,
7 generic k, A(R) `transl_inv_full_from_chk`, %:

| run | internal | cross | external | total |
|---|---|---|---|---|
| old PBE 192344, tol 1e-3 | 1.06 (1.41) | 18.74 (25.37) | 44.18 (57.50) | 2.02 (2.53) |
| PBE 210600 (edge on a band), tol 1e-3 | 0.61 (0.76) | 8.90 (13.28) | 85.94 (103.90) | 4.04 (5.00) |
| PBE 212616, tol 1e-3 | 0.81 (0.97) | 3.50 (3.63) | 16.39 (17.18) | 0.67 (0.84) |
| PBE 212616, tol 5e-3 | 0.71 (0.77) | 1.29 (1.47) | 0.22 (0.39) | 0.68 (0.75) |
| LDA 212615, tol 1e-3 | 0.58 (0.82) | 2.37 (2.77) | 17.80 (18.14) | 0.68 (0.84) |
| LDA 212615, tol 5e-3 | 0.54 (0.79) | 1.06 (1.20) | 0.13 (0.18) | 0.52 (0.82) |
| LDA 212615, pinned | 0.56 (0.84) | 0.84 (1.28) | 0.14 (0.17) | 0.52 (0.85) |
| LDA 212615, pinned, H from .chk | 0.70 (0.90) | 1.04 (1.19) | 0.16 (0.22) | 0.67 (0.92) |

**Per class:** the 16 z-preserving operations (C2z, C4z, S4z, C2x/y, C2[110], i, mirrors) are ≤ 0.01-0.04% in
every part on the clean LDA runs (≤ 0.1-0.2% on PBE). The 32 spin-mixing ones (C3, S6 about [111], C4/S4 about x/y,
C2[011], m[011]) carry the whole residual: internal 0.78-0.90%, cross ~1.2%, external ≤ 0.2% at 5e-3.

**Reading it** (as first written; superseded by §3.1). The internal part depends only on H(R). It does not move with
the tolerance, with pinning, or with a full-precision H from the `.chk`, and it is the same for LDA and PBE although
LDA's NSCF eigenvalues are 100× more symmetric. **Recommendation: accept 0.5-0.9%.** It is below the old (pre-fix) 2% and far below the 4-14% floor, the
dθ endpoint/linearity and k-mesh errors are of the same order, and a residual of this size in Ω changes dθ by a
similar relative amount at most. The candidate cause (NSCF eigenvector accuracy, the spinor start along z) would cost
a ~9 h NSCF + Wannier rerun to test.

### 3.1 The internal floor is numerical truncation of H (bug #22, found later the same day)

The `.chk`-H row above was built to test whether `_hr.dat`'s 1e-6 eV print caused the internal floor; it didn't move,
because every model was also built with PythTB's `min_hopping_norm=1e-5`, which drops elements one by one in the
spin-z basis and so cannot commute with the spin-mixing operations. Without the cut (`hopping_cut_test.py`, pinned,
op #02, internal/total): `_hr.dat` 0.10/0.08%, `.chk` **0.000/0.000%**. `_hr.dat` rounding alone breaks the
spin-mixing bands by 63 μeV (`.chk`: 0.6 μeV; `hr_rounding_eigenvalues.py`). So the YIO DFT, NSCF and Wannierization
are symmetric to the precision of the test; the earlier hypothesis (NSCF eigenvector accuracy) is ruled out. Fix
proposed: H from the `.chk` and no cut. See `notes/bugs/22_hopping_cut_and_hr_rounding.md`.

**Full group with the fix, production fold (2026-10-12):** total median (max) LDA 212615 **0.000% (0.000)**, PBE 212616
0.013% (0.016), both only on the spin-mixing operations (`base_full_group_hchk_nocut.log`; bug #22 has the parts).

## 4. MnBi2Te4

**Per-WF centres are not C3-covariant** (bug #21; `archive/diagnostics/2026-10-11_ws_positions/`). In 210651 the
Bi and Te p_x/p_y functions sit 0.06-0.15 Å off their atoms in the plane (p_z and s on the axis, Mn d on the atom).
Inversion and C2′ map centres onto centres (3e-8 Å); C3 does not (0.144 Å). Interpolated eigenvalues at 7 generic k
under the 12 operations of −3m′:

| fold positions, tolerance | max \|E(gk) − E(k)\| | \|Tr Ω\| / Σₙ\|Ωₙ\| |
|---|---|---|
| WF centres, 1e-3 | 61.6 meV | 3.2e-6 |
| WF centres, 5e-3 | 64.0 meV | 2.8e-6 |
| atom-pinned, 5e-3 | 0.000 meV | 1.3e-6 |

A first full-group curvature table gave 8-45% "residuals". They are noise: the base keeps P·{T|½,½,½}, every band is
Kramers-degenerate and the traced curvature is zero (1e-6 of the per-band scale, ~4e2 Å²); the table divides by
max|Tr Ω| ≈ 1e-4 Å². On the pinned model, H is exactly inversion-covariant: |H(−k)| = |H(k)| on the same orbital of
the partner atom, phases k-independent to 2e-6, and a Kubo curvature from H alone agrees with `omega_internal`.
**Consequence:** MBT needs atom-pinned positions in `run_axion.py` (not yet an option); symmetry checks on the MBT
base must use eigenvalues (or the modes, where PT is broken).

**Mesh: cost of 8×8×8.** Timings of the 4×4×4 runs (64 k):
- NSCF on 64 cores, no pools (base 209653): 5h49m, 5.5 min per k. With `-nk 8` on 48 cores (the mode runs) it took
  16.7 h, so pools are slower here.
- pw2wannier90 (`-nk 8`, 48 cores): 3h25m.
- `tmp/`: 29 GB at 300 bands; tarball 13.9 GB.

| mesh | k | NSCF | pw2wannier90 | tmp/ + tarball per structure | base + 4 modes |
|---|---|---|---|---|---|
| 4×4×4 | 64 | 5.8 h | 3.4 h | 43 GB | 0.2 TB |
| 6×6×6 | 216 | ~20 h | ~11.5 h | ~145 GB | ~0.73 TB |
| 8×8×8 | 512 | ~47 h | ~27 h | ~340 GB | ~1.7 TB |

**Fewer bands** (corrected the same day: an earlier "nbnd ~180" forgot the 56 semicore bands below −7 eV, which the
NSCF must compute). From 210651's `.eig` and `.chk` (300 bands, window −7 to 20 eV, frozen to 10 eV; occupied 138,
gap 7.27 → 7.47 eV):

| outer window top | bands in window (min, max over k) | nbnd needed (+ ~10 buffer) | Wannier-subspace weight above it now (max k, of 116) |
|---|---|---|---|
| 20 eV (now) | 210, 232 | 288 → 300 | — |
| 17 eV | 160, 198 | ~265 | 0.57 |
| 15 eV | 140, 166 | ~235 | 1.11 |
| 13 eV | 128, 136 | ~205 | 1.47 |
| 12 eV | 118, 124 | ~190 | 1.84 |

Lowering the top only changes the model above the frozen window (10 eV, 2.7 eV above the gap), where it moves about one
band's worth of weight. At ~205 bands: NSCF and `tmp/` ~⅓ less, the overlaps ~½ less. Not tested yet: rerun wannier90 on
210651's `.amn`/`.mmn` with `dis_win_max = 13` (minutes) and compare bands and spreads. Decision: the user's.

## 5. YIO mode 2 at Q = 0.01 from the Q = 0.005 density

`Y2Ir2O7/phonon/mode_2/Q_0p01_from_Q0p005_2026-10-11/`, job 212707, pw.x `db4b99b4`. Input = the Q = 0.005
production input (local-TF 0.05, ndim 20, conv_thr 1e-7, tprnfor, `startingpot = 'file'`) with the Q = 0.01 Ir
positions. `tmp.tar.gz` holds only the Q = 0.005 SCF charge density, `paw.txt`, `occup.txt` and the xml (36 MB). The
NSCF (512 k, `-nk 8`) runs only if the SCF converges. First check (19:55): |m|abs 3.39, total 0, accuracy 7e-5 Ry.

**Result (2026-10-12): not converged.** 300 iterations, stopped at 3.0e-5 Ry. It touched 1.2e-7 at iteration 51
(conv_thr 1e-7), then wandered between 1e-6 and 2e-4 Ry with E oscillating ±0.2 mRy; total magnetization 0, |m|abs
3.39 throughout. The same plateau as bug #17 (PBE + segni fix); the NSCF did not run. LDA SCFs on the same binary
converge in ~22 iterations.

### 5.1 MBT T₂⁻ reruns on `db4b99b4` (2026-10-12)

All four converged (11-16 iterations), 16 k, total magnetization 0, |m|abs 9.99. Energy vs force against base 209653
(E = −6947.39550993 Ry; F_base·Δr = 0 for a non-totally-symmetric mode), Q = 0.01 Å:

| mode (job) | atoms moved | E − E_base | −½ F·Δr | ratio |
|---|---|---|---|---|
| 1 (212708), Mn | 1, 2 | +0.0245 mRy | +0.0197 mRy | 1.24 |
| 2 (212709) | 3-6 | +0.1138 mRy | +0.1099 mRy | 1.04 |
| 3 (212710) | 7-10 | +0.0757 mRy | +0.0749 mRy | 1.01 |
| 4 (212711) | 11-14 | +0.0470 mRy | +0.0467 mRy | 1.01 |

All positive (the bug #19 runs were 0.2-0.3 mRy *below* the base). Mode 1's 1.24 is within the force precision: its
gap is 4.8e-6 Ry, and the SCF force correction (3.9e-4 Ry/bohr) × |Δr| (0.027 bohr) is ~1e-5 Ry.

## 6. Cluster storage

`~/mnt/w/axion` uses **6.5 TB** (target < 2 TB): MnBi2Te4 1.9, SrTiO3 0.44, Y2Ir2O7 4.1. QE scratch directories
alone hold 4.7 TB. `scripts/cluster_cleanup_2026-10-11.py` deletes `tmp/`, `tmp_*`, `*.save`, tarballs, `.mmn`,
`.amn`, `UNK` everywhere except the keep list:
- **YIO:** YP and YL2 bases and NSCF `tmp/` + W90 folders; modes 1 and 3 from YP; the new mode-2 run; only the wave
  functions of 210654 (Q = 0.005) and 210623.
- **MBT:** base 209653 `tmp/` (start for a denser NSCF) and its W90 tarballs; the kfix reruns; the base-state tarballs.
- **SrTiO3:** one restart per DFT run.

It deletes the pre-segni YIO runs (192344, 207546, 178365, 192345, 207547, 174937, 178377, 192347, 207549, nk tests),
all MBT PBE runs (#18), the MBT LDA T₂⁻ runs with #19, killed/test NSCFs (210587, 210588), every `tmp_bands` copy, and
nested tarball copies. Outputs (`*.out`, `.wout`, `.chk`, `_hr.dat`, `_r.dat`, NOTE) stay. A dry run listed 1420 items,
~6.0 TB before the SrTiO3 and frozfix keeps were added. Running it from this session was blocked by the permission
check; the user runs it (dry run first). It writes `~/cleanup_2026-10-11_manifest.tsv` and a NOTE line in every
touched folder.

## 7. Upstream drafts

`notes/upstream/`: QE MRs for new_ns_nc (bug #7 says an MR was sent in 2026-09; check), segni (#16) and
irreducible_BZ (#19), a QE discussion issue for the noncollinear GGA state (#18), Wannier90 issues for the centre wrap
(#20) and PR 702 (#5), and a PythTB issue for the `lattice` copy (#2). The three QE patches apply cleanly to develop
c5a3cf95b.
