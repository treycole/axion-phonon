# Which QE and Wannier90 builds to use

The cluster has several QE and Wannier90 builds, and some of them are known to give wrong answers. Run
`md5sum pw.x` before debugging anything at source level.

**Frozen copies for job scripts (2026-10-19):** `/mnt/wkdk/bc889/axion/bin/qe_db4b99b4/` (pw.x, pw2wannier90.x,
bands.x, projwfc.x, ld1.x) and `bin/wannier90_26c8b506/`, with `MD5SUMS` and a `README` mapping every md5 seen in run
folders to its build. Pseudopotentials: one copy each in `/mnt/wkdk/bc889/axion/pseudo/{rel-pz,rel-pbe,pbe}/`
(`MD5SUMS`, `README`). Use absolute paths (`pseudo_dir`, `$BIN/pw.x`, `md5sum $BIN/pw.x` in the job log) and
`-l hostname=!n138`: jobs run in `/scratch/<job>` with only the submission folder's top-level files copied in, so
relative paths and subfolders don't reach them, and n138 has no `/mnt/wkdk` (test jobs 214399 failed there, 214400 passed).
Old run folders: `BINARIES.md5` records the md5 of the executable copies before they were deleted.

## Quantum ESPRESSO

| build | where | pw.x md5 | use for |
|---|---|---|---|
| **production** (QE develop, 8.0dev) | cluster `~/q-e-kfix`, branch `production/axion-phonon`, commit `fa9b19c8f` | **`db4b99b4`** | **all production runs since 2026-10-11**, MnBi₂Te₄ and Y₂Ir₂O₇ |
| segni-fixed | cluster `~/q-e-segnifix`, branch `fix/paw-segni-lsign`, commit `4b480ac68` | `e270b05b` | superseded: has the `irreducible_BZ` bug (#19). Same k lists as `db4b99b4` for YIO and the MBT base, so those runs stay valid |
| `~/q-e` | cluster | `c41c409` | superseded: has the segni bug (#16); `v_rad` fix in |
| any build with pw.x md5 `1541022` | | `1541022` | never: `v_rad`/`g_rad` bug (#12) |
| QE 7.5 / 7.6 stock | modules | — | never: `new_ns_nc` bug (#7) |

The production build carries four fixes on top of upstream (three of them not yet upstream; drafts in
[`../upstream/`](../upstream/README.md)):

1. **`new_ns_nc`** symmetrizer pairs the rotation with the right atom ([bug #7](../bugs/07_qe_new_ns_nc_symmetrizer.md)).
   Patch: [`qe_patches/new_ns_nc_preimage.patch`](../../qe_patches/new_ns_nc_preimage.patch).
2. **PAW `v_rad`/`g_rad`** initialized in noncollinear GGA ([bug #12](../bugs/12_qe_paw_vrad_grad_uninitialized.md)).
3. **PAW `segni_rad` `lsign` guard** ([bug #16](../bugs/16_qe_paw_segni_stale_ux.md)). Patch:
   `archive/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/paw_segni_lsign.patch`.
4. **`irreducible_BZ`** unfolds with plain lattice rotations, not with k → −S k for time-reversed ones
   ([bug #19](../bugs/19_qe_sym_ftrans_trev_energy.md)). Patch:
   [`qe_patches/irreducible_BZ_unfold_trev.patch`](../../qe_patches/irreducible_BZ_unfold_trev.patch).

**Source control.** `~/Repos/q-e` (local) holds `production/axion-phonon` = `fix/new_ns_nc-preimage` (934f8cdd3) +
4b480ac68 (segni) + fa9b19c8f (irreducible_BZ), the exact source of `db4b99b4`; the cluster tree `~/q-e-kfix` is on
the same branch. All are on the GitLab fork (`origin`, pushed 2026-10-11): `production/axion-phonon`,
`fix/new_ns_nc-preimage`, and the one-commit merge-request branches on `upstream/develop`,
`fix/paw-segni-lsign-develop` and `fix/irreducible_BZ-unfold-trev`.

**Checked 2026-10-11** (dry runs at low cutoff, `~/kcheck_2026-10-11`): `db4b99b4` and `e270b05b` print byte-identical
k lists and weights for the YIO base and modes 1-3 and the MBT base; only the MBT T2⁻ modes change (13 → 16 k).

`ld1.x` (for generating pseudopotentials) is built in the same tree; it needs `module load intel/2024 intel/ompi`.

Noncollinear +U runs as `kind = 0` (Dudarev) in this build, against `kind = 1` in QE 7.2. At J = 0 they are
equivalent.

## Wannier90

- **Wannierization (`.chk`):** `~/wannier90-dev` wannier90.x `26c8b506`, which reports **Release 4.0.3** (develop v4.0.1-160-gc4a924a5) in every production `.wout` (checked 214043, 214125 on 2026-10-19; this line said 3.1.0 before). Projection only: `num_iter = 0`, `dis_num_iter = 0`.
- **Plot stage (`_hr.dat`/`_r.dat`):** Jae-Mo Lihm's fork, branch `plan5-write-ndegen-applied`, cluster
  `~/wannier90-jaemo/wannier90-jaemo.x` (called by full path). It is used with `restart = plot` from the 3.1.0 `.chk`
  and the three flags `use_ws_distance`, `transl_inv_full`, `write_ndegen_applied`. Recipe and reasons:
  [`wannier/README.md`](../../wannier/README.md).
- **Production A(R)** is `transl_inv_full_from_chk` with `SHARED_WS = True`: Jae-Mo's formula rebuilt in Python from
  the `.chk` (CLAUDE.md pitfall 2). State the A(R) source with every number.
- **Run `wannier90.x` with `mpirun -np 1`** on a fully reserved node. It replicates its arrays on every rank, so
  `-np $NSLOTS` multiplies the memory (job 202714 was killed this way). `pw.x` and `pw2wannier90.x` keep
  `-np $NSLOTS`.

## Per-run files to keep

From CLAUDE.md's cluster practices: one good DFT run per structure, with `tmp.tar.gz` (for QE restarts and
`pw2wannier90.x`), an untarred `tmp/` (for the Julia rigid-shift `.amn`), and one `.mmn`. Delete extra copies.
The cluster total should stay under 2 TB.
