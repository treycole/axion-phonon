# Which QE and Wannier90 builds to use

The cluster has several QE and Wannier90 builds, and some of them are known to give wrong answers. Run
`md5sum pw.x` before debugging anything at source level.

## Quantum ESPRESSO

| build | where | pw.x md5 | use for |
|---|---|---|---|
| **segni-fixed** (QE develop, 8.0dev) | cluster `~/q-e-segnifix`, branch `fix/paw-segni-lsign`, commit `4b480ac68` | **`e270b05b`** | **all production runs**, MnBi₂Te₄ and Y₂Ir₂O₇ |
| `~/q-e` | cluster | `c41c409` | superseded: has the segni bug (#16); `v_rad` fix in |
| any build with pw.x md5 `1541022` | | `1541022` | never: `v_rad`/`g_rad` bug (#12) |
| QE 7.5 / 7.6 stock | modules | — | never: `new_ns_nc` bug (#7) |

The production build carries three fixes on top of upstream:

1. **`new_ns_nc`** symmetrizer pairs the rotation with the right atom ([bug #7](../bugs/07_qe_new_ns_nc_symmetrizer.md)).
   Patch: [`calculations/qe_patches/new_ns_nc_preimage.patch`](../../calculations/qe_patches/new_ns_nc_preimage.patch).
2. **PAW `v_rad`/`g_rad`** initialized in noncollinear GGA ([bug #12](../bugs/12_qe_paw_vrad_grad_uninitialized.md)).
3. **PAW `segni_rad` `lsign` guard** ([bug #16](../bugs/16_qe_paw_segni_stale_ux.md)). Patch:
   `calculations/diagnostics/2026-10-05_ws_tolerance_window/qe_segni_patch/paw_segni_lsign.patch`.

`ld1.x` (for generating pseudopotentials) is built in the same tree; it needs `module load intel/2024 intel/ompi`.

Noncollinear +U runs as `kind = 0` (Dudarev) in this build, against `kind = 1` in QE 7.2. At J = 0 they are
equivalent.

## Wannier90

- **Wannierization (`.chk`):** Wannier90 3.1.0. Projection only: `num_iter = 0`, `dis_num_iter = 0`.
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
