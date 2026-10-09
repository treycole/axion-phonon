# Stage 1: everything up to `_hr.dat` and `_r.dat`

This folder is the only place that knows about Quantum ESPRESSO, `pw2wannier90`
and Wannier90. Its job ends when two files exist for each structure:

```
seed_hr.dat     H(R)  = <0s|H|Rt>
seed_r.dat      A(R)  = <0s|r|Rt>
(seed_wsvec.dat, seed_centres.xyz, seed.win  must sit beside them)
```

After that, PythTB reads them and `modules/curvature.py` + `modules/axion.py` do
all the physics. Nothing in `modules/` reads a `.chk`, `.mmn` or `.amn`.

The Python in this folder:

| file | what it is |
|---|---|
| `run_plot_stage.sh`, `check_outputs.py` | the current path: write the files, check them |
| `wannier_io.py`, `position_matrix.py` | the older in-house route: read `.chk`/`.mmn`/`.nnkp`, load `Lambda(R)`, build `A(R)` in Python. Not needed by the current path; kept for the `.mmn` route and the validation notebooks |

```
QE scf -> nscf -> wannier90 -pp -> pw2wannier90 -> .amn .mmn .eig
      -> wannier90.x (Wannierise, gives .chk)                 [both structures]
      -> wannier90.x restart=plot  with the three flags      -> _hr.dat _r.dat
         (Jae-Mo's fork)                                        _wsvec.dat _centres.xyz
```

## The recipe

1. **Base structure.** SCF, NSCF on a Gamma-centred `mp_grid` (8x8x8 for YIO),
   `nosym = .true.`, then `pw2wannier90` and `wannier90.x`. For YIO the base uses
   the projection gauge (`num_iter = 0`, `dis_num_iter = 0`).
2. **Displaced structure.** Same SCF/NSCF settings and the same `mp_grid`. The
   two Wannierisations must share a gauge or the finite difference in beta is
   noise: generate the displaced `.amn` from the *undisplaced* Wannier functions
   with Jae-Mo's rigid-shift tool, then run Wannier90 in the projection gauge
   only. See [`../notes/theory/rigid_shift.md`](../notes/theory/rigid_shift.md).
3. **Write H(R) and A(R)** with `run_plot_stage.sh` (below). This is a
   `restart = plot` run from the existing `.chk`, so it needs no `.mmn` or `.amn`.

## The final step, and why each flag

`run_plot_stage.sh RUN_DIR` copies the `.win` into `RUN_DIR/jaemo/`, symlinks
`.chk` and `.eig`, sets the keywords below and runs Jae-Mo's `wannier90.x`
(`~/Repos/wannier90-JaeMo`, branch `plan5-write-ndegen-applied`). The production
files in `RUN_DIR` are never touched.

| keyword | why |
|---|---|
| `restart = plot` | read the finished `.chk`; skip `.mmn`, `.amn` and re-Wannierisation (verified to read a Wannier90 3.1.0 `.chk`) |
| `use_ws_distance = true` | pair-dependent Wigner-Seitz rule; the centre-independent `min\|R\|` breaks crystal symmetry |
| `transl_inv_full = true` | the translation-equivariant position-matrix formula (symmetric link phase, `exp(-i b.R/2)` at the final image) |
| `write_ndegen_applied = true` | required with the two above; divides every Wigner-Seitz weight out before writing, so PythTB reads a plain sum over R |
| `write_hr = true`, `write_rmn = true` | the two files |
| `bands_plot`, `write_tb` = false | not needed, and they cost time and disk |

Leave any one of the first four out and the result is a different, worse
position matrix (see `notes/log/2026-09-16_jaemo_rdat_production_closure.md`).
Disk: about 6 GB per structure for 176 Wannier functions.

## Check what you wrote

```bash
python wannier/check_outputs.py RUN_DIR/jaemo                     # one structure
python wannier/check_outputs.py BASE/jaemo MODE/jaemo             # a base/mode pair
python wannier/check_outputs.py RUN_DIR/jaemo --quick             # headers only, no PythTB
```

It checks the flags in `_wsvec.dat` and the `.win`, that the `_hr.dat` degeneracy
block is all ones, and (loading through PythTB) that the `R = 0` diagonal of
`A(R)` equals the Wannier centres and that `A(R) = A(-R)^dag`. For a pair it
also reports how much of `H_mode - H_base` is a change of Wigner-Seitz image
assignment (`ws_pair_consistency`); that number is information, not a failure.

## Known open issue

Wannier90 picks the Wigner-Seitz images per run from that run's own centres, so a
displaced structure can spread a tied hopping differently from the base. This is
measured, present in the older `.mmn` route too, and its effect on dtheta is not
yet known: `../notes/log/2026-09-20_jaemo_ndegen_yio_mode1.md`.
