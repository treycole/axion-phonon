# pipeline/

One script per step of the calculation, the same for every structure (`notes/restructure/design.md` §2, §5). Each
`<n>_<step>.sh` runs on the cluster login node: set the settings block at its top (the structure folder first), run
it, and wait for `OK` in the step folder before running the next step. A step refuses to run unless the steps it reads
are `OK` and their jobs have left the queue. Nothing is ever queued behind another job.

| step | script | folder it writes | reads | job |
|---|---|---|---|---|
| 1 | `1_scf.sh` | `dft/<structure>/scf/` | `scf/<prefix>.scf.in` (you write it) | pw.x |
| 2 | `2_nscf.sh` | `dft/<structure>/nscf_k<N>/` | step 1; the base's `nscf_k<N>/<prefix>.nscf.in` (you write it once) | pw.x -nk |
| 3 | `3_pw2wannier.sh` | `dft/<structure>/pw2wannier_k<N>/` | step 2; the base's `pw2wannier_k<N>/<prefix>.win` (you write it once) | pw2wannier90.x |
| 4 | `4_overlap.sh` | `dft/<structure>/overlap_k<N>/` | step 2 of the structure and of the base (displaced structures only) | `overlap.py` |
| 5 | `5_wannier.sh` | `windows/<window>_k<N>/<structure>/` | step 3; displaced: step 4 and the base's step 5 | `mode_amn.py`, wannier90.x |
| 6 | `6_c2.sh` | `windows/<window>_k<N>/dtheta_<structure>/c2/` | step 5 of the structure and of the base | `c2.py` × N_TASKS |
| 7 | `7_dtheta.sh` | `windows/<window>_k<N>/dtheta_<structure>/` | every step-6 task | `dtheta.py`, on the login node |
| – | `cs_theta.sh`, `cs_theta_sum.sh` | `windows/<window>_k<N>/dtheta_<structure>/cs/` | step 5 of the structure and of the base | `cs_theta.py` × N_TASKS; `cs_theta_sum.py` on the login node |
| – | `dfpt.sh` | `dft/<structure>/ph/` | step 1; `ph/<prefix>.ph.in` (you write it) | pw.x (the SCF again, with wavefunctions), ph.x |

Paths are relative to the material folder (`/mnt/wkdk/bc889/axion/<material>/`; the same tree in `data/<material>/`):
`dft/` holds what depends only on a structure, `windows/` what depends also on the Wannier window (layout of
`data/MnBi2Te4/README.md`). The base runs steps 1, 2, 3, 5. A displaced structure runs 1 to 7, after the base's step 5
for the same window. Two checks sit outside the chain: `dfpt.sh` (SrTiO₃: the Born-charge reference) and
`cs_theta.sh` (θ from the Chern–Simons 3-form of each structure, so [θ(structure) − θ(base)]/Q checks ∂θ/∂Q without
the frozen gauge; θ = 0 materials only, log 2026-10-14). Step 6 on a uniform mesh (Y₂Ir₂O₇): `NK_T=10 M=1 HALF=-1
N3=10` is the plain 10³ grid, and step 7's sum is then its mean.

**What you write by hand:** each structure's `dft/<structure>/scf/<prefix>.scf.in`, and for the base once per mesh its
`nscf_k<N>/<prefix>.nscf.in` and `pw2wannier_k<N>/<prefix>.win`. A displaced structure's NSCF input and `.win` are the
base's with its own atomic positions (`inputs.py`); a displaced structure differs from the base only in its positions.

**Other materials and trial sets** (added 2026-10-31 for the SrTiO₃ Born-charge validation; empty or default, every
step behaves as before):
- `PSEUDO` (step 1): the pseudopotential set, a folder of `/mnt/wkdk/bc889/axion/pseudo/` with its `MD5SUMS`
  (default `rel-pz`; SrTiO₃ uses `oncv-pbe-rel`, norm-conserving, which step 4 handles with no augmentation term;
  `dfpt.sh` uses the set of step 1).
- `TRIAL` (steps 3, 5): another Wannier trial set of the same structures. Step 3 with `TRIAL` writes only that set's
  base `.amn`, from `dft/base/pw2wannier_k<N>_<TRIAL>/<prefix>.win` (written by you: its projections and lower window
  edges); the `.mmn` and `.eig` of `pw2wannier_k<N>/` serve every set. Step 5 with `TRIAL` reads that `.win` and
  `.amn` (a displaced structure's `.win` is made from it) and writes `windows/<TRIAL>_f<froz>_w<win>_k<N>/`.

**Markers** in each step folder: `SUBMITTED` (job id, at submission), then `OK` or `FAIL` (job id, written by the
job), and a `NOTE` line per submission. The job's files come back in `<step folder>/<job id>bc889/`. To redo a step
that is `OK`, move its folder aside (`archive/`); a `FAIL` step is resubmitted by running the script again.

**Every job** runs in node scratch, copies its own inputs there and checks them (it does not rely on the prolog),
prints the md5 of its binary, writes `OK` or `FAIL`, and exits 1 on failure. Every step script checks the md5s of
the binaries, the pseudopotentials and the code first (`common.sh`).

**Code on the cluster:** `/mnt/wkdk/bc889/axion/code/2026-10-29/` holds `axion_phonon/`, `pipeline/` and `pythtb/`
(external 5776792), with `MD5SUMS`. The jobs run `pipeline/*.py` from a copy in the step folder, with the settings
filled in by the step script, and import `axion_phonon` from there.

**Dry run:** `SUBMIT=false` writes the job script and prints the `qsub` line without submitting.

Files: `common.sh` (paths of the binaries, pseudopotentials and code; the helpers), `inputs.py` (displaced
structure's NSCF input and `.win`; the window check), `overlap.py`, `mode_amn.py`, `c2.py`, `dtheta.py` (the Python run
by the jobs; each reads its inputs from the job's folder and calls `axion_phonon`).
