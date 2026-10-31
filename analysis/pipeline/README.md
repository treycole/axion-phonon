# pipeline/

One script per step of the calculation, the same for every structure (`notes/restructure/design.md` §2, §5). Each
`<n>_<step>.sh` runs on the cluster login node: set the settings block at its top (the structure folder first), run
it, and wait for `OK` in the step folder before running the next step. A step refuses to run unless the steps it reads
are `OK` and their jobs have left the queue. Nothing is ever queued behind another job.

| step | script | folder it writes (nested under what it was made from) | reads | job |
|---|---|---|---|---|
| 1 | `1_scf.sh` | `<structure>/scf/` | `scf/<prefix>.scf.in` (you write it) | pw.x |
| 2 | `2_nscf.sh` | `scf/nscf_k<N>/` | step 1; the base's `scf/nscf_k<N>/<prefix>.nscf.in` (you write it once) | pw.x -nk |
| 3 | `3_pw2wannier.sh` | `scf/nscf_k<N>/pw2wannier/` | step 2; the base's `pw2wannier/<prefix>.win` (you write it once) | pw2wannier90.x |
| 4 | `4_overlap.sh` | `scf/nscf_k<N>/overlap/` | step 2 of the structure and of the base (displaced structures only) | `overlap.py` |
| 5 | `5_wannier.sh` | `scf/nscf_k<N>/wannier/` (or `studies/windows/<window>/`; setting RUN) | step 3; displaced: step 4 and the base's step 5 | `mode_amn.py`, wannier90.x |
| 6 | `6_c2.sh` | `scf/nscf_k<N>/wannier/dtheta/c2/` | step 5 of the structure and of the base | `c2.py` × N_TASKS |
| 7 | `7_dtheta.sh` | `scf/nscf_k<N>/wannier/dtheta/` | every step-6 task | `dtheta.py`, on the login node |

The base runs steps 1, 2, 3, 5. A displaced structure runs 1 to 7, after the base's step 5 for the same window.

**What you write by hand:** each structure's `scf/<prefix>.scf.in`, and for the base once per mesh its
`scf/nscf_k<N>/<prefix>.nscf.in` and `scf/nscf_k<N>/pw2wannier/<prefix>.win`. A displaced structure's NSCF input and `.win` are the
base's with its own atomic positions (`inputs.py`); a displaced structure differs from the base only in its positions.

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
