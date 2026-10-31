## Restart mode

To use restart mode in Quantum ESPRESSO, we need the temporary save directory from the previous run. This is specified in the `&control` section of the input file with the `restart_mode` and `outdir` parameters.

The cluster does not copy over directories, so we need to first tar the temporary save directory 

```bash
tar -czvf tmp.tar.gz tmp/
```

## What is where

`/mnt/wkdk/bc889/axion/INDEX.md` lists every run folder with its status; each folder has a `_STATUS` file
(PRODUCTION / COMPARISON / VALIDATION / MIXED / SUPERSEDED). Regenerate with `scripts/cluster_mark_status_2026-10-17.sh`
after editing its list. Julia for compute-node jobs: `/mnt/wkdk/bc889/julia/` (`JULIA_DEPOT_PATH=.../depot`).

## Jobs run in node scratch, never in /mnt/wkdk (admin warning, 2026-10-23)

The cluster admin sent a serious warning about job 217301, which ran with `-cwd` and `cd /mnt/wkdk/...`. Every job
must instead work in the node's scratch:
- **What the prolog and epilog do** (`/opt/sge/default/prolog_i36.sh`):
  - The prolog creates `/scratch/<JOB_ID>bc889` on the node and copies the submit folder's *files* into it, but not its
    subfolders.
  - `source $TMPDIR/sge_init.sh` cds there.
  - At the end the epilog pigzips files over 100 MB and copies the folder back to `<submit folder>/<JOB_ID>bc889/`.
- **In the job script:**
  - No `-cwd` and no `cd /mnt/wkdk`.
  - Copy or zcat the large inputs into scratch first. Reading from `/mnt/wkdk` is fine.
  - Delete the inputs before the job exits, so that only outputs come back.
  - Write only small markers or logs directly to `/mnt/wkdk`.
- **No SGE array jobs.** All tasks of an array share one `JOB_ID`, so tasks on the same node share one scratch folder,
  and the first task to finish deletes it. Submit separate jobs with `qsub -v TASK=i` instead.
- **Fixed scripts:** `archive/setup_templates/MnBi2Te4/modes/setup_mbt_dtheta_tube.sh` and
  `setup_mbt_window_study.sh`.
  - The tube chain was stage, c2 array and assemble, run in place. It is now 10 separate c2 jobs in scratch and a sum
    job.
  - The chunks land in `<pair>/c2/<job>bc889/chunks/`.
- **Jobs that ran in place (before the fix):**
  - the tube chains of 10-19 to 10-22, whose staged `<pair>/{base,mode}/` copies are about 42 GB of uncompressed `.chk`
    and `_hr.dat` duplicates;
  - mode_symmetry (217301);
  - the 10-19 refined-BZ and 10-21 tube-convergence diagnostics.

## Nodes

- **dko48m4:** 48 cores. **n138 has no `/mnt/wkdk` mount: exclude it** (`qsub -l hostname='!n138' ...`). The
  submit directory is staged to the node's `$TMPDIR` and copied back by the epilog, so a job that never `cd`s into
  `/mnt/wkdk` itself would also run there, but excluding n138 is simpler.
  **Any node can lose the mount after an outage:** on 2026-10-26 n113 (dko64m4) came back without `/mnt/wkdk`, and a
  job reading its binaries from there failed at `mpirun` ("could not access or execute"). If a job dies at launch,
  check `ls /mnt/wkdk` on its node and add that node to the exclusion (`hostname='!n138&!n113'`) until it is fixed.
- **dko64m4:** n108, n109, n113, 64 cores, 251 GB, h_rt up to 248 h; often idle (`#$ -q dko64m4`, `#$ -pe orte_one 64`).

### Faster nodes with 48 cores (NVMe `/scratch`)

1. n126
2. n246
3. n247

Use:

```bash
qsub -l hostname=n126 send_job.sh
```

## Copying files from cluster 

Use `scripts/sync_from_cluster.sh` (set `CLUSTER`, `REMOTE` and `LOCAL` at the top, then
`bash scripts/sync_from_cluster.sh`). It copies only the file types in the allow list
`scripts/sync_filter.txt`, so executables (`*.x`), QE scratch (`tmp/`, `*.save/`), pseudopotential copies
and `.mmn`/`.amn` stay on the cluster. It shows what would come over (sizes, total, files over 100 MB) and
asks before copying. Fetch a one-off big file (a `.mmn`) by name with plain `rsync`.