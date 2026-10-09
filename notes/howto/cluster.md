## Restart mode

To use restart mode in Quantum ESPRESSO, we need the temporary save directory from the previous run. This is specified in the `&control` section of the input file with the `restart_mode` and `outdir` parameters.

The cluster does not copy over directories, so we need to first tar the temporary save directory 

```bash
tar -czvf tmp.tar.gz tmp/
```

## Faster nodes with 48 cores

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