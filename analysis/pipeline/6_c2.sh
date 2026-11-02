#!/bin/bash
# Step 6: c2 density of a base/structure pair on a k list, as N_TASKS jobs (qsub -v TASK=i, no array job). Run on the
# login node; wait for windows/<window>/dtheta_<structure>/c2/OK_1 .. OK_<N_TASKS>, then run 7_dtheta.sh. Running it
# again submits only the tasks that have neither OK_i nor a job in flight (a failed task is resubmitted).
#   k list  the tube list of axion_phonon.tube (written once per pair folder): an NK_T x NK_T grid in each of N3
#           planes, refined M x M within HALF cells of the small-gap line. HALF = -1 and M = 1, N3 = NK_T: the plain
#           NK_T^3 grid.
#   job     copies both Wannier runs (.chk, .eig, .win, .nnkp) into scratch, runs c2.py.
#   output  windows/<window>/dtheta_<structure>/: klist.npz, c2/out/chunks/chunk_<i-1>.npz
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
NK=8                     # Wannier mesh
WINDOW=f10.962_w12.810_k8   # a folder of $MAT/windows/: f<frozen top>_w<outer top>_k<mesh>
Q=0.01                   # displacement amplitude of STRUCTURE, Angstrom
N_OCC=82                 # occupied Wannier bands (MBT 82; YIO 156)
NK_T=16                  # tube k list (production MBT: 16, 4, 3, 8)
M=4
HALF=3
N3=8
N_TASKS=10
BATCH=50
QUEUE=dko48m4,dko64m4
SLOTS=8                  # memory = SLOTS x 3.6 GB: 8 for an 8^3 mesh; 24 for 10^3 (the fold over 27 x 1000 R needs
                         # more than 28.8 GB: the 10^3 tubes of 10-30 ran out of memory at 8)
H_RT=24:00:00
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
D=$MAT/windows/$WINDOW/dtheta_$STRUCTURE
NAME=$STRUCTURE
mkdir -p $D/c2
BW=$(need_ok $MAT/windows/$WINDOW/base)
MW=$(need_ok $MAT/windows/$WINDOW/$STRUCTURE)
(cd $D && PYTHONPATH=$CODE $PY - <<PY
from pathlib import Path
import numpy as np
from axion_phonon import tube
kl = tube.klist($NK_T, $M, $HALF, $N3)
if Path("klist.npz").exists():                       # a rerun: the list must be the one the chunks were made on
    old = np.load("klist.npz")
    assert all(np.array_equal(old[k], kl[k]) for k in kl), "klist.npz differs from these settings"
else:
    np.savez("klist.npz", **kl)
print(f"k list: {len(kl['k'])} points")
PY
)
cp -p $D/klist.npz $D/c2/
fill $CODE/pipeline/c2.py $D/c2/c2.py "PREFIX=\"$PREFIX\"" "N_OCC=$N_OCC" "Q=$Q" "N_TASKS=$N_TASKS" "BATCH=$BATCH"

{ header c2_$NAME $SLOTS 3.6G $H_RT $QUEUE '_$TASK'; printf "%s\n" \
    "export OMP_NUM_THREADS=$SLOTS OPENBLAS_NUM_THREADS=$SLOTS MKL_NUM_THREADS=$SLOTS" \
    "for f in c2.py klist.npz; do cp -p \$SGE_O_WORKDIR/\$f . && [ -s \$f ] || fail \"input \$f\"; done" \
    "mkdir base mode" \
    "for f in $PREFIX.chk $PREFIX.eig $PREFIX.win $PREFIX.nnkp; do (cd base && get $BW \$f) && (cd mode && get $MW \$f) || fail \"\$f\"; done" \
    "date; $PY c2.py > c2_task\$TASK.log 2>&1; rc=\$?" \
    "rm -rf base mode klist.npz; cat c2_task\$TASK.log" \
    "[ \$rc = 0 ] && [ -s chunks/chunk_\$(printf %04d \$((TASK - 1))).npz ] || fail \"c2.py\"" \
    "ok"; } > $D/c2/send_job.sh

for i in $(seq 1 $N_TASKS); do
  if [ -e $D/c2/OK_$i ]; then echo "task $i: OK (job $(cat $D/c2/OK_$i))"; continue; fi
  if [ -e $D/c2/SUBMITTED_$i ] && [ ! -e $D/c2/FAIL_$i ]; then
    echo "task $i: job $(cat $D/c2/SUBMITTED_$i) has no OK_$i or FAIL_$i yet (qstat; if it died, delete SUBMITTED_$i)"
    continue
  fi
  submit $D/c2 send_job.sh _$i "-v TASK=$i"
  [ "$SUBMIT" = true ] && sleep 5
done
note $D "c2 of $NAME vs base, $WINDOW, Q $Q, n_occ $N_OCC, tube ($NK_T, $M, $HALF, $N3), $N_TASKS tasks: base $BW, structure $MW."
cp -p $0 $D/
