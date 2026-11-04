#!/bin/bash
# Chern-Simons 3-form check of a base/structure pair (theta = 0 materials: Y2Ir2O7), as N_TASKS jobs (qsub -v TASK=i,
# no array job), one k slice each. Not a step of the dθ chain; it reads step 5 of the structure and of the base. Run on
# the login node; wait for windows/<window>/dtheta_<structure>/cs/OK_1 .. OK_<N_TASKS>, then run cs_theta_sum.sh.
# Running it again submits only the tasks that have neither OK_i nor a job in flight.
#   job     copies both Wannier runs (.chk, .eig, .win, .nnkp) into scratch, runs cs_theta.py.
#   output  windows/<window>/dtheta_<structure>/cs/out_<i>/chunks/chunk_<i-1>.npz
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/Y2Ir2O7
STRUCTURE=GM2-_p1_Q0.01  # a folder of $MAT/dft/
PREFIX=Y2Ir2O7
WINDOW=f12.422_w22.644_k8   # a folder of $MAT/windows/: f<frozen top>_w<outer top>_k<mesh>
N_OCC=156
NKS="10, 12"             # nk^3 meshes for the zone sum (YIO: equal to 1e-6 at 10 and 12, log 10-17)
NK_TRIAL=6
N_TASKS=2
BATCH=32
QUEUE=dko48m4,dko64m4
SLOTS=16                 # memory = SLOTS x 3.6 GB (the 10-26 series used 60 GB per task)
H_RT=24:00:00
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
D=$MAT/windows/$WINDOW/dtheta_$STRUCTURE/cs
NAME=$STRUCTURE
mkdir -p $D
BW=$(need_ok $MAT/windows/$WINDOW/base)
MW=$(need_ok $MAT/windows/$WINDOW/$STRUCTURE)
fill $CODE/pipeline/cs_theta.py $D/cs_theta.py "PREFIX=\"$PREFIX\"" "N_OCC=$N_OCC" "NKS=[$NKS]" "NK_TRIAL=$NK_TRIAL" \
  "N_TASKS=$N_TASKS" "BATCH=$BATCH"

{ header cs_$NAME $SLOTS 3.6G $H_RT $QUEUE '_$TASK'; printf "%s\n" \
    "export OMP_NUM_THREADS=$SLOTS OPENBLAS_NUM_THREADS=$SLOTS MKL_NUM_THREADS=$SLOTS" \
    "cp -p \$SGE_O_WORKDIR/cs_theta.py . && [ -s cs_theta.py ] || fail \"input cs_theta.py\"" \
    "mkdir base mode" \
    "for f in $PREFIX.chk $PREFIX.eig $PREFIX.win $PREFIX.nnkp; do (cd base && get $BW \$f) && (cd mode && get $MW \$f) || fail \"\$f\"; done" \
    "date; $PY cs_theta.py > cs_task\$TASK.log 2>&1; rc=\$?" \
    "rm -rf base mode; cat cs_task\$TASK.log" \
    "[ \$rc = 0 ] && [ -s chunks/chunk_\$(printf %04d \$((TASK - 1))).npz ] || fail \"cs_theta.py\"" \
    "ok"; } > $D/send_job.sh

for i in $(seq 1 $N_TASKS); do
  if [ -e $D/OK_$i ]; then echo "task $i: OK (job $(cat $D/OK_$i))"; continue; fi
  if [ -e $D/SUBMITTED_$i ] && [ ! -e $D/FAIL_$i ]; then
    echo "task $i: job $(cat $D/SUBMITTED_$i) has no OK_$i or FAIL_$i yet (qstat; if it died, delete SUBMITTED_$i)"
    continue
  fi
  submit $D send_job.sh _$i "-v TASK=$i"
  [ "$SUBMIT" = true ] && sleep 5
done
note $D "3-form theta of $NAME and the base, $WINDOW, n_occ $N_OCC, nk $NKS, SCDM trials of the base on $NK_TRIAL^3, $N_TASKS tasks: base $BW, structure $MW."
cp -p $0 $D/
