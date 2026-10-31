#!/bin/bash
# Step 4 (displaced structures only): PAW overlaps O(k) = <psi_structure,mk| S |psi_base,nk> between the structure's
# and the base's NSCF states, every k (bug #23: replaces generate_amn.jl). Run on the login node; wait for
# dft/<structure>/overlap_k<NK>/OK. The window does not enter: step 5 turns O into the .amn for any base Wannier run (A = O V_base).
#   job     copies both saves into scratch, runs N_PROCESSES overlap.py processes over k slices.
#   output  dft/<structure>/overlap_k<NK>/out/chunks/chunk_*.npz
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
NK=8
N_PROCESSES=8            # overlap.py processes, THREADS threads each: N_PROCESSES x THREADS = SLOTS
THREADS=6
QUEUE=dko48m4,dko64m4
SLOTS=48
H_RT=48:00:00
SCRATCH_GB=300           # both saves (MBT 2 x 87 GB at 8^3)
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
[ $STRUCTURE != base ] || die "step 4 is for displaced structures only"
[ $((N_PROCESSES * THREADS)) = $SLOTS ] || die "N_PROCESSES x THREADS != SLOTS"
D=$MAT/dft/$STRUCTURE/overlap_k$NK
NAME=$STRUCTURE
ready $D
NB=$(need_ok $MAT/dft/base/nscf_k$NK)
NM=$(need_ok $MAT/dft/$STRUCTURE/nscf_k$NK)
NKPT=$((NK * NK * NK))
cmp -s <(sed -n '/^K_POINTS/,$p' $MAT/dft/base/nscf_k$NK/$PREFIX.nscf.in) <(sed -n '/^K_POINTS/,$p' $MAT/dft/$STRUCTURE/nscf_k$NK/$PREFIX.nscf.in) \
  || die "the two NSCF inputs have different k lists"
fill $CODE/pipeline/overlap.py $D/overlap.py "PREFIX=\"$PREFIX\"" "N_TASKS=$N_PROCESSES"

{ header ovl${NK}_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "cp -p \$SGE_O_WORKDIR/overlap.py . && [ -s overlap.py ] || fail \"input\"" \
    "scratch $SCRATCH_GB" \
    "date; get_save $NB/tmp/$PREFIX.save base $NKPT" \
    "date; get_save $NM/tmp/$PREFIX.save mode $NKPT" \
    "date; for t in \$(seq 1 $N_PROCESSES); do TASK=\$t OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS $PY overlap.py > overlap_task\$t.log 2>&1 & done; wait" \
    "rm -rf base mode; tail -n 2 overlap_task*.log; date" \
    "[ \$(ls chunks/chunk_*.npz 2> /dev/null | wc -l) = $N_PROCESSES ] || fail \"\$(ls chunks 2> /dev/null | wc -l) of $N_PROCESSES chunks\"" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "PAW overlaps k$NK of $NAME onto the base: saves $NM and $NB; $N_PROCESSES processes. Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
