#!/bin/bash
# Symmetry table of Y2Ir2O7 for the paper (tab:symmetry-yio) and Fig. 6(b), 2026-11-03: run_yio_symmetry.py on the
# base and one displaced structure at the production window, with the models of the production dtheta. One job per
# structure; run on the login node with run_yio_symmetry.py next to this script; wait for OK (every threshold met) or
# FAIL (a crash, or a threshold: the end of symmetry.log says which) in windows/<window>/symmetry_<structure>/.
# The table is the job's symmetry.log (the base rows of every job agree); copy the folder to data/.
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/Y2Ir2O7
STRUCTURE=GM2-_p1_Q0.01  # a folder of $MAT/windows/$WINDOW/: GM2-_p1_Q0.01, GM2-_p2_Q0.01, GM2-_p3_Q0.01
PREFIX=Y2Ir2O7
WINDOW=f12.422_w22.644_k8
N_OCC=156
EIG_BANDS=260            # lowest DFT bands in test 1
N_GENERIC=7              # random k
QUEUE=dko48m4,dko64m4
SLOTS=24                 # memory: the A(R) fold of a 176-WF 8^3 run needs > 28.8 GB (as step 6)
H_RT=12:00:00
SUBMIT=true
# ======================================================================================
source /mnt/wkdk/bc889/axion/code/2026-10-29/pipeline/common.sh
check_tools
HERE=$(cd $(dirname $0) && pwd)
D=$MAT/windows/$WINDOW/symmetry_$STRUCTURE
NAME=sym_$STRUCTURE
ready $D
BW=$(need_ok $MAT/windows/$WINDOW/base)
MW=$(need_ok $MAT/windows/$WINDOW/$STRUCTURE)
BS=$(need_ok $MAT/dft/base/scf)          # the pw.x outputs: each structure's symmetry group
MS=$(need_ok $MAT/dft/$STRUCTURE/scf)
[ -s $HERE/run_yio_symmetry.py ] || die "run_yio_symmetry.py must sit next to $0"
fill $HERE/run_yio_symmetry.py $D/run_yio_symmetry.py "PREFIX=\"$PREFIX\"" "N_OCC=$N_OCC" "EIG_BANDS=$EIG_BANDS" \
    "N_GENERIC=$N_GENERIC"

{ header $NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "export OMP_NUM_THREADS=$SLOTS OPENBLAS_NUM_THREADS=$SLOTS MKL_NUM_THREADS=$SLOTS" \
    "cp -p \$SGE_O_WORKDIR/run_yio_symmetry.py . && [ -s run_yio_symmetry.py ] || fail \"input run_yio_symmetry.py\"" \
    "mkdir base mode" \
    "for f in $PREFIX.chk $PREFIX.eig $PREFIX.win $PREFIX.nnkp; do (cd base && get $BW \$f) && (cd mode && get $MW \$f) || fail \"\$f\"; done" \
    "(cd base && get $BS $PREFIX.scf.out) && (cd mode && get $MS $PREFIX.scf.out) || fail \"$PREFIX.scf.out\"" \
    "date; $PY run_yio_symmetry.py > symmetry.log 2>&1; rc=\$?" \
    "rm -rf base mode; cat symmetry.log" \
    "[ \$rc = 0 ] && [ -s symmetry_base.npz ] && [ -s symmetry_mode.npz ] || fail \"run_yio_symmetry.py (crash or threshold: tail symmetry.log)\"" \
    "ok"; } > $D/send_job.sh

submit $D send_job.sh
note $D "symmetry of $STRUCTURE and base, $WINDOW, n_occ $N_OCC, lowest $EIG_BANDS DFT bands, $N_GENERIC random k: base $BW, structure $MW; groups from $BS, $MS."
cp -p $0 $D/
