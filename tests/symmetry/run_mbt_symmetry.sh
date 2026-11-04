#!/bin/bash
# Symmetry table of MnBi2Te4 for the paper (tab:symmetry-mbt), 2026-10-31: run_mbt_symmetry.py on the base and one
# displaced structure at the production window, with the models of the production dθ. One job per structure; run on
# the login node with run_mbt_symmetry.py and symmetry_lib.py next to this script; wait for OK (every threshold of
# symmetry_lib.THRESHOLDS met) or FAIL (a crash, or a threshold: the end of symmetry.log says which) in
# windows/<window>/symmetry_<structure>/. The table is the job's symmetry.log (the base rows of every job should agree);
# copy the folder to data/ and tests/symmetry/test_mbt_symmetry.py checks it.
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/windows/$WINDOW/: T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
WINDOW=f10.962_w12.810_k8
FROZEN_TOP=10.962        # eV, dis_froz_max of WINDOW
N_OCC=82
N_GENERIC=7              # random k
N_TUBE=5                 # k within 0.02 1/A of the k_z axis
QUEUE=dko48m4,dko64m4
SLOTS=8                  # memory: the A(R) fold needs > 14.4 GB (4 slots ran out, 11-03; as step 6)
H_RT=04:00:00
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
for f in run_mbt_symmetry.py symmetry_lib.py; do [ -s $HERE/$f ] || die "$f must sit next to $0"; done
cp -p $HERE/symmetry_lib.py $D/
fill $HERE/run_mbt_symmetry.py $D/run_mbt_symmetry.py "PREFIX=\"$PREFIX\"" "N_OCC=$N_OCC" "FROZEN_TOP=$FROZEN_TOP" "N_GENERIC=$N_GENERIC" \
    "N_TUBE=$N_TUBE"

{ header $NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "export OMP_NUM_THREADS=$SLOTS OPENBLAS_NUM_THREADS=$SLOTS MKL_NUM_THREADS=$SLOTS" \
    "for f in run_mbt_symmetry.py symmetry_lib.py; do cp -p \$SGE_O_WORKDIR/\$f . && [ -s \$f ] || fail \"input \$f\"; done" \
    "mkdir base mode" \
    "for f in $PREFIX.chk $PREFIX.eig $PREFIX.win $PREFIX.nnkp; do (cd base && get $BW \$f) && (cd mode && get $MW \$f) || fail \"\$f\"; done" \
    "(cd base && get $BS $PREFIX.scf.out) && (cd mode && get $MS $PREFIX.scf.out) || fail \"$PREFIX.scf.out\"" \
    "date; $PY run_mbt_symmetry.py > symmetry.log 2>&1; rc=\$?" \
    "rm -rf base mode; cat symmetry.log" \
    "[ \$rc = 0 ] && [ -s symmetry_base.npz ] && [ -s symmetry_mode.npz ] || fail \"run_mbt_symmetry.py (crash or threshold: tail symmetry.log)\"" \
    "ok"; } > $D/send_job.sh

submit $D send_job.sh
note $D "symmetry of $STRUCTURE and base, $WINDOW (frozen top $FROZEN_TOP eV), n_occ $N_OCC, $N_GENERIC random + $N_TUBE tube k: base $BW, structure $MW; groups from $BS, $MS."
cp -p $0 $D/
