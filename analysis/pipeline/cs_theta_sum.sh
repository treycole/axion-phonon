#!/bin/bash
# The Chern-Simons 3-form check of a pair, once every cs_theta.sh task is OK. Runs here on the login node (seconds).
#   output  windows/<window>/dtheta_<structure>/cs/: theta_cs.log (the result), theta_cs.npz, OK
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/Y2Ir2O7
STRUCTURE=GM2-_p1_Q0.01  # a folder of $MAT/dft/
WINDOW=f12.422_w22.644_k8   # a folder of $MAT/windows/
Q=0.01                   # displacement amplitude of STRUCTURE, Angstrom
# ======================================================================================
source $(dirname $0)/common.sh
D=$MAT/windows/$WINDOW/dtheta_$STRUCTURE/cs
[ -e $D/OK ] && die "$D: already OK"
N_TASKS=$(grep "^N_TASKS = " $D/cs_theta.py | awk '{print $3}')
CHUNKS=""
for i in $(seq 1 $N_TASKS); do
  J=$(need_ok $D _$i)
  CHUNKS="$CHUNKS${CHUNKS:+, }\"$J/chunks/chunk_$(printf %04d $((i - 1))).npz\""
done
fill $CODE/pipeline/cs_theta_sum.py $D/cs_theta_sum.py "CHUNKS=[$CHUNKS]" "Q=$Q"
(cd $D && PYTHONPATH=$CODE $PY cs_theta_sum.py > theta_cs.log 2>&1) || { cat $D/theta_cs.log; die "cs_theta_sum.py"; }
cat $D/theta_cs.log
echo "cs_theta_sum.sh $(date +%F\ %H:%M)" > $D/OK
echo "$(date +%F\ %H:%M) theta from $N_TASKS chunks: $(tail -1 $D/theta_cs.log)" >> $D/NOTE
