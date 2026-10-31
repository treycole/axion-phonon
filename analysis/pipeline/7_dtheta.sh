#!/bin/bash
# Step 7: dtheta/dQ of a pair, once every step-6 task is OK. Runs here on the login node (seconds, no job).
#   output  dft/<structure>/windows/<window>/dtheta_<structure>/: dtheta.log (the result), c2.npz (k-resolved c2 at both ends), dtheta.npz, OK
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
NK=8
WINDOW=f10.962_w12.810_k8   # a folder of $MAT/windows/: f<frozen top>_w<outer top>_k<mesh>
# ======================================================================================
source $(dirname $0)/common.sh
D=$MAT/windows/$WINDOW/dtheta_$STRUCTURE
[ -e $D/OK ] && die "$D: already OK"
N_TASKS=$(grep "^N_TASKS = " $D/c2/c2.py | awk '{print $3}')
CHUNKS=""
for i in $(seq 1 $N_TASKS); do
  J=$(need_ok $D/c2 _$i)
  CHUNKS="$CHUNKS${CHUNKS:+, }\"$J/chunks/chunk_$(printf %04d $((i - 1))).npz\""
done
fill $CODE/pipeline/dtheta.py $D/dtheta.py "CHUNKS=[$CHUNKS]"
(cd $D && PYTHONPATH=$CODE $PY dtheta.py > dtheta.log 2>&1) || { cat $D/dtheta.log; die "dtheta.py"; }
cat $D/dtheta.log
echo "7_dtheta.sh $(date +%F\ %H:%M)" > $D/OK
echo "$(date +%F\ %H:%M) dtheta from $N_TASKS chunks: $(tail -1 $D/dtheta.log)" >> $D/NOTE
