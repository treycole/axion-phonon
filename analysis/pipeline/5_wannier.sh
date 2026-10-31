#!/bin/bash
# Step 5: Wannier functions (wannier90.x, projection only) for one window. Run on the login node; wait for
# windows/<window>/<structure>/OK. The base runs first: a displaced structure's trial functions are the base's Wannier
# functions for the same window.
#   input   the step-3 .win with dis_froz_max = FROZ_MAX and dis_win_max = WIN_MAX. The window is checked here, on the
#           login node, against the structure's NSCF eigenvalues (and, for a displaced structure, the base's).
#   job     base: the pw2wannier .amn (atomic trials). Displaced: mode_amn.py, A = O V_base from step 4 and the base
#           .chk. Then wannier90.x -np 1 with the step-3 .mmn and .eig.
#   output  windows/<window>/<structure>/out/: <PREFIX>.chk(.gz), .wout, .eig, .win, .nnkp (everything step 6 reads);
#           displaced structures also keep their .amn.gz
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
NK=8
FROZ_MAX=10.278          # eV (recorded in the folder's NOTE)
WIN_MAX=15.169
TOL=0.005                # eV: every window edge at least this far from every NSCF eigenvalue
QUEUE=dko48m4,dko64m4
SLOTS=8                  # wannier90.x runs on one rank (-np 1); the slots reserve its memory (8 x 3.6 GB)
H_RT=24:00:00
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
WINDOW=f${FROZ_MAX}_w${WIN_MAX}_k$NK        # the folder of $MAT/windows/ for this window and mesh
D=$MAT/windows/$WINDOW/$STRUCTURE
NAME=$STRUCTURE
ready $D
P=$(need_ok $MAT/dft/$STRUCTURE/pw2wannier_k$NK)
NS=$(need_ok $MAT/dft/$STRUCTURE/nscf_k$NK)
sed -E -e "s/^(\s*dis_froz_max\s*=\s*).*/\1$FROZ_MAX/" -e "s/^(\s*dis_win_max\s*=\s*).*/\1$WIN_MAX/" \
  $MAT/dft/$STRUCTURE/pw2wannier_k$NK/$PREFIX.win > $D/$PREFIX.win
grep -Eq "^\s*dis_froz_max\s*=\s*$FROZ_MAX$" $D/$PREFIX.win && grep -Eq "^\s*dis_win_max\s*=\s*$WIN_MAX$" $D/$PREFIX.win \
  || die "could not set the window in $D/$PREFIX.win"
[ "$(diff $MAT/dft/$STRUCTURE/pw2wannier_k$NK/$PREFIX.win $D/$PREFIX.win | grep -c '^>' || true)" -le 2 ] || die "the .win changed in more than the window"
cp -p $MAT/dft/$STRUCTURE/pw2wannier_k$NK/$PREFIX.nnkp $D/
if [ $STRUCTURE = base ]; then
  $PY $CODE/pipeline/inputs.py window $NS/$PREFIX.nscf.out $D/$PREFIX.win $TOL > $D/window_check.txt \
    || { cat $D/window_check.txt; die "window"; }
  AMN_LINES=("get $P $PREFIX.amn")
  KEEP_AMN="rm -f $PREFIX.amn                              # the base .amn stays in pw2wannier"
else
  NSB=$(need_ok $MAT/dft/base/nscf_k$NK)
  O=$(need_ok $MAT/dft/$STRUCTURE/overlap_k$NK)
  BW=$(need_ok $MAT/windows/$WINDOW/base)
  $PY $CODE/pipeline/inputs.py window $NS/$PREFIX.nscf.out $D/$PREFIX.win $TOL $NSB/$PREFIX.nscf.out > $D/window_check.txt \
    || { cat $D/window_check.txt; die "window"; }
  fill $CODE/pipeline/mode_amn.py $D/mode_amn.py "PREFIX=\"$PREFIX\""
  AMN_LINES=("cp -p \$SGE_O_WORKDIR/mode_amn.py . && [ -s mode_amn.py ] || fail \"input mode_amn.py\""
             "mkdir overlap && cp -p $O/chunks/chunk_*.npz* overlap/ || fail \"overlap chunks\""
             "find overlap -name '*.gz' -print0 | xargs -0 -r -P 8 gunzip -f   # the epilog gzips files > 100 MB"
             "get $BW $PREFIX.chk && mv $PREFIX.chk base.chk"
             "export OMP_NUM_THREADS=$SLOTS OPENBLAS_NUM_THREADS=$SLOTS MKL_NUM_THREADS=$SLOTS"
             "date; $PY mode_amn.py; rc=\$?"
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"
             "rm -rf overlap base.chk                  # before any fail: inputs must not be copied back"
             "[ \$rc = 0 ] || { rm -f $PREFIX.mmn; fail \"mode_amn.py\"; }")
  KEEP_AMN="pigz -p \$NSLOTS $PREFIX.amn                     # the computed input: kept"
fi
cat $D/window_check.txt

{ header w90${NK}_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "for f in $PREFIX.win $PREFIX.nnkp; do cp -p \$SGE_O_WORKDIR/\$f . && [ -s \$f ] || fail \"input \$f\"; done" \
    "md5sum $W90" "get $P $PREFIX.mmn" "get $P $PREFIX.eig" "${AMN_LINES[@]}" \
    "date; mpirun -np 1 $W90 $PREFIX > $PREFIX.wout" \
    "rm -f $PREFIX.mmn; $KEEP_AMN" \
    "grep -q 'All done' $PREFIX.wout || fail \"wannier90.x\"" \
    "[ -s $PREFIX.chk ] || fail \"no .chk\"" \
    "grep 'Omega I ' $PREFIX.wout | tail -1" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "wannier90 k$NK of $NAME, dis_froz_max $FROZ_MAX eV, dis_win_max $WIN_MAX eV: .mmn/.eig $P$( [ $STRUCTURE = base ] && echo ", atomic .amn" || echo ", .amn = O V_base from $O and $BW"). Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
