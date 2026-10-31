#!/bin/bash
# Step 2: NSCF on the Wannier mesh (pw.x), from the step-1 density. Run on the login node; wait for dft/<structure>/nscf_k<NK>/OK.
#   input   the base's dft/<structure>/nscf_k<NK>/<PREFIX>.nscf.in (written by you, once per material and mesh: explicit k list,
#           nbnd, nosym). For a displaced structure this script writes its own copy with the structure's positions.
#   job     copies the input and the SCF save, runs pw.x -nk POOLS, checks it finished on every k; keeps the save.
#   output  dft/<structure>/nscf_k<NK>/out/: <PREFIX>.nscf.out, tmp/<PREFIX>.save (the Bloch states for steps 3 and 4)
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
NK=8                     # Wannier mesh NK x NK x NK
POOLS=8
QUEUE=dko64m4
SLOTS=64
H_RT=248:00:00
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
D=$MAT/dft/$STRUCTURE/nscf_k$NK
TEMPLATE=$MAT/dft/base/nscf_k$NK/$PREFIX.nscf.in
NAME=$STRUCTURE
ready $D
SCF=$(need_ok $MAT/dft/$STRUCTURE/scf)
[ -s $TEMPLATE ] || die "no $TEMPLATE"
$PY $CODE/pipeline/inputs.py nscf $TEMPLATE $MAT/dft/$STRUCTURE/scf/$PREFIX.scf.in > $D/$PREFIX.nscf.in.new
if [ $STRUCTURE = base ]; then
  cmp -s $D/$PREFIX.nscf.in.new $TEMPLATE || die "$TEMPLATE: positions differ from the base scf.in"
  rm $D/$PREFIX.nscf.in.new
else
  mv $D/$PREFIX.nscf.in.new $D/$PREFIX.nscf.in
fi
IN=$D/$PREFIX.nscf.in
grep -Eq "^\s*calculation\s*=\s*'nscf'" $IN || die "$IN: calculation is not 'nscf'"
for key in prefix outdir pseudo_dir; do
  [ "$(grep -E "^\s*$key\s*=" $IN)" = "$(grep -E "^\s*$key\s*=" $MAT/dft/$STRUCTURE/scf/$PREFIX.scf.in)" ] || die "$key differs from the scf.in"
done
NKPT=$(grep -A1 -E "^K_POINTS" $IN | tail -1 | awk '{print $1}')
[ "$NKPT" = $((NK * NK * NK)) ] || die "$IN: K_POINTS has $NKPT points, not $NK^3"

{ header nscf${NK}_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "cp -p \$SGE_O_WORKDIR/$PREFIX.nscf.in . && [ -s $PREFIX.nscf.in ] || fail \"input\"" \
    "md5sum $QE/pw.x" \
    "get_save $SCF/tmp/$PREFIX.save tmp 0" \
    "date; mpirun -np \$NSLOTS $QE/pw.x -nk $POOLS -input $PREFIX.nscf.in > $PREFIX.nscf.out" \
    "grep -q 'JOB DONE' $PREFIX.nscf.out || fail \"pw.x did not finish\"" \
    "nk=\$(grep 'number of k points=' $PREFIX.nscf.out | awk '{print \$5}'); [ \"\$nk\" = $NKPT ] || fail \"\$nk k points\"" \
    "[ \$(ls tmp/$PREFIX.save | grep -c '^wfc[0-9]*\\.dat\$') = $NKPT ] || fail \"wfc files\"" \
    "find tmp -mindepth 1 -maxdepth 1 ! -name $PREFIX.save -exec rm -rf {} +   # per-rank files: copies of the save" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "nscf k$NK of $NAME: $IN ($NKPT k), density from $SCF; -nk $POOLS, $QUEUE $SLOTS slots. Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
