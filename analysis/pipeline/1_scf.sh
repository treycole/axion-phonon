#!/bin/bash
# Step 1: SCF of one structure (pw.x). Run on the cluster login node; wait for dft/<structure>/scf/OK.
#   input   dft/<structure>/scf/<PREFIX>.scf.in, written by you: pseudo_dir = the shared pseudo folder, outdir './tmp/'.
#   job     copies the input (and START, the save to start from when the input has startingpot = 'file'; its
#           wavefunctions are not copied), runs pw.x,
#           checks convergence, the number of k points and the total and absolute magnetization; keeps the save
#           without wavefunctions (the NSCF needs only the density and the PAW and Hubbard occupations).
#   output  scf/out/: <PREFIX>.scf.out, tmp/<PREFIX>.save
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
START=/mnt/wkdk/bc889/axion/MnBi2Te4/dft/base/scf/out/tmp/MnBi2Te4.save   # or "" (input without startingpot)
N_KPOINTS=16             # k points pw.x must find: the orbit count (MBT base 13, T2- patterns 16; bug #19)
M_TOTAL=0.00             # expected |total magnetization| (Bohr mag/cell)
M_ABS=9.99               # expected absolute magnetization (Bohr mag/cell)
M_TOL=0.05               # allowed deviation of both
QUEUE=dko64m4
SLOTS=64
H_RT=248:00:00
SUBMIT=true              # false: write the job script and print the qsub line, submit nothing
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
D=$MAT/dft/$STRUCTURE/scf
IN=$D/$PREFIX.scf.in
NAME=$STRUCTURE
ready $D
[ -s $IN ] || die "no $IN"
grep -Eq "^\s*calculation\s*=\s*'scf'" $IN || die "$IN: calculation is not 'scf'"
grep -Eq "^\s*prefix\s*=\s*'$PREFIX'" $IN || die "$IN: prefix is not '$PREFIX'"
grep -Eq "^\s*outdir\s*=\s*'\./tmp/?'" $IN || die "$IN: outdir is not './tmp/'"
grep -Eq "^\s*pseudo_dir\s*=\s*'$PSEUDO/?'" $IN || die "$IN: pseudo_dir is not '$PSEUDO/'"
if grep -Eq "^\s*startingpot\s*=\s*'file'" $IN; then
  [ "$(basename "$START")" = $PREFIX.save ] && [ -s $START/charge-density.dat ] || die "START is not a $PREFIX.save with a density"
else
  [ -z "$START" ] || die "START is set but $IN does not read it (no startingpot = 'file')"
fi

{ header scf_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "cp -p \$SGE_O_WORKDIR/$PREFIX.scf.in . && [ -s $PREFIX.scf.in ] || fail \"input\"" \
    "md5sum $QE/pw.x" \
    "mkdir -p tmp/$PREFIX.save; [ -z \"$START\" ] || find $START -maxdepth 1 -type f ! -name 'wfc*' -exec cp -p {} tmp/$PREFIX.save/ \\;" \
    "find tmp -name '*.gz' -print0 | xargs -0 -r gunzip -f; [ -z \"$START\" ] || [ -s tmp/$PREFIX.save/charge-density.dat ] || fail \"start density\"" \
    "date; mpirun -np \$NSLOTS $QE/pw.x -input $PREFIX.scf.in > $PREFIX.scf.out" \
    "grep -q 'convergence has been achieved' $PREFIX.scf.out || fail \"not converged\"" \
    "nk=\$(grep 'number of k points=' $PREFIX.scf.out | awk '{print \$5}'); [ \"\$nk\" = $N_KPOINTS ] || fail \"\$nk k points, expected $N_KPOINTS\"" \
    "grep -E 'total magnetization|absolute magnetization' $PREFIX.scf.out | tail -2" \
    "mt=\$(grep 'total magnetization' $PREFIX.scf.out | tail -1 | awk '{s = 0; for (i = 4; i <= NF - 2; i++) s += \$i * \$i; print sqrt(s)}')" \
    "ma=\$(grep 'absolute magnetization' $PREFIX.scf.out | tail -1 | awk '{print \$4}')" \
    "awk -v a=\$mt -v b=$M_TOTAL -v t=$M_TOL 'BEGIN {exit !((a - b) ^ 2 <= t ^ 2)}' || fail \"|total magnetization| \$mt, expected $M_TOTAL\"" \
    "awk -v a=\$ma -v b=$M_ABS -v t=$M_TOL 'BEGIN {exit !((a - b) ^ 2 <= t ^ 2)}' || fail \"absolute magnetization \$ma, expected $M_ABS\"" \
    "rm -f tmp/$PREFIX.save/wfc*.dat; find tmp -mindepth 1 -maxdepth 1 ! -name $PREFIX.save -exec rm -rf {} +" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "scf of $NAME: $IN, start ${START:-none}, expect $N_KPOINTS k, |m| $M_TOTAL, |m|abs $M_ABS; $QUEUE $SLOTS slots. Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
