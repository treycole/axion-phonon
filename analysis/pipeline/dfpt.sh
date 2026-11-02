#!/bin/bash
# DFPT at q = 0 (ph.x): the reference Born effective charges and dielectric tensor of one structure. Not a step of the
# dθ chain; it reads step 1. Run on the cluster login node; wait for dft/<structure>/ph/OK.
#   input   dft/<structure>/ph/<PREFIX>.ph.in, written by you: outdir './tmp/', prefix as the SCF.
#   job     copies step 1's input, reruns that SCF (step 1 keeps no wavefunctions, which ph.x needs), checks the total
#           energy against step 1's, runs ph.x, checks it finished and printed the effective charges; keeps no save.
#   output  ph/out/: <PREFIX>.scf.out, <PREFIX>.ph.out
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/SrTiO3
STRUCTURE=base           # a folder of $MAT/dft/
PREFIX=SrTiO3
PSEUDO=/mnt/wkdk/bc889/axion/pseudo/oncv-pbe-rel   # as in step 1
E_TOL=1e-6               # Ry: allowed difference of the rerun SCF's total energy from step 1's
QUEUE=dko48m4
SLOTS=16
H_RT=24:00:00
SUBMIT=true              # false: write the job script and print the qsub line, submit nothing
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
[ "$(md5sum < $QE/ph.x | cut -c1-8)" = $PH_MD5 ] || die "ph.x md5"
D=$MAT/dft/$STRUCTURE/ph
IN=$D/$PREFIX.ph.in
NAME=$STRUCTURE
ready $D
[ -s $IN ] || die "no $IN"
SCF=$(need_ok $MAT/dft/$STRUCTURE/scf)
grep -Eq "^\s*prefix\s*=\s*'$PREFIX'" $IN || die "$IN: prefix is not '$PREFIX'"
grep -Eq "^\s*outdir\s*=\s*'\./tmp/?'" $IN || die "$IN: outdir is not './tmp/'"
E1=$(grep '^!' $SCF/$PREFIX.scf.out | tail -1 | awk '{print $5}')
[ -n "$E1" ] || die "$SCF/$PREFIX.scf.out: no total energy"
N_KPOINTS=$(grep 'number of k points=' $SCF/$PREFIX.scf.out | awk '{print $5}')

{ header ph_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "cp -p $MAT/dft/$STRUCTURE/scf/$PREFIX.scf.in \$SGE_O_WORKDIR/$PREFIX.ph.in . && [ -s $PREFIX.scf.in ] && [ -s $PREFIX.ph.in ] || fail \"inputs\"" \
    "md5sum $QE/pw.x $QE/ph.x" \
    "date; mpirun -np \$NSLOTS $QE/pw.x -input $PREFIX.scf.in > $PREFIX.scf.out" \
    "grep -q 'convergence has been achieved' $PREFIX.scf.out || fail \"SCF not converged\"" \
    "nk=\$(grep 'number of k points=' $PREFIX.scf.out | awk '{print \$5}'); [ \"\$nk\" = $N_KPOINTS ] || fail \"\$nk k points, step 1 had $N_KPOINTS\"" \
    "e=\$(grep '^!' $PREFIX.scf.out | tail -1 | awk '{print \$5}'); echo \"total energy \$e Ry, step 1 $E1 Ry\"" \
    "awk -v a=\$e -v b=$E1 -v t=$E_TOL 'BEGIN {exit !((a - b) ^ 2 <= t ^ 2)}' || fail \"total energy \$e, step 1 $E1\"" \
    "date; mpirun -np \$NSLOTS $QE/ph.x -input $PREFIX.ph.in > $PREFIX.ph.out" \
    "grep -q 'JOB DONE' $PREFIX.ph.out || fail \"ph.x did not finish\"" \
    "grep -q 'Effective charges (d Force / dE) in cartesian axis with asr applied' $PREFIX.ph.out || fail \"no effective charges\"" \
    "rm -rf tmp" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "dfpt of $NAME: $IN, SCF rerun from $MAT/dft/$STRUCTURE/scf/$PREFIX.scf.in (step 1 E = $E1 Ry, $N_KPOINTS k); $QUEUE $SLOTS slots. Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
