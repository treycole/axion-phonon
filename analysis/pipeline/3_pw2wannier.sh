#!/bin/bash
# Step 3: overlaps between neighbouring k (pw2wannier90.x): the .mmn and .eig, one per structure and mesh; for the
# base also the .amn onto its atomic trial orbitals. Run on the login node; wait for dft/<structure>/pw2wannier_k<NK>/OK.
#   input   the base's dft/<structure>/pw2wannier_k<NK>/<PREFIX>.win (written by you, once per material and mesh: projections,
#           windows, num_iter = 0, dis_num_iter = 0). For a displaced structure this script writes its .win with each
#           atom moved by its displacement. Here: wannier90.x -pp on the login node (the .nnkp).
#   job     copies the NSCF save into scratch, runs pw2wannier90.x -nk POOLS -pd true, gzips .mmn (and .amn).
#   output  dft/<structure>/pw2wannier_k<NK>/out/: <PREFIX>.mmn.gz, .eig (, .amn.gz); <PREFIX>.win and .nnkp in the step folder
set -euo pipefail
# ============================================================================== settings
MAT=/mnt/wkdk/bc889/axion/MnBi2Te4
STRUCTURE=T2-_p3_Q0.01   # a folder of $MAT/dft/: base, T2-_p1_Q0.01, ...
PREFIX=MnBi2Te4
NK=8
POOLS=8
QUEUE=dko48m4,dko64m4
SLOTS=48
H_RT=96:00:00
SCRATCH_GB=200           # free scratch the job needs: the save (MBT 87 GB at 8^3) and the .mmn
SUBMIT=true
# ======================================================================================
source $(dirname $0)/common.sh
check_tools
D=$MAT/dft/$STRUCTURE/pw2wannier_k$NK
TEMPLATE=$MAT/dft/base/pw2wannier_k$NK/$PREFIX.win
NAME=$STRUCTURE
ready $D
NSCF=$(need_ok $MAT/dft/$STRUCTURE/nscf_k$NK)
NKPT=$((NK * NK * NK))
[ -s $TEMPLATE ] || die "no $TEMPLATE"
$PY $CODE/pipeline/inputs.py win $TEMPLATE $MAT/dft/base/scf/$PREFIX.scf.in $MAT/dft/$STRUCTURE/scf/$PREFIX.scf.in > $D/$PREFIX.win.new
if [ $STRUCTURE = base ]; then
  diff <(sed -n '/begin atoms_cart/,/end atoms_cart/p' $TEMPLATE | awk 'NF == 4 {printf "%s %.10f %.10f %.10f\n", $1, $2, $3, $4}') \
       <(sed -n '/begin atoms_cart/,/end atoms_cart/p' $D/$PREFIX.win.new | awk 'NF == 4 {printf "%s %.10f %.10f %.10f\n", $1, $2, $3, $4}') \
    > /dev/null || die "inputs.py win changed the base's own atoms"
  rm $D/$PREFIX.win.new
  WRITE_AMN=.true.
else
  mv $D/$PREFIX.win.new $D/$PREFIX.win
  WRITE_AMN=.false.                                  # the displaced structure's .amn comes from steps 4 and 5
fi
WIN=$D/$PREFIX.win
grep -Eq "^\s*num_iter\s*=\s*0\s*$" $WIN && grep -Eq "^\s*dis_num_iter\s*=\s*0\s*$" $WIN || die "$WIN: not num_iter = 0, dis_num_iter = 0"
grep -Eiq "^\s*mp_grid\s*=?\s*$NK\s+$NK\s+$NK\s*$" $WIN || die "$WIN: mp_grid is not $NK $NK $NK"
printf "%s\n" "&inputpp" "    outdir = './tmp/'" "    prefix = '$PREFIX'" "    seedname = '$PREFIX'" "    write_amn = $WRITE_AMN" \
  "    write_mmn = .true." "    write_unk = .false." "    write_eig = .true." "/" > $D/${PREFIX}_pw2wan.in
(cd $D && $W90 -pp $PREFIX > /dev/null 2>&1 && [ -s $PREFIX.nnkp ]) || die "wannier90.x -pp failed in $D"
for blk in kpoints nnkpts; do                        # projection centres move with the atoms; these must not
  [ $STRUCTURE = base ] || cmp -s <(sed -n "/begin $blk/,/end $blk/p" $D/$PREFIX.nnkp) \
                                   <(sed -n "/begin $blk/,/end $blk/p" $MAT/dft/base/pw2wannier_k$NK/$PREFIX.nnkp) \
    || die "the .nnkp $blk differ from the base's"
done

{ header pw2w${NK}_$NAME $SLOTS 3.6G $H_RT $QUEUE; printf "%s\n" \
    "for f in $PREFIX.win $PREFIX.nnkp ${PREFIX}_pw2wan.in; do cp -p \$SGE_O_WORKDIR/\$f . && [ -s \$f ] || fail \"input \$f\"; done" \
    "scratch $SCRATCH_GB" "md5sum $QE/pw2wannier90.x" \
    "date; get_save $NSCF/tmp/$PREFIX.save tmp $NKPT" \
    "date; mpirun -np \$NSLOTS $QE/pw2wannier90.x -nk $POOLS -pd true -input ${PREFIX}_pw2wan.in > ${PREFIX}_pw2wan.out" \
    "rm -rf tmp; grep -q 'JOB DONE' ${PREFIX}_pw2wan.out || fail \"pw2wannier90.x did not finish\"" \
    "[ -s $PREFIX.mmn ] && [ -s $PREFIX.eig ] || fail \"no .mmn or .eig\"" \
    "[ $WRITE_AMN = .false. ] || [ -s $PREFIX.amn ] || fail \"no .amn\"" \
    "date; pigz -p \$NSLOTS $PREFIX.mmn; [ ! -e $PREFIX.amn ] || pigz -p \$NSLOTS $PREFIX.amn" \
    "ok"; } > $D/send_job.sh
submit $D send_job.sh
note $D "pw2wannier90 k$NK of $NAME: save $NSCF, write_amn $WRITE_AMN; .win from $TEMPLATE. Job $(cat $D/SUBMITTED 2>/dev/null)."
cp -p $0 $D/
