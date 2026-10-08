#!/bin/bash
# Build pw.x with the |m| regularization (paw_mag_regularize.patch) on top of the segni fix, in a COPY
# of ~/q-e-segnifix, then run base SCF stall tests R1-R3 (delta = 1e-5, 1e-4, 1e-3 e/bohr^3) from the
# same stalled 207546 state and settings as T1/T1c. Does not modify ~/q-e-segnifix. Stops on the first error.
set -euo pipefail
SRC=$HOME/q-e-segnifix
DST=$HOME/q-e-segnireg
PATCH=$HOME/paw_mag_regularize.patch
B=$HOME/mnt/w/axion/Y2Ir2O7/phonon/base/u_3.0
S=$B/stall_test_2026-10-07
[ -e "$DST" ] && { echo "$DST exists"; exit 1; }
[ "$(md5sum < $SRC/bin/pw.x | cut -c1-8)" = e270b05b ] || { echo "unexpected segnifix pw.x"; exit 1; }
cp -a "$SRC" "$DST"; cd "$DST"
sed -i "s#$SRC\([/ ]\|\$\)#$DST\1#g" make.inc
grep -c "$SRC\b" make.inc || true
patch -p1 --dry-run < "$PATCH"; patch -p1 < "$PATCH"
git checkout -q -b fix/paw-mag-regularize
git -c user.name="Trey Cole" -c user.email="treco713@gmail.com" commit -q -am \
  "PAW noncollinear GGA (local frame): optional smooth |m| (QE_PAW_MAG_DELTA) to remove the kink at m = 0"
git log -1 --format='%h %s'
module load intel/2024 intel/ompi
nice make -j8 pw > build_reg.log 2>&1 || { tail -40 build_reg.log; exit 1; }
md5sum bin/pw.x
for t in "R1 1e-5" "R2 1e-4" "R3 1e-3"; do
  set -- $t; R=$S/${1}_delta$2
  mkdir $R
  cp -p $S/T1c_patched/{*.UPF,tmp.tar.gz,Y2Ir2O7.scf.in} $R/
  cp -p $DST/bin/pw.x $R/
  sed -e "s/^#\$ -N .*/#\$ -N YIOst_$1/" -e "s|^mpirun -np \$NSLOTS ./pw.x|export QE_PAW_MAG_DELTA=$2\nmpirun -x QE_PAW_MAG_DELTA -np \$NSLOTS ./pw.x|" $S/T1c_patched/send_job.sh > $R/send_job.sh
  grep -n -E "QE_PAW|pw.x -input" $R/send_job.sh
  echo "Stall test (2026-10-07): as T1c (segni-fixed pw.x, plain beta 0.03, CG, from the stalled 207546 state) but with the |m| regularization patch (~/q-e-segnireg), QE_PAW_MAG_DELTA=$2 e/bohr^3. T1 (unpatched) converged in 153 iterations; T1c (patched) stalled at 2e-5..1e-4 Ry for 338 iterations." > $R/NOTE
  (cd $R && qsub -l hostname='!n126&!n246&!n247' send_job.sh)
done
