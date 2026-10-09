# Rigid-shift `.amn` for a distorted structure: cluster commands

Why this is done: [`../theory/rigid_shift.md`](../theory/rigid_shift.md). A scripted version is [`../../scripts/generate_rigid_shift_amn.sh`](../../scripts/generate_rigid_shift_amn.sh).

To generate the `*.amn` file in the cluster, we run

```bash
GEN=~/mnt/w/QEWavefunctions.jl/scripts/generate_amn.jl
ls "$GEN" || echo "!! GEN path wrong"

BASE=/home/bc889/mnt/w/Y2Ir2O7/phonon/base/u_3.0/*bc889
BASE_SAVE=$BASE/tmp/Y2Ir2O7.save 

gunzip -k $BASE_SAVE/wfc*.dat.gz # if necessary, unzip the undistorted structure's wavefunctions
```

for each mode we do 

```bash
M=/home/bc889/mnt/w/Y2Ir2O7/phonon/mode_*/*bc889     # the mode's W90 dir (has .nnkp)
M_SAVE=$M/tmp/Y2Ir2O7.save
gunzip -k $M_SAVE/wfc*.dat.gz # if necessary, unzip the distorted structure's wavefunctions
mkdir -p $M/rigid_shift
echo "check: [$BASE_SAVE] [$M_SAVE]"                       # both must be non-empty
$GEN Y2Ir2O7  $BASE  $M  $BASE_SAVE  $M_SAVE  $M/rigid_shift/Y2Ir2O7.amn
```