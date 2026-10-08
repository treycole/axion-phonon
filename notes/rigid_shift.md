# Rigidly shifted Wannier functions

If we independently Wannierize the undistorted and distorted structures, we have different gauges for the Wannier functions and the finite difference calculation of the Wannier Hamiltonian will be unpredictable or unphysical. We expect that the combined internal and external contributions to the Berry curvature will lead to a gauge-independent result for the axion response. 


## Jae-Mo's code

We use Jae-Mo's Julia code that takes the MLWFs of the undistorted structure as trial orbitals for the distorted structure. These MLWFs are then rigidly shifted according to the difference in atomic positions between the undistorted and distorted structures. This ensures that the Wannier functions are in the same gauge and that the finite difference calculation of the Wannier Hamiltonian is predictable. 

The result is that we get a new `*.amn` file for the distorted structure. We then Wannierize, using only the projection gauge, so that we don't introduce a new gauge transformation.

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