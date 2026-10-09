# Rigidly shifted Wannier functions

If we independently Wannierize the undistorted and distorted structures, we have different gauges for the Wannier functions and the finite difference calculation of the Wannier Hamiltonian will be unpredictable or unphysical. We expect that the combined internal and external contributions to the Berry curvature will lead to a gauge-independent result for the axion response. 


## Jae-Mo's code

We use Jae-Mo's Julia code that takes the MLWFs of the undistorted structure as trial orbitals for the distorted structure. These MLWFs are then rigidly shifted according to the difference in atomic positions between the undistorted and distorted structures. This ensures that the Wannier functions are in the same gauge and that the finite difference calculation of the Wannier Hamiltonian is predictable. 

The result is that we get a new `*.amn` file for the distorted structure. We then Wannierize, using only the projection gauge, so that we don't introduce a new gauge transformation.

How to run it on the cluster: [`../howto/rigid_shift_commands.md`](../howto/rigid_shift_commands.md).
