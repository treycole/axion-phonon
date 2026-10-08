# Roadblocks


- The Hubbard U value is chosen to be 3.0 eV, which is a reasonable value for Ir 5d orbitals. We want it large enough so that the gap is sizable, and so that we are in the trivial insulator phase. 

- There are several local minima in the Wannierization procedure that is due to the DFT+U procedure. Different combinations of Hubbard occupations lead to different ground state energies. Sometimes, depending on the starting magnetization, the scf will converge to a different AIAO configuration with quenched moments. We need to make sure that the undistorted and distorted scf ground states are in the same AIAO configuration and in the same local minimum branch. We need small step sizes as well for this same reason.


# Fixes

