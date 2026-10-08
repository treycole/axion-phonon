# YIO Wannierization

## Atoms in unit cell

There are 

- 4 Y atoms
- 4 Ir atoms
- 14 O atoms

## Character of bands w/o SOC

- Y s: form deeper valence states around -30 eV
- Y p: flat manifold around -10 eV
- Y d: predominantly conduction states above E_F
- O p: predominantly valence states around [0, 15] eV
- O s: predominantly valence states around [-12, -5] eV
- Ir s, p, d: span valence and conduction from [-12, 15] eV

## Trial wavefunctions

We want all dispersive occupied bands and some conduction bands to be captured by the trial wavefunctions. The following are the trial wavefunctions we will use:

- O p: large DOS in primary valence bands from [0, E_F] eV
- O s: large DOS in semi-core valence states from [-12, 0] eV
- Ir p: span valence (semi-core and primary) and conduction states. Important for effective spin-1/2 around Fermi level.
- Ir s: similar DOS spread as Ir p.
- Ir d: largest DOS that spans primary valence and conduction bands.
- Y d: predominantly conduction states above E_F. 

Total number of trial wavefunctions when spin-orbit coupling is included: 

- O p = 14 atoms * 3 orbitals * 2 spin = 84
- O s = 14 atoms * 1 orbital * 2 spin = 28
- Ir p = 4 atoms * 3 orbitals * 2 spin = 24
- Ir s = 4 atoms * 1 orbital * 2 spin = 8
- Ir d = 4 atoms * 5 orbitals * 2 spin = 40
- Y d = 4 atoms * 5 orbitals * 2 spin = 40

Total = 224 trial wavefunctions. 

This is a large number of trial wavefunctions, but we want to capture all the relevant physics in the valence and conduction bands.

## Disentanglement and frozen windows

We want enough conduction bands (`nbnds`) in the nscf calculation so that the frozen window can be large enough to include most of the conduction bands, minus a few Wannier bands, and the remaining disentanglement window is large enough so that the iterative subspace selection can find localized Wannier functions. 

We have been using `nbnds = 300` but the remaining disentanglement window is somewhat small. We are now increasing to `nbnds = 350` so that the disentanglement window is larger.

### Frozen window

The frozen window is chosen to include all occupied bands and most of the conduction bands so that most of the Wannier bands exactly match the DFT bands.

Previously, we chose the inner window to be from [0, 23] eV, which makes most of the Wannier bands exactly match the DFT bands. The disentanglement window was chosen to be [0, 28] eV. This is small enough so that the unconverged upper most conduction bands are excluded. However, now that we are including Os trial orbitals, the windows need to be adjusted to incorporate the semi-core valence states. 

Looking at the the bands of the undistorted structure (with SOC), the semi-core valence states are within [-14, 0] eV and there are 52 of them. If we extend the bottom of the inner window to there, and leave the upper end of the inner window as is, then there are at most 244 bands across the BZ which is too many. Changing the upper end of frozen window to be 19 eV means there are at most 216 states between [-14, 19] eV across the BZ, which is less than the number of trial wavefunctions (224). This is a reasonable choice for the frozen window.

**Frozen window** = [-14, 19] eV

### Disentanglement window

The remaining states above the frozen window are included in the disentanglement window. The disentanglement window is chosen to be so that the unconverged upper most conduction bands are excluded. The upper end of the disentanglement window was previously chosen to be 28 eV. Now that we are including more bands, we can increase the upper end of the disentanglement window to 30 eV for a first try. This is small enough so that the unconverged upper most conduction bands are excluded.

**Disentanglement window** = [-14, 30] eV


