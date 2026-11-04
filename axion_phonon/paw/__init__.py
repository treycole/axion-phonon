"""Step 4: the displaced structure's .amn from overlaps with the PAW metric (bug #23); theory and sources in
``notes/theory/paw_overlap.md``.

- ``qe_port.py``: S of one structure, ported from QE (projectors, q_ij, spin-orbit, ``s_psi_nc``). Nothing new.
- ``cross_structure.py``: the overlap of two structures, O = <psi_Q|psi_0>, and A = O V_0. Its plane-wave part
  and A = O V_0 are what QEWavefunctions.jl's ``generate_amn.jl --no-displacement`` does; the projectors on each
  side's own atoms and the d . D_ij partial-wave term (exact to first order in d) are new.
"""
from axion_phonon.paw.cross_structure import mode_amn, overlap, pw_overlap
from axion_phonon.paw.qe_port import PAWSystem

__all__ = ["PAWSystem", "overlap", "pw_overlap", "mode_amn"]
