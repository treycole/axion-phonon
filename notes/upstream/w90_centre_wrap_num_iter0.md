# Wannier90 issue draft: centres wrap for sites beyond 2/3 of a long cell (num_iter = 0), and the WS images follow

**Target:** github.com/wannier-developers/wannier90, issue. No patch; asks for a warning or a documented
recommendation.

---

**Title:** With num_iter = 0 on a coarse mesh, centres of functions far from the origin of a long cell come out
wrapped by a reciprocal-lattice period; use_ws_distance then picks wrong images silently

**Summary.** A centre is −Σ_b w_b b Im ln M_nn(k,b), and Im ln is defined modulo 2π. In a long cell on a coarse
mesh, the shortest b-vector along the long axis has a short wrap period, and functions far from the origin get an
inconsistent centre.
- **Example:** rhombohedral MnBi2Te4 (hexagonal c = 81.85 Å), 4×4×4 mesh, 116 spinor functions (Bi s,p; Te s,p;
  Mn d), projection only. The shortest axial b-vector is (b₁+b₂+b₃)/4, wrap period 4c/3 = 109.1 Å. Every function
  at crystal coordinate x > 2/3 along the axis is affected: Bi at z = 58.32 Å is reported at −50.4 Å
  (58.32 − 109.1 = −50.8).
- **Symptoms:** Ω_D = 1597 of 1978 Å² total spread; Bi spreads up to 56 Å²; exact bands on the mesh but an
  interpolated gap of −346 meV (DFT +198 meV). `use_ws_distance` chooses the images from these centres; with it off
  the gap is +83 meV.
- **Workaround:** place every projection site at |crystal x| ≤ ½ along the long axis (shift by a lattice vector in
  the .win). Total spread 382 Å², Ω_D 0.46 Å², gap +104 meV.

**Request.** A warning when Ω_D is large compared with Ω_OD, or when a centre lies far from its projection site;
and/or a note in the documentation that projection sites should be placed near the origin of the home cell in long
cells.

**Attachments to prepare:** the two .win files (original and recentred) and the .wout of both.
