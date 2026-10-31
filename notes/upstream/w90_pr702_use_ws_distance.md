# Wannier90 PR comment draft: PR 702's `_r.dat` ignores `use_ws_distance`

**Target:** github.com/wannier-developers/wannier90, comment on PR 702 (`transl_inv_full`). Jae-Mo Lihm's
`write_ndegen_applied` branch (jaemolihm/wannier90, plan5-write-ndegen-applied, 88134317) fixes it; coordinate with
him before posting.

---

**Summary.** With `transl_inv_full = .true.`, the `_r.dat` written by PR 702 equals postw90's A(R) to machine
precision, but it is frozen at `use_ws_distance = .false.`: the file is byte-identical whether the flag is true or
false, because the flag never reaches `plot_write_rmn`. `_hr.dat`/`_tb.dat` users who set `use_ws_distance = .true.`
then get H(R) and A(R) on different Wigner-Seitz representatives.

**Effect.** On a symmetry-forced zero of the Berry curvature (Y2Ir2O7 base, P·T), PR 702's A(R) is 15× better than
stock `_r.dat` but still 3×10⁵× worse than the `.mmn` route. Printing at higher precision changes this by ~1%; the
representative is the problem.

**Fix (Jae-Mo's branch).** Apply the WS distance and degeneracy in the writer (`write_ndegen_applied = .true.`). With
all three flags, its `_r.dat` matches postw90's `write_aa_r` at `use_ws_distance = .true.` and an independent export
to ~1e-8, and is Hermitian with every −R partner present.

**Open point to raise.** Under `transl_inv_full`, images of one aliasing class differ by (−1)ⁿ factors, so a tie rule
applied after the fact cannot reproduce Wannier90's; the tie tolerance and tie-breaking would need to be applied
inside Wannier90 (`plot.F90:198`).
