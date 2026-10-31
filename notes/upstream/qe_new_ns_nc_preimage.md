# QE MR draft: `new_ns_nc` symmetrizes noncollinear Hubbard occupations with the wrong atom

**Target:** gitlab.com/QEF/q-e, merge request into `develop`. **Patch:** `patches/qe_new_ns_nc_preimage.patch`
(commit 74bc68eb3, GitLab `treycole/q-e` branch `fix/new_ns_nc-preimage`, already pushed).
**Status:** `notes/bugs/07_qe_new_ns_nc_symmetrizer.md` records an MR as submitted in 2026-09. Check that it is
open on QEF; if not, send this. Upstream develop (c5a3cf95b, 2026-10-11) still has the bug.

---

**Title:** new_ns_nc: build the symmetrized occupation of atom na from the pre-image irt(invs(isym), na)

**Summary.** `PW/src/new_ns.f90::new_ns_nc` (noncollinear DFT+U) symmetrizes the occupation of atom `na` from
`nb = irt(isym, na)`, the **image** of `na`, while applying the forward rotation (`D nr Dᵀ` with
`conj(d_spin_ldau)`). Consistency requires the **pre-image**, `irt(invs(isym), na)`. The collinear `new_ns`
contracts `Dᵀ nr D` and is correct as it stands.

**When it matters.** Image and pre-image coincide for self-inverse operations (E, i, C2, mirrors). They differ for
operations of order 3, 4 or 6 that permute Hubbard atoms. These appear e.g. with `symmetry_with_labels = .true.`
and several Hubbard species of the same element.

**Example.** Y2Ir2O7 all-in-all-out (four Ir species, 48 operations, 28 of them permuting Ir with order > 2). QE's
analytic all-in-all-out ns is invariant under all 48 forward actions to 6e-15, yet one call of `new_ns_nc` as
written takes |m| from 3.0 to 1.0. In a real SCF the four Ir moments differ by ~2.5% and tilt up to 0.7° off ⟨111⟩;
patched: 0.00% and 0.00°, and the symmetrizer is an exact projector (1e-15).

**Not checked:** `new_nsg.f90` (DFT+U+V) uses `irt` and `d2` in a similar pattern.

**Attachments:** `qe_patches/verification/` (unpatched vs patched outputs and the moment comparison),
`archive/diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py`.
