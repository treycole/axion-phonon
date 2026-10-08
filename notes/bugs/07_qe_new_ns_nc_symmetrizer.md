# #7 — QE `new_ns_nc` symmetrizes noncollinear Hubbard occupations with the wrong atom

> **Status: FIXED** (2026-09-21). A two-line patch (`calculations/qe_patches/`) is in every pw.x used since
> 2026-09-22: the `q-e-nsfix` build, `~/q-e` (`c41c409`) and `~/q-e-segnifix` (`e270b05b`). An upstream MR was
> submitted from GitLab `fix/new_ns_nc-preimage`.

**Found:** 2026-09-21. **Affects:** QE 7.5 and 7.6 noncollinear DFT+U with symmetry, when symmetry operations
permute Hubbard atoms. In this project: every Y2Ir2O7 DFT run through 2026-09-21, and MnBi2Te4 T2⁻ with
`symmetry_with_labels` (#9).

## Symptom
- The four YIO Ir moments differed by ~9% and tilted off their ⟨111⟩ axes.
- A Hubbard spin channel came out too weak (the minimum gap was 0.72 eV, 1.04 eV after the fix).
- There was a Berry-curvature symmetry floor of 4-14%. Part of that floor was this bug; the rest was #16.

## Root cause
`PW/src/new_ns.f90::new_ns_nc` builds the symmetrized ns of atom `na` from `nb = irt(isym, na)`, the **image**
of `na`, while applying the forward rotation (`D nr Dᵀ`, `conj(d_spin_ldau)`). Pairing the forward rotation with
the image is inconsistent: it must use the **pre-image**, `irt(invs(isym), na)`. The collinear `new_ns`
contracts `Dᵀ nr D` and is correct.

## Why it was a bug
- **Harmless cases:** for self-inverse operations (E, i, C2, mirrors), image = pre-image.
- **Wrong cases:** every operation of order 3, 4 or 6 that **permutes** Hubbard atoms, which is 28 of 48 in
  YIO (8 C3, 8 S6, 6 C4T, 6 S4T).
- **What lets them through:** `symmetry_with_labels = .true.` with four Ir species (atoms matched by chemical
  symbol) is what lets such operations cross the Ir species.
- **Effect:** the "symmetrized" ns isn't symmetric. QE's own analytic AIAO ns is invariant under all 48 forward
  actions (6e-15), but `new_ns_nc` as written doesn't return it (|m| 3.0 → 1.0).

## Fix
Use the pre-image, `irt(invs(isym), na)`: two lines (`calculations/qe_patches/`). With it, the symmetrizer is an
exact projector (1e-15).

## Verification
- **Ported routines:** QE's own `irt`/`sr`/`t_rev`/`find_u`/`d_matrix` ported to Python
  (`calculations/diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py`).
- **Real pw.x build:** identical raw data, patched vs unpatched. Unpatched: Ir |m| spread 2.5%, tilts up to 0.7°.
  Patched: 0.00% / 0.00°.
- **Saved ns exactly covariant:** on 2026-10-05 the saved ns of the patched 192344 and 193145 SCFs was exactly
  covariant under all 48 operations, phases included (3e-10).

## Remaining
- **The floor remained after the fix** (H(R) residuals fell 2-4×). Its cause turned out to be #16.
- **Workarounds at the time:** `nosym = .true.`, or dropping `symmetry_with_labels`. The user prefers the patch,
  keeping symmetry on.

**Sources:** [`../progress/2026-09-21_qe_new_ns_nc_bug.md`](../progress/2026-09-21_qe_new_ns_nc_bug.md),
[`../progress/2026-09-29_yio_qe76_new_ns_nc_pair.md`](../progress/2026-09-29_yio_qe76_new_ns_nc_pair.md);
memory `qe-new-ns-nc-bug`; CLAUDE.md pitfall 7.
