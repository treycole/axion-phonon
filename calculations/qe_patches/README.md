# QE 7.5 patch: noncollinear DFT+U occupation symmetriser

`new_ns_nc_preimage.patch` fixes `PW/src/new_ns.f90::new_ns_nc`. Two lines: import `invs`, and build the
symmetrised occupation of atom `na` from the atom that `isym` maps *into* `na`, `irt(invs(isym), na)`, not from
`irt(isym, na)`.

    cd /path/to/qe-7.5
    patch -p0 < .../new_ns_nc_preimage.patch
    make pw

The full explanation, the evidence and the verification are in
[`notes/log/2026-09-21_qe_new_ns_nc_bug.md`](../../notes/log/2026-09-21_qe_new_ns_nc_bug.md);
the check script that reproduces it from any `verbosity = 'high'` `pw.x` output is
[`diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py`](../diagnostics/2026-09-21_qe_new_ns_nc/qe_ns_symmetriser_check.py).

`verification/compare_hubbard_moments.py` prints the Ir Hubbard moments of two `pw.x` outputs. The unpatched and
patched runs it was checked on are in `data/Y2Ir2O7/qe_patches/verification/`, and the production-scale
pair is the SCF runs `178382bc889` (unpatched) and `192271bc889` (patched) in `data/`:

    python calculations/qe_patches/verification/compare_hubbard_moments.py \
        unpatched data/Y2Ir2O7/qe_patches/verification/scf_unpatched.out \
        patched   data/Y2Ir2O7/qe_patches/verification/scf_patched.out

Only `new_ns_nc` is changed. The collinear `new_ns` contracts the rotation on the other index and is correct
as it stands. Not checked: `new_nsg.f90` (DFT+U+V), which uses `irt` and `d2` in a similar pattern.
