# 22. PythTB's hopping cut and `_hr.dat` rounding break the spin-mixing symmetries (the YIO 0.5-0.9% floor)

**Status: FIXED in production 2026-10-12** (user's decision; found 2026-10-11). `load_wannier_pair(...,
hamiltonian_source="chk")` and `min_hopping_norm=None` are the defaults; `run_axion.py` settings `H_SOURCE`,
`MIN_HOPPING_NORM`. Evidence: `archive/diagnostics/2026-10-05_ws_tolerance_window/` `hopping_cut_test.py`,
`hr_rounding_eigenvalues.py` and their logs.

## Symptom

After the DFT fixes (#7, #16) and the WS tolerance (#13), the YIO base curvature still broke the 32 spin-mixing
operations of m-3m′ by 0.5-0.9% in the H-only (internal) part, the same in LDA and PBE, unchanged by the fold rule,
pinning or the tolerance. The 16 operations that keep z were exact (≤ 0.01%).

## Cause

Two numerical truncations of H(R), both made element by element in the spin-z spinor basis:
1. **`min_hopping_norm=1e-5`** in `W90.model(...)` (PythTB drops every element with |H_ij(R)| < 1e-5 eV), used by
   `archive/code/run_axion.py`, `modules/curvature.py`, `modules/axion.py`, `wannier/wannier_io.py` and the symmetry
   scripts.
2. **`_hr.dat` prints H to 1e-6 eV.**

A z-preserving operation maps each basis function to one basis function times a phase (±1, ±i), so it permutes the
|H_ij| exactly and both truncations commute with it. A spin-mixing operation maps a function to a combination of its
same-site spin partners, so an element just above the cut maps to combinations of elements just below it (and
rounding errors don't map onto rounding errors). The truncated H is then not covariant.

## Evidence (LDA base 212615, pinned fold, 7 generic k, %, internal / cross / external / total)

| H(R) | cut | #01 (keeps z) | #02 (C3-type, spin-mixing) | #05 (spin-mixing, primed) |
|---|---|---|---|---|
| `_hr.dat` | 1e-5 (production) | 0.008 / 0.035 / 0.007 / 0.008 | 0.828 / 0.822 / 0.142 / 0.844 | 0.825 / 0.836 / 0.145 / 0.841 |
| `.chk` | 1e-5 | 0.007 / 0.011 / 0.005 / 0.007 | 0.899 / 1.039 / 0.161 / 0.918 | — |
| `_hr.dat` | none | 0.007 / 0.029 / 0.006 / 0.007 | 0.101 / 0.540 / 0.025 / 0.083 | 0.097 / 0.535 / 0.028 / 0.079 |
| **`.chk`** | **none** | **0.000 / 0.005 / 0.005 / 0.000** | **0.000 / 0.003 / 0.005 / 0.000** | **0.000 / 0.004 / 0.006 / 0.000** |

- **Bands** (`hr_rounding_eigenvalues.py`): `_hr.dat` breaks the spin-mixing operations by 63 μeV (2.6 μeV for #01);
  the `.chk` H by 0.6 μeV. The folded H of the two routes differ by ≤ 7e-7 eV per element; that is enough.
- **Independent check:** a Kubo curvature from the folded H alone (no PythTB, no cut) gives 0.318% on #02 from
  `_hr.dat` and 0.000% from the `.chk`.
- **The `.chk` H is the same Hamiltonian:** its class sums agree with `_hr.dat` to 3.4e-6 eV (the print rounding).

**Full group, production fold** (`base_full_group_hchk_nocut.py`: Wannier centres, tolerance 5e-3, H from `.chk`, no cut,
47 non-identity operations, median (max), %; rows also in `symmetry_comparison_table.md`):

| base | internal | cross | external | total |
|---|---|---|---|---|
| LDA 212615 | 0.008 (0.008) | 0.220 (0.224) | 0.012 (0.020) | **0.000 (0.000)** |
| PBE 212616 | 0.012 (0.019) | 0.401 (0.441) | 0.252 (0.306) | **0.013 (0.016)** |

- **LDA is exact** in the total. The nonzero parts sit on the 32 spin-mixing operations only: the split into
  internal/cross/external depends on the embedding τ, the total does not, and τ = Wannier centres is not exactly
  covariant (spin-partner centres differ by ≤ 8.6e-4 Å). With pinned τ each part is exact too (≤ 0.007%).
- **PBE keeps 0.013-0.016%,** again only on the spin-mixing operations (the 16 z-preserving ones are 0.00). Likely the
  PBE DFT's own residual (its NSCF eigenvalues are ~100× less symmetric than LDA's); not traced.

## Fix (applied 2026-10-12)

Production H(R) from the `.chk` + `.eig` (V†εV, as A(R) already is with `transl_inv_full_from_chk`), and no hopping
cut (`min_hopping_norm=None`). The `.eig` must be kept next to every `.chk` (the cleanup script keeps it).
- **Code:** `wannier/wannier_io.py` `hamiltonian_classes_from_checkpoint` gives one class sum per aliasing class.
  `load_wannier_pair(hamiltonian_source="chk")`, the default, puts it on one representative before the
  orbital-dependent fold.
- **Consistency check:** it raises if a class sum differs from `_hr.dat` by more than 1e-4 eV, i.e. the files are
  from different runs. On base 212615 the difference is 3.4e-6 eV, which is print rounding.
- **`run_axion.py`:**
  - `H_SOURCE = "chk"` and `MIN_HOPPING_NORM = None`.
  - The output folder gets `_hrdat`/`_cut` suffixes for the old settings.
  - Both settings are stored in the results metadata.

## Remaining

- Check dθ before and after on the YIO LDA mode-1 pair (`H_SOURCE = "hr_dat"`, `MIN_HOPPING_NORM = 1e-5` for the
  old value).
- MBT: same change, together with the pinned option (#21).
