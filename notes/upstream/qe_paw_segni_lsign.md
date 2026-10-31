# QE MR draft: PAW noncollinear GGA uses a stale `ux` when `lsign` is false

**Target:** gitlab.com/QEF/q-e, merge request into `develop` + issue. **Patch:** `patches/qe_paw_segni_lsign.patch`
(branch `fix/paw-segni-lsign-develop`, one commit on develop c5a3cf95b).

---

**Title:** PAW one-centre GGA (noncollinear) applies the fixed-axis spin sign with a leftover axis when
`lsign = .false.`, breaking the symmetry of noncollinear magnets

**Summary.** For noncollinear GGA, QE builds n↑,↓ = ½(n ± s|m|). `compute_ux` sets `lsign = .true.` and a common
axis `ux` only if all starting moments are (anti)parallel; then s = sign(m·ux). Otherwise `lsign = .false.` and the
local frame (s = 1) should be used.
- `PW/src/compute_ux.f90` sets `ux` to the first magnetic atom's starting moment and **never resets it** when it
  later finds a non-parallel moment.
- `PW/src/compute_rho.f90` (smooth grid) uses `segni = SIGN(1, m·ux)` only `IF (lsign)` — correct.
- `PW/src/paw_onecenter.f90::compute_rho_spin_lm` (line ~2340) and `compute_drho_spin_lm` (~2554) use
  `segni_rad = SIGN(1, m·ux)` **unconditionally** (present since at least QE 6.5).

So for a noncollinear order, the PAW spheres use the first magnetic atom's axis. Inside each sphere m rotates, so
s flips where m ⊥ ux at finite |m|, and n↑, n↓ jump there. The GGA turns the jumps into a spurious one-centre
potential D_ij. Because `ux` is one fixed direction, symmetry-equivalent atoms get inequivalent spurious terms.

**Why it is easy to miss.** The SCF symmetrizes ρ, m, becsum and `ddd_paw` (`PAW_symmetrize_ddd`) every iteration,
so every printed and stored quantity is exactly symmetric. A non-self-consistent run builds `ddd_paw` once from the
stored becsum and never symmetrizes it, so its eigenvalues show the breaking.

**Example.** Y2Ir2O7, all-in-all-out order (magnetic group m-3m′, 48 operations), PBE+U, PAW, SOC. NSCF on the full
8³ mesh with `nosym`: eigenvalues of symmetry-related k points differ by up to **4.8 meV**; only the 12 operations
of −3m′ about the first Ir atom's [111] axis survive (to ~1 meV). With the patch the worst difference is **46 μeV**
(all 47 operations). LDA shows no breaking (no gradients); moving `ux` to z with a tiny seed moment on another atom
moves the pattern accordingly.

**Fix.** In both routines, `segni_rad = SIGN(1, m·ux)` only if `lsign`, and 1 otherwise, the rule of
`compute_rho`. (`use_ux = lsign` is passed into the OpenACC loop as `firstprivate`.) Clearing `ux` in `compute_ux`
alone is not enough: with signed zeros (gfortran default) `SIGN(1, m·0)` is −1 in part of the sphere.

**Who is affected.** Noncollinear magnetic + GGA + PAW runs whose starting moments are not all parallel or
antiparallel (all-in-all-out pyrochlores, canted and spiral magnets). Collinear orders run noncollinearly
(`lsign = .true.`) are unchanged by the patch.

**Attachments to prepare:** the Y2Ir2O7 input (22 atoms), unpatched vs patched NSCF eigenvalue-symmetry tables
(`archive/diagnostics/2026-10-05_ws_tolerance_window/eig_symmetry_xml_full.py`), the step-by-step code
walk-through in `notes/bugs/16_qe_paw_segni_stale_ux.md`.
