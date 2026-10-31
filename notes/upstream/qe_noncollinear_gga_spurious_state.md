# QE issue draft: noncollinear GGA with `lsign` finds a spurious lower-energy state in a collinear antiferromagnet

**Target:** gitlab.com/QEF/q-e, issue (discussion; not a code regression, no patch). Optional; the project worked
around it by using LDA for MnBi2Te4.

---

**Title:** Noncollinear PBE (+SOC) for a collinear AFM (MnBi2Te4): spurious layer-asymmetric state 55 mRy below the
symmetric one; SCF stalls

**Summary.** MnBi2Te4 with SOC must run noncollinear. The moments are collinear (±z), so `lsign = .true.` and
n↑,↓ = ½(n ± sign(m_z)|m|). From a random start (`startingwfc = 'atomic+random'`) the SCF drifts into a state with a
transverse magnetization texture (∫|m⊥| = 3.5 μB vs 0.1 in the symmetric state), |m|abs from 10.3 to 12.6-13.5 μB,
and inequivalent Mn layers. It lies **55 mRy below** the symmetric state, the gain entirely in the grid XC energy,
and never converges (stalls at 1e-5 to 1e-3 Ry).

**Mechanism (our reading).** With a transverse texture, m_z changes sign where |m| is finite, so the fixed-axis spin
density jumps by |m| (84,258 grid links in one run). PBE's gradient term rewards these jumps; the potential omits
∂s/∂m, so the SCF cannot settle. The symmetry that would forbid the state, the anti-translation {1′|½,½,½}, is not
enforced: QE treats the cell as a supercell and drops it.

**Controls.** Collinear PBE without SOC converges symmetric in 13 iterations. Noncollinear LDA (+SOC) from a random
start converges in 12 to the symmetric state. Noncollinear PBE from `startingwfc = 'atomic'` converges (the start is
exactly in the symmetric subspace), but the symmetric state is a saddle.

**Questions for the developers.** Is this a known limitation of the fixed-axis GGA? Would enforcing anti-translations
(magnetic supercell symmetries) be possible in `find_sym`?

**Attachments to prepare:** the base input, the outputs of the random-start, atomic-start and LDA runs, and
`archive/diagnostics/2026-10-08_mbt_noncollinear_gga/magnetization_texture.py`.
