# Deriving the axion selection rule by hand (Bilbao / ISOTROPY)

Written 2026-10-25 for David's request: get the Wyckoff site of the magnetic atom, the magnetic groups before and after
the distortion, the axion-active irrep, and check that the spin cannot tilt, for YIO and MBT, without trusting the
script in [`../theory/axion_selection_rules.md`](../theory/axion_selection_rules.md).

Independent check from our own QE inputs (spglib, not Bilbao):
`archive/diagnostics/2026-10-25_magnetic_sites/magnetic_sites.py` (output `magnetic_sites.out` next to it).
Every number in the tables below comes from that output. Labels marked *(check)* are my derivation of a
Bilbao/ISOTROPY label and still need comparing with the server.

## 0. The whole argument in one line

θ is a time-odd pseudoscalar: under an operation g with rotation R it picks up χ_θ(g) = det R · (−1 if g is primed).
In the **paramagnetic parent** (grey group G1′), θ belongs to **mΓ1⁻** (A₁ᵤ, odd under T). The magnetic order L
belongs to some irrep Γ_L. A phonon Q changes θ at first order only if θ ∝ L·Q is allowed, i.e.

$$\Gamma_Q = \Gamma_L \otimes m\Gamma_1^-\qquad(\text{both 1D}).$$

Q must sit at the same k as L, because θ is uniform (k = 0). Everything else is bookkeeping and checks.

## 1. Recipe

| step | what | Bilbao tool | alternative |
|---|---|---|---|
| 1 | parent space group + Wyckoff letter of each atom | **WPASSIGN** (Structure Utilities) on the parent structure; or read off **WYCKPOS** by matching coordinates | FINDSYM (iso.byu.edu) |
| 2 | magnetic order: its k, its irrep Γ_L, the magnetic space group (BNS) | **MAGNDATA** entry (lists both); **MAXMAGN** or **k-SUBGROUPSMAG** (parent group + k → possible MSGs) | ISODISTORT, magnetic irreps |
| 3 | θ's irrep | always mΓ1⁻ for a centrosymmetric parent (A₁ᵤ, time-odd) | |
| 4 | Γ_Q = Γ_L ⊗ mΓ1⁻; displacement patterns of Γ_Q | **SMODES** (iso.byu.edu) for the patterns; Bilbao **REPRES** / **Get_irreps** for the character table at k | ISODISTORT |
| 5 | MSG of the distorted magnetic structure; no θ-odd op left | **MGENPOS** of the MSG: list each op, compute det R · (±1); none may give −1 | ISODISTORT with L and Q together gives the MSG directly |
| 6 | magnetic site symmetry of the magnetic atom after distortion; allowed moment | **MWYCKPOS** of the distorted MSG: the allowed moment components at that site | |

- **SMODES does not report the Wyckoff position; you type it in.** Each atom line of its input is
  `Mn a 0 0 0` (label, Wyckoff letter, coordinates), and the output only echoes it. That's why step 1 is a separate tool.
- **Wyckoff letters depend on the origin choice.** In Fd-3m, 16c and 16d are both .−3m sites, and moving the origin by
  (½,½,½) swaps them. Our YIO input puts Ir on 16d (spglib); the Materials Project table in
  [`../materials/Y2Ir2O7.md`](../materials/Y2Ir2O7.md) has Ir on 16c. Same physics. Compare site symmetries, not letters.
- **Irrep labels at k ≠ 0 depend on the origin too** (see §3, the half translation). Use one origin in every tool.

## 2. Y₂Ir₂O₇ (all-in-all-out)

| | value | from |
|---|---|---|
| parent space group | Fd-3m (227) | spglib |
| Ir Wyckoff site | 16c or 16d (origin), site symmetry .−3m (D₃d) | spglib: 16d in our cell |
| magnetic k | Γ (the 4 Ir of the primitive cell already carry the whole AIAO pattern) | |
| Γ_L (AIAO) | mΓ2⁺ (A₂g, time-odd) *(check: MAGNDATA 0.862, Eu₂Ir₂O₇)* | characters |
| magnetic space group | Fd-3m′ (BNS 227.131), point group m-3m′ (48 ops, 24 primed) | spglib, MAGNDATA |
| Ir magnetic site symmetry | −3m′: E, 2C₃, 3C₂′, I, 2S₆, 3m′ | script |
| θ-odd ops (θ quantized) | I, S₆, m (improper, unprimed); C₄′, C₂′ (proper, primed) | |
| Γ_Q | mΓ2⁺ ⊗ mΓ1⁻ = **Γ2⁻ (A₂ᵤ)** | §0 |
| distorted MSG | F-4′3m′ (BNS 216.76), point group −4′3m′ (24 ops: E, 8C₃, 3C₂, 6S₄′, 6m′) | spglib |
| θ-odd ops left | none (det R · s = +1 for all 24) | |
| Ir site after distortion | 16e (x,x,x), magnetic site symmetry **3m′** (E, 2C₃, 3m′) | spglib |
| allowed Ir moment | along the local ⟨111⟩ only → **cannot tilt**; only its length can change | script |

## 3. MnBi₂Te₄ (A-type AFM, moments ∥ c)

| | value | from |
|---|---|---|
| parent space group | R-3m (166), nonmagnetic hex cell a = 4.3336, c = 40.926 Å (3 septuple layers) | spglib, SMODES input |
| Wyckoff sites | Mn **3a** (0,0,0), site −3m (D₃d); Bi, Te1, Te2 all **6c** (0,0,z), site 3m | spglib, `data/MnBi2Te4/phonons/SMODES.in` |
| magnetic k | **T = (0,0,3/2)** hex = (½,½,½) rhombohedral (also called Z): neighbouring layers opposite | |
| Γ_L (A-type, m ∥ c) | **mT2⁺** with Mn at the origin *(check: MAGNDATA 1.309)* | characters, below |
| magnetic space group | R_I-3c (BNS 167.108), c doubled to 81.852 Å; point group −3m1′ (grey, because of {1′\|t½}) | spglib, MAGNDATA 1.309 |
| Mn magnetic site symmetry | **−3m′**: E, I, 2C₃, 2S₆, 3C₂′, 3m′ (Wyckoff letter in R_I-3c: read from MWYCKPOS 167.108) | script |
| θ-odd ops (θ quantized) | I (at Mn), S₆, the m that carry t½; Θ = {1′\|t½}, C₃·Θ, the C₂′ | |
| Γ_Q | mT2⁺ ⊗ mΓ1⁻ = **T2⁻** (4 copies; SMODES T2⁻) | §0 |
| distorted MSG | **R-3′m′ (BNS 166.100)**, same doubled cell; point group −3′m′ (12 ops: E, 2C₃, 3C₂, I′, 2S₆′, 3m′) | spglib |
| θ-odd ops left | none: I′, S₆′, m′ are improper *and* primed (+1); C₂ is proper and unprimed (+1) | |
| Mn site after distortion | 6c (0,0,z), magnetic site symmetry **3m′** (E, 2C₃, 3m′) | spglib |
| allowed Mn moment | along c only → **cannot tilt** | script |

### The half translation

This is the part David flagged. Call t the translation of the nonmagnetic lattice that takes the Mn1 layer to the Mn2
layer (t½ in the doubled cell). In the magnetic structure t alone is not a symmetry, because it flips the spins;
{1′|t} is (the anti-translation).

- **Why Q must be at T.** L is odd under t (k = T, phase e^{ik·t} = −1); θ is even under t. So Q must also be odd
  under t: k = T, the layers moving in **antiphase**. The in-phase partner (A₂ᵤ at Γ, the IR-active Γ2⁻ in SMODES) has the same
  pattern inside each layer but cannot couple.
- **Why the point group looks grey but θ is still allowed to move.** The magnetic point group −3m1′ contains 1′, which
  would seem to force θ = 0 or π forever. But 1′ only appears as {1′|t}. On Q at T, {1′|t} acts as t, giving −1, and
  χ_θ({1′|t}) = −1 too. So Q and θ transform alike; the point group loses that sign because it forgets t.
  In the magnetic cell, T2⁻ is a Γ mode that is odd under the anti-translation.
- **Which inversion survives.** In R-3m there are inversion centres on Mn (3a) and in the van der Waals gap (3b,
  (0,0,½) hex). In the magnetic group the gap one is I′: it equals t∘I_Mn, and t comes with 1′. T2⁻ breaks I at
  Mn (Mn1 moves +z, inversion sends it to −z) but keeps I′ at the gap (it sends Mn1 at +δ onto Mn2 at −δ, and the
  1′ fixes the spin). The script finds exactly one inversion in each distorted structure, primed, at z = 20.463 Å = the gap.
- **Labels flip with the origin.** I at the gap = t∘I at Mn, and χ(t) = −1 at T, so moving the origin from 3a to 3b
  flips every ± at T: T2⁻ ↔ T2⁺, and mT2⁺ ↔ mT2⁻. The product rule still holds, but SMODES (our input has Mn at the
  origin) and Bilbao must use the same origin, or the labels won't match.

### Hand check: L, Q, θ under the generators

Q = T2⁻ mode 1 (Mn1 +z, Mn2 −z); L = m_z(Mn1) = −m_z(Mn2). An axial z-vector is unchanged by I, flipped by an in-plane
C₂ and by a vertical mirror.

| operation (parent grey group) | L | Q | θ | L·Q = θ? |
|---|---|---|---|---|
| C₃ about c | + | + | + | ✓ |
| I at Mn1 | + | − | − | ✓ |
| C₂ in-plane through Mn1 | − | − | + | ✓ |
| vertical mirror through Mn1 | − | + | − | ✓ |
| t (Mn1 layer → Mn2 layer) | − | − | + | ✓ |
| T (time reversal) | − | + | − | ✓ |

These six generate the grey group; for YIO do the same with C₃ about [111], I at Ir, C₂ about [110], C₄ about z and T.
The script runs this check for all 48 (MBT) and 96 (YIO) grey-group operations and every mode: 0 failures.

### David's character-table method (his pyrochlore PDF, `emails/symm.pdf`) applied to MBT

For YIO he lists the classes of m-3m, marks each one antiunitary (A) and/or improper (I), and gets θ's characters
(+ for unmarked or for both marks, − for exactly one mark). The matching m-3m irrep is A₂ᵤ. For MBT the magnetic
point group is −3m1′, and 1′ alone is θ-odd, so the same steps on the point group would say no phonon ever couples.
Use the 24 operations of R_I-3c modulo the magnetic lattice instead: each −3m operation about Mn1 appears alone
({g|0}) or with t ({g|t}), and t swaps primed ↔ unprimed because it flips the spins.

| | E | 2C₃ | 3C₂ | I | 2S₆ | 3σ_d |
|---|---|---|---|---|---|---|
| {g\|0} marks | | | A | I | I | I A |
| θ on {g\|0} | + | + | − | − | − | + |
| {g\|t} marks | A | A | | I A | I A | I |
| θ on {g\|t} | − | − | + | + | + | − |

- The {g|0} row is A₂ᵤ of −3m (characters 1, 1, −1, −1, −1, 1), and the {g|t} row is −A₂ᵤ.
- So the mode must be A₂ᵤ-like in each layer **and odd under t**: k = T, i.e. T2⁻.
- The {g|0} row alone is the **ferromagnetic** group −3m′. That's why a skeleton built on −3m′ (R-3m′, 166.101) gives
  A₂ᵤ = Γ2⁻ and a distorted group R3m′ (160.67). spglib confirms both for FM order, but that is FM MBT, not A-type.
  Filled-in reply in David's format: `emails/yiro-mbt-symm_filled.txt`.

### Reading MGENPOS 167.108

- MGENPOS prints 72 lines = 24 operations × the 3 R-centring translations (0,0,0), (⅔,⅓,⅓), (⅓,⅔,⅔) of hexagonal
  axes. Read only the 24 in the (0,0,0) block.
- Coordinates are in the **doubled** hex cell (c = 81.852 Å), so (0,0,½) is the layer-to-layer translation t. The last entry
  of each triplet (±1) is time reversal.
- The 24, each checked as a symmetry of our structure with its moments (`2026-10-25_magnetic_sites/`, scratch check):

| block | operations (x,y,z form) | improper? | primed? | θ |
|---|---|---|---|---|
| 1-6 | x,y,z · −y,x−y,z · −x+y,−x,z · y,x,−z+½ · x−y,−y,−z+½ · −x,−x+y,−z+½ | no | no | + |
| 7-12 | −x,−y,−z · y,−x+y,−z · x−y,x,−z · −y,−x,z+½ · −x+y,y,z+½ · x,x−y,z+½ | yes | no | − |
| 13-18 | blocks 1-6 with z → z+½: x,y,z+½ · … · y,x,−z · x−y,−y,−z · −x,−x+y,−z | no | yes | − |
| 19-24 | blocks 7-12 with z → z+½: −x,−y,−z+½ · … · −y,−x,z · −x+y,y,z · x,x−y,z | yes | yes | + |

- The T2⁻ distortion keeps exactly the 12 operations with θ = + (blocks 1-6 and 19-24). That is R-3′m′ (166.100).
  This is general: a mode in θ's 1D irrep keeps exactly the θ-even operations. −x,−y,−z+½ (op 19) is the
  inversion centre in the van der Waals gap, at z = ¼ of the doubled cell = 20.46 Å.

## 4. Bilbao path for MBT (to follow, like David's pyrochlore PDF)

Bilbao's programs sit behind a Cloudflare check (2026-10-25), so they work from a browser but not from scripts.

1. **MAGNDATA** → search MnBi₂Te₄ → entry 1.309. Note: parent R-3m, k = (0,0,3/2), irrep, MSG R_I-3c (167.108). Its cell and
   z parameters are exactly our input.
2. **WPASSIGN** (or WYCKPOS 166 by hand) with the nonmagnetic cell from `SMODES.in` → Mn 3a (−3m), the rest 6c.
3. **MWYCKPOS 167.108** → the Mn site, its magnetic site symmetry (expect −3m′) and allowed moment (expect m ∥ c only).
4. **REPRES / Get_irreps** for 166 at T → character table; confirm A-type m_z order is mT2⁺ (characters in §3's table, L column)
   and that T2⁻ ⊗ mT2⁺ has the characters of mΓ1⁻.
5. **SMODES** (already done, `data/MnBi2Te4/phonons/SMODES.out`) → the four T2⁻ patterns.
6. **ISODISTORT** with mT2⁺ (magnetic) + T2⁻ (displacive) → resulting MSG: expect R-3′m′ (166.100) in the doubled
   cell. Or **k-SUBGROUPSMAG** of 166 at T, then pick the subgroup without θ-odd ops.
7. **MGENPOS 166.100** → list the 12 ops; check det R · (±1) = +1 for each.
8. **MWYCKPOS 166.100** → Mn on 6c, site 3m′, moment (0,0,m_z) only: no tilt.

## 5. What could go wrong (if Bilbao disagrees)

- The *(check)* labels: different tools name the k point T or Z, and the ± depends on the origin (§3).
- BNS vs OG numbering: spglib gives BNS 166.100 (OG 166.4.1330); Bilbao lists both.
- A Wyckoff letter that differs only by an origin shift (YIO 16c/16d).
