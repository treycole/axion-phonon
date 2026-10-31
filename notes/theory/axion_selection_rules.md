# Axion selection rules for phonons and magnons

Computed, not just argued: `archive/diagnostics/2026-10-14_symmetry_selection_rules/selection_rules.py`.
- It finds each magnetic group by brute force from the DFT input structures and moments: YIO 210527, MBT 210628, and
  the mode SCF inputs.
- It then applies the character formulas below.
- It is pure group theory on the inputs; no wavefunctions are used. Runtime ~1 s.

## 1. Which operations quantize θ

θ multiplies E·B in the Lagrangian:
- **E** is a polar vector, time-even. **B** is an axial vector, time-odd.
- So E·B, and with it θ, is unchanged by proper rotations, flips under improper ones (inversion, mirrors, S_n), and
  flips under time reversal.
- For an operation g with rotation R, primed (with T) or not:

$$\chi_\theta(g) = \det R \cdot s(g),\qquad s = +1 \text{ unprimed},\; -1 \text{ primed}.$$

If the crystal has any g with χ_θ(g) = −1, then θ = −θ mod 2π, so **θ ∈ {0, π}**. Such operations are θ-odd.

| | θ-odd operations found (script) |
|---|---|
| YIO, AIAO, m-3m′ (48 ops, 24 primed) | i, S₆, σ_h (improper, unprimed); C₄′, C₂′ (proper, primed) |
| MBT, A-type AFM (24 ops, 12 primed) | i, S₆, σ_d·S (improper, unprimed); Θ = E·S′, C₃·S′, C₂′′ (proper, primed) |

- **The paper's "θ quantized solely by inversion" for YIO is incomplete.** C₄′ and C₂′ also quantize θ.
- In MBT, the primed proper ops are Θ = T·S_{1/2} and its products. The in-plane C₂ is primed on its own, because it
  flips m_z.

## 2. When does a phonon change θ at first order?

Expand θ(Q) = θ₀ + θ′Q + O(Q²) for a mode amplitude Q (time-even).
- Apply any g in the group: θ must transform by χ_θ(g), and Q goes to (gQ). Hence θ′·(gQ) = χ_θ(g) θ′Q for every g.
- **So θ′ ≠ 0 only if the mode transforms exactly like θ**: gQ = χ_θ(g)Q for all g. In other words, Q belongs to
  "θ's irrep", which is one-dimensional.
- A multi-dimensional irrep can never couple linearly.
- "The mode breaks every θ-quantizing symmetry" is **necessary, not sufficient**. If a mode breaks them all but sits
  in another irrep, θ appears only at a higher order.
- The order is the lowest k with θ's irrep inside Sym^k of the mode's irrep.
- The number of linearly active modes in a representation D is n_θ(D) = (1/|G|) Σ_g χ_θ(g) tr D(g).

**YIO** (22-atom cell, all Γ phonons):

| irrep | copies | θ linear | θ first at | IR | Raman |
|---|---|---|---|---|---|
| A₁g | 1 | no | — | — | Raman |
| **A₂ᵤ (Γ₂⁻)** | 3 | **yes** | Q | — | — |
| E_g | 1 | no | — | — | Raman |
| Eᵤ | 3 | no | Q³ | — | — |
| T₁g | 2 | no | — | — | — |
| T₁ᵤ | 8 | no | Q³ | IR | — |
| T₂g | 4 | no | — | — | Raman |
| T₂ᵤ | 4 | no | > Q⁴ | — | — |

**Only A₂ᵤ is linearly axion-active.** The paper's YIO table, which marks Eᵤ, T₁ᵤ and T₂ᵤ as axion-active, and the
sentence about T₁ᵤ being "both axion-active and IR-active", are wrong at linear order. Each displaced Γ₂⁻ structure
has 24 ops (12 primed), with no θ-odd operation left: −4′3m′.

**MBT** (14-atom AFM cell; labels are D₃d irreps of the point ops about Mn1, with (Γ) or (Z) for whether the
half-translation S acts as +1 or −1, i.e. in-phase or antiphase between the two septuple layers):

| irrep | copies | θ linear | θ first at | IR | Raman |
|---|---|---|---|---|---|
| A₁g(Γ) / A₁g(Z) | 3 / 3 | no | — | — / — | Raman / — |
| **A₂ᵤ(Z) = T₂⁻** | 4 | **yes** | Q | — | — |
| A₂ᵤ(Γ) | 4 (one acoustic) | no | — | IR | — |
| E_g(Γ) / E_g(Z) | 3 / 3 | no | — | — | Raman / — |
| Eᵤ(Γ) / Eᵤ(Z) | 4 / 4 | no / no | — / Q³ | IR / — | — |

- The four axion-active modes are exactly the SMODES T₂⁻ set. Each displaced structure has 12 ops (6 primed), with
  θ allowed.
- **A₂ᵤ(Z) and A₂ᵤ(Γ) have the same pattern within each septuple layer.** They differ only in whether neighbouring
  layers move in phase (Γ: polar, IR-active) or in antiphase (Z: axion-active, dark).
  - So each axion phonon has an IR-active twin, split only by the interlayer coupling.
  - For the rigid-layer mode, the twin is the acoustic translation.

## 3. Why axion-active phonons are optically dark: a general rule

- **IR activity:** the mode must transform like a polar vector component r_i (it carries a dipole).
- **Axion activity:** it must transform like θ.
- **Both:** θ·r_i must be invariant. θ·r_i is time-odd and axial, so it transforms like a magnetization component
  M_i.
- The character sums are literally identical: Σ χ_θ(g) tr R = Σ s det R tr R. Hence:

> **An axion-active phonon is IR-active only if the magnetic group allows a net magnetization along that direction.**
> With inversion (unprimed) in the group, θ is P-odd and Raman tensors are P-even, so it is never Raman-active.

| | net M allowed | θ in a vector (IR and axion) | θ in a Raman tensor |
|---|---|---|---|
| YIO AIAO | 0 | 0 | 0 |
| MBT AFM | 0 | 0 | 0 |
| control: MBT with FM order | 1 (M_z) | 1 | 0 |

- **The control confirms the rule.** With FM order, the axion-active set moves from A₂ᵤ(Z) to A₂ᵤ(Γ): the in-phase,
  c-polarized, IR-active modes. Only the 3 optical ones are physical; a uniform translation cannot change θ.
- **So a light-drivable axion phonon needs a centrosymmetric ferromagnetic axion insulator.**
  - Bulk FM MnBi₂Te₄ itself is, I believe, predicted to be a Weyl semimetal, with θ undefined (*to verify*).
  - MnBi₈Te₁₃ is reported as an FM axion insulator (*to verify*).
- **In AFM MBT, a magnetic field ∥ c induces M_z.** The dark A₂ᵤ(Z) mode then gains an IR dipole, so its absorption
  grows ∝ B², with E ∥ c, next to its IR-active A₂ᵤ(Γ) twin.
  - This happens through every mechanism that rides on M, not only the axion E·B term.

## 4. Why magnons cannot change θ in bulk MBT or YIO

- In both materials the magnetic ion sits on an inversion centre that is unprimed in the magnetic group. The script
  finds one such inversion fixing every magnetic site (mod the magnetic lattice) in each.
- Spins are axial: inversion leaves them unchanged.
- So **every q = 0 spin configuration keeps P, and θ stays quantized to all orders**. Uniform magnons, canting, spin
  flop and the paramagnetic-to-ordered transition all leave it at 0 or π.
- The script also gives n_θ(spin fluctuations) = 0 for both.
- In YIO this is the known statement that all q = 0 orders of the pyrochlore iridates preserve inversion (Wan et al.
  2011).
- **For MBT it means the dynamical axion of the bulk must be a lattice mode.**
  - The 2025 observation (Qiu et al., Nature 641, 62) used a six-layer film. An even-N film has PT instead of P, so
    magnons can couple there; the effect comes from the finite thickness.
  - Whether Qiu et al. state the bulk no-go: *to check* in their SI.

## 5. Caveats

- Linear order and q = 0 of the magnetic cell only. Two-phonon and phonon-magnon (Q·δm) terms are not analysed.
- "θ first at Q³" means allowed by symmetry, not that it is large.
- Labels come from characters on the input structures, so they inherit any symmetry the inputs lack. Both bases have
  their full groups (48 and 24 ops), and all mode structures have the expected subgroups.

## 6. Cross-checks: MAGNDATA and PythMPG (2026-10-18)

- **Magnetic space groups (MAGNDATA).** All-in-all-out pyrochlores are Fd-3m′ (BNS 227.131): entries 0.2 and 0.916
  (Cd₂Os₂O₇), 0.862 (Eu₂Ir₂O₇), 0.954 (Nd₂Ir₂O₇). There is no Y₂Ir₂O₇ entry. MnBi₂Te₄ is R_I-3c (BNS 167.108): entries 1.308
  (Yan 2019) and 1.309 (Ding 2020). The 1.309 cell (a = 4.33360, c = 81.85200 Å) and its z parameters (Bi 0.21245,
  Te 0.06668 / 0.14703) are exactly our MBT input.
- **Te labels.** Te1 (z = 0.06668, mode 3) is the outer Te, 5.46 Å from Mn and facing the 2.72 Å vdW gap. Te2 (z = 0.14703,
  mode 4) is the inner Te, 1.61 Å from Mn along c. The old `main.tex` had these reversed.
- **PythMPG** (conda env `pythmpg`, `mpg_tools.get_num_indep`). Counts of independent components:

  | Jahn | m-3m′ | -3m1′ | -4′3m′ | -3′m′ |
  |---|---|---|---|---|
  | `ae` (θ) | 0 | 0 | 1 | 1 |
  | `aeV` (M) | 0 | 0 | 0 | 0 |
  | `V` (polar vector) | 0 | 0 | 0 | 0 |
  | `aeV2` (linear ME tensor α) | 0 | 0 | 1 | 2 |

  - θ is quantized in both bases and allowed in both distorted groups.
  - No net M in any of the four, so the axion-active modes are IR-dark (§3).
  - In -4′3m′ the only allowed α is isotropic, so symmetry cannot separate the Chern-Simons part of the phonon-induced α
    from the rest.
  - PythMPG works with point groups only. The mode counting (§2, which irrep and how many copies) still needs the
    displacement representation, i.e. the script.
