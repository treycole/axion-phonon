# Results at a glance

The settings, trial sets and responses of every axion-active pattern, and every convergence test with its status. Each
number here is in the paper draft (`paper/draft/main_draft.tex`) with a `%` comment naming its data; where a number
comes from and how to remake it is in [`REPRODUCE.md`](REPRODUCE.md); the runs are listed in
`data/<material>/README.md`. **Update this file when a run finishes** (one row per result; status column).

Units: ∂θ/∂Q in rad/Å, Q the displacement per Cartesian component of the SMODES basis vector **v**ᵢ. "Undistorted end"
and "displaced end" are the two estimates of ∂θ/∂Q from one base/displaced pair: the momentum-space quantities (H(k)
eigenstates, k-curvatures) of the undistorted or of the displaced model, the finite difference in Q the same. The
undistorted end is the estimate of the derivative at Q = 0 and the one the paper uses; the displaced end converges to
the same limit as Q → 0 (difference ∝ Q²) and serves only as a check of the pair construction.

## Settings

| | MnBi₂Te₄ | Y₂Ir₂O₇ |
|---|---|---|
| DFT | LDA+U, U = 4 eV on Mn 3d, SOC, A-type AFM, exact hex-CIF cell (a = 4.33360, c = 81.85200 Å), 14 atoms | LDA+U, U = 3 eV on Ir 5d, SOC, all-in-all-out, 22 atoms |
| PAW (`rel-pz`) | Mn `spn` 0.3.1, Bi `dn` 1.0.0, Te `n` 1.0.0 (generated, `MnBi2Te4/pseudo_gen/`) | Y `spn` 1.0.0, Ir `n` 0.2.3, O `n` 0.1 |
| code | pw.x md5 `db4b99b4`, conv_thr 1e-10, ecutwfc 60 Ry | same, ecutwfc 50 Ry |
| base runs | SCF + NSCF 214290 (8³, 240 bands), pw2wannier90 214300 | SCF 213998, NSCF 214124 (8³, 300 bands), pw2wannier90 220421 |
| VBM (8³ NSCF) | 7.268 eV | 10.4437 eV |
| trial set (production) | Bi s,p; Te s,p; Mn d (spinor): **116 WFs, 82 occupied** | Y p; Ir d; O s,p (spinor): **176 WFs, 156 occupied** |
| window (production) | `f10.962_w12.810_k8`: frozen −7 … 10.962 eV (VBM + 3.69, the highest frozen top 116 WFs allow), outer −7 … 12.810 eV | `f12.422_w22.644_k8`: frozen −14 … 12.422 eV (VBM + 1.98), outer −14 … 22.644 eV |
| Wannier90 | 4.0.3, projection only (`num_iter = 0`, `dis_num_iter = 0`); base 219737 | same; base 220455 (Ω_I 117.10 Å²) |
| displaced structures | trial functions = the base's WFs through the PAW overlaps with first-order augmentation (pipeline step 4) | same |
| model | H(R) from `.chk` + `.eig`, no hopping cut; A(R) translationally equivariant from the `.chk`; shared WS images; τ pinned to the base's atoms | same |
| c₂ k list | tube (16, 4, 3, 8): 7,928 k around the k_z axis, m = 4, 8 planes, half-width 3 | uniform 10³ |

## Response per pattern

**Y₂Ir₂O₇, Γ₂⁻** (11-02 pipeline run; window `f12.422_w22.644_k8`; Q → 0 from a + bQ² + cQ⁴ over Q = 0.01, 0.02,
0.03, 0.05 Å):

| pattern | atoms (\|v_κ\|) | undistorted end, Q = 0.01 / 0.02 / 0.03 / 0.05 | displaced end, same Q | Chern–Simons 3-form [θ(Q) − θ(0)]/Q, 12³ | **Q → 0** |
|---|---|---|---|---|---|
| 1 | 4 Y (√3) | −0.017114 / −0.017093 / −0.017058 / −0.016906 | −0.01721 / −0.01748 / −0.01793 / −0.01939 | −0.017127 / −0.017143 / −0.017170 / −0.017212 | **−0.017120** |
| 2 | 4 Ir (√3) | 0.055021 / 0.054984 / 0.054915 / 0.054782 | 0.05460 / 0.05338 / 0.05163 / 0.04853 | 0.054924 / 0.054617 / 0.054164 / 0.053346 | **+0.055037** |
| 3 | 12 O(48f) (1) | 0.063867 / 0.063757 / 0.063680 / 0.063420 | 0.06383 / 0.06361 / 0.06334 / 0.06245 | 0.063844 / 0.063667 / 0.063475 / 0.062846 | **+0.063880** |

Every estimator (either end, their mean, the 3-form) gives the same Q → 0 limit to ≤ 2×10⁻⁶.

**MnBi₂Te₄, T₂⁻** (10-30 production; window `f10.962_w12.810_k8`; tube list):

| pattern | atoms | undistorted end, Q = 0.01 | displaced end, Q = 0.01 | Q = 0.03 (undistorted / displaced) | **used** |
|---|---|---|---|---|---|
| 1 | 2 Mn | −1.95890 | −1.95500 | – | **−1.9589** |
| 2 | 4 Bi | 3.81959 | 3.87149 | – | **+3.8196** |
| 3 | 4 Te1 (vdW gap) | 6.38695 | 6.53961 | 6.36436 / 7.82703 | **+6.3870** (Q → 0, a + bQ²: 6.3898 / 6.3787) |
| 4 | 4 Te2 (bonded to Mn) | −4.54479 | −4.53670 | – | **−4.5448** |

**Axion charges** z_κ = (∂θ/∂Qᵢ)/(Nᵢ\|v_κ\|) (rad/Å): MnBi₂Te₄ Mn −0.98, Bi +0.95, Te1 +1.60, Te2 −1.14;
Y₂Ir₂O₇ Y −0.0025, Ir +0.0079, O(48f) +0.0053.

**Normal modes** (frozen-phonon force constants; ∂θ per Å of total displacement): MnBi₂Te₄ 17.8 / 82.7 / 128.7 /
143.3 cm⁻¹: −1.45 / +0.10 / +3.63 / −0.62; Y₂Ir₂O₇ 228.6 / 328.3 / 471.6 cm⁻¹: −0.024 / +0.010 / −0.017.

## Convergence and checks

| test | MnBi₂Te₄ | Y₂Ir₂O₇ | status |
|---|---|---|---|
| k integration | dense tube (16, 16, 5, 32), 995,552 k: production within +0.003 / +0.13 / +0.04 / +0.17% (patterns 1–4) | uniform 4³ / 6³ / 8³ / 10³ / 12³: 8³–12³ agree to 4×10⁻⁶ | done, current |
| uniform mesh instead of the tube | 4³ / 8³ / 16³: 54.7 / 13.6 / 3.20 (pattern 3), against 6.39 | (uniform is production) | done |
| amplitude Q | pattern 3 undistorted end: 6.387 (0.01) → 6.364 (0.03), −0.4%; displaced end +2.4% (0.01), +23% (0.03) of the undistorted, ∝ Q² | undistorted end within 1.2% of its limit at 0.05 Å (pattern 1; ≤ 0.15% at 0.02 Å); displaced end up to 13% at 0.05 Å; 3-form chord ≤ 3.1% | done, current |
| energy window | pattern 3, 15 windows on 8³ (10-28/30, current settings; table below): frozen top at VBM + 3.69 eV, 5 outer tops: 6.352–6.398 (0.7%); lower frozen tops scatter 5.78–7.56 | 6 windows (10-19, earlier trial functions without the PAW metric and the Q = 0 self-test base; table below): ≤ 0.75% on the plateau, frozen top 14.396 eV fails (Ω_I +35%, pattern 1 +5.7%) | MBT current; YIO not rerun (the PAW-metric correction moved the production values by ≤ 0.6%) |
| Wannier mesh | 8³ → 10³ at windows `f10.698_w13.915`: 6.29562 → 6.23783 (−0.9%), `f8.395_w13.915`: 7.07453 → 6.92349 (−2.1%); **production window on 10³** (`f10.951_w12.797_k10`, edges moved 11/13 meV into holes of the 10³ NSCFs; wannier90 221248, 221296; c2 221324-34): **6.44024 / 6.59839, +0.83% from 8³** (6.38695 / 6.53961); smallest direct gap 196.3 meV (8³: 195.7) | – | MBT done |
| larger trial set | probe set (Bi p,d; Te p,d; Mn s,d; s at the vdW gaps; 220 functions; pw2wannier90 221258, `dft/base/pw2wannier_k8_probe/projectability/`): the bands above VBM + 4 eV are **Te d** (0.69–0.81 of their weight; Bi d, Mn s, gap s ≤ 0.06), equally on Te1 and Te2 and over all five d orbitals; even the probe set describes every band only to VBM + 4.0 eV (lowest band 0.74 there). **The 240-band NSCF caps any set at ~144 WFs** (the outer window must stay below band 240, VBM + 8.42 eV, and hold num_wann + 56 bands at every k), so Te d (+80) needs NSCFs with more bands for the base and the 4 patterns: decision pending | `Yd`: Y p,d; Ir d; O s,p (216 WFs), pw2wannier90 221257; projectability ≥ 0.97 for every band to VBM + 3 eV, ≥ 0.82 to VBM + 4.25 eV (production set ≥ 0.93 only to VBM + 2.5). Base Ω_I 169.2 / 182.1 Å², largest spread 2.46 / 2.47 Å² (no diffuse functions). Q = 0.01 Å, undistorted end, change from production: window `Yd_f13.444_w22.644_k8` (VBM + 3.0): −0.016838 (−1.62%), 0.053515 (−2.74%), 0.063600 (−0.42%); window `Yd_f14.554_w22.644_k8` (VBM + 4.11): −0.017125 (+0.06%), 0.054086 (−1.70%), 0.062915 (−1.49%); displaced ends move alike; smallest direct gap 0.955 eV in all. **Model dependence up to 2.7%, not the 0.75% of the window-only study** | YIO done (11-04); MBT decision pending |
| projectability of the production sets (weight of each band in the span of the trial functions, from the `.amn`; `dft/base/pw2wannier_k8/projectability/`) | mean 0.84–0.87 per 1-eV bin up to VBM + 3 eV, 0.78 at 3–4 eV, **0.23 at 4–5 eV**, ≤ 0.16 above: the production frozen top (3.69) is at the edge | ≥ 0.93 for every band up to VBM + 2.5 eV (the t₂g manifold; production frozen top 1.98), 0.85 at 3–4 eV (Ir e_g + O p), **0.38 at 4–5 eV** (Y d), ≤ 0.26 above | done (221259/60) |
| magnetic symmetry of the models | all patterns (and pattern 3 at 0.03 Å) pass: total ≤ 6×10⁻⁵, DFT ≤ 0.063 μeV; broken inversion 0.4–1.8% (5.8% at 0.03 Å) | base m-3m′ and patterns −4′3m′ pass: total ≤ 5×10⁻⁶, DFT ≤ 1.35 μeV; controls ~200%, broken C4z·T 1.8–12% | done, current |
| exchange–correlation | PBE not usable (bug #18) | PBE+U vs LDA+U at the 10-12 settings: Y −3%, O −11% | YIO old settings |
| Hubbard U | – | pattern 3, Q = 0.01 Å, undistorted end (`data/Y2Ir2O7/studies/u_scan/`): U = 3.00: 0.063867 (gap 0.69 eV indirect); **U = 2.75: 0.06100 (−4.5%; gap 0.587 eV indirect, smallest direct on 10³ 0.846 eV)**; **U = 2.50: 0.05490 (−14.0%; 0.491 / 0.732 eV)**: the response falls as the gap closes; U = 2.25 (|m|abs 2.98, Ir 0.266 μB) and 2.00 to come | running |

### MnBi₂Te₄ window study (pattern 3, Q = 0.01 Å, undistorted end, tube list, 8³)

Band errors: interpolated vs DFT on the 104-point path. "gap edges": the valence top and conduction bottom; "≤ 3 eV":
every DFT band between the VBM and VBM + 3 eV to its nearest Wannier band (max / median over k).

| window | frozen top − VBM (eV) | ∂θ/∂Q | gap-edge error, max on the path (meV) | ≤ 3 eV error, max / median (meV) |
|---|---|---|---|---|
| **f10.962_w12.810 (production)** | 3.69 | **6.3869** | 7.6 | 47.7 / 8.4 |
| f10.962_w13.046, w13.697, w14.020, w15.171 | 3.69 | 6.3970 / 6.3515 / 6.3645 / 6.3976 | 7.1–8.7 | 47–53 / 8.5–9.5 |
| f10.790_w15.169 | 3.52 | 6.4672 | 9.4 | 51.8 / 8.2 |
| f10.698_w13.915 | 3.43 | 6.2956 | 8.2 | 54.7 / 8.8 |
| f10.278_w12.814, w13.696, w15.169, w15.412 | 3.01 | 6.4877 / 6.4364 / 6.4323 / 6.4163 | 11.7–14.0 | 116–197 / 10–15 |
| f9.564_w15.169 | 2.30 | 5.7773 | 44.0 | 168.6 / 119.1 |
| f8.395_w13.915 | 1.13 | 7.0745 | 43.9 | 355.3 / 264.9 |
| f8.344_w15.169 | 1.08 | 7.5605 | 28.3 | 468.6 / 301.9 |
| f7.809_w15.169 | 0.54 | 6.2628 | 580.2 | 575.2 / 379.8 |

The smallest direct gap on the tube list is 191–197 meV in every window (DFT 195 meV on the path), and the gap at the
smallest-gap path point is within 4 meV of DFT in every window: the band energies at the gap do not decide the
response. The θ-odd mass m₅ a T₂⁻ displacement opens enters the gap only as √(m² + m₅²), at second order, while θ
responds at first order; the curvature also involves the interband matrix elements to every conduction state of the
model, so conduction bands misplaced by 0.1–0.6 eV (the low frozen tops) change it. The scatter shrinks as the frozen
top rises to the 116-WF ceiling. Data: `data/MnBi2Te4/studies/window_tests/assemble*.log`, per-window band and gap
data in `windows/<window>_k8/bands_and_parts_T2-_p3_Q0.01/parts.npz`.

### Y₂Ir₂O₇ window study (Q = 0.01 Å, 10³; 10-19 settings, see the table above)

| frozen / outer top (eV) | Ω_I (Å²) | pattern 1 | pattern 2 | pattern 3 |
|---|---|---|---|---|
| **12.422 / 22.644** (production window) | 117.10 | reference | reference | reference |
| 10.786 / 22.644 | 116.84 | −0.75% | +0.25% | −0.56% |
| 13.120 / 22.644 | 117.37 | −0.24% | +0.02% | +0.06% |
| 12.422 / 19.993 | 119.47 | +0.44% | +0.10% | −0.17% |
| 12.422 / 24.522 | 117.34 | −0.14% | −0.00% | −0.03% |
| 14.396 / 22.644 | 158.55 | +5.67% | −1.07% | +0.32% |
| 14.396 / 24.522 | 159.25 | +5.56% | −1.20% | +0.34% |

Log: `notes/log/2026-10-17_yio_amplitude_series.md`, "YIO window study: result".

## PDOS: where the shells are (states per cell, base 8³ NSCF, projwfc.x)

MnBi₂Te₄ (VBM 7.268 eV): Bi p and Te p carry the bands to VBM + 4 eV, Mn d sits at VBM + 2–3 eV (6.1 states).
**Above VBM + 4 eV most of the DOS projects onto no pseudo-atomic orbital** (VBM + 4–5 eV: 5.5 of 8.0 states; 6–7 eV:
9.2 of 9.8): higher shells (Bi 6d, Te 5d/6s) or nearly-free states in the van der Waals gaps. This is why 116 WFs
cannot freeze above VBM + 3.69 eV.

Y₂Ir₂O₇ (VBM 10.4437 eV): Ir t₂g to VBM + 2 eV; Ir e_g (with O p) at VBM + 3–5 eV; **Y 4d from VBM + 4 eV** (6.6,
5.0, 6.7, 7.3 states in the 1-eV bins from 4 to 8 eV); Ir 6s/6p, Y 5s/5p above VBM + 6 eV. The production set freezes
only the t₂g manifold; Ir e_g, in the set, lies in the outer window among Y 4d states that are not.
