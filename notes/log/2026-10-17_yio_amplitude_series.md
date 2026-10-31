# 2026-10-17 — YIO Γ₂⁻ amplitude series: E(Q), F(Q) and θ(Q)

**Purpose (user, 10-17).** A paper figure showing E quadratic and F linear in Q, and θ linear in Q, along the three
LDA basis modes. MBT gets the same E/F series later (4 modes × 0.01 / 0.02 / 0.05 Å, plus the existing 0.03 Å set);
no θ(Q) for MBT (4.7 days of NSCF per amplitude).

## Choices

- **Basis modes, not normal modes.** Basis-mode E(Q)/F(Q) is the harmonic test along the coordinates dθ/dQ uses;
  normal-mode runs would repeat the 4×4 / 3×3 eigen-decomposition (about 40 SCFs instead of 26). For MBT, the
  rigid-layer runs (log 10-15) already sample the 17.8 cm⁻¹ eigenvector (overlap 0.9995).
- **Q > 0 only.** Γ₂⁻ is odd under inversion, which is a unitary symmetry of the m-3m′ base, so the −Q structure
  is the inversion image of +Q: E(−Q) = E(Q), F reflects, θ(−Q) = −θ(Q) exactly. The figure reflects the data.
- **conv_thr 1e-10** everywhere, including a new base SCF (with `tprnfor`; the 210527 base did not print forces).
  The 0.01 Å runs repeat the 2026-10-12 ones at 1e-10, which is also the direct test of whether conv_thr 1e-7
  limits dθ (log 10-15). The Wannier base stays 212615 (1e-7), so that test covers the mode side only.
- **Unrelaxed base:** its forces are inversion-even, so they have no Γ₂⁻ projection: no linear term in E(Q) and no
  offset in F(Q).

## Runs (submitted 2026-10-17 13:36, pw.x `db4b99b4`)

Script: `archive/setup_templates/Y2Ir2O7/modes/setup_amplitudes_YL2.sh` (cluster copy in `Y2Ir2O7/phonon/`).
Positions: base + (Q/0.01)·(Q = 0.01 Å positions − base), so Q is the displacement per Cartesian component along
the raw SMODES vector (the 2026-10-12 convention); the 0.01 Å inputs reproduce the 10-12 positions exactly. Each
mode job is SCF (1e-10, `tprnfor`, `startingpot = 'file'` from the YL2 density) then, if converged, the 512-k NSCF
with the base NSCF input (`nbnd`, k list unchanged).

| Q (Å) | mode 1 (Y) | mode 2 (Ir) | mode 3 (O) |
|---|---|---|---|
| 0.01 | 213999 | 214003 | 214016 |
| 0.02 | 214000 | 214004 | 214011 |
| 0.03 | 214001 | 214005 | 214012 |
| 0.05 | 214002 | 214015 | 214013 |

Base at 1e-10: 213998 (SCF only). Mode 3 at 0.02 / 0.03 / 0.05 Å were moved from the full dko48m4 queue to the idle
dko64m4 nodes (64 slots, `send_job_64.sh`; first ids 214008-214010 deleted while queued). The last two queued jobs were resubmitted excluding only n138
(mode 2 0.05 Å 214006 → 214015 on n126, mode 3 0.01 Å 214007 → 214016 on n246): n126/n246/n247 had been excluded with
n138 without a reason, but they mount `/mnt/wkdk`. All 13 running from 13:50. Cluster folders: `Y2Ir2O7/phonon/mode_k/Q_0p0X_conv1e-10_2026-10-17/`,
base `.../YL2_lda_relpz/conv1e-10_2026-10-17/`.

## Next

1. Per SCF: 24 operations (base 48), total magnetization 0, |m|abs ≈ 3.21, E − E_base against −½F·Δr.
2. Per NSCF: the base windows (14.396 / 24.522 eV and the lower edges) in gaps of every k; the 0.05 Å runs move
   bands ~5× more than 0.01 Å.
3. Mode `.amn` on the login node (`--no-displacement`, base 212615 `.chk`), wannier90, then the c₂ sweep (shared
   WS, `.chk` H and A(R)) to nk 12 for each run.
4. Figure runs: delete each run's `tmp/` after its Wannier step (about 108 GB each untarred); no `tmp.tar.gz`.

## Wannier chains as held SGE jobs; julia on compute nodes (14:00)

- **Julia now runs on compute nodes** from a copy on the shared filesystem: `/mnt/wkdk/bc889/julia/{depot,QEWavefunctions.jl}`
  (`JULIA_DEPOT_PATH`; test job 214017 on n247 precompiled once in 80 s, then loaded). The login-node-only rule
  (compute nodes mount a different `/mnt/home`) no longer forces the `.amn` step onto the login node.
- **YIO amplitude series:** `setup_amp_w90.sh` (template `modes/`, cluster `Y2Ir2O7/phonon/amp_pipeline/`). Per run, a
  w90 folder with the base 212615 `.win` (this run's `atoms_cart`, `-pp` done) and two held jobs:
  - amn (8 slots, held on the SCF+NSCF job): NSCF check, gunzip the save, `window_check.py` (every base window edge
    ≥ 5 meV from every eigenvalue), `generate_amn.jl --no-displacement`;
  - w90 (48 slots, held on amn): pw2wannier90 `-nk 8` (mmn, eig), wannier90 `-np 1`; then deletes the run's `tmp/`.

  | Q (Å) | mode 1 amn / w90 | mode 2 | mode 3 |
  |---|---|---|---|
  | 0.01 | 214018 / 214019 | 214026 / 214027 | 214034 / 214035 |
  | 0.02 | 214020 / 214021 | 214028 / 214029 | 214036 / 214037 |
  | 0.03 | 214022 / 214023 | 214030 / 214031 | 214038 / 214039 |
  | 0.05 | 214024 / 214025 | 214032 / 214033 | 214040 / 214041 |

  Markers in each w90 folder: `AMN_OK`, `W90_OK`, or `NSCF_FAIL` / `GUNZIP_FAIL` / `WINDOW_FAIL` / `AMN_FAIL` / `W90_FAIL`.
- **Base SCF at 1e-10 (213998):** E = −4499.85691367 Ry (5e-8 Ry below the 1e-7 run), 11 iterations, 48 operations,
  |m|abs 3.21, total force 0.0068 Ry/bohr (the O x parameter; no Γ₂⁻ part).
- **MBT 8×8×8 chain** (`archive/setup_templates/MnBi2Te4/modes/setup_mbt_k8_w90.sh`, cluster
  `MnBi2Te4/w90_k8_pipeline/`): base prep 214042 (held on NSCF 212997: gunzip, window check including band 240 above
  15.047 eV, `tmp.tar.gz`) → base w90 214043 (pw2wannier90 amn/mmn/eig, wannier90; h_rt 248 h); per mode k = 1-4
  an amn job 214044/46/48/50 (held on 214043 and the mode NSCF; same checks, then the `.amn` from the base `.chk`)
  → w90 214045/47/49/51. The `.win` is the recentred 4×4×4 `.win` (window_o15_f10) with num_bands 240, mp_grid 8 8 8,
  the NSCF's 512 k in order; mode sites get the same lattice shift (0 or −(a₁+a₂+a₃)) as the base site. Production
  runs: `tmp/` and `tmp.tar.gz` are kept. pw2wannier90 at 8×8×8 is expected to take ~8× the 4×4×4 3h25m.
- **Cleanup script fixed before anyone ran it:** `scripts/cluster_cleanup_2026-10-11.py` would have deleted the MBT base
  8×8×8 NSCF (`209653bc889/nscf_k8_2026-10-12`) and the amplitude-series saves; both are now in its keep list.

## Cluster tidy: every run folder marked (user, 14:15)

`scripts/cluster_mark_status_2026-10-17.sh` (run on the cluster): a `_STATUS` file in each run folder (PRODUCTION /
COMPARISON / VALIDATION / MIXED / SUPERSEDED, with what it is and where production lives) and an index
`/mnt/wkdk/bc889/axion/INDEX.md` (75 folders, sizes). The 44 `.bak` files and 2 stopped SCF outputs moved into `_bak/`
next to where they were; 11 superseded YIO setup scripts into `Y2Ir2O7/phonon/_old_scripts/`; 3 old cleanup scripts
into `axion/_old_scripts/`. Nothing deleted, and no run folder moved: held jobs and the cleanup script use their
paths. **SUPERSEDED totals ~5.2 TB** (MBT PBE base U scan 805 G and PBE phonons 935 G; YIO pre-fix runs ~3.3 TB):
deleting those folders would bring the cluster from 7.0 TB to ~1.8 TB (user's decision).

## E(Q) and F(Q) results (SCFs, 15:00)

All 12 SCFs converged: 20-28 iterations, 24 operations, total magnetization 0, |m|abs 3.21 (mode 2: 3.20 at 0.03 Å,
3.18 at 0.05 Å, the Ir moment following the Ir displacement), SCF force correction ≤ 9e-4 eV/Å, v·F(0) ≤ 2e-14.
Notebook: `analysis/amplitude_series.ipynb` (data in `data/Y2Ir2O7/phonon/GM2-/modeN/Q_0.0XA/lda/scf_nscf_conv1e-10_<job>/`,
base `base/lda/scf_conv1e-10_213998/`).

| | ΔE at 0.01 / 0.02 / 0.03 / 0.05 Å (meV) | 2ΔE/Q² ÷ (−f/Q) | k from E with Q⁴ / from F with Q³ (eV/Å²) | a₄ from E / b/4 from F (eV/Å⁴) |
|---|---|---|---|---|
| mode 1 (Y) | 21.95 / 88.07 / 198.46 / 553.95 | 0.997 / 0.999 / 0.997 / 0.993 | 439.8 / 440.0 | +669 / +645 |
| mode 2 (Ir) | 24.62 / 98.52 / 221.25 / 611.51 | 1.000 / 1.002 / 1.003 / 1.008 | 493.1 / 492.7 | −777 / −719 |
| mode 3 (O) | 7.52 / 30.01 / 67.44 / 187.10 | 1.004 / 1.002 / 1.001 / 1.001 | 150.0 / 149.8 | −67 / −19 |

- **Harmonic to ≤ 1.5% at 0.05 Å** and ≤ 0.3% at 0.01 Å, in both E and F. The deviation is a quartic term: mode 1
  stiffens, mode 2 softens, mode 3 is harmonic. Energy and force agree on it (b ≈ 4a₄ for modes 1 and 2), and their
  separate fits agree on k to 0.04 / 0.08 / 0.17%. No Q³ term, as inversion requires.
- **The scatter at 0.01 Å** (E/F up to 0.4% off for modes 1 and 3) is a few tens of μeV in ΔE, the smallest energy
  differences of the set (7.5-25 meV); it shrinks as 1/Q² for mode 3, as an energy offset would.
- **Phonon frequencies** from the full 3×3 Φ at each amplitude: 228.6 / 328.3 / 471.6 cm⁻¹ at 0.01 Å, 227.8 / 330.3 /
  471.0 at 0.05 Å (log 10-12 §16: 229 / 328 / 471 at 1e-7). Φ is symmetric to 0.03% at 0.01 Å (0.3% at conv_thr 1e-7)
  and 0.7% at 0.05 Å, where the cubic force terms enter.

## Window edges move with the mesh and with Q (16:00-18:15)

**The window check now also compares with the base.** `window_check.py` takes the base NSCF as a fourth argument and
requires the same number of bands below each edge at every k, not only a 5 meV distance. A band that crosses an edge
between base and mode freezes (or admits) different states on a whole k-star, so the two Wannier models are built
differently even though both edges sit in gaps. `window_search.py` (cluster `amp_pipeline/`) scans a 1 meV grid for
edges that pass in every run of a set.

**YIO amplitude series.** With the production edges (dis_froz_max 14.396, dis_win_max 24.522):
- mode 3, band 177 on a 48-k star (k 93, 95, 100, ...): 14.4267 eV (base) → 14.4209 (0.01) → 14.4092 (0.02) → 14.3959
  (0.03, 0.1 meV from the edge) → 14.370 (0.05, below it: one more frozen state on the star); mode 2 at 0.05 Å crosses
  the same edge (0.2 meV); the outer edge 24.522 also changes counts within the set.
- The production 0.01 Å pairs (212755-7 with 212615) have identical counts at both edges (margins 9-20 meV): the
  published dθ are not affected.
- **Series edges:** dis_froz_max **12.422** (centre of a gap 11.78-13.06 eV open in all 13 NSCFs, ≥ 645 meV from
  every eigenvalue; 160 frozen states = 156 occupied + 4), dis_win_max **22.644** (20 meV margin; 238-244 states in the
  outer window). Base re-Wannierized with wannier90 only on 212615's `.amn`/`.mmn`/`.eig`:
  `nscf_from_210527/w90_trial02_base_ampwin_2026-10-17/` (job 214064; `.chk`/`.nnkp` in `ref/`). Chain v2:
  `setup_amp_w90_v2.sh`, per run `w90_trial02_ampwin_2026-10-17/`, amn/w90 jobs 214065-214088 (mode 1: 214065-72,
  mode 2: 214073-80, mode 3: 214081-88, in Q order). The first chain (214018-214041) was cancelled before it changed
  anything; its `w90_trial02_no-displacement/` folders in the runs hold only setup files and failure markers.
- The series' 0.01 Å point against the production 0.01 Å pair is then a test of dθ against the window choice.
- **Series base (214064, 18:44):** Ω_total 147.47 Å² (Ω_I 117.10, Ω_D 0.003, Ω_OD 30.37) against 193.76 for 212615 (Ω_I
  159.26): the narrower window gives more localized functions; centres unchanged (sum identical).

**Julia cache race.** Two `.amn` jobs that started in the same minute on different nodes both re-precompiled into the
shared depot; one died with "Stale file handle" on a replaced `.so` (mode 3 0.01 Å, 214034). Every job now uses
`JULIA_DEPOT_PATH=$TMPDIR/jdepot:/mnt/wkdk/bc889/julia/depot` (compile cache on the node, 80 s; tested with two
simultaneous jobs 214054/214055).

**MBT 8×8×8 windows.** The 4×4×4 edges are 0.2 meV (10.009) and 0.5 meV (15.047) from bands on the 8×8×8 mesh (the
denser mesh fills those small gaps). Edges open in the base and all four mode NSCFs, with base-consistent counts:
**dis_froz_max 10.278** (8.6 meV margin; 162-170 bands below, i.e. 106-114 frozen states for 116 WF) and
**dis_win_max 15.169** (7.8 meV; 196-226 bands below, under num_bands 240); dis_froz_min = dis_win_min = −7 is in a
1.1 eV gap. Set in the five `.win` files before any job ran (base prep 214042 held, edited, released 18:15). Mode
jobs were resubmitted with the julia fix and the base-consistency check: amn 214056/58/60/62, w90 214057/59/61/63.
8×8×8 NSCFs: JOB DONE in all five (4.6 days each), gap 7.268 → 7.466 eV (198 meV, as on 4×4×4), 174 GB save each.

## MBT amplitude series submitted; which SCFs the dθ runs use (22:50)

**MBT T₂⁻ amplitude series** (`archive/setup_templates/MnBi2Te4/modes/setup_mbt_amplitudes.sh`, cluster
`MnBi2Te4/phonon_lda/`): base SCF at conv_thr 1e-10 with tprnfor (`base/soc_lda/u_4.0/forces_conv1e-10_2026-10-17/`,
214107) and 4 modes × Q = 0.01 / 0.02 / 0.05 Å at 1e-10 (`T2-/mode_k/Q_0p0X_conv1e-10_2026-10-17/`, 214108-214120;
mode 1: 214108-10, mode 2: 214111/13/14, mode 3: 214115-17, mode 4: 214118-20); with Q_0p03_2026-10-15 (213555-8) that
is 4 amplitudes per mode. Inputs: Q_0p01_kfix scf.in with scaled positions and conv_thr 1e-10, base 209653 start.
SCF only (save deleted in the job).

**Convergence of the SCFs behind each number:**

| | dθ (production) | force constants / frequencies |
|---|---|---|
| YIO | base 210527 and modes 212755-7, conv_thr 1e-7 (→ NSCF → W90) | amplitude series, 1e-10 (base 213998) |
| YIO θ(Q) series | modes 1e-10, but the base Wannier model is from the 1e-7 base NSCF 210590 | — |
| MBT | base 209653 and modes 212708-11, conv_thr 1e-7 (→ 8×8×8 NSCFs) | Q = 0.03 at 1e-10 (213555-8), base 210628 (1e-7, reached 3.8e-9); series at 1e-10 submitted |

A base error δ enters the finite difference as δ/Q, so a 1e-7 base matters most at the smallest amplitude. The series
saves are kept as hard links (`<job>bc889/tmp_keep/`, no extra space) so the `.amn` can be redone from a 1e-10 base
if wanted.

**YIO base at 1e-10 through Wannier (user, 23:00):** `archive/setup_templates/Y2Ir2O7/base/setup_base_conv10_w90.sh`.
NSCF 214124 from SCF 213998 (`YL2_lda_relpz/conv1e-10_2026-10-17/nscf_from_213998/`), then 214125: pw2wannier90 once,
wannier90 with the series windows (`w90_trial02_base/ref_ampwin/`) and the production windows (`ref_prodwin/`), both
window sets checked against this NSCF and the 1e-7 base NSCF 210590. Next, after it and the current series jobs:
the 12 series `.amn` from ref_ampwin (wannier90 only, reusing tonight's `.mmn`/`.eig`) and the three Q = 0.01 runs
from ref_prodwin. Comparisons: series on the 1e-7 base vs on the 1e-10 base = base-threshold effect at every Q; the
1e-10 production-window pairs vs the production pairs = the full conv_thr effect on dθ; the 1e-10 pairs with the two
window sets = the window effect. The MBT amplitude SCFs got priority −100 (`qalter -p`) so this NSCF runs first.

**Shared-filesystem bottleneck (10-18 00:45):** the 11 series w90 jobs started together, each `cp -r` of an uncompressed
108 GB save from `/mnt/wkdk`: after ~2 h they held 32-72 GB each (aggregate ~100 MB/s), so pw2wannier90 starts
~2-3 h late. For the next chain: stage the compressed `tmp.tar.gz` (~48 GB, `pigz -d` on the node), or stagger the
jobs (`-hold_jid` in pairs); pw2wannier90 reading the save over NFS directly would move the same data.

**MBT base 8×8×8 Wannier done (214043, 10-18 ~02:00):** windows −7 / 10.278 / 15.169 eV; Ω_total 432.70 Å² (Ω_I 325.53,
Ω_D 0.72, Ω_OD 106.45), spreads 0.60-7.51 Å² (4×4×4, windows 10.009/15.047: 345.5, largest 5.67; the coarse mesh
underestimates spreads). All 116 centres within 0.24 Å of their projection sites (median 0.05 Å): no wrap (#20).
Base SCF at 1e-10 (214107): E = −6947.39550995 Ry, 1e-8 Ry from 210628, 7 iterations, |m|abs 9.99, total force
0.0326 Ry/bohr as before. Mode amn jobs 214056/58/60/62 released.

**MBT 8×8×8 base checks (10-18, `archive/diagnostics/2026-10-18_mbt_k8_base_symmetry/`, pinned positions, WS 5e-3):**
- Mesh: Wannier = DFT to 0.000 meV for the 82 occupied bands. Interpolated gap on a shifted 12³ grid 199.1 meV (DFT on
  the mesh 198.1); the 4×4×4 model gave 105-119 meV there, so the 8×8×8 mesh removes the off-mesh error.
- Eigenvalues: E(gk) = E(k) to ≤ 1e-4 meV for all 12 operations of −3m′.
- Curvature: Tr F is PT-forced to zero in the base, so the test uses T_ab = Re Tr[F_a F_b] (gauge invariant,
  T(gk) = R T Rᵀ for every operation). Residuals: internal 0.0000%, cross ≤ 0.019%, external ≤ 0.13%, total
  ≤ 0.0065%, all in the C₃-related operations and in the A(R) parts (A(R) uses the off-site p_x/p_y centres, #21);
  the external part is 0.05% of T. 3.7 h on the laptop (H(R)/A(R) with 941 R).

**1e-10 base and rebase (10-18 morning):** NSCF 214124 (gap 10.4437 → 11.1286 eV; 1e-7 base 10.4433 → 11.1283) and
Wannier 214125 (pw2wannier90 once, wannier90 with both window sets): Ω_total 147.466 Å² (series windows; 1e-7 base
147.472) and 193.744 Å² (production; 1e-7 base 193.757); both window sets pass against this NSCF with the same bands
below each edge as 210590. `ampwin_full/` (214162) reruns the series-window wannier90 for a complete output set
(`_centres.xyz` from the same run as the `.chk`; the base job copied only `.chk`/`.nnkp`/`_hr.dat` to `ref_*/`).
Rebase jobs (`setup_amp_rebase10.sh`): `w90_trial02_ampwin10_2026-10-18/` for all 12 runs (214147-214158, mode 1
Q order then mode 2, mode 3) and `w90_trial02_prodwin10_2026-10-18/` for the three 0.01 Å runs (214159-214161); each
an 8-slot job held on the run's series w90 job: `.amn` from the 1e-10 base, wannier90 with the series job's
`.mmn`/`.eig`, mode save from `tmp_keep/`.

**dθ/dQ so far (series windows, 1e-7 base, nk 10; base end / mode end / mean):**

| | 0.01 Å | 0.03 Å | 0.05 Å |
|---|---|---|---|
| mode 1 | | −0.016971 / −0.018027 / −0.017499 | |
| mode 2 | | 0.054803 / 0.050763 / 0.052783 | 0.053957 / 0.047325 / 0.050641 |
| mode 3 | 0.063491 / 0.064077 / 0.063784 | | 0.063297 / 0.063605 / 0.063451 |

nk 10 = nk 12 to 1e-7 (mode 2, 0.03 Å), so the series runs at nk 10 only. Mode 3 at 0.01 Å is −2.6% from the
production pair (0.065488): windows (and the base threshold, which barely moves the Wannier model); the 1e-10
production-window pairs separate the two.

## MBT T₂⁻ amplitude series results (10-18 ~10:00)

All 13 SCFs JOB DONE (15-21 iterations, 16 k, 12 operations, |m|abs 9.99, E above the base). Base at 1e-10 (214107).
Data: `data/MnBi2Te4/archive/phonon/T2-/modeN/Q_0.0XA/lda/scf_conv1e-10_<job>/` (0.03 Å: `scf_21355x`), base
`data/MnBi2Te4/phonons/amplitude_series_conv1e-10/base/214107bc889/`.

| basis mode | ΔE 0.01 / 0.02 / 0.03 / 0.05 Å (meV) | 2ΔE/Q² ÷ (−f/Q) | k (eV/Å², E fit / F fit) |
|---|---|---|---|
| 1 (Mn) | 0.331 / 1.330 / 2.986 / 8.275 | 0.993 / 0.998 / 1.001 / 1.002 | 6.646 / 6.654 |
| 2 (Bi) | 1.549 / 6.195 / 13.931 / 38.685 | 1.001 / 0.999 / 1.000 / 1.002 | 30.966 / 30.994 |
| 3 (Te1) | 1.029 / 4.133 / 9.312 / 25.828 | 0.992 / 0.997 / 1.000 / 1.004 | 20.696 / 20.753 |
| 4 (Te2) | 0.641 / 2.550 / 5.735 / 15.911 | 1.007 / 1.001 / 1.001 / 1.000 | 12.755 / 12.733 |

- The Mn E/F of 1.24 (Q = 0.01, conv_thr 1e-7, log 10-15) was SCF noise: 0.993 at 1e-10.
- Harmonic within 0.5% at 0.05 Å; v·F(0) ≤ 3e-13 eV/Å.
- Frequencies (4×4 Φ, Phi c = ω² M c): 17.95 / 82.80 / 128.89 / 143.50 cm⁻¹ at 0.01 Å; 17.62 / 82.74 / 128.41 / 143.20 at
  0.05 Å; **Q → 0 (Φ_ij + c_ij Q² fit over 4 amplitudes): 18.02 / 82.78 / 128.95 / 143.52 cm⁻¹**, Φ symmetric to
  0.073% (force noise; YIO 0.007%). Soft mode 18.0 cm⁻¹ as the rigid-layer energy estimate (18.0-18.2).

## conv_thr matters for dθ through the MODE SCF, not the windows (10-18 ~14:30)

| mode 3, Q = 0.01 Å, nk 10 | windows | base / mode conv_thr | base end | mode end | mean |
|---|---|---|---|---|---|
| production pair (212615 + 212995) | production | 1e-7 / 1e-7 | 0.065017 | 0.065959 | 0.065488 |
| prodwin10 (214125 + 214161) | production | 1e-10 / 1e-10 | 0.063530 | 0.064464 | 0.063997 |
| series rebased (214162 + 214155) | series | 1e-10 / 1e-10 | 0.063493 | 0.064076 | 0.063784 |
| series (212615 inputs + 214082) | series | 1e-7 / 1e-10 | 0.063491 | 0.064077 | 0.063784 |

- **Base threshold: negligible** (series vs rebased, mode 3 at all four Q: ≤ 2e-5 relative).
- **Windows: +0.3%** (prodwin10 vs series rebased, both 1e-10).
- **Mode SCF threshold: −2.3%** (production vs prodwin10, same windows): the displaced structure's SCF at 1e-7 shifts
  dθ by ~2%, amplified ~1/Q by the finite difference. This corrects the earlier reading (above) that the −2.6% was
  the window choice, and the 10-15 estimate (≲ 0.2% from force bookkeeping).
- **Consequences:** the production YIO dθ (1e-7) carry ~2% (mode 3; modes 1-2 to follow); the MBT 8×8×8 chain is
  built on 1e-7 mode SCFs (212708-11). Use 1e-10 mode SCFs for every dθ.

## MBT mode 3 Wannier: Mn d spreads 30% larger than the base (10-18 ~17:00)

Mode 3 (Te1, 214061): Ω_total 453.70 Å², Ω_I 344.50 (base 214043: 432.70 / 325.53), +21 Å² concentrated in the Mn d
WFs (100-106 on Mn1, 109-115 on Mn2: 4.016 → 5.362 and 3.234 → 4.128 Å²); centres unchanged (≤ 0.02 Å). YIO mode
runs match their base to 0.5%. The DFT is not the cause: base vs mode NSCF eigenvalues differ by ≤ 8-27 meV per band
(median 3-6 meV), no offset. Q = 0 self-test submitted (214245, `w90_sp_d_k8_2026-10-17/selftest_Q0_2026-10-18/`):
the mode pipeline (`.amn` from the base WFs, `--no-displacement`) applied to the base's own Bloch states must give
back Ω_I = 325.53; if not, the pipeline itself changes the MBT subspace and the mode dθ would carry that artifact.

**YIO mode 2 (Q = 0.01) at 1e-10:** production windows 0.053541 / 0.054106 (mean 0.053824; production 1e-7 0.054165:
−0.6%); series windows (1e-7 base, 1e-10 mode) 0.055258 / 0.053728: the window choice moves the base-end estimate by
3.2% for the Ir mode (mode 3: 0.06%).

## The mode pipeline does not reproduce the base at Q = 0 (10-18 ~18:15)

**MBT Q = 0 self-test (214245):** `.amn` from the base WFs (`--no-displacement`) onto the BASE Bloch states, wannier90
with the base `.mmn`/`.eig`/`.win`: Ω_I 344.45, Ω_total 453.64 Å², vs 325.53 / 432.70 for the base Wannierized from
atomic trials (214043), and 344.50 / 453.70 for mode 3 at 0.01 Å. The 0.01 Å displacement changes almost nothing; the
6% jump is the construction (base-WF projections vs atomic projections through the same projection-only
disentanglement). The cause in the `.amn`/disentanglement path is not identified (in exact arithmetic, base WFs as
trials should return the base subspace).

**Consequence:** every mode model tends, as Q → 0, to the self-test model, not to the base it is paired with: the
finite difference carries a fake difference δ/Q. δ is inversion-even (both base models are symmetric), so it does not
change θ at first order at the symmetric base end, and enters the mode end as δ·O(Q)/Q = O(1). This matches the YIO
series exactly: base end a + bQ² (clean) in every mode, mode end with a Q-independent 0.5-1% offset. YIO has the
same mismatch at a smaller scale (mode 3 at 0.01 Å: Ω 146.77 vs base 147.47).
*(2026-10-18 late: the test below does not support this as stated; see "Self-test base: results" at the end.)*

**Fix:** pair the modes with the self-test model as the base (same construction). MBT: 214245 is that base. YIO: self-
test on the 1e-10 base with the series windows submitted (214249, `.../w90_trial02_base/selftest_ampwin_2026-10-18/`);
then recompute a few amplitudes against it and check that the two ends meet as Q → 0.

**Mode-SCF threshold, all three YIO modes (Q = 0.01, production windows, 1e-10 vs 1e-7 means):** mode 1 −0.018140 vs
−0.018051 (+0.5%), mode 2 0.053824 vs 0.054165 (−0.6%), mode 3 0.063997 vs 0.065488 (−2.3%).

**Window dependence at fixed 1e-10 (Q = 0.01 Å; series windows 12.422/22.644 → production 14.396/24.522):**

| | base end | mode end |
|---|---|---|
| mode 1 | −0.017045 → −0.018082 (+6.1% in magnitude) | −0.017327 → −0.018198 (+5.0%) |
| mode 2 | 0.055297 → 0.053541 (−3.2%) | 0.053689 → 0.054106 (+0.8%) |
| mode 3 | 0.063493 → 0.063530 (+0.06%) | 0.064076 → 0.064464 (+0.6%) |

The base end is immune to the base/mode construction mismatch, so this is a genuine dependence on which conduction
states (12.4-14.4 eV) are frozen; the occupied states are frozen and identical in both. Largest systematic found so
far (up to 6%). Needs a window convergence study (frozen top at several gaps open at every amplitude).
Chain complete 18:41: 12 series + 12 rebased (series windows) + 3 rebased (production windows) pairs.

## Self-test base: results (2026-10-18 late)

Same mode model, two bases; `transl_inv_full_from_chk`, shared WS, H from `.chk`, no cut.

**YIO mode 1, Q = 0.01 Å, nk 10, series windows** (mode 214147; `run_selftestbase.py`):

| base | Ω_I (Å²) | base end | mode end | (mode − base)/base |
|---|---|---|---|---|
| atomic trials, 1e-10 (214162) | 117.10 | −0.017045 | −0.017327 | 1.65% |
| self-test (214249) | 116.42 | −0.017134 | −0.017234 | 0.58% |

**MBT mode 3, Q = 0.01 Å** (mode 214061, mode SCF 1e-7, atom-pinned; `../../archive/diagnostics/2026-10-18_mbt_k8_dtheta/`):

| nk | self-test base 214245 (Ω_I 344.45): base / mode end | atomic base 214043 (Ω_I 325.53): base / mode end |
|---|---|---|
| 4 | 48.333 / 49.892 | 48.084 / 49.881 |
| 6 | 21.105 / 21.795 | 20.966 / 21.794 |
| 8 | 12.260 / 12.645 | 12.184 / 12.647 |

- MBT: changing the base moves the mode end ≤ 0.02% and the base end 0.5-0.7%; the end gap (3.1-3.7%) does not depend
  on the base. The predicted O(1) mode-end offset from δ/Q is absent despite the 6% Ω_I difference.
- YIO: both ends move ~0.5%; the end gap drops from 1.65% to 0.58%. (I first called the 0.58% "far above O(Q²)",
  comparing it with the base end's own Q² coefficient; wrong comparison, see the Q dependence below.)
- The base end is **not** immune to the base choice: it is evaluated at H₀, which changes by δ (an O(δ) shift, like the
  window dependence). The "base end immune" statements above (consequence paragraph, window section) apply only to
  the δ/Q part of ∂H.
- So the construction mismatch is at most part of the YIO end gap and none of MBT's. Open: whether the remaining YIO
  gap is Q-independent (pairs at 0.02-0.05 Å running) and what sets MBT's 3% (only Q = 0.01 Å exists for MBT Wannier).
  MBT dθ still falls fast with nk (48 → 21 → 12 rad/Å); nk 12-24 running for the self-test pair.

**YIO mode 1, end gap vs Q (nk 10, series windows):** gap = mode end − base end, fit c₀ + c₂Q².

| Q (Å) | atomic base 214162: base / mode end | gap | self-test base 214249: base / mode end | gap |
|---|---|---|---|---|
| 0.01 | −0.017045 / −0.017327 | −0.000282 | −0.017134 / −0.017234 | −0.000100 |
| 0.02 | −0.017011 / −0.017591 | −0.000580 | −0.017101 / −0.017500 | −0.000400 |
| 0.03 | −0.016957 / −0.018041 | −0.001084 | −0.017046 / −0.017952 | −0.000905 |
| 0.05 | −0.016742 / −0.019495 | −0.002753 | −0.016832 / −0.019411 | −0.002580 |

Atomic base: c₀ = −1.68e-4 rad/Å (1.0% of dθ), c₂ = −1.03 rad/Å³. Self-test base (two points): c₀ = −2.7e-7 (0.002%),
c₂ = −1.00. With the self-test base the two ends meet as Q → 0 and the gap is the ordinary c₂Q² with the same c₂:
the YIO Q-independent offset was the base construction mismatch. MBT has one Q, so its 3% gap (independent of the
base) can't be split into c₀ and c₂Q² yet.

All four Q (fit a + bQ², each end; gap c₀ + c₂Q², fit residuals ≤ 1.3e-5): self-test base: base end −0.017152,
mode end −0.017139 (Q → 0 limits agree to 0.08%), gap c₀ = +1.3e-5 (zero within the residuals), c₂ = −1.04.
Atomic base: base end −0.017063, mode end −0.017231, c₀ = −1.7e-4 (1.0%). The self-test base also moves the Q → 0
value itself by 0.5% (−0.017063 → −0.017152): an O(δ) model dependence of the same kind as the windows.


## MBT c2 integrand: a small-gap tube along k_z (2026-10-19)

MBT mode 3 (Q = 0.01 Å, mode SCF 1e-7, self-test base 214245) on uniform meshes: dθ base end 48.33 / 21.10 / 12.26 /
5.80 / 2.72 rad/Å at nk 4 / 6 / 8 / 12 / 16, falling ~1/nk². Not a code error: c2 = 1.09e3 at (½,½,½) and 5.1e2 at Γ on
every even mesh, against a BZ mean of ~3, so coarse meshes over-weight a narrow feature (`c2_peak_cut.py`: FWHM ≈ 0.05-0.075
reduced units across it; smooth; symmetric where symmetry requires).

Box refinement around Γ and (½,½,½) only (`../../archive/diagnostics/2026-10-19_mbt_refined_bz/`, cluster job 214345):
inside part converged (h/4 → h/8: 0.17%), total 4.357 / 4.450 rad/Å, with the outside on nk 16. (Pilot folder
`MnBi2Te4/dtheta_refined_2026-10-19/mode3_Q0p01_kfix_selftest/` deleted 10-19 as superseded by the tube run 214347:
100 tasks × ~41 min on 8 cores; dθ(m) base end 2.7217 / 4.4680 / 4.3647 / 4.3573 at m = 1/2/4/8 (mode end 2.8214 /
4.5563 / 4.4578 / 4.4502); outside the boxes 0.9855 / 1.0185 at nk 16; c2 at (½,½,½) 1.09e3, at Γ 5.1e2; min gap
in the boxes 0.192 / 0.188 eV.) A uniform nk 32 run
(214346) showed the outside is NOT converged (0.99 → 3.09 rad/Å; uniform nk 32 total 6.58 / 6.67). Cause: the direct gap
is ≈ 0.20 eV along the whole k_z line through Γ (reduced direction (1,1,1); MBT bands barely disperse along k_z,
c = 81.85 Å), and c2 ~ 1e3 all along it. The integrand lives on a tube, not two points.

Refinement in the plane across the tube instead (`make_klist_tube.py`: basis (b1, b2, b1+b2+b3), tube at u1 = u2 = 0;
32 k_z planes × in-plane nk 16 with an M = 8 subgrid in a 7 × 7-cell box; nested estimates in m, number of planes and box
size): job 214347.

**Tube result (job 214347, 106,976 k, `assemble_tube.py`).** dθ base end / mode end (rad/Å), box half-width 3 cells:

| in-plane level m (spacing h/m, h = 1/16) | 32 k_z planes | 16 | 8 | 4 |
|---|---|---|---|---|
| 1 (= uniform nk 16) | 2.7217 / 2.8214 | 2.7217 / 2.8214 | 2.7219 / 2.8216 | 2.8361 / 2.9374 |
| 2 (= uniform nk 32 in-plane) | 6.5812 / 6.6779 | 6.5812 / 6.6779 | 6.5824 / 6.6792 | 6.9582 / 7.0619 |
| 4 | 6.3733 / 6.4779 | 6.3733 / 6.4779 | 6.3743 / 6.4790 | 6.7248 / 6.8356 |
| 8 | 6.3707 / 6.4753 | 6.3707 / 6.4753 | 6.3717 / 6.4763 | 6.7222 / 6.8329 |

- k_z: 8 planes agree with 32 to 0.015%; 4 planes are 5.5% high. The integrand barely depends on k_z.
- In-plane: m 4 → 8 moves 0.04%.
- Box size (`outside_tube.py`: the outside region from the uniform nk 32 run 214346 at spacing h/2 instead of the
  coarse h): inside(m = 8, half-width H) + outside(nk 32) = 6.365 / 6.383 / 6.366 base end for H = 1 / 2 / 3. With the
  coarse outside, H = 2 is 1.0% low (6.3225) and H = 3 is 0.08% high (6.3707 vs 6.3656).
- **MBT mode 3, Q = 0.01 Å (1e-7 mode SCF, self-test base 214245): dθ = 6.37 / 6.47 rad/Å (base / mode end),
  converged to ~0.1%.** The end gap is 1.6% (3.1% at nk 8 was partly a mesh artifact). Min direct gap on the list
  0.192 / 0.188 eV.
- Production scheme for the rebuilt pairs: H = 3, m = 4, 8 planes = 8 × (49 × 16 + 207) = 7,928 k per pair (13× fewer
  than 214347), error budget ≲ 0.15%.

## YIO self-test base: all modes (2026-10-19 04:30)

`run_selftestbase.py`, self-test base 214249, nk 10, series windows. dθ (rad/Å), base end / mode end:

| Q (Å) | mode 1 | mode 2 | mode 3 |
|---|---|---|---|
| 0.01 | −0.017134 / −0.017234 | 0.054677 / 0.054297 | 0.063787 / 0.063777 |
| 0.02 | −0.017101 / −0.017500 | 0.054507 / 0.053066 | 0.063710 / 0.063668 |
| 0.03 | −0.017046 / −0.017952 | 0.054215 / 0.051262 | 0.063685 / 0.063588 |
| 0.05 | −0.016832 / −0.019411 | 0.053356 / 0.047690 | 0.063596 / 0.063298 |

Q → 0 limits (each end fit separately; gap = mode end − base end):

| mode | fit | base end | mode end | gap c₀ | gap c₂ (rad/Å³) | max residual |
|---|---|---|---|---|---|---|
| 1 | a + bQ² | −0.017152 | −0.017139 | +1.3e-5 (0.08%) | −1.04 | 1.3e-5 |
| 2 | a + bQ² + cQ⁴ | 0.054738 | 0.054737 | −9e-7 (0.002%) | −3.85 (Q⁴: 634) | 4.6e-7 |
| 3 | a + bQ² | 0.063763 | 0.063769 | +6e-6 (0.009%) | −0.12 | 5.0e-6 |

Mode 2's mode end is not quadratic up to 0.05 Å (a + bQ² leaves a 5e-4 residual, c₀ = −5.4e-4); with Q⁴ (one degree of
freedom left) the ends meet. With the self-test base the two ends agree at Q → 0 for all three modes; the remaining
gap is the ordinary O(Q²) curvature of dθ(Q). Window dependence: window study (18 pairs) running.

## YIO window study: result (2026-10-19 13:36)

`archive/diagnostics/2026-10-19_yio_window_study/run_window_study.py`, all 18 pairs. Q = 0.01 Å, conv_thr 1e-10
(base 213998 NSCF 214124, modes 1e-10), nk 10, each window set with its own Q = 0 self-test base, production settings
otherwise. dθ (rad/Å), base end / mode end, and the mean relative to the series windows (12.422 / 22.644):

| froz max / win max (eV) | Ω_I (Å²) | mode 1 | mode 2 | mode 3 |
|---|---|---|---|---|
| 10.786 / 22.644 | 116.84 | −0.017053 / −0.017056 (−0.75%) | 0.054894 / 0.054349 (+0.25%) | 0.063440 / 0.063405 (−0.56%) |
| **12.422 / 22.644** (series) | 117.10 | −0.017134 / −0.017234 | 0.054677 / 0.054297 | 0.063787 / 0.063777 |
| 13.120 / 22.644 | 117.37 | −0.017093 / −0.017193 (−0.24%) | 0.054688 / 0.054309 (+0.02%) | 0.063825 / 0.063815 (+0.06%) |
| 14.396 / 22.644 | 158.55 | −0.018144 / −0.018172 (+5.67%) | 0.054109 / 0.053696 (−1.07%) | 0.063982 / 0.063994 (+0.32%) |
| 12.422 / 19.993 | 119.47 | −0.017215 / −0.017304 (+0.44%) | 0.054742 / 0.054340 (+0.10%) | 0.063679 / 0.063671 (−0.17%) |
| 12.422 / 24.522 | 117.34 | −0.017110 / −0.017209 (−0.14%) | 0.054675 / 0.054297 (−0.00%) | 0.063770 / 0.063758 (−0.03%) |
| 14.396 / 24.522 (production) | 159.25 | −0.018126 / −0.018153 (+5.56%) | 0.054030 / 0.053633 (−1.20%) | 0.063990 / 0.064003 (+0.34%) |

Ω_I from each window's atomic-trial base Wannier run (`window_study_2026-10-19/<set>/base/2*/`; production and series
from `ref_prodwin/`, `ref_ampwin/`).

- **Plateau:** frozen top 10.786-13.120 eV and outer top 19.993-24.522 eV all agree with the series windows to
  ≤ 0.75% (mode 1), ≤ 0.25% (mode 2), ≤ 0.56% (mode 3), with Ω_I 116.8-119.5 Å².
- **The outer window doesn't matter** (19.993 → 24.522 at frozen 12.422: ≤ 0.6%; at frozen 14.396, 22.644 → 24.522:
  ≤ 0.2%).
- **Frozen top 14.396 eV is the outlier**, and it alone causes the production/series difference: Ω_I jumps 117 → 159 Å²
  (+35%) and mode 1 moves +5.7%. Freezing up to 14.396 forces states into the 176-WF set (Y p, Ir d, O s p) that it
  does not describe well; the gauge-invariant spread says the model changed, not the physics.
- **Use the series windows (12.422 / 22.644)** for YIO dθ; the window uncertainty on the plateau is ≲ 0.75%. The
  self-test Q → 0 values above (−0.017145 / 0.054737 / 0.063766 rad/Å) are on the plateau.

## YIO θ(Q) from the Chern–Simons 3-form (2026-10-20 03:00)

`archive/diagnostics/2026-10-19_yio_cs_theta_series/` (`cs_series.py` on the cluster, array job 214446, 24 tasks
= 12 pairs × 2 k slices on n126/n246/n247; first try 214445 ran out of memory at 30 GB per task, the A(R) build peaks
above that; `assemble.py`, output `assemble.log`, data `data/Y2Ir2O7/phonon/GM2-/cs_theta_series_214446/`). Same pairs
as the c2 series: base = self-test 214249, modes = series-window rebased runs 214147-214158, conv_thr 1e-10. Method as
log 10-14 §2: one SCDM trial set from the base (σ_min 0.541), full connection, θ(Q) = θ_CS(mode) − θ_CS(base).
θ_CS(base) = −1.2e-9 (nk 10), +9e-11 (nk 12): zero. **nk 10 = nk 12 to 1e-6 relative in every θ(Q).**

θ(Q) (rad, nk 12) and θ(Q)/Q (rad/Å):

| Q (Å) | mode 1 | mode 2 | mode 3 |
|---|---|---|---|
| 0.01 | −1.7151e-4 (−0.017151) | 5.4624e-4 (0.054624) | 6.3754e-4 (0.063754) |
| 0.02 | −3.4335e-4 (−0.017168) | 1.08634e-3 (0.054317) | 1.27154e-3 (0.063577) |
| 0.03 | −5.1588e-4 (−0.017196) | 1.61589e-3 (0.053863) | 1.90158e-3 (0.063386) |
| 0.05 | −8.6202e-4 (−0.017240) | 2.65206e-3 (0.053041) | 3.13788e-3 (0.062758) |

**The comparison with c2 is the Q → 0 intercept.** θ(Q)/Q is the secant, dθ/dQ|₀ + O(Q²); fit a + bQ² + cQ⁴ over
the four Q (a + bQ² leaves 1.5e-4 for mode 2, which is not quadratic to 0.05 Å, as its c2 mode end):

| mode | a, 3-form (rad/Å) | c2 Q → 0 limit | difference | θ/Q at 0.05 Å vs a |
|---|---|---|---|---|
| 1 | −0.017144 | −0.017145 | −0.006% | +0.6% |
| 2 | 0.054737 | 0.054737 | 0.000% | −3.1% |
| 3 | 0.063791 | 0.063766 | +0.04% (fit scatter 2e-5) | −1.6% |

- The two routes give the same dθ/dQ to ≤ 0.04%. They share H(R) and A(R), so this checks the c2 code, the gauge
  (A_β = 0 is exact here: projection gauge, Λ = 0 at O(λ)) and the k sums, not the Wannier model.
- **θ is linear in Q to 0.6% (mode 1), 3.1% (mode 2), 1.6% (mode 3) at 0.05 Å.** For the amplitude figure plot θ_CS(Q)
  with the line Q·dθ/dQ|₀.
- **Correction:** the c2 two-end mean at finite Q is not θ(Q)/Q; it differs by O(Q²) (mode 1: −0.19 / −0.77 / −1.7 /
  −4.9% at 0.01-0.05 Å, scaling as Q²). Its +5.5% / −7.3% / −0.5% change over 0.01-0.05 Å (quoted on 10-19 as θ's
  nonlinearity) is wrong; the 3-form numbers above replace it. Comparing the two at finite Q, I also first blamed the
  frozen gauge; that was wrong (A_β = 0 is exact at O(λ)).
- Cost: load 12 min per pair (two 2.7 GB .chk, A(R) on 13,824 R; peak > 30 GB), 2.2 s per k point per structure on
  one EPYC core (the gauge construction is a per-k Python loop, so BLAS threads don't help).

## Band figures: k conventions and the two interpolants (2026-10-20)

DFT bands of both 1e-10 bases on each base's `.win` kpoint_path (`archive/setup_templates/MnBi2Te4/base/
setup_dft_bands.sh`): MBT 214640 (104 k, SCF 214290), YIO 214641 (208 k, SCF 213998). `analysis/bands.ipynb`
rewritten for the overlay.
- **Coordinates.** The path is the bands.in `K_POINTS crystal_b`, expanded by `read_kpath` exactly as pw.x does:
  reduced coordinates of the `CELL_PARAMETERS` basis. PythTB `solve_ham` takes reduced k of the `.win`
  `unit_cell_cart` basis; the two cells are identical row for row for both materials (checked in the notebook).
  pw.x prints bands.out k in Cartesian 2π/alat; the notebook converts back (k·a_i / alat) and checks point for point.
  x axis = cumulative |Δk| in Å⁻¹ with 2π (Wannier90's `_band.dat` convention; YIO Γ–X = 0.6166 Å⁻¹ = 2π/10.19 Å);
  Wannier90's node positions (`_band.labelinfo.dat`) equal ours to 1e-6 Å⁻¹.
- **Our models vs Wannier90's own `_band.dat`, at its k points exactly** (`archive/diagnostics/2026-10-20_band_interpolants/w90_band_check.py`):

  | run | `_hr.dat` model | production `.chk` model |
  |---|---|---|
  | YIO self-test 214249 | 0.09 meV | 0.16 meV (E_F ± 3 eV), 0.6 meV all bands |
  | MBT self-test 214245 (1e-7 base), pinned | 0.07 meV | **16.8 meV** |

  MBT: identical at the 12 path points on the 8³ mesh (5e-7 meV); off-mesh up to 16 meV on T–P1 and F–Γ (bands 2-3 eV
  below E_F), 11 meV on Γ–L (at −0.67 eV), 0.1 meV on Γ–T. Two interpolants of the same mesh data: the atom-pinned
  fold (production, C3-exact, bug #21) vs Wannier90's fold on the Wannier centres. The DFT bands 214640 decide which is
  closer to DFT between mesh points. YIO has no pinning, and the two agree.
- **YIO DFT bands 214641 (done 10-20 12:53, 19 min on n126, pw.x db4b99b4).** k match the path to 7e-5 (reduced). DFT
  bands inside each run's frozen window vs the nearest Wannier band (production `.chk` model), 208 path k:
  self-test 214249 (series windows, −14 to 12.422 eV): **max 2.3 meV, median 0.22 meV**; 214125 (old production
  windows, −14 to 14.396 eV): max 235 meV, median 15 meV, the same window as the dθ outlier of the window study.
  Figure `data/Y2Ir2O7/base/lda/dft_bands_214641/bands_dft_wannier_214249_214125.png`.
- **YIO conduction bands above the frozen window (user, 10-20: "poorly aligned").** Of the 20 empty Wannier bands, the
  4 Ir j=1/2 bands (157-160, 0.35-1.0 eV above E_F) are inside the frozen window: |W − DFT| ≤ 0.10 meV (214249). The
  other 16 (161-176, Ir e_g-like) lie above it (frozen top E_F + 1.65 eV), where the DFT manifold at 2.3-8 eV mixes Ir
  e_g with Y 4d, which the set leaves out (Y p only; log 09 "faithful set"). With projection-only disentanglement
  (dis_num_iter = 0) these are the projected Hamiltonian of the e_g trials: band by band vs DFT +466 meV on average,
  up to 1.4 eV; the lowest at Γ is 2.85 vs DFT 2.32 eV, even at mesh k. 214125 (frozen top E_F + 3.62 eV) reproduces
  the 2.3-3.4 eV bands but interpolates them with wiggles up to 235 meV near its frozen edge (Ω_I 117 → 159 Å²).
  dθ moves ≤ 0.75% across frozen tops 10.8-13.1 eV (window study), so the e_g bands' placement is not in the result
  at that level.
- `analysis/bands.ipynb` restructured (10-20, user): per material an inline compute cell (settings at the top; saves
  each run's eigenvalues, windows (`W90.dis_window`, `dis_froz_window`), mesh (`.chk` `kpt_red`) and `.eig` in
  `<BANDS_IN folder>/wannier_eigs/`, read back while the path, `.chk` and `_hr.dat` are unchanged) and an inline plot
  cell (plot settings line by line; only draws). Whole notebook 2 s from the saved eigenvalues.
- **PythTB `external`: hopping table vectorised (10-20, user; uncommitted in the worktree).** `W90.model` already
  collected the hoppings as arrays, but `TBModel._append_hops` converted them one by one (`_val_to_block`, `.tolist()`)
  and `HoppingTable.extend` did again (`_coerce_amplitude`) and built a `(i, j, R)` → row dict with one Python tuple
  per hopping (17 M for YIO). Now the amplitudes go in as arrays and the dict is built on the first lookup (`find`,
  `update`, `set_hop`). `W90.model()` YIO 214249: ~59 s → 0.4 s; production load 232 → 57 s (left: two `_hr.dat`
  reads 2 × 13.6 s, already one `np.loadtxt`; our WS fold 24 s). Eigenvalues bit-identical (YIO 214249, MBT 214245);
  pythtb tests 152 passed.

## Band figures: MBT (2026-10-21 00:30)

DFT bands 214640 (104 k on the `.win` kpoint_path, SCF 214290 at conv_thr 1e-10, pw.x db4b99b4, 8h52m) against the
rebuild's Q = 0 self-test 214301 (the dθ base; production model: H(R) from `.chk`, atom-pinned WS fold, bug #21), in
`analysis/bands.ipynb` (eigenvalues cached in `dft_bands_214640/wannier_eigs/`; 27 s with the PythTB
R-block storage). Base 214300 (atomic trials) and self-test 214301: Ω_I 325.53 / 344.46 Å², Ω_total 432.71 / 453.65 Å²
(Release 4.0.3). Window check of NSCF 214290: dis_froz_max 10.278 is 8.7 meV from the nearest eigenvalue, OK.

DFT band in the frozen window (−7 to 10.278 eV) minus the nearest Wannier band, by E − E_F (E_F 7.4058 eV):

| E − E_F (eV) | max (meV) | median (meV) |
|---|---|---|
| < −3 | 10.1 | 0.26 |
| −3 to −1 | 10.6 | 0.80 |
| −1 to 0 | 13.4 | 0.76 |
| 0 to 1 | 10.3 | 1.08 |
| 1 to 2 | 23.2 | 2.25 |
| 2 to 2.87 | 249.1 | 3.17 |

At the 12 path points on the 8×8×8 mesh the maximum is 0.05 meV (exact, as the frozen window requires). The large
errors are all above E_F + 2.75 eV, on bands within ~0.1 eV of the frozen top, off the mesh (worst 249 meV at
(0.232, 0, 0.232) on F–Γ, 172 meV at (0.232, 0, 0) on Γ–L). The frozen window only pins the Wannier bands at the mesh
points; off the mesh, bands just under the frozen top are presumably mixed with the disentangled states above it (not
checked further). Smallest gap on the path: DFT 194.8 meV, Wannier 191.8 meV, both at k 11 (next to T on T–P1).

## Does the conduction-band misfit matter? Occupied-band errors (2026-10-21 01:00)

Question (user): YIO's e_g bands above the frozen window (+0.47 eV) and MBT's bands under its frozen top (up to 0.25-0.3
eV) are off DFT; does that make dθ wrong?

- **Theory.** ∂θ/∂λ (Eq. 4curv) is the second Chern form of the occupied bundle: a functional of the occupied projector
  only. The unoccupied part of the Wannier space only completes the basis; it decides the internal/cross/external
  split, not the sum (paper: "only their sum is independent of the size and composition of the Wannier space"),
  provided the Wannier space contains the occupied states at every k and the gap stays open. So the misfit itself does
  not enter. What enters is how well the occupied states are interpolated between mesh points; the band errors near
  the gap measure that.
- **YIO evidence.** Window plateau (frozen top VBM + 0.34 to 2.68 eV: ≤ 0.75%; outer top VBM + 9.5 to 14.1 eV: ≤ 0.6%).
  The 14.396 outlier (VBM + 3.95) is the case where the occupied interpolation got worse: occupied bands vs DFT on the
  214641 path, index-matched: 214249 (series) max 2.29 meV, median 0.033, top 2 eV of valence 0.23 meV; 214125
  (14.396) max 6.50 meV, median 0.21, top 2 eV 4.50 meV; Ω_I 117 → 159; dθ mode 1 +5.7%. Bands bordering the gap in
  214249: 0.06 meV (top valence), 0.09 meV (bottom conduction).
- **MBT 214301 (pinned, .chk), index-matched on the 214640 path:** occupied max 13.4 meV (median 0.38); conduction up to
  VBM + 2 eV 17.6 meV; VBM + 2-2.5 eV 38 meV; VBM + 2.5-2.8 eV 170 meV; VBM + 2.8 eV to the frozen top 319 meV. Top
  valence band max 11.5 meV (on Z–P1), bottom conduction 10.0 meV. By segment (4 top valence + 4 bottom conduction):
  Γ–Z 0.05 meV (exact along the k_z line), Z–P1 11.5, P1–F 1.0, F–Γ 10.3, Γ–L 10.5 meV. **Near Γ off the Γ–Z line**
  (|k| ≤ 0.09 reduced, inside the dθ tube): near-gap errors up to 10.5 meV peaking midway between mesh points (0.054
  between the mesh points 0 and 0.125), the direct gap 3-4 meV too small (226.6 vs 230.7 meV at (0.018,0,0.018); 305.4
  vs 309.1 at 0.071). DFT VBM is off the mesh, at (0.036, 0, 0) on Γ–L, 7.3032 eV (mesh VBM 7.2681); Wannier VBM 10 meV
  higher. Gap: 163 meV indirect on the path (CBM 7.4663 at Z), 198 meV on the mesh (direct, at Z); Wannier 153 meV.
- **Verdict.** YIO: the misfit is above the frozen window and the window plateau shows dθ insensitive to it; occupied
  and near-gap bands are good to 2.3 / 0.1 meV. Fine. MBT: the large misfit is at VBM + 2.5-3 eV and does not enter by
  itself, but the near-gap interpolation error near Γ (10 meV, gap −2%) sits where the c2 integrand is concentrated, and
  is ~3x the occupied error of the YIO outlier on a 5x smaller gap. It comes from the 8×8×8 mesh (interpolation midway
  between mesh points), not from the windows. Its effect on dθ is NOT measured. Tests: (1) MBT window test (frozen top
  lowered into another gap; reuse .amn/.mmn/.eig, rerun wannier90 + self-test/mode .amn), (2) dθ vs Wannier mesh
  (a denser NSCF, e.g. 10×10×10; ~0.5 h per k per pool on dko64m4), the direct test of the interpolation error.

## Paper band figure (2026-10-21 01:40)

`analysis/bands.ipynb` last cell "Paper figure": (a) YIO 214249, (b) MBT 214301, DFT solid / Wannier dashed, zero at the
DFT VBM on the path, shared −3 to 2.5 eV, dotted frozen-window top (YIO + 1.98 eV; MBT + 2.98 eV off scale), MBT Γ–T
labelled Γ–Z. Writes `paper/figures/bands.pdf` (+ png). `paper/main2.tex`: new Sec. IV E "Accuracy of the interpolated
bands" + Fig. 3 (the paragraph above in paper form, MBT test as \Red TBD); windows paragraph and Table (details) now the
series windows for YIO (12.422 / 22.644; 14.396 superseded 10-19); MBT gap 163 meV indirect (198 on the mesh) in the
structure paragraph and the proximity paragraph; Wannier90 3.1.0 → 4.0.3; placeholder figure numbers +1; MBT Fig.
placeholder (d) bands removed (now Fig. 3b). Backup of the previous main2.tex in the session scratchpad only.

## Paper: YIO production numbers, band physics, MBT integrand (2026-10-21 02:30)

- **YIO couplings at production settings** (`archive/diagnostics/2026-10-21_yio_mode_couplings/mode_couplings.py`,
  log next to it): Φ from the 1e-10 SCFs at Q = 0.01 Å (asymmetry 0.03%), ω 228.6 / 328.3 / 471.6 cm⁻¹, characters
  Y/Ir/O 0.04/0.67/0.29, 0.80/0.02/0.19, 0.01/0.00/0.99; with dθ/dQ = −0.017145 / 0.054737 / 0.063766 rad/Å (self-test,
  series windows, Q → 0): dθ/dq −0.00204 / +0.00118 / −0.00415 rad Å⁻¹ amu⁻½, per Å −0.024 / +0.010 / −0.017,
  δθ_zp 5.5e-4 / 2.7e-4 / 7.8e-4 rad; axion charges z −0.00247 / +0.00790 / +0.00531 rad/Å. Check: the old dθ/dQ
  with these eigenvectors gives the old −0.00205. All in main2.tex (tables V, VI, appendix D, abstract, summary).
- **YIO band edges** (`2026-10-21_yio_band_edges/band_edges.py`, self-test 214249): VBM at L (0,½,0), CBM at Γ, gap
  684.9 meV indirect; smallest direct 954.9 meV at (0.6,0,0) on the Γ–L line. Edges: Ir d 0.47-0.54 / O p 0.46-0.52
  (valence top), Ir d 0.70-0.77 (conduction bottom) at Γ, X, L, W: the j_eff = 1/2 manifold, same composition
  everywhere (no character exchange across the gap).
- **MBT band inversion** (`2026-10-21_mbt_band_inversion/band_edges.py`, self-test 214301, pinned; inversion check
  4.5e-7 eV): at Γ the top valence doublet (7.2548 eV) is even, p_z 0.72, Bi 0.37 / Te1 (outer) 0.53 / Te2 0.10; the
  bottom conduction doublet (7.4730) is odd, Bi 0.43. At F (gap 1.04 eV) and L (0.98 eV) the valence top is Te2-heavy
  and the conduction bottom is even Bi p_z (Bi 0.56-0.59): the atomic order. Occupied parities 48/34 at Γ (atomic
  "ionic" set 46/36, log 10-14 §1), 46/36 at F, 41/41 at Z and L. At Z every doublet holds one even and one odd state.
- **MBT integrand paragraph** in main2 Sec. VII C uses the mode-3 numbers of "MBT c2 integrand: a small-gap tube" (1e-7
  mode SCF, labelled "a first calculation"); final values come from the 1e-10 chain. The per-k c2 of every chain pair
  is saved (`chunks/*.part.npz` + `klist.npz`; README) for Fig. 7c.
- **Figure:** per-panel energy ranges (YIO −23.5 to 2.5, MBT −14 to 2.5 eV) showing every occupied Wannier band; the
  user had restyled the cell (dict literals, comments) before this edit, which was applied on top of that version.

## c₂ in Cartesian k-space, both materials (2026-10-21 03:00, draft figure)

`analysis/figures/c2_distribution.ipynb` (draft, `paper/figures/c2_distribution.pdf`). Data: MBT mode 3, Q = 0.01 Å,
1e-7 pair, tube run 214347 (106,976 k; chunks + klist copied to
`data/MnBi2Te4/archive/phonon/T2-/mode3/c2_distribution_mode3_1e-7_214347/`, also the uniform nk 32 run); YIO production
pairs, `c2_density` on the 10³ mesh.
- **MBT:** Σw c₂ = 6.4230 rad/Å (= mean of 6.3707 / 6.4753, check). Fraction of ∂θ/∂Q within k⊥ of the k_z axis
  through Γ: 47.7% within 0.05 Å⁻¹ (0.34% of the zone), 84.1% within 0.1 (1.29%), 97.1% within 0.2 (5.2%); in-plane
  zone width 1.67 Å⁻¹. On the k_z ≈ 0 plane: a hexagonally warped core with c₂ > 10² inside the 0.25 eV direct-gap
  contour, and six negative lobes (|c₂| ~ 10-10²) at k⊥ ≈ 0.1 Å⁻¹ along the six in-plane directions; negative parts
  cancel 7% of the total. log|c₂| vs log gap correlation −0.79 (nk 32).
- **YIO:** no hot spot; c₂ ~ 0.1-0.5 spread over the zone with sign changes. Ir pattern on the (110) plane (contains Γ,
  X, L, K, U): negative cross centred on Γ along Γ–L, positive near X/U/K; its negative regions cancel 79% of the
  positive (Y: mostly negative; O: −5%). Half of Σ|c₂| comes from ~30% of the zone (MBT: 84% of the integral from
  1.3%).
- Tr[𝔅𝔅] is pointwise gauge invariant (unlike the 3-form density), so these maps are physical for the given interpolant.
- To finish: MBT from the 1e-10 chain chunks (all four modes); YIO dense (110) plane per pattern (cluster, the pair
  load needs > 18 GB); decide per pattern vs physical mode (c₂ is linear in ∂_λ: mode map = Σ cᵢ c₂ᵢ) and placement
  (replaces Figs. 6d and 7c?).

## Correction: c₂ along k_z, and the tube-mesh geometry (2026-10-21 03:30)

- **c₂ is not constant along k_z.** On the axis (u₁ = u₂ = 0, all 32 planes of 214347) it runs from 519 (Γ) to
  1117 (Z). "The integrand barely depends on k_z" (section "MBT c2 integrand") meant only that 8 planes integrate it to
  0.015%: the k_z dependence is smooth (periodic trapezoid rule), not small. main2 Sec. VII C now says "varies smoothly,
  by about a factor of two between Γ and Z".
- **Geometry of `make_klist_tube.py`:** in-plane parts of b₁, b₂ at 120°, |b⊥| = 1.674 Å⁻¹; coarse spacing h|b⊥| =
  0.105 Å⁻¹, fine (M = 8) 0.013 Å⁻¹. Refined region = 7×7 cells, a 120° rhombus in Cartesian k, half-side 0.366 Å⁻¹,
  inradius 0.317 Å⁻¹ (not a cylinder). The constant-u₃ planes are tilted (b₁, b₂ have k_z component 0.077 Å⁻¹ each):
  ±0.03 Å⁻¹ across the refined box. Integration cell = primitive cell in (b₁, b₂, b₁+b₂+b₃), not the Wigner–Seitz zone;
  max distance from the axis after folding 0.912 Å⁻¹ (hexagon edge 0.837, corner 0.967).
- **Azimuthal structure:** c₂ in rings on the k_z ≈ 0 plane: 793-815 at k⊥ 0.03-0.04, −4 to 166 at 0.07-0.09 (the six
  negative lobes), 16-83 at 0.15-0.17 Å⁻¹. Cumulative fraction: 1.0077 at 0.3, 1.0072 at 0.6 Å⁻¹; k beyond 0.6 Å⁻¹ = 53%
  of the zone, −0.7% of ∂θ/∂Q.
- `c2_distribution.ipynb` panel (b) now shows k⊥ shells (weighted mean, range) out to 0.91 Å⁻¹; the first version's
  per-point scatter formed columns because each in-plane point recurs in all 32 planes at the same k⊥ (403 distinct
  values, ~192 points each).

## MBT tube convergence at 1e-10 and the faster code (2026-10-21 22:10)

- **Why:** the 1e-10 production pairs moved more than mode 3 at 1e-7 between the nested estimates: m 2 → 4 (8 planes)
  8.7% (mode 1) / 5.8% (mode 2); 4 → 8 planes (m 4) 5.5% / 8.4%; box half-width 2 → 3 0.6% / 0.9%. The m 2 → 4 change
  sits in the ring of 8 cells around the axis cell (mode 1 +0.213 of +0.173 rad/Å; mode 2 −0.259 of −0.221), where the
  sign-changing lobes at k⊥ ≈ 0.1 Å⁻¹ are; rings 2-3 move < 0.6%.
- **Why 16 planes:** the k_z sum is a periodic trapezoid rule (error falls exponentially with the number of planes), and
  the 16-plane list contains the 8 production planes, so 8 vs 16 measures the production error directly while
  reusing every production point (mode 3 at 1e-7: 8 vs 32 planes 0.015%).
- **Cost:** the lists nest exactly (M 4 subgrid = even points of M 8; coordinates multiples of 1/256), so only the new
  45,560 of 53,488 k per pair are computed (`new_points.py`), and `assemble_convergence.py` merges them with the
  production chunks (checked locally: the (m 4, 8 planes) entry reproduces mode 1's production table exactly).
- **Faster code:** with pythtb 9db1288 (cluster code 2026-10-19), 89% of each 50-k batch is `np.unique` re-binning the
  hopping table in every `hamiltonian()`/`velocity()` call (19.9 of 22.4 s, laptop, independent of batch size).
  pythtb 5776792 (per-R blocks) does the batch in 1.93 s (11.5×); model load 143 → 112 s (dominated by the .chk). Same
  models: H to 7e-12 eV, dH to 2e-10, E identical, Tr[ΩΩ] of the occupied curvature to 5e-7 relative (the occupied
  blocks themselves differ by ~1e-2: Kramers-pair gauge). Cluster code `code/2026-10-21`; validated on the cluster
  against production c2 before the convergence tasks run (`validate_code.py`, VALIDATED marker).
- **SGE:** a held job is released when its predecessor exits with any status; exit 100 did not keep the mode 3 Q 0.03
  chain held (it ran on no model and wrote NaN). Guards must sit in the dependent jobs themselves.
- **Result (22:40):** validation 214981: 100 production c2 points of mode 1 recomputed with code/2026-10-21 to 5.8e-7
  relative (max |c2| 983). Convergence tasks: 45,560 k in 386-1055 s per task of 5,695 k, load 190-430 s (production
  7,928 k took 870-1685 s per task of 793 k). dθ base / mode end (rad/Å), box half-width 3:

  | | m 4, 8 planes (production) | m 4, 16 | m 8, 8 | m 8, 16 |
  |---|---|---|---|---|
  | mode 1 | −1.98812 / −1.98312 | −1.98776 / −1.98277 | −1.98618 / −1.98119 | −1.98583 / −1.98083 |
  | mode 2 | 3.84110 / 3.89289 | 3.84032 / 3.89209 | 3.83837 / 3.89010 | 3.83760 / 3.88930 |

  Planes and in-plane level are converged (≤ 0.12% together). Box half-width 2 → 3 at m 8, 16 planes: −1.99706 →
  −1.98583 (0.57%), 3.80674 → 3.83760 (0.80%); the next ring (half-width 4) is not computed, so the box error is
  probably a few tenths of a percent, not measured.
- **Box size (23:20, `box_H5`: half-width 5, m 4, 8 planes; 8,640 new k per pair):** base end, rad/Å, half-width
  2 / 3 / 4 / 5: mode 1 −1.99953 / −1.98812 / −1.98765 / −1.98701; mode 2 3.80747 / 3.84110 / 3.83800 / 3.83671; mode 3
  6.38442 / 6.43441 / 6.43143 / 6.43119. The 2 → 3 step (0.6-0.8%) was the ring-3 cells' own refinement, not a slowly
  converging tail: 3 → 5 is 0.05-0.11%. Summary: production (half-width 3, m 4, 8 planes) is within 0.2% of the
  converged tube integral for modes 1-3 (box and m/planes errors have the same sign: mode 1 +0.17%, mode 2 +0.20%,
  mode 3 +0.11% relative to the estimate with all three refinements added).

- **All pairs (10-22 ~01:00).** Converged (box H5 + m8/16-plane correction at H3; base / mode end, rad/Å): mode 1
  −1.9847 / −1.9797, mode 2 3.8332 / 3.8849, mode 3 6.4274 / 6.5956, mode 4 −4.5426 / −4.5306 (Q 0.01), mode 3 Q 0.03
  6.4058 / 8.0205; all 0.11-0.21% below the H3/m4/8-plane values. Mode 3's end gap 0.168 → 1.615 rad/Å for Q 0.01 →
  0.03 (×9.6, ∝ Q²); Q² extrapolation 6.430 (base end) / 6.418 (mode end). Three assembles failed on OpenBLAS memory
  even with one thread in two of them; rerun on the login node.
- **Mode 3 end difference (10-22).** Binned by base direct gap: k with gap < 220 / 220-260 meV carry 45 / 36% of dθ and
  61-63 / 33% of the end difference at both Q; their gap falls 3.2 / 2.6 meV (Q 0.01) and 23.7 / 19.7 meV (Q 0.03).
  Inversion (T₂⁻ odd) makes θ − π odd in Q, so gaps and dθ/dQ are even: the change is Q². The base end is the
  Q → 0 derivative (6.427 → 6.406); the mode end is the integrand on the displaced structure.
- **Which estimate (10-22, user asked).** For an axion-active pattern θ − θ₀ is odd in Q, so both estimates are
  O(Q²) estimates of ∂θ/∂Q|₀ (end e = F(Q_e)[secant]; F(0)[∂²H] = 0 because ∂²H is fully symmetric; θ''(0) = 0 then
  forces F'(0)[∂H] = 0). Neither end is "the physical one": the Q → 0 limit is. YIO: ends agree to ≤ 0.08% at Q → 0,
  so the mean was harmless. MBT: the undistorted-end Q² coefficient is ~65× smaller than the displaced end's (pattern 3),
  so the paper uses the undistorted end at Q = 0.01 Å (main2 methods paragraph rewritten; old text claimed O(Q)).
