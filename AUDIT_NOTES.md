# CMSTG audit — working notes

Running log of the independent audit, bug fixes and reruns of the CMSTG simulation suite.
Started 2026-10-22. Auditor: Claude (Opus 5.5), Claude Code session.

- Working copy: `/home/aion/CMSTG-fixes` (fresh clone of `cisomorph/Curvature-Memory_Scalar-Tensor_Gravity`, branch `audit-fixes`, base commit `cd03006`).
- Committed on branch `audit-fixes` and pushed to `origin` on 2026-10-22 at the user's request (commit `7b5f22c`). The SIM124 results were added in a follow-up commit, also pushed at the user's request. `main` is untouched and no PR has been opened.
- Not available on this machine: CLASS, classy, CAMB, the Planck plikHM likelihood and clipy. Sims that need them (SIM88, SIM97, SIM98) can't be rerun here.
- Available: Python 3.10, NumPy 2.2.6, SciPy 1.15.3, emcee 3.1.6, matplotlib.

---

## 0. Executive summary (as of the end of the rerun)

### What was wrong (shared across the suite)
- **Data (B5):** 22 sims compared models to a hardcoded "DESI H(z)" table that is not DESI data. There are two variants of it. They are now replaced by the official DESI DR1 (2024) and DR2 (2025) Gaussian BAO likelihoods.
- **Statistics (B2):** "tension" = √(χ²/N) is not a significance. It is replaced by χ², p-value, GOF σ, and Δχ² against ΛCDM run through the same pipeline.
- **Sound horizon / θ* (B3, B4):**
  - r_s was truncated at z = 1e4 or 1e5, hardcoded at 144.7 Mpc, or integrated over 0 → z*.
  - θ* used r_s(z_drag) instead of r_s(z*), and θ_MC in place of θ*.
- **Friedmann equation (B9):**
  - Phase 3 / Phase 2 family: `H²(3F + 6F'Ψ' − ½Ψ'²) = ω…` has no GR limit (H² = ω/1.5, i.e. Ω_m,eff = 0.21).
  - Every family: F'Ψ' and the kinetic term were 2× too large.
  - Model-specific errors: Gauss–Bonnet term omitted (SIM133); Horndeski energy 3× too small and field equation unchanged (SIM136); G3 term entered linearly and the 1/K was missing (SIM142); U force/energy factor 2 (SIM143); matter-coupling sign (SIM144); slow roll with R = 12H² (SIM123); no F'R source (SIM124).
- **ε_H (B10):** taken at fixed field values, dropping the field's own evolution.
- **Fixed H0:** the Phase 3/4 sims held H0 = 67.4 or 67.59 and r_s fixed, then treated θ* as a pass/fail test. Any model that changes H(z) "fails θ*" that way. This is the basis of Paper I's Theorem 2.
- **Solar System (B1):** SIM107 never evaluated the canonical Ψ̄ = 2.62 and used a wrong coupling formula, a hand-chosen mass sign, and two unit errors.

### Before → after (DR1 Δχ² is against ΛCDM in the same pipeline; > 0 means better than ΛCDM)

| sim | model | published | audited (faithful / θ*-matched) | key numbers (θ*-matched) | local gravity |
|---|---|---|---|---|---|
| SIM107 | Solar System | PASS | **FAIL** | \|γ−1\| = 9.5e-4 at Ψ̄ = 2.62 | 41× Cassini, no screening |
| SIM121C | frozen F0 + CPL, MCMC | PARTIAL ("2.63 / 2.77σ") | PARTIAL | F0 adds Δχ² = 0.00 (exact degeneracy); w0wa gives +10.5 DR1 / +23.8 DR2 | — |
| SIM122 | Λ0 + CPL | FAIL | PARTIAL | Λ0 unconstrained (profile spread 1.5e-8) | — |
| SIM121B | H0 tension | FAIL | FAIL | CMSTG lowers the true H0 (65.9–66.6), so the tension worsens to 5.6–6.2σ | — |
| SIM130 | **Phase 1 unlocked** | FAIL | FAIL / **PASS** | +5.4 / +13.9 at Ψ_ini = 2.62; +7.9 / +19.7 at Ψ_ini = 3.5 | **49–85× Cassini** |
| SIM131 | ξΨR | FAIL | FAIL / FAIL | best −49 | — |
| SIM132 | subtractive F (attractor) | FAIL | PARTIAL / **PASS** | +9.3 / +21.2 (Λ0 = 0.2; robust ±0.3) | **Cassini and LLR pass** (F' ≈ 0 today); BBN ΔG ≈ +7% (marginal) |
| SIM133 | Gauss–Bonnet | FAIL | FAIL (no solution) | constraint breaks down | — |
| SIM134 | exponential F (attractor) | FAIL | FAIL / **PASS** | +8.9 / +19.2 (Λ0 = 0.2; robust) | **Cassini and LLR pass**; BBN ΔG ≈ +7% |
| SIM135 | bi-scalar ξφR | FAIL | FAIL / PASS | +9.0 / +19.6 (ξ = 0.03) | 145× Cassini |
| SIM136 | Horndeski G^{μν}∂Ψ∂Ψ | FAIL | **FAIL (c_T)** | +4.6 at best | c_T² − 1 ≥ 2.4e-11 ≫ GW170817 |
| SIM123 | hilltop DE, no Λ_b | FAIL | PASS / PASS | +9.3 / +21.5 | 162× Cassini |
| SIM124 | K(Ψ) = e^{−αΨ²} (hilltop) | FAIL | **FAIL** (fine-tuned: needs Ψ_ini ≈ 2.25–2.40) | the two solutions that reach Ψ0 = 2.62: +4.2 / +11.0 (α = 0), +6.1 / +15.6 (α = 0.1) | 41× / 82× Cassini |
| SIM141 | running Λ0(a) | PARTIAL ("θ* 63σ") | PARTIAL / PASS | +9.5 / +21.4 (tanh) | LLR 14σ or Cassini 13× |
| SIM142 | G3(Ψ)□Ψ (kinetic) | FAIL | FAIL / PASS | ≈ the Phase 1 gain (+5.4 … +6.5); quadratic c3 = 0.1 is a **ghost** | ~50× Cassini |
| SIM143 | Ψ + quintessence | FAIL | TRIVIAL | φ adds nothing beyond Phase 1 | ~50× Cassini |
| SIM144 | mixed-source φ | FAIL | FAIL / PASS | up to +9.1 / +19.9 | 66–135× Cassini |
| SIM138 | per-bin diagnostic | DISTRIBUTED ("DESI prefers H0 ≈ 74") | diagnostic | ΛCDM misfit is in the LRG bins (z = 0.51 / 0.71 in DR1; 0.71 / 0.93 in DR2); "H0 ≈ 74" was a table artifact | — |
| SIM126 / 129 | memory-field toys | N/A | N/A (observables only) | SIM126 ε = 2e-10: +13 vs its own baseline, θ* 4σ low; SIM129: γ > 0 always worse | — |
| SIM150 | Λ0 sweep | INCOMPATIBLE | not rerun (robust) | its own docstring finds \|γ−1\| ≈ 9.5e-4 at the locked Λ0 ("Something is wrong") | — |
| SIM88 / 97 / 98 | CLASS / Planck | — | **not rerun** (CLASS and plikHM unavailable here) | SIM88 ran with Ψ ≈ 0.01 (≈ ΛCDM) | — |
| SIM120 | spec only | — | not run | — | — |

Summary figure: `audit_logs/audit_summary.png`, made by `audit_checks/summary_figure.py`. It shows the best θ*-matched Δχ² vs ΛCDM per model (DR1 and DR2), coloured by local-gravity status.

### What this means for the papers (the user's call; nothing in Papers/ was edited)
1. **Paper II** (Phase 1 passes; "2.77σ structural floor"):
   - The canonical model fails Cassini by 41× (no screening).
   - The 2.77σ floor is not a CMSTG property: frozen F0 is exactly degenerate with H0, and ΛCDM itself sits at 2.1σ (DR1) / 3.1σ (DR2) in this setup.
   - The DE is a bare Λ (or CPL), contrary to the README.
   - The README's "G_eff/G = 1 − 15 ppm" contradicts F0 = 0.521.
2. **Paper I** (no-go theorems):
   - The empirical basis of Theorem 2 (SIM141's "θ* 63σ") is an artifact of holding H0 and r_s fixed. With H0 re-solved from θ*, late-time modifications improve DESI + θ*: w0waCDM, and several CMSTG variants with Δχ² ≈ +5 … +9.5 (DR1) and +14 … +22 (DR2).
   - Theorem 1's "wrong direction" argument also assumed a fixed H0. Growing F (SIM130, SIM123) improves the fit once H0 floats.
   - SIM144 does not establish completeness of a no-go.
   - **What does survive is a different obstacle: local gravity.** Every variant with an unscreened scalar at Ψ ~ O(1) fails Cassini, LLR or GW170817.
3. **The interesting survivors are the attractor models SIM132 and SIM134.** Ψ settles where F'(Ψ*) = 0, so the fifth force switches off. They improve DESI + θ* (≈ +9 DR1, ≈ +20 DR2) and pass Cassini and LLR. Still untested:
   - BBN (G about 7% higher than today).
   - The full CMB spectrum (G changes ~6% between z ~ 3000 and z*), which needs hi_class or an equivalent.
   - The Λ0 = 0.2 optimum lies at the edge of the published grid.
   - The DE is still a bare Λ.
4. DESI DR2 (three years of data, published March 2025) is the "Y3" test the README says is awaited for ~2027. It is included in every rerun.

---

## 1. Session log (chronological)

1. (Unrelated earlier task) Merged `Robin Hood Prince of Thieves (1991).mkv` and `... Part 2.mkv` in `/data/Media/MKV2` with `mkvmerge` (lossless append, Part 2's extra subtitle tracks dropped). Originals deleted at the user's request; result renamed to `Robin Hood Prince of Thieves (1991).mkv`.
2. (Unrelated) Checked GitHub account `cisomorph`: no notifications, open PRs or issues; 7 repos; last activity 2026-05-18.
3. First CMSTG audit pass on a shallow clone in the session scratchpad. Read README, RESEARCH_RULES, SIM88, SIM107 and SIM121C. Wrote independent check scripts `check.py` and `check2.py`; these and the later check scripts are now copied to `audit_checks/`.
4. User asked to fix the bugs and rerun the sims. Made a full clone at `/home/aion/CMSTG-fixes` on branch `audit-fixes`.
5. Grepped all sims for the bug patterns (section 3). Downloaded the official DESI BAO data (section 4).
6. Read SIM131 (Phase 3 template) and SIM141 (Phase 4 template) in full. Found the Friedmann-equation and ε_H bugs (B9 and B10 below).
7. User asked for notes of all work: this file.
8. Built and validated `common/cmstg_obs.py`. It reproduces DESI's own ΛCDM fits and Planck θ*.
9. Rewrote SIM121C around the shared module; MCMC reruns for DR1 and DR2. Found the F0–H0 degeneracy (B11).
10. Rewrote SIM107 (Cassini). Verdict PASS → FAIL.
11. Derived the corrected Friedmann / Raychaudhuri / scalar equations. Validated: GR limit exact; constraint-ε equals closed-form ε to 2e-10. Built `common/st_solver.py` and `common/audit_driver.py`.
12. Ported and reran SIM130–134. SIM130, 132 and 134 flip to θ*-matched PASS.
    - Robustness checked independently (closed-form ε, three integrators, three z_ini).
    - Local-gravity check (`audit_checks/viability.py`): SIM132 and SIM134 pass Cassini and LLR via their attractor.
13. Added a blow-up guard to the solver; reran SIM130–134 with identical numbers.
14. Built `common/mf_solver.py` (multi-field, K(φ), growth with G_eff) and validated it. Ported and reran SIM141–144, SIM135, SIM136, SIM138.
15. Ported and reran SIM122 (DR1 and DR2), SIM123 and SIM121B. Evaluated the SIM126/129 toys with corrected observables.
16. SIM124 hung on one case. Added a per-case 60 s watchdog (unresolved cases recorded in its output) and reran.
17. Wrote the executive summary (section 0).
18. Committed (`7b5f22c`) and pushed branch `audit-fixes` at the user's request, while SIM124 was still running.
19. SIM124 finished: verdict FAIL (fine-tuned; the solutions reaching Ψ0 = 2.62 fail Cassini by 41× and 82×). Added it to section 0, section 6, its `results.md` and the summary figure, then committed and pushed as a follow-up at the user's request.

---

## 2. Bug list

Status key: FOUND, FIXED (code changed), RERUN (fixed and re-executed).

### B1. Cassini test evaluated at the wrong field value; wrong coupling formula (SIM107). FOUND
- SIM107 evaluates PPN γ only at Ψ0 = 1.0. Its own output gives `Psi0_critical_Cassini = 1.139`, but the canonical value is Ψ̄ = 2.62, so the published PASS does not apply to the canonical model.
- The coupling formula is also wrong. SIM107 uses α = Λ0Ψ/(1+2Λ0Ψ²) from "A² = 1+2Λ0Ψ²". The correct Einstein-frame coupling for S = ∫√-g[F R − ½(∂Ψ)²] is derived next.
  - Write the action in DEF form, (1/16πG*)∫√-g[F̃R − Z(∂Ψ)²], with 8πG* = 1. Then F̃ = 2F = 1 + 2Λ0Ψ² and Z = 1.
  - α² = F̃'²/(3F̃'² + 2ZF̃).
  - γ − 1 = −2α²/(1+α²) = −F̃'²/(ZF̃ + 2F̃'²).
  - Checks: Brans–Dicke (F̃ = Φ, Z = ω/Φ) gives −1/(ω+2). The Jordan-frame BD mapping gives ω_BD = F/(2F'²), the same result.
- At Ψ = 2.62, Λ0 = 0.003: F̃' = 4Λ0Ψ = 0.03144, and |γ−1| = 9.48e-4. That is 41× the Cassini bound of 2.3e-5.
  - CORRECTION: in my first chat report I said 21×. That figure used Z = 1 with the README normalization F = (1+2Λ0Ψ²)/2, where Z = ½ is correct. 41× is the right number.
  - SIM107's own (wrong) formula gives 1.14e-4, which is 5× the bound. Cassini fails either way.
- The effective-mass sign was chosen by hand. The code comment reads "taking the magnitude; the chameleon direction adds mass in dense regions": m_eff² = m0² + 2Λ0|R|.
  - Derived result (trace equation substituted into the scalar equation): □Ψ(1 + 3F'²/F) = V' − (F'/2F)ρ + …, i.e. V_eff = V − (ρ/2) ln F.
  - Matter contribution to the mass: m² ⊃ −(ρ/2)(F''/F − F'²/F²) < 0. That is anti-chameleon, not chameleon.
  - The Cassini verdict does not depend on this sign. With either sign the field's range at 1 AU is ≫ 1 AU (λ ~ 30–300 Mpc), and SIM107 itself finds no thin shell (χ_sun = 1).

### B2. "DESI tension" is not a significance (≥ 10 sims). FOUND
- `tension = sqrt(chi2/N)` is used as if it were a number of σ. It is an RMS pull.
- SIM121C example: χ² = 41.49 for N = 6 gives "2.63σ". The real goodness of fit is p = 2.3e-7, about 5.2σ with dof = 6, or 5.8σ with dof = 3.
- Files: SIM120, SIM121C, SIM122, SIM123, SIM130, SIM131, SIM138, SIM141, SIM142, SIM143, SIM144 (see the grep in section 3).

### B3. Sound-horizon integral truncated (13+ sims). FOUND
- r_s = ∫_{z_d}^{∞} c_s/H dz is cut off at z = 1e4 in SIM121C, 122, 123, 130, 131, 132, 133, 134, 135, 136, and at z = 1e5 in SIM143 and 144.
- The missing tail in radiation domination is ≈ (c/√3)/(100·√ω_r) / Z = 267,700 Mpc / Z:
  - Z = 1e4 loses ≈ 27 Mpc (r_s = 122.6 instead of 147.0 Mpc, −17%).
  - Z = 1e5 loses ≈ 2.7 Mpc (−1.8%).
  - The θ* error bar is 0.03%, so both are large.
- With the 1e4 cutoff, ΛCDM cannot reproduce Planck θ* for any H0 in [50, 90]. Verified by running SIM121C's own functions.
- SIM141 and SIM142 instead hardcode r_s = 144.7 Mpc, independent of the model.
- SIM129 integrates r_s from 0 to z* (wrong limits). To check.

### B4. θ* uses r_s(z_drag) instead of r_s(z*). FOUND
- θ* = r_s(z*)/D_M(z*). The code uses r_s(z_drag = 1059.6), about 2% too large.
- θ_obs = 1.04101 ± 0.00029 is Planck 100θ_MC from TT,TE,EE+lowE+lensing+BAO, not 100θ*. Using θ_MC with an exact θ* calculation is inconsistent, and the BAO-combined value double-counts BAO.

### B5. DESI data in the code are not DESI measurements (22 sims). FOUND
- Hardcoded table: z = [0.295, 0.510, 0.706, 0.930, 1.317, 2.330], H = [81.7, 97.9, 110.7, 128.1, 156.4, 240.8] ± [4.5, 4.4, 6.2, 5.6, 8.6, 11.0].
- DESI publishes D_M/r_d, D_H/r_d and D_V/r_d, not H(z). Official DR1 values, converted with r_d = 147.09, give H = 97.2, 101.5, 114.0, 147.5, 239.2 at z = 0.51, 0.706, 0.93, 1.317, 2.33.
  - The code's values at 0.706 and 0.93 are off by 3.0σ and 6.3σ.
  - z = 0.295 (BGS) has only D_V, so there is no H(z) measurement there at all.
- Files: SIM120, 121B, 121C, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 135, 136, 138, 141, 142, 143, 144, 150.

### B6. H0 label is not H(z=0) in SIM121C. FOUND
- E²(z) = (0.5/F0)(Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_DE f_DE) with Ω_DE = 1 − Ω_m − Ω_r, so E(0)² = 0.5/F0 ≠ 1.
- True H0 = H0_param·√(0.5/F0). The published MAP "H0 = 75.96" is really 71.77.
- The MAP F0 = 0.55998 sits at the prior edge (0.56).

### B7. Hardcoded, inconsistent baselines. FOUND
- SIM131, SIM141 and others hardcode `SIM121C_chi2 = 18.26` and `tension = 2.77`. SIM121C's own published output is χ²_DESI = 41.49 and "tension" = 2.63. The 2.77σ in the README is not reproduced by SIM121C's output.

### B8. Bare cosmological constant in "no bare Λ" theory; Ψ ≈ 0 in "Phase 1" CMB runs. FOUND (model-level, documented, not "fixed")
- Ω_L = 1 − Ω_m − Ω_r (or Λ_bare) is used in SIM87, 88, 90, 92, 93, 97, 98, 108, 109, 127, 139 and in all Phase 3/4 sims (`Lambda_bare`). The README says "Dark energy emerges from the Ψ–geometry coupling without a bare cosmological constant."
- SIM88 starts at Ψ_ini = 0.01, so G_eff/G = 0.99998; that is effectively ΛCDM, not the Ψ̄ = 2.62 model. Its CMB difference from ΛCDM comes only from the different H0/Ω_m inputs.
- SIM88's `results.md` calls it "G_eff Solar System / GW Speed", but the code is a CLASS CMB run.
- The README says "G_eff/G = 1 − 15 ppm at z = 0". For F_eff = 0.5206, G_eff/G = 0.5/F = 0.960 (−4%).

### B9. Friedmann equation coefficients wrong (Phase 3/4 families). FOUND, verifying per sim
- Derivation (minisuperspace, lapse N; signature −+++): for S = ∫√-g[F(Ψ)R − ½(∂Ψ)²] + S_m,
  6FH² + 6HḞ − ½Ψ̇² = ρ (GR check: F = ½ gives 3H² = ρ).
  In e-folds with ρ = 3·ω·H100²: H²(2F + 2F'Ψ' − Ψ'²/6) = ω_m a⁻³ + ω_r a⁻⁴ + ω_Λ.
- SIM131 uses H²(3F + 6F'Ψ' − ½Ψ'²) = ω_m a⁻³ + …. In the GR limit (F = ½, Ψ' = 0) that gives H² = ω/1.5, not ω, so the matter era runs with G_eff = (2/3)G. The F'Ψ' and kinetic terms are also 2× too large relative to the F term.
- SIM141 uses coef = 3F + 6F'Ψ' + 6Λ0'Ψ² − ½Ψ'² with rhs = 3F0·ω + Λ_b.
  - The rhs normalization is a deliberate choice: early-time H matches GR with Planck ω.
  - Relative to 3F, the correct coefficients are 3F'Ψ' and −¼Ψ'². The F'Ψ' and kinetic terms are 2× too large.

### B10. ε_H ignores the scalar's own evolution (Phase 3/4 families). FOUND
- ε_H is computed by finite-differencing H²(Ψ, Ψ', N) in N with Ψ and Ψ' held fixed. That drops the dΨ/dN and dΨ'/dN contributions.
- Closed form from the ij equation plus the scalar EOM, with F' = dF/dΨ and Ψ' = dΨ/dN:
  ε ≡ −Ḣ/H² = [(ρ+p)/H² + (1 + 2F'')Ψ'² + 24F'² − 8F'Ψ'] / (4F + 12F'²).
  GR checks: matter gives 3/2, radiation gives 2.
- Source: DEF/BEP-form equations with F_BEP = 2F and Z = 1:
  −4FḢ = ρ + p + Ψ̇² + 2F̈ − 2HḞ, and Ψ̈ + 3HΨ̇ = 6F'(Ḣ + 2H²).

---

## 3. Grep results (where each bug lives)

- Fake DESI table (B5): SIM120 (2 files), 121B, 121C, 122, 123, 124, 126, 129, 130, 131, 132, 133, 134, 135, 136, 138, 141, 142, 143, 144, 150.
- sqrt(χ²/N) metric (B2): SIM120, 121C, 122, 123, 130, 131, 138, 141, 142, 143, 144.
- r_s cut at 1e4 (B3): SIM121C, 122, 123, 130, 131, 132, 133, 134, 135, 136. Cut at 1e5: SIM143, 144. Hardcoded r_s = 144.7: SIM141, 142. Odd limits (0 → z*): SIM129. Correct (np.inf): SIM121B.

---

## 4. External data used

- DESI DR1 (2024) BAO, official Gaussian likelihood: `CobayaSampler/bao_data`, files `desi_2024_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`. Saved as `Ordered_Simulations/common/data/desi_dr1_bao_{mean,cov}.txt`.
- DESI DR2 (2025, three years of data) BAO: `desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`. Saved as `common/data/desi_dr2_bao_{mean,cov}.txt`.
  - Note: the README says the "next decisive test" is "DESI Y3 (~2027)". DESI DR2, based on the first three years of data, was released in March 2025.

---

## 4b. Shared module `Ordered_Simulations/common/cmstg_obs.py` (new) and its validation

Contents:
- Official DESI DR1/DR2 Gaussian BAO likelihood (`DesiBAO`).
- `sound_horizon` integrated to z = ∞ (in ln a from a = 1e-9), plus `theta_star_100` = r_s(z*)/D_M(z*) and `h_from_theta`.
- `H_from_table`: wraps a tabulated background and extends it above its last redshift with matter+radiation scaling, i.e. the gravity modification is held frozen up there.
- Statistics: `chi2_to_sigma` (p-value to Gaussian σ, safe for huge χ²), `delta_chi2_to_sigma` (Wilks), `summarize_fit`.
- `lcdm_reference`: ΛCDM computed in the same pipeline as the sims.
- Planck reference values: 100θ* = 1.04110 ± 0.00031, z* = 1089.92, z_drag = 1059.94 (Planck 2018 VI, TT,TE,EE+lowE+lensing).

Calibration:
- The simple r_s integral in the sims' simplified background (ω_m = 0.1430 including the massive ν; ω_r = 4.18e-5 counting all 3.046 species as radiation) is 0.08% below Planck's Boltzmann-code values. Fixed factors are applied: K_STAR = 1.00081, K_DRAG = 1.00083.
- Independent check: with a proper massive-neutrino background (1 × 0.06 eV with the Fermi–Dirac integral, plus 2 massless species), the raw integral gives 144.438 / 147.098 Mpc against Planck's 144.43 / 147.09. So the offset is the neutrino double-count, with a −0.01% residual.

Validation (`common/validate_cmstg_obs.py`, output in `common/validation_output.json`):
- 100θ* at h = 0.6736 is 1.04134 (+0.76σ). h solved from θ* is 0.6728 (Planck 0.6736 ± 0.0054).
- DESI BAO-only flat ΛCDM fit, compared with the published DESI papers:
  - DR1: Ω_m = 0.2939 (+0.015/−0.014), h·r_d = 101.94. DESI 2024: 0.295 ± 0.015, 101.8 ± 1.3.
  - DR2: Ω_m = 0.2975 ± 0.0087, h·r_d = 101.54. DESI 2025: 0.2975 ± 0.0086, 101.54 ± 0.73. The likelihood is reproduced.
- The old truncated r_s(z_drag) is reproduced: 122.57 Mpc (cut at 1e4) and 144.33 Mpc (cut at 1e5), against 147.09.
- ΛCDM in the sims' setup (Planck ω's fixed, H0 from θ*):
  - DR1: χ² = 22.04/12, p = 0.037 (2.08σ GOF).
  - DR2: χ² = 32.30/13, p = 0.002 (3.07σ GOF).
  - So ΛCDM itself sits about 2–3σ from DESI in this fixed-ω setup. That is the real baseline any CMSTG mechanism has to beat.
- w0waCDM in the same setup is preferred over ΛCDM by 2.87σ (DR1) and 4.61σ (DR2).
  - Best fits: w0 = −0.35, wa = −2.0 (DR1); w0 = −0.39, wa = −1.84 (DR2).
  - DESI's own DESI+CMB numbers are 2.6σ and 3.1σ. The pipeline overstates the preference because ω_m is fixed rather than marginalized. This is a sanity check only, not a claim.

## 5. Plan / status

- [x] `common/cmstg_obs.py`: DESI likelihood (DR1/DR2), r_s to ∞, θ* at z*, significance helpers, ΛCDM reference. Validated (section 4b).
- [x] SIM121C: fixed B2–B6 (script rewritten around `cmstg_obs`; `--release DR1|DR2`; H0 kept as an emcee blob; prior w0+wa<0). MCMC reruns for DR1 and DR2 are done (section 6).
- [x] SIM107: fixed B1 plus two unit bugs; rerun. Verdict FAIL (section 6).
- [x] Phase 3 (SIM130–136) and Phase 4 (SIM141–144): fixed B2–B5, B7, B9, B10 and the model-specific bugs; rerun in faithful and θ*-matched modes; GR limit and ε consistency validated.
- [x] Other DESI-table sims: SIM121B, 122, 123, 124, 138 rerun; SIM126/129 observables corrected; SIM150 and SIM120 documented, not rerun.
- [x] Before/after table (section 0); summary figure (`audit_logs/audit_summary.png`).
- [ ] Not possible here: SIM88/97/98 need CLASS and the Planck plikHM likelihood.
- [ ] Suggested next steps (not done): a full Boltzmann (hi_class) plus BBN check of the SIM132/SIM134 attractor models; extend their Λ0 grid below 0.2.

### New finding B11. F0 is exactly degenerate with H0 in the frozen-Ψ model (SIM121C)
- With Ψ frozen, 6F0H² = ρ is GR with every ω multiplied by g = 0.5/F0, while ω_b/ω_γ in the sound speed is unchanged.
- Then H(z) → √g·H(z), r_s → r_s/√g and D → D/√g, so θ* and every D/r_d are unchanged once H0 is re-solved from θ*.
- Checked numerically: F0 = 0.5 and F0 = 0.5206 (with w = −1) both give χ²_DESI(DR1) = 22.0391 exactly. H0 moves from 67.28 to 65.94 and r_d from 147.09 to 150.09.
- Consequences:
  - F0 has no effect on the DESI/θ* fit; all of SIM121C's improvement over ΛCDM comes from the phenomenological CPL (w0, wa).
  - Raising F0 only lowers the CMB-inferred H0, which worsens the Hubble tension (the original figure already noted the direction).
  - The published MAP sat at the F0 prior edge (0.56) only because of bugs B3–B6.

### SIM107 extra bugs found while fixing
- G·M_sun/c² was written as 1.476e-12 Mpc; the correct value is 4.78e-20 Mpc. That made Φ_N(Sun) = 65.5 instead of 2.1e-6, so the ad-hoc thin-shell factor clipped to 1.
- 1 AU was written as 4.848e-9 Mpc; the correct value is 4.848e-12 Mpc.
- The no-screening conclusion survives a proper Khoury–Weltman test: m_in·R_sun = 2.6e-4 ≪ 1 (with either mass sign), and ΔR/R = 1.3e7 ≫ 1. The published code reached it by accident.

## 6. Results (before → after)

### SIM107 — Solar System / Cassini (rerun 2026-10-22, log `audit_logs/SIM107.log`)
| quantity | published | audited |
|---|---|---|
| Ψ0 evaluated | 1.0 | 2.62 (canonical) |
| \|γ−1\| | 1.78e-5 | 9.48e-4 (41× the 2.3e-5 bound) |
| Ψ0 allowed by Cassini | < 1.139 | < 0.400, i.e. F0 − ½ < 4.8e-4 |
| screening | χ_sun = 1 (unit bug) | none: m_in·R_sun = 2.6e-4, ΔR/R = 1.3e7 |
| m_eff² from matter | +2Λ0ρ (hand-chosen sign) | −(ρ/2)(F''/F − F'²/F²) ≈ −0.0053ρ |
| **verdict** | **PASS** | **FAIL** |

### SIM121C — Joint DESI + θ* MCMC (rerun 2026-10-22; logs `audit_logs/SIM121C_DR{1,2}.log`; DR2 outputs in `SIM121C/Outputs/DR2/`)
Setup: 32 walkers × 3000 steps after 800 burn-in. Acceptance 0.56. τ ≈ 46–49 steps, so the chains are about 60τ long and converged.

All models below were run through the same pipeline (Planck ω's fixed, H0 solved from θ*):

| model | free params | DR1 χ²/N | DR1 GOF | DR1 Δχ² vs ΛCDM | DR2 χ²/N | DR2 GOF | DR2 Δχ² vs ΛCDM |
|---|---|---|---|---|---|---|---|
| ΛCDM | 0 | 22.04/12 | 2.08σ | — | 32.30/13 | 3.07σ | — |
| Phase 1 frozen F0 = 0.5206, w = −1 | 0 | 22.04/12 | 2.08σ | 0.00 | 32.30/13 | 3.07σ | 0.00 |
| F0 free, w = −1 | 1 | 22.04 | 2.26σ | 0.00 | 32.30 | 3.23σ | 0.00 |
| SIM113 point (F0 = 0.5206, w0 = −0.973, wa = −0.41) | 0 | 14.17 | 1.06σ | +7.87 | 16.71 | 1.25σ | +15.59 |
| w0waCDM (F0 = ½) | 2 | 11.51 | 1.00σ | +10.53 | 8.51 | 0.43σ | +23.78 |
| CMSTG full (F0, w0, wa) | 3 | 11.51 | 1.17σ | +10.53 | 8.51 | 0.56σ | +23.78 |

- F0 adds Δχ² = 0.00 beyond w0waCDM in both releases (B11). The F0 posterior is flat over its prior; the MAP F0 (0.493 for DR1, 0.509 for DR2) means nothing.
- The improvement over ΛCDM is entirely the phenomenological CPL (w0, wa) dark energy, i.e. DESI's own published result. w0 runs to the original prior edge (−0.60).
- MAP (DR1): w0 = −0.600, wa = −1.357, H0 = 66.19 (CMB-derived); χ² = 11.51 for 12 points, dof 9, GOF 1.17σ.
- MAP (DR2): w0 = −0.600, wa = −1.28, H0 = 64.73; χ² = 8.51 for 13 points, GOF 0.56σ.
- Verdict PARTIAL (w0 outside the sim's own (−1.3, −0.7) window). Published verdict: PARTIAL, "2.63σ" from the non-DESI table.
- Published MAP (F0 = 0.560, "H0" = 75.96) does not survive. In the corrected pipeline that point has no θ* solution with H0 in [55, 85].

### Infrastructure for the Phase 3/4 reruns (2026-10-22)
- `common/cmstg_bg.py`: constraint-consistent ε, plus the closed-form Raychaudhuri ε for cross-checks.
- `common/st_solver.py`: corrected single-field solver for F(Ψ,N)R + V + f(Ψ)𝒢, with Λ_bare calibration and θ*-matching.
  - The Gauss–Bonnet Friedmann term was derived via minisuperspace: H²(6F + 6F'Ψ' − ½Ψ'²) + 24f'Ψ'H⁴ = ρ + V. It agrees with Guo & Schwarz 2010 after mapping f = −ξ/2.
- `common/audit_driver.py`: runs every case in two modes:
  - faithful: H0 fixed as published.
  - θ*-matched: H0 solved so 100θ* matches Planck. This is the fair test; a fixed H0 makes θ* fail for any change in H(z).
- `common/validate_cmstg_bg.py` (output `validation_bg_output.json`):
  - GR limit is exact (6.7e-16). The published equations' GR limit is off by 18%: Ω_m,eff = 0.21, 100θ* = 1.0848, DESI DR1 χ² = 464.
  - Constraint-ε equals the Raychaudhuri closed form to 2.3e-10 along the SIM131 ξ = 1/6 trajectory, where ε runs from 0.36 to 1.99. So the corrected Friedmann, scalar and Raychaudhuri equations are mutually consistent.
  - Solver checks: it reproduces the validation run exactly; its GR limit with θ* matched reproduces the ΛCDM reference to Δχ² = 0.0000; f_GB → 0 agrees to 2e-13.

### SIM131 (rerun; log `audit_logs/SIM131.log`) — FAIL, unchanged
- Reference ξ = 1/6, Λ0 = 0.003: Ψ0 = 2.668, F0 = 0.966 (G_eff roughly halved).
  - Faithful: 100θ* = 0.892, χ²(DR1) = 554.
  - θ*-matched: no H0 in [40, 100] reproduces θ*.
- ξ ≥ 1/6: no θ*-matched solution.
- Best case ξ = 1/20: θ*-matched H0 = 74.5, χ² = 71.0 (DR1), Δχ² = −49 vs ΛCDM; DR2 Δχ² = −165.
- Verdict FAIL in both modes. The published FAIL holds (its numbers were wrong: 100θ* = 0.797, "5.97σ").

### SIM130 (rerun; log `audit_logs/SIM130.log`) — VERDICT CHANGES: published FAIL → θ*-matched PASS
The model is the Phase 1 canonical action with Ψ unlocked: F = ½ + 0.003Ψ², no potential, Λ_bare, scanned over Ψ_ini.
- Ψ_ini = 2.62: Ψ grows to 2.882 today (+10%); F runs from 0.5212 at z* to 0.5249 today.
  - θ*-matched: H0 = 66.64, χ²(DR1) = 16.61 vs ΛCDM 22.04, so **Δχ² = +5.43**. DR2: χ² = 18.36 vs 32.30, so **Δχ² = +13.95**.
  - Faithful (H0 = 67.4): χ² = 13.15 but 100θ* = 1.0434, which is 7.5σ off. So FAIL in that mode.
- Best of the scan, Ψ_ini = 3.5 (θ*-matched): DR1 Δχ² = +7.89 (GOF 1.21σ), DR2 Δχ² = +19.68 (GOF 0.85σ). Verdict PASS by the sim's own thresholds.
- Ψ_ini ≤ 0.5 gives Δχ² ≈ 0, consistent with the GR limit.
- Robustness, from an independent script (`audit_checks/robust130.py`: closed-form ε; DOP853 / LSODA / Radau; z_ini = 1e5 / 1e6 / 1e7): Δχ² = +5.429 / +5.433 / +5.434 (DR1) and +13.94 / +13.95 / +13.95 (DR2). Not a numerical artifact.
- **Not viable as is:**
  - Ψ0 ≈ 2.88 gives |γ−1| ≈ 1.15e-3, about 50× Cassini (SIM107 corrected), with no screening in the locked action.
  - Ġ/G ≈ −Ḟ/F = −1.8e-3·H0 ≈ −1.2e-13 /yr, about 2.5σ from LLR (Hofmann & Müller 2018: (7.1 ± 7.6)e-14).
  - Λ_bare is still present.
- Implications:
  - The published SIM130 conclusion ("curvature sourcing negligible / the 2.77σ is structural") is reversed for cosmology.
  - Paper I's "wrong direction" (Theorem 1) argument assumed a fixed H0. With H0 re-solved from θ*, growing F with Ψ_ini ≈ 2.6 improves the DESI fit.
  - The real obstacle is the Solar System, not DESI.

### SIM132 / SIM133 / SIM134 (rerun; logs `audit_logs/SIM13{2,3,4}.log`)

**SIM132** (subtractive F = ½(1 − 2Λ0Ψ²) + ξΨ, ξ = 1/6, Ψ_ini = 0): **published FAIL → θ*-matched PASS** (faithful: PARTIAL).

| Λ0 | θ*-matched H0 | χ² (DR1) | Δχ² vs ΛCDM (DR1) |
|---|---|---|---|
| 0.2 | 67.64 | 12.77 | +9.27 (GOF 1.02σ); DR2 +21.24 (GOF 0.64σ) |
| 0.33 | — | — | +4.25 |
| 0.5 | — | — | +1.57 |
| 1.0 | — | — | +0.39 |

- The best value sits at the edge of the published grid, so the true optimum may lie at smaller Λ0.

**SIM134** (F = ½e^{−2Λ0Ψ²} + ξΨ): **published FAIL → θ*-matched PASS**.
- Λ0 = 0.2: H0 = 67.93, Δχ² = +8.88 (DR1) / +19.20 (DR2).
- Cosmetic bug: the published `Psi_attractor` bracket [1e-6, 10] has the same sign at both ends, so it always reported no attractor. Fixed to [1e-6, 1/(2√Λ0)].

**SIM132 and SIM134 robustness** (`audit_checks/robust132.py`: closed-form Raychaudhuri ε; DOP853 / LSODA / Radau; z_ini = 1e5 / 1e6 / 1e7):

| sim | Δχ² DR1 | Δχ² DR2 | H0 | Ψ0 |
|---|---|---|---|---|
| SIM132 | +9.27 / +9.30 / +9.30 | +21.24 / +21.75 / +21.78 | 67.64 → 67.44 | 0.4163 |
| SIM134 | +8.88 / +9.12 / +9.14 | +19.20 / +20.34 / +20.43 | — | — |

- The small drift with z_ini comes from Ψ creeping during the radiation era (F'(0) = ξ ≠ 0). The conclusion is unchanged.

**SIM133** (Gauss–Bonnet α̃Ψ𝒢): **FAIL (unchanged)**, but for a different reason.
- With the derived GB Friedmann term (8α_GBΨ'H⁴) the constraint loses its physical root almost immediately for every α̃ ≥ 1e-8.
- In radiation domination 𝒢 = −24H⁴ ≠ 0, so the source is about 2e9 at z = 1e5 and Ψ runs away.
- The published run omitted the GB Friedmann term and reported 100θ* = 63.

**Solar-System viability of the improved models** (θ*-matched best cases; script `audit_checks/viability.py`):

| case | Ψ0 | F_Ψ(Ψ0) | \|γ−1\| (× Cassini) | Ġ/G [1/yr] | G_BBN/G0 − 1 |
|---|---|---|---|---|---|
| SIM130 Ψ_ini = 2.62 | 2.882 | 0.0173 | 1.1e-3 (49×) | −1.2e-13 | +0.8% |
| SIM130 Ψ_ini = 3.5 | 3.850 | 0.0231 | 2.0e-3 (85×) | −2.1e-13 | +1.4% |
| SIM132 Λ0 = 0.2 | 0.4163 | 1.5e-4 | 7.9e-8 (0.003×) | −4e-18 | +6.9% |
| SIM134 Λ0 = 0.2 | 0.4527 | −1.7e-4 | 1.1e-7 (0.005×) | −2e-17 | +7.2% |

- SIM132 and SIM134 sit on their attractors today (F'(Ψ*) = 0), so the Einstein-frame coupling α ∝ F' vanishes and Cassini and LLR pass easily.
- Open issues for SIM132/SIM134:
  - G at BBN is about 7% higher than today (marginal; needs a proper BBN bound, e.g. Alvey et al. 2020).
  - G changes by about 6% between z ~ 3000 and z*. θ* + BAO does not test that; a full Boltzmann code (hi_class / CLASS-EFT) is needed.
  - Dark energy is still a bare Λ.

The blow-up guard (exception when the constraint has no physical root, plus |Ψ| > 1e3) was added to `st_solver.py` afterwards. SIM130–134 were rerun and all numbers are unchanged.

### Multi-field solver `common/mf_solver.py` (validation `common/validate_mf_solver.py`, output `validation_mf_output.json`)
- Supports F(φ_i, N), K_i(φ_i), V, a phenomenological V_int, and growth with the quasi-static G_eff = (1/2F)(2F + 4S)/(2F + 3S), where S = Σ F_i²/K_i.
- Single field reproduces `st_solver` exactly: SIM130 θ*-matched H0 = 66.6395, χ² = 16.6102.
- Two-field GR limit: 6.7e-16.
- GR growth: σ8 = 0.811, fσ8(0, 0.5, 1) = 0.428 / 0.475 / 0.431, the standard Planck-ΛCDM values.

### SIM141 — running Λ0(a) (phenomenological), rerun; log `audit_logs/SIM141.log`
- Published: PARTIAL, best DESI case "tanh a_tr = 0.50, γ = 1" with "θ* 63σ". That is the basis of Paper I Theorem 2 ("DESI fix and CMB mutually exclusive").
  - Cause: H0 and r_s were both held fixed (r_s hardcoded at 144.7 Mpc), plus the coefficient bugs.
- Consistency check: the audited constant-Λ0 reference, θ*-matched, gives χ²(DR1) = 16.61, Δχ² = +5.43, H0 = 68.00. That is exactly SIM130 with H0 × √(2F0), as the density-normalisation degeneracy predicts.
- Audited, θ*-matched:

| case | H0 | Δχ² DR1 | Δχ² DR2 | RSD χ²/N | verdict |
|---|---|---|---|---|---|
| tanh a_tr = 0.67, γ = 1.0 | 65.93 | +9.52 | +21.39 | 0.78 | PASS |
| tanh a_tr = 0.67, γ = 0.5 | — | +7.65 | +18.18 | — | PASS |
| linear β = 0.1 / exp β = 0.1 | — | +4.4 | — | — | PARTIAL_PHENOM |
| tanh a_tr = 0.50, γ = 1.0 (published "best") | 64.69 | −10.49 | — | — | FAIL_DESI |

  - Overall: faithful PARTIAL, θ*-matched PASS. The cosmology-level anti-correlation is gone once H0 floats and r_s is computed.
- **Local gravity rules these out** (`audit_checks/viability141.py`):

| case | \|γ−1\| (× Cassini) | Ġ/G [1/yr] | LLR |
|---|---|---|---|
| tanh 0.67 / 1.0 | 1.5e-6 (0.1×) | +1.14e-12 | 14σ (Λ0 switching off today makes G run fast) |
| tanh 0.67 / 0.5 | 13.5× | +5.3e-13 | 6σ |
| constant (Phase 1) | 49× | −1.2e-13 | −2.5σ |

### SIM135 — bi-scalar (Ψ frozen, F = F0 + ξφ, φ_ini = 0); log `audit_logs/SIM135.log`
- Published: FAIL ("3.73σ"). The published Friedmann had the group-A bug, rhs = ω with coef 3F.
- Audited, θ*-matched:
  - ξ = 0.03: H0 = 68.26, Δχ² = +9.00 (DR1) / +19.60 (DR2). Cosmology PASS.
  - ξ = 0.01: Δχ² = +2.3 / +6.1.
  - ξ ≥ 0.05: FAIL; ξ = 0.15 has no θ* solution.
- **Local gravity:** φ's coupling ξ is unscreened, giving Cassini 17× (ξ = 0.01) and 145× (ξ = 0.03). Not viable.

### SIM136 — Horndeski α G^{μν}∂Ψ∂Ψ; log `audit_logs/SIM136.log`
Derived equations (minisuperspace; cross-checked against the Horndeski G4 = F + αX background equations):
- Friedmann: 6FH² + 6HḞ = ρ + ½Ψ̇²(1 + 18αH²). The published sim used (1 + 6αH²): the coupling scales as N⁻³ under the lapse, so its energy is 3× its Lagrangian value.
- Field: (1 + 6αH²)Ψ̈ + [3H + 6α(3H³ + 2HḢ)]Ψ̇ = F'R. Published: "unchanged at leading order".
- Tensor speed: c_T² = (F + αX)/(F − αX). This was not checked in the published sim, although RESEARCH_RULES 2.3.2 requires c_T = c.

Results:
- Test A (Ψ_ini = 0): F'(0) = 0, so Ψ stays 0 and the result is trivially ΛCDM.
- Test B (Ψ_ini = 2.62), θ*-matched: Δχ² = +4.6 → +2.5 for α = 1e-8 → 1e-5 (DR1). This is the Phase 1 unlocked gain, slightly degraded.
- c_T² − 1 = 2.4e-11 … 2.3e-8 today, against GW170817's < ~6e-15. **FAIL** for every α in the published scan.
- α ≥ 1e-4: the constraint has no physical root.

### SIM142 — "Galileon" G3(Ψ)□Ψ ≡ kinetic rescaling K = 1 − 2G3_Ψ; log `audit_logs/SIM142.log`
- Published bugs: the G3 term entered Friedmann linearly in Ψ' (it must be quadratic); the F'Ψ' and kinetic terms were doubled; the field equation lacked the 1/K.
- Audited, θ*-matched: cosmology PASS for all non-ghost cases, Δχ² = +5.4 … +6.5 (DR1) and +13.9 … +16.5 (DR2); the best case is linear c3 = 0.1. The gain is inherited from the Phase 1 unlocked background (SIM130); c3 adds at most ~1.1.
- Quadratic form with c3 = 0.1: K(Ψ̄) = −0.048, a **ghost**. The published "α_K > 0, no ghost" was wrong for that case.
- Cassini 49–57× in every case. **Not viable.**

### SIM138 — DESI per-bin diagnostic, official likelihood; log `audit_logs/SIM138.log`
Per-bin χ² is exact, because DESI's covariance is block-diagonal by tracer bin.
- DR1:
  - ΛCDM total 22.04, driven by z = 0.51 (9.4) and z = 0.706 (9.5).
  - Phase 1 unlocked (Ψ_ini = 2.62): 16.61 (z = 0.51: 7.6, z = 0.706: 6.8).
  - w0waCDM: 11.05.
- DR2:
  - ΛCDM total 32.30, driven by z = 0.706 (14.6), z = 0.934 (7.7), z = 0.51 (4.7), z = 0.295 (2.6).
  - Phase 1 unlocked: 18.36.
  - w0waCDM: 7.48.
- The published conclusion "DESI prefers H0 ≈ 74" (ΛCDM best fit to the table: H0 = 74.0, Ω_m = 0.290) was an artifact of the non-DESI table.

### SIM143 — Ψ + quintessence φ (two fields); log `audit_logs/SIM143.log`
- Published bugs:
  - Friedmann coefficients (as in SIM141).
  - U's energy weight vs its force were inconsistent by a factor of 2. The energy weight was kept (V = 2U in ρ-units) and the force corrected.
  - r_s was anchored to 144.7 + Δr_s. Removed.
- θ*-matched: every non-failing case is **TRIVIAL**. φ adds ≤ 0 beyond the Phase 1 unlocked gain (Δχ² = +5.43 DR1 / +13.94 DR2), and large U0 makes things worse.
- Faithful mode reports PASS_PHENOMENOLOGICAL for exp λ = 1 / U0 = 0.2 and power n = 1 / U0 = 0.2. There φ just stands in for the H0 shift needed to match θ* at H0 = 67.59; it is not a real improvement.
- Cassini ~50× (from Ψ). Not viable.

### SIM144 — mixed-source φ (ξ_R, β_m) completeness probe; log `audit_logs/SIM144.log`
- Published bugs:
  - The Friedmann term +2β_mφρ_m had the opposite sign to its own φ equation; the energy must be −L_int.
  - The φ equation used ω_m a⁻³ instead of the physical ρ_m.
  - Coefficient bugs as in SIM141.
- θ*-matched results:
  - ξ_R = 0.01, β_m = 0.01: H0 = 69.65, Δχ² = +9.07 (DR1) / +19.94 (DR2).
  - ξ_R = 0, β_m = 0.01: +7.93 / +19.76.
  - ξ_R = 0.01, β_m = 0: +6.84 / +17.34.
  - Any β_m ≥ 0.1 or ξ_R ≥ 1: no solution. ξ_R = 0.1: huge FAIL.
- **Cassini 66–135× for every cosmology pass** (the direct matter coupling adds α² = 8β_m²). Not viable.
- The published claim "no case satisfies DESI and θ*" (the completeness argument for Paper I's no-go) is false in the corrected pipeline. Small-coupling cases pass the cosmology tests; local gravity is what excludes them.

### SIM124 — non-canonical kinetic term K(Ψ) = e^{−αΨ²} with the hilltop potential; log `audit_logs/SIM124.log`
- Published bugs:
  - Friedmann with ρ = ω divided by 3F (no GR limit).
  - **The F'R curvature source was missing from the field equation**, although F(Ψ)R is in the action.
  - θ* integrated r_s from z = 0 to 1090; it must be z* → ∞.
  - A second variant of the non-DESI table.
- Audited: corrected equations with the source restored; λ calibrated to H0 (no bare Λ); z_ini = 1500; Ψ'_ini = 0.
  - Each case gets a 20 s watchdog. For α ≥ 0.3 and large Ψ_ini, K drops (e.g. 0.30) and the field runs away; the integration crawls and the case is recorded as unresolved. That is 41 cases, listed in `Outputs/sim124_results.json`.
- Part A, Ψ(z=0) over (α, Ψ_ini): no roll from the hilltop. Ψ(z=0) ≈ 1.09–1.2 × Ψ_ini for small Ψ_ini. For α ≥ 0.3 the field runs away from Ψ_ini ≈ 2 (α = 0.3–1.0), ≈ 1 (α = 1.5–3) or ≈ 0.5 (α = 5).
- Part B, Ψ_ini(α) with Ψ(z=0) = 2.62:
  - α = 0: Ψ_ini = 2.403 (ΔΨ = 0.22).
  - α = 0.1: Ψ_ini = 2.246 (ΔΨ = 0.37).
  - α ≥ 0.2: no bracketed crossing. The upper bracket runs away, so the Ψ0 = 2.62 locus lies between the largest finite trial and the runaway boundary and is not resolved here.
- Part C (θ*-matched):
  - α = 0: H0 = 66.66, Δχ² = +4.24 (DR1) / +10.99 (DR2), Cassini 41×.
  - α = 0.1: H0 = 67.09, Δχ² = +6.09 / +15.60, Cassini 82× (K ≈ 0.5 strengthens the coupling).
- **Verdict FAIL**, same as published but on a corrected basis. Reaching Ψ0 = 2.62 needs Ψ_ini ≈ 2.25–2.40, i.e. fine-tuning, not a roll from the hilltop. The published version needed Ψ_ini ∈ [1.84, 2.57]. Every solution that reaches Ψ0 = 2.62 fails Cassini.

### SIM123 — running F with hilltop dark energy, V = λ(Ψ² − v²)², v = 13.16, no bare Λ; log `audit_logs/SIM123.log`
- Published bugs:
  - 3F H² = ρ + V with ρ = ω (no GR limit).
  - First-order slow roll with R = 12H² at all times (the de Sitter value). The true R/H² = 6(2 − ε) is 3 in matter and ≈ 0 in radiation.
  - Non-DESI table, truncated r_s.
- Audited (full second-order equations; λ calibrated to H0, no Λ_b): **published FAIL → θ*-matched PASS**.
  - Best: Λ0 = 0.0064, Ψ_ini = 2.0 → Ψ0 = 2.47, F0 = 0.539, H0 = 67.57; Δχ² = +9.30 (DR1) / +21.51 (DR2).
  - All 100 grid points have Δχ² ≥ 0. Small Ψ_ini gives ≈ 0, because the field stays on the hilltop.
- **Cassini 162×** at the best point. Not viable.
- The gain is curvature-sourced growth of Ψ (G_eff falling at late times), not hilltop roll. With v ≫ Ψ the potential is effectively a cosmological constant λv⁴, so "no bare Λ" holds only formally.

### SIM121B — H0 tension; log `audit_logs/SIM121B.log`
- The published background had the 3F normalisation (no GR limit), used the non-DESI table, and compared "ΛCDM-inferred H0" against SH0ES without computing the CMSTG universe's true H0.
- Audited, all models θ*-matched:

| model | true H0 | tension with SH0ES |
|---|---|---|
| ΛCDM | 67.28 | 5.03σ |
| Phase 1 frozen | 65.94 | 6.21σ |
| Phase 1 unlocked | 66.64 | 5.60σ |

- **FAIL**: F > ½ means G_eff < G, which lowers the true H0, so CMSTG worsens the Hubble tension.

### SIM122 — unlocked Λ0 + CPL (frozen Ψ); logs `audit_logs/SIM122_DR{1,2}.log`
- Profile of min χ² over Λ0 ∈ [1e-4, 0.01]: total spread Δχ² = 1.5e-8. **Λ0 is exactly unconstrained** (B11).
- MCMC, DR1: MAP w0 = −0.600 (prior edge), wa = −1.355, χ² = 11.51 (GOF 1.17σ), Δχ² = +10.53 vs ΛCDM. That is exactly the w0waCDM value from SIM121C.
- Verdict PARTIAL, only because the unconstrained Λ0 wandered to 0.009, putting F0 = 0.562 outside the sim's [0.5, 0.56] window. Not evidence for CMSTG.
- Published: FAIL ("F0 forces H(z) ~10% too high"; that statement came from bugs B3–B6).

### SIM126 / SIM129 — memory-field toy models (phenomenological, no action); `common/audit_memory_field.py`, log `audit_logs/SIM126_SIM129.log`
- The models are kept unchanged; only the observables were corrected (official DESI, r_s, θ*).
- Published bugs:
  - A second variant of the non-DESI table (errors [2.9, 3.8, 4.5, 5.3, 7.3, 15.7]).
  - SIM129's r_s was integrated from 0 to z* (wrong), with Ω_b in place of ω_b (R off by ≈ 2.2).
  - E(1) ≠ 1 when ε or γ ≠ 0, so H(z=0) ≠ 67.4.
- The baseline (ε = γ = 0) has 100θ* = 1.0360, which is 16σ low, because the toys use Ω_r = 9.4e-5 instead of about 9.2e-5. DR1 χ² = 27.8.
- SIM126: ε = 2e-10 improves on its own baseline (Δχ² = +13.3 DR1, +33.2 DR2), partly via H(z=0) drifting to 69.6. θ* is 1.0398 (4σ low). Not θ*-consistent; a proper normalisation is needed before any claim. Larger ε fails badly.
- SIM129: every γ > 0 worsens the fit monotonically (Δχ² = −16 … −560 DR1). Published "N/A".
- Original scripts annotated with a header pointing to the audit evaluator.

### SIM150 — Phase 5 Λ0 sweep (not rerun)
- The verdict INCOMPATIBLE (for Λ0_required = 0.154) is robust on several independent grounds: Cassini 2e4×, fσ8, UV.
- The BAO component used the non-DESI table but does not change the verdict.
- **Important: SIM150's own docstring computes, at the LOCKED Λ0 = 0.003, ω_BD ≈ 2109 and |γ−1| ≈ 9.5e-4, writes "this exceeds 2.3×10⁻⁵! Something is wrong", then moves on.** That is the same violation the audited SIM107 finds (9.48e-4, 41×). The project's own later code had caught it.

### SIM120 — "Not run (spec only)" per its results.md
- Not rerun. It contains the non-DESI table and would need the same treatment if it is ever run.

What this means for the "2.77σ structural floor":
- The floor is not a property of CMSTG. The Phase 1 frozen model is observationally identical to ΛCDM for BAO + θ*.
- ΛCDM itself sits at 2.08σ (DR1) / 3.07σ (DR2) from DESI in this fixed-ω setup.


## 7. Stale published figures

68 published figure files (PNG and PDF) in rerun sims were not regenerated; they show pre-audit results.
Each affected sim folder has a `STALE_FIGURES.md` listing them. Nothing was deleted (they are still in git); removing or regenerating them is the user's call.

Affected sims (file count): SIM121B (4), SIM122 (8), SIM123 (8), SIM124 (6), SIM126 (2), SIM129 (2), SIM130 (6), SIM131 (2), SIM132 (4), SIM133 (4), SIM134 (4), SIM135 (4), SIM136 (4), SIM138 (2), SIM142 (2), SIM143 (2), SIM144 (4).
