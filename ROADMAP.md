# CMSTG roadmap (locked 2026-10-22)

The detailed, step-by-step record of this work is `RESEARCH_LOG.md`.

This file is the plan agreed after the audit (`AUDIT_NOTES.md`). Pass/fail criteria are fixed here
**before** the runs, following the project's own SPEC-before-run rule. Changing a criterion afterwards
requires a dated note explaining why.

Status key: `[ ]` not started, `[~]` in progress, `[x]` done, `[!]` blocked.

---

## Phase 1 — Lock in the corrected pipeline
Goal: no result gets a verdict unless the code first reproduces known answers.

- [x] 1.1 Automated tests (`tests/`, pytest): 16 tests, all passing locally in about 8 s.
  - ΛCDM: DESI DR1/DR2 BAO-only fits reproduce the DESI papers (Ω_m within 0.25σ).
  - Planck: H0 solved from θ* matches Planck within 1σ.
  - GR limit: `st_solver` and `mf_solver` reproduce ΛCDM to 1e-10.
  - ε consistency: constraint-ε equals Raychaudhuri ε (1e-8).
  - Frozen-F0 degeneracy: χ²(F0 = 0.5206) = χ²(ΛCDM).
  - Cassini formula: Brans–Dicke limit γ − 1 = −1/(ω + 2).
  - Statistics helpers.
- [x] 1.2 Continuous integration (green on GitHub) (`.github/workflows/tests.yml`): GitHub Actions runs the tests on every push to `audit-fixes`.
- [x] 1.3 Shared code (policy; all audited sims already comply): all new sims import `Ordered_Simulations/common/`; constants and data live there only.

## Phase 2 — Attractor models (SIM132 / SIM134 families), background level → **SIM154**
Goal: map where the Solar-System-safe attractor models live, and whether they survive BBN.

- [x] 2.1 Wide scan, θ*-matched, DESI DR1 + DR2:
  - subtractive F = ½(1 − 2Λ0Ψ²) + ξΨ
  - exponential F = ½e^{−2Λ0Ψ²} + ξΨ
  - Λ0 from 0.02 to 0.5, ξ from 1/24 to 1/4, Ψ_ini = 0
- [x] 2.2 Local gravity at z = 0: Cassini γ and LLR Ġ/G.
- [x] 2.3 BBN: speed-up S = √(G_BBN/G0) → Y_p and D/H (Kneller & Steigman 2004 fits), against observations.
- [x] 2.4 Verdict per point; list the survivors. **Result (SIM154): 2 survivors**, family A and B at Λ0 = 0.15, ξ = 1/12. Δχ² vs ΛCDM = +7.0 / +17.7 (A) and +7.1 / +18.1 (B) for DR1 / DR2; BBN is the binding constraint (Y_p +1.8σ). Exploratory follow-up: map the survivor band finely and confirm BBN with a BBN code.

**Pre-registered survival criteria for a parameter point (all must hold):**

| # | test | criterion |
|---|---|---|
| C1 | cosmology (θ*-matched) | DESI DR2 GOF < 2σ **and** Δχ²(DR2) vs ΛCDM ≥ +4 |
| C2 | Solar System | \|γ − 1\| < 2.3e-5 (Cassini) **and** \|Ġ/G − 7.1e-14\| ≤ 3 × 7.6e-14 /yr (LLR, Hofmann & Müller 2018) |
| C3 | BBN | predicted Y_p and D/H each within 2σ of Y_p = 0.245 ± 0.003 and D/H = (2.547 ± 0.025)×10⁻⁵, with fit uncertainties added in quadrature (±0.0006 for Y_p, ±3% for D/H) |
| C4 | stability | F > 0 throughout; a canonical scalar is then ghost-free |

Points passing C1–C4 go to Phase 3. If none pass, the attractor route is closed at background level and Phase 3 is not run for it.

## Phase 3 — Full CMB + BAO + supernovae (professional tools)
- [x] 3.1 Environment (venv; Cobaya 3.6.2; hi_class v3.4 + Python wrapper; native Planck, DESI DR1/DR2, Pantheon+/Union3/DES-Y5 data): Python venv; Cobaya; hi_class built from source; likelihoods:
  - Planck 2018 plik-lite TTTEEE, low-ℓ, lensing
  - DESI DR2 BAO
  - Pantheon+ / Union3 / DES-Y5 SNe
- [x] 3.2 Validation (Planck-point χ²: low-ℓ TT 23.45, low-ℓ EE 396.12, plik-lite 583.16, lensing 8.86, matching Planck's published values): hi_class ΛCDM reproduces Planck 2018 best fit; Cobaya reproduces published ΛCDM χ² and parameters.
- [x] 3.3 Implement the surviving model(s) in hi_class as a covariant Horndeski model (G2 = X, G4 = F(φ), with bare Λ); validate its background against `st_solver`.
- [x] 3.4 Best-fit and MCMC comparison: ΛCDM vs w0waCDM vs model(s). Report Δχ², AIC/BIC and posteriors.
  - **Result 2026-10-24 (SIM155/RESULT.md): CMSTG-A FAILS the pre-registered criteria.**
    - Δχ² vs ΛCDM = +4.65 (D1), +4.50 (D2), +3.19 (D3), all < 6.
    - At the best fit C3 fails (Y_p +3.25σ).
    - w0waCDM gains +9.2 … +19.6.
    - MCMC not run (pre-registered: only if a model passes).
  - **Pre-registered 2026-10-23, before any Phase 3 fit was run:**
    - Data, baseline D1: Planck 2018 (low-ℓ TT + low-ℓ EE + plik-lite TTTEEE + lensing, native Cobaya likelihoods)
      + DESI DR2 BAO + Pantheon+. Variants D2 and D3 swap Pantheon+ for Union3 and DES-Y5.
    - Models:
      - ΛCDM (6 parameters)
      - w0waCDM (8)
      - CMSTG attractor family A (8: + Λ0, ξ; φ_ini = 0 at a = 1e-10, on the radiation-era attractor, with bare Λ)
    - Same code for all: hi_class + our adapter. Same helium treatment: YHe = 0.2454 + 0.438(ω_b − 0.02238) + 0.16(S − 1),
      with S = √(G_BBN/G0) = √(2F0); S = 1 for ΛCDM and w0wa.
    - Method: best fit (maximum likelihood) from at least 3 starting points per model; MCMC only if a model passes.
    - **CMSTG passes Phase 3** if, on D1:
      - (a) its best fit beats ΛCDM by Δχ² ≥ 6 (2 extra parameters);
      - (b) the best-fit point still passes Phase 2's C2 (Cassini, LLR) and C3 (BBN);
      - (c) Δχ² vs ΛCDM ≥ 6 also on at least one of D2/D3.
    - Reported regardless: Δχ² and ΔAIC (= −Δχ² + 2Δk) vs ΛCDM and vs w0waCDM.
    - "Competitive with w0waCDM" means a χ² within 3 of w0waCDM's.

## Phase 4 — Theory
- [~] 4.1 (notes in THEORY_NOTES.md) Can the attractor form F'(Ψ*) = 0 be derived from the curvature-memory idea rather than chosen?
- [x] 4.2 (THEORY_NOTES.md: CMSTG attractors are a DN-class model; |γ−1| ≈ 1e-7 as DN predicts) Relation to the Damour–Nordtvedt mechanism (GR as a cosmological attractor, PRL 70, 2217, 1993).
- [~] 4.3 (THEORY_NOTES.md: not by any audited model so far) Can the bare Λ be removed?
- [ ] 4.4 Status of Paper III's UV claims (not audited yet).

## Phase 5 — Distinctive predictions
- [x] 5.1 (survivor: |d_GW/d_EM − 1| < 1e-6 for z < 10, so sirens cannot test it) Gravitational-wave standard sirens: d_GW/d_EM = √(F(0)/F(z)).
- [x] 5.2 G(z) history as seen by BBN and the CMB. At the SIM155 best fit: G_BBN/G0 = 1.080, G_rec/G0 = 1.0025, G = G0 for z ≲ 10. At the best BBN-allowed point: 1.023 and 1.005. The modification is early (`SIM155/growth.json`).
- [x] 5.3 Growth of structure: fσ8(z), S8.
  - Survivor (background level): fσ8 −2.5% (z = 0), −1.3% (z = 0.5) vs ΛCDM; σ8 −0.55%.
  - **Superseded by the full Boltzmann result (SIM155):**
    - best fit: S8 = 0.871 (+7.5%), σ8 +3.9%, fσ8 +3.3 … +3.9% (z = 0–1);
    - BBN-allowed point: S8 +2.5%.
    - The background-level estimate held the other parameters at fixed Planck values. With every parameter refit, growth is enhanced, not suppressed.
    - At ΔF ≈ 0.011, close to the old survivor: σ8 +1.9%, fσ8 +1.2% (z = 0) and +1.6% (z = 0.5).

## Phase 6 — Papers
- [~] 6.1 (draft in PAPER_CORRECTIONS.md) Paper II: Cassini and "2.77σ" claims.
- [~] 6.2 (draft in PAPER_CORRECTIONS.md) Paper I: basis of Theorems 1 and 2.
- [~] 6.3 (draft in PAPER_CORRECTIONS.md) README: "no bare Λ", "G_eff/G = 1 − 15 ppm", the DESI "Y3 ~2027" wording.
- [x] 6.4 If a model survives Phase 3: one focused paper. **Not applicable**: no model survived (SIM155). PAPER_CORRECTIONS.md records what the papers should say instead.

## Phase 7 — Independent review
- [ ] 7.1 Share the audit, branch and code with a scalar-tensor cosmologist; my fixes need checking too.

## Phase 8 — Screened curvature memory (exploratory; added 2026-10-24 at the user's request)
- [x] 8.1 Theory map (SCREENING_NOTES.md):
  - **Thin-shell screening.** The CMSTG attractor is a symmetron in the curvature. Screening the Sun and Milky Way needs Λ0 ≳ 8e4, after which the cosmological ΔF ≲ 3e-6 (WHK no-go): no-go.
  - **Vainshtein screening** (genuine cubic Galileon). Cassini is screened, but LLR caps |α_M(0)| ≲ 0.002, since Ġ is not screened. Late-time effects driven by the coupling are then ≲ 0.1% unless tuned or reduced to quintessence.
  - Ways out need a different theory (dark-sector-only coupling; screening of Ġ; tuned bands).
- [x] 8.2 SIM156, the quantitative check of the Vainshtein branch: "pinned-thawing" CMSTG vs ΛCDM, w0wa and thawing quintessence on compressed Planck + DESI DR2 + Pantheon+.
  - **Result 2026-10-24: FAIL, as predicted** (SIM156/RESULT.md). Validation V1–V3 passed (w0wa +7.68 vs SIM155's 9.19).
  - Within |α_M(0)| ≤ 0.002, PT = thawing quintessence = +2.77; the fit switches the coupling off.
  - Unconstrained, PT reaches +4.95, but at α_M(0) = −0.11 (56× LLR).
  - **Pre-registered (SIM156/SPEC.md, committed before any run):**
    - Validation V1–V3 (distance priors reproduced; Pantheon+ χ² matches SIM155 to 1.0; w0wa's Δχ² within ±3 of SIM155's 9.19).
    - PT passes if, with |α_M(0)| ≤ 0.002, Δχ² vs ΛCDM ≥ 6 and BBN (C3) passes.
    - Prediction: PT ≈ TQ < w0wa.
- [x] 8.3 Only if 8.2 passes: full hi_class fit (cubic Galileon + F(Ψ)R) with SIM155's machinery. **Not applicable** (8.2 failed).
- [~] 8.4 Open directions, each a different theory (SCREENING_NOTES §4). The user chose option 1 (dark-sector-only coupling) and option 2 (screening of Ġ) on 2026-10-24, so they continue in Phase 9.

## Phase 9 — Dark curvature memory and hiding the drift of G (added 2026-10-24 at the user's request)
- [x] 9.1 Option 2, theory/literature (SCREENING_NOTES §6). No GW-consistent theory hides Ġ for a universally coupled F(Ψ)R:
  - luminal Horndeski has G_N = 1/(16πG4(φ0(t))) inside Vainshtein (Kobayashi 2019);
  - the surviving DHOST theories keep a cosmologically set coupling (Hirano, Kobayashi & Yamauchi 2019);
  - thin-shell hides Ġ but leaves no cosmology (WHK).
  - Option 2 merges into option 1: coupling only the dark sector means baryons never feel a drift.
- [x] 9.2 Option 1, SIM157 "dark curvature memory": CMSTG's A = (2F)^{-1/2} acting on dark matter only, with a pinned-thawing trigger (DPT). Background (compressed Planck + DESI DR2 + Pantheon+) plus linear growth (eBOSS DR16 fσ8).
  - **Result 2026-10-24: FAIL (a)** (SIM157/RESULT.md).
    - Δχ²_bg = +2.88 (TQ +2.77, w0wa +7.68); growth fine; coupling off early.
    - The apparent phantom crossing is in the right direction but has the wrong shape (w_eff −0.96 → −1.22 at z = 1).
    - S8 +6%, from a dark fifth force ≈ 2× gravity.
  - **Pre-registered (SIM157/SPEC.md, committed before any run):** DPT passes if Δχ²_bg vs ΛCDM ≥ 6, growth no worse than ΛCDM by more than 4 in χ²_RSD, and β ≈ 0 before recombination.
  - Expectation: TQ < DPT < w0wa.
- [x] 9.3 Only if 9.2 passes: a Boltzmann implementation of the coupled perturbations (CLASS) and a full fit with SIM155's machinery. **Not applicable** (9.2 failed).



---

## Log
- 2026-10-22: roadmap locked at the user's request. Starting Phase 1, then Phase 2 (SIM154).
- 2026-10-22: Phase 1 done: 16 tests pass locally; CI workflow added.
- 2026-10-22: Phase 2 done: SIM154 run per its pre-registered SPEC. 2 of 70 points survive C1–C4. Next: Phase 3 environment (hi_class + Cobaya).
- 2026-10-23: Phase 3.3 done: `cmstg_attractor` implemented in hi_class and validated against ΛCDM (GR limit), regulator convergence, and the project's own solver (H to 2e-5, F to 1e-5). Patch and docs in `hi_class_model/`. Next: build the hi_class Python wrapper, install the likelihoods, validate Cobaya on ΛCDM.
- 2026-10-23: Phase 3.1–3.2 done (Cobaya + hi_class adapter `common/cobaya_hiclass.py`; Planck-point χ² match). SIM155 launched: 27 best fits (3 models × D1–D3 × 3 starts).
- 2026-10-23: Phase 4 and 6 drafts (THEORY_NOTES.md, PAPER_CORRECTIONS.md); Phase 5.1 and 5.3 predictions for the survivor.
- 2026-10-23: SIM155 v1 (Cobaya `minimize`/BOBYQA defaults) was stopped after D1.
  - Starting points disagreed by up to 2.3 (ΛCDM) and 5.5 (w0wa) in χ², i.e. early termination; the χ² noise floor is only about 0.005.
  - All CMSTG starts crashed: hi_class's c_s² test fails on turning-point 0/0 noise, and Cobaya passes −inf to BOBYQA.
  - Fixes (method only; criteria unchanged): `cs2_safe_smg` 1e-5 → 1e-2 (justified in hi_class_model/README.md), and a new two-stage driver `SIM155/fit.py` (BOBYQA with a finite penalty for failures, then repeated MIGRAD polish).
  - All models are rerun with v2 for a uniform comparison. The v1 results stay in `SIM155/fits/` for the record.
- 2026-10-23: v2 was too slow (a ΛCDM fit had not finished after 50 min; no v2 result exists) and was replaced by v3 (method only; criteria unchanged; details in RESEARCH_LOG.md):
  - The optimiser is `common/fitcore.py`: coarse BOBYQA, then rounds of Hessian whitening plus BOBYQA, with failure retries and a continuous prior wall. It is unit-tested in CI on synthetic problems at the measured noise.
  - Measured on the real likelihood: χ² noise 0.07 RMS (hi_class sampling at default precision), so best-fit χ² is good to about ±0.1. No CMSTG failures in 40 random points with the new c_s² tolerance.
  - All 27 fits relaunched 19:22 EDT.
- 2026-10-23 (evening): interim SIM155 results for D1. ΛCDM: seeds agree to 0.025 in −log(post); H0 = 68.30 ± 0.28. w0waCDM: Δχ² = +9.2 vs ΛCDM (w0 = −0.85, wa = −0.55). CMSTG-A: about +3.3 … +4.0 so far, but its seeds stop at different points of a flat, curved (Λ0, ξ, H0) valley along which the BBN-relevant ΔF = ξ²/(4Λ0) varies from 0.009 to 0.024. Method addition (criteria unchanged): a profile likelihood in ΔF fixes the best fit (`SIM155/profile.py`; RESEARCH_LOG decision 6).
- 2026-10-24: **Phase 3 done; SIM155 verdict FAIL** (pre-registered criteria).
  - The CMSTG attractor improves on ΛCDM by at most Δχ² = +4.65 (D1), only where BBN excludes it (Y_p +3.25σ). Where BBN allows it the gain is +3.25 (ΔAIC +0.75).
  - w0waCDM gains +9.2 … +19.6.
  - The model's G change happens before recombination, so it cannot mimic late-time dark energy, and it raises S8.
  - Phase 5.2/5.3 done at the full-CMB best fits.
  - THEORY_NOTES 4.5: structural tension between local gravity and late-time evolution; screening is the open question.
  - Remaining: 4.1, 4.3, 4.4 (theory), 7.1 (independent review).
- 2026-10-24: Phase 8 opened at the user's request ("can we explore that idea?" — screening).
  - 8.1 theory map done (SCREENING_NOTES.md; literature verified: WHK 2012; Babichev et al. 2011; Burrage & Dombrowski 2020; Creminelli et al. 2020; Chen, Huang & Wang 2019).
  - The project's SIM142 "Galileon" was a rescaled kinetic term with no screening; a genuine cubic Galileon was never tested.
  - SIM156 pre-registered.
- 2026-10-24: **SIM156 done: FAIL, as predicted.** Screened (Vainshtein) curvature memory within the LLR bound adds nothing to thawing quintessence (+2.77 vs ΛCDM; w0wa +7.68). Even ignoring LLR, the maximum is +4.95 at 56× the bound. Conclusion (SCREENING_NOTES): universally coupled F(Ψ)R curvature memory cannot produce the DESI-preferred dynamics, early (BBN) or late (LLR / WHK).
- 2026-10-24: Phase 9 opened ("i like option 1 and 2").
  - 9.1: option 2 is closed in theory; it merges into option 1 (SCREENING_NOTES §6).
  - DARK_COUPLING_NOTES.md written (model, equations, literature: Khoury–Lin–Trodden 2025 and 2026 coupled-DE papers).
  - eBOSS DR16 fσ8 read from the paper's Table III.
  - SIM157 pre-registered.
- 2026-10-24: **SIM157 done: FAIL (a).**
  - Dark curvature memory (coupling on dark matter only) adds +0.11 over thawing quintessence (background +2.88 vs w0wa +7.68).
  - It produces a right-direction but wrong-shape phantom crossing and raises S8 by 6%.
  - Across SIM155–157, every tested form of curvature memory tops out near Δχ² ≈ +3 (limits: BBN, LLR, growth), against +8…+9 for w0waCDM.
