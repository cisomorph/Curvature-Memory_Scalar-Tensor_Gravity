# CMSTG roadmap (locked 2026-10-22)

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
- [~] 1.2 Continuous integration (`.github/workflows/tests.yml`): GitHub Actions runs the tests on every push to `audit-fixes`.
- [x] 1.3 Shared code (policy; all audited sims already comply): all new sims import `Ordered_Simulations/common/`; constants and data live there only.

## Phase 2 — Attractor models (SIM132 / SIM134 families), background level → **SIM154**
Goal: map where the Solar-System-safe attractor models live, and whether they survive BBN.

- [ ] 2.1 Wide scan, θ*-matched, DESI DR1 + DR2:
  - subtractive F = ½(1 − 2Λ0Ψ²) + ξΨ
  - exponential F = ½e^{−2Λ0Ψ²} + ξΨ
  - Λ0 from 0.02 to 0.5, ξ from 1/24 to 1/4, Ψ_ini = 0
- [ ] 2.2 Local gravity at z = 0: Cassini γ and LLR Ġ/G.
- [ ] 2.3 BBN: speed-up S = √(G_BBN/G0) → Y_p and D/H (Kneller & Steigman 2004 fits), against observations.
- [ ] 2.4 Verdict per point; list the survivors.

**Pre-registered survival criteria for a parameter point (all must hold):**

| # | test | criterion |
|---|---|---|
| C1 | cosmology (θ*-matched) | DESI DR2 GOF < 2σ **and** Δχ²(DR2) vs ΛCDM ≥ +4 |
| C2 | Solar System | \|γ − 1\| < 2.3e-5 (Cassini) **and** \|Ġ/G − 7.1e-14\| ≤ 3 × 7.6e-14 /yr (LLR, Hofmann & Müller 2018) |
| C3 | BBN | predicted Y_p and D/H each within 2σ of Y_p = 0.245 ± 0.003 and D/H = (2.547 ± 0.025)×10⁻⁵, with fit uncertainties added in quadrature (±0.0006 for Y_p, ±3% for D/H) |
| C4 | stability | F > 0 throughout; a canonical scalar is then ghost-free |

Points passing C1–C4 go to Phase 3. If none pass, the attractor route is closed at background level and Phase 3 is not run for it.

## Phase 3 — Full CMB + BAO + supernovae (professional tools)
- [ ] 3.1 Environment: Python venv; Cobaya; hi_class built from source; likelihoods:
  - Planck 2018 plik-lite TTTEEE, low-ℓ, lensing
  - DESI DR2 BAO
  - Pantheon+ / Union3 / DES-Y5 SNe
- [ ] 3.2 Validation: hi_class ΛCDM reproduces Planck 2018 best fit; Cobaya reproduces published ΛCDM χ² and parameters.
- [ ] 3.3 Implement the surviving model(s) in hi_class as a covariant Horndeski model (G2 = X, G4 = F(φ), with bare Λ); validate its background against `st_solver`.
- [ ] 3.4 Best-fit and MCMC comparison: ΛCDM vs w0waCDM vs model(s). Report Δχ², AIC/BIC and posteriors.
  - Criteria to be fixed in this file before 3.4 runs.

## Phase 4 — Theory
- [ ] 4.1 Can the attractor form F'(Ψ*) = 0 be derived from the curvature-memory idea rather than chosen?
- [ ] 4.2 Relation to the Damour–Nordtvedt mechanism (GR as a cosmological attractor, PRL 70, 2217, 1993).
- [ ] 4.3 Can the bare Λ be removed?
- [ ] 4.4 Status of Paper III's UV claims (not audited yet).

## Phase 5 — Distinctive predictions
- [ ] 5.1 Gravitational-wave standard sirens: d_GW/d_EM = √(F(0)/F(z)).
- [ ] 5.2 G(z) history as seen by BBN and the CMB.
- [ ] 5.3 Growth of structure: fσ8(z), S8.

## Phase 6 — Papers
- [ ] 6.1 Paper II: Cassini and "2.77σ" claims.
- [ ] 6.2 Paper I: basis of Theorems 1 and 2.
- [ ] 6.3 README: "no bare Λ", "G_eff/G = 1 − 15 ppm", the DESI "Y3 ~2027" wording.
- [ ] 6.4 If a model survives Phase 3: one focused paper.

## Phase 7 — Independent review
- [ ] 7.1 Share the audit, branch and code with a scalar-tensor cosmologist; my fixes need checking too.

---

## Log
- 2026-10-22: roadmap locked at the user's request. Starting Phase 1, then Phase 2 (SIM154).
