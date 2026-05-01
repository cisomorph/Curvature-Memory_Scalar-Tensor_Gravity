# SIM146 — Distinctive Prediction Extraction

**Tier:** 3 (Validation)
**Depends on:** SIM144 complete, SIM145 at least partial PASS
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, SIM144 RESULT, winner's RESULT

---

## Question

What observational predictions does Phase 4 CMSTG (Phase 1 + winning mechanism) make that differ from ΛCDM by more than current precision, and which can be tested decisively by ongoing or imminent surveys (DESI Y3, Euclid, CMB-S4, LiteBIRD)?

## Motivation

SIM144 will likely show CMSTG-Phase 4 at a competitive but not dominant ln B. The Bayesian evidence alone does not move forward the scientific status of a theory; decisive observation does.

A theory must make predictions that future data can confirm or refute. Until Phase 4 CMSTG has such predictions, its status is "consistent with data" — not "a real physical theory." SIM146 extracts those predictions.

## Target observables

For each, compute prediction and compare to:
- (i) ΛCDM best-fit
- (ii) Current best measurement with error bar
- (iii) Forecast sensitivity of the next-generation survey

### Cosmological

1. **w(z) shape at z = 0.5 − 2.5** (DESI Y3, DESI Y5)
   - Compute w(z) from the winning mechanism's Ω_DE(z)
   - Compare to ΛCDM constant w = −1
   - Is the difference larger than DESI Y3 projected σ(w_a) ≈ 0.1?

2. **H(z) at z = 1 − 2** (Euclid)
   - Absolute H(z) values
   - Euclid will measure H(z) to ~1% at z ~ 1.5
   - Is the CMSTG-ΛCDM difference above this precision?

3. **σ₈(z) evolution at z = 0 − 2** (Euclid, LSST)
   - Growth rate
   - KiDS and DES measure this at ~1% today; Euclid will reach 0.3%
   - Does the winning mechanism predict a distinct σ₈(z) slope?

4. **CMB spectral distortions (μ and y)** (PIXIE or successor, probably not near-term but worth noting)
   - Phase 1 predicted <2.93% RMS deviation in TT
   - Does Phase 4 predict distinctive μ or y distortions?

5. **CMB B-mode tensor-to-scalar ratio r** (LiteBIRD, CMB-S4)
   - Phase 1 inflation sim (SIM110) failed slow-roll without ξ ~ 10⁴ non-minimal coupling
   - Does Phase 4 mechanism change this? If winner was a late-time mechanism, probably not — but verify.

### Astrophysical (if SIM137 found a structural pattern for χ-DM)

6. **SPARC χ-DM predictions** — for the galaxy subclass identified in SIM137, make specific m₂₂ predictions that can be tested by deeper rotation-curve surveys or by sub-galactic probes (MW dwarf dynamics).

7. **Small-scale structure cutoff** — if χ-DM is real at m₂₂ ≈ 0.28, the de Broglie wavelength implies a characteristic cutoff in the matter power spectrum. This is testable via Lyα forest (already partially constrains SIM121), 21cm surveys, or satellite galaxy counts.

### Gravitational

8. **GW propagation at z > 1** (LISA, ET, CE)
   - Phase 1 passed c_T = c at GW170817 precision (z ≈ 0.009)
   - At higher z, does the winning mechanism predict distinguishable luminosity-distance offsets?
   - Standard siren cosmology sensitivity is the comparison point.

9. **Modified gravity growth parameter μ and slip η** (Euclid, LSST joint analyses)
   - If F_eff(z) evolves in Phase 4 CMSTG, μ(z, k) and η(z, k) differ from ΛCDM
   - Euclid sensitivity: σ(μ) ~ 0.05 today, improving

## Procedure

1. For each observable above, compute the Phase 4 CMSTG prediction using the SIM144 best-fit (or posterior mean with uncertainty).
2. Compute the ΛCDM prediction for the same observable.
3. Compute Δ(observable) = CMSTG − ΛCDM.
4. Compare to:
   - Current best σ(observable) — tells you if already testable
   - Projected σ(observable) from upcoming surveys — tells you if near-term testable
5. Compute detection S/N = |Δ(observable)| / σ_projected for each.
6. Rank predictions by S/N. Predictions with S/N > 3 are the falsification targets.

## Success criteria

- **STRONG_PASS:** at least 2 observables with S/N > 5 testable within 3 years. → CMSTG Phase 4 is a genuine predictive theory. Paper VIII can make publication-grade claims.
- **MODERATE_PASS:** at least 1 observable with S/N > 3 within 5 years. → Paper VIII proceeds with specific falsifiability claim; continued monitoring of survey data is the exit strategy.
- **WEAK_PASS:** predictions exist but all are below current or projected precision. → Paper VIII states the theory is consistent but not uniquely testable; CMSTG's status remains "compatible with data" rather than "predictively confirmed."

Document verdict honestly. A WEAK_PASS is not a failure — it is the correct scientific statement about where Phase 4 CMSTG stands.

## Failure modes to watch

- **Optimism bias:** it's easy to cherry-pick parameter values where CMSTG differs maximally from ΛCDM. Use the MCMC posterior, not the best-fit, to get realistic prediction bands.
- **Apparent vs real:** some "differences" are artifacts of fitting to the same data. Only count genuinely independent predictions.
- **Survey timeline:** DESI Y3 ~ 2026-27, Y5 ~ 2028, Euclid first data ~ 2025-26 with full survey ~ 2030, CMB-S4 ~ 2030, LiteBIRD ~ 2032. Check current timelines via web search before publishing, they shift.
- **Systematic floor:** some observables are systematics-limited, not statistics-limited. σ_projected in such cases should reflect the systematic, not Poisson.

## Deliverables

- `output.json` with per-observable Δ, σ, S/N, timeline
- `RESULT.md` with ranked prediction table, 2-3 highest-S/N predictions discussed in detail, detection strategy
- Figures for top 2-3 predictions: CMSTG vs ΛCDM with current data + projected error bars
- Paper VIII "Falsifiable Predictions" section draft

## Estimated time

Medium. Mostly standard cosmology forecasting against published Fisher matrices and survey specifications.
