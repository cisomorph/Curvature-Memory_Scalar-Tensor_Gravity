# SIM137 — SPARC Failure-Mode Analysis

**Tier:** 1 (Diagnostic)
**Depends on:** SIM119 outputs (65/161 PASS, χ-field DM fit)
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md

---

## Question

Why did 96 of 161 SPARC galaxies fail the χ-DM fit in SIM119? Is the failure pattern structural (physical signal) or random (incomplete mechanism)?

## Motivation

SIM119 is one of only two genuine positive results from the Phase 2 programme (alongside SIM128 S₈ suppression). A 65/161 pass rate is informative but not universal. Before committing Phase 4 effort to mechanism design, we need to know whether χ-DM is a real mechanism applied to the wrong galaxies, or a lucky fit on a subset with no structural basis.

## Action spec

No new physics. Re-analyze SIM119 outputs.

## Inputs

- SIM119 per-galaxy fit results (χ², pass/fail flag, best-fit m₂₂, residuals)
- SPARC catalog metadata (Lelli et al. 2016): morphological type, HI mass, stellar mass, effective radius, mean surface brightness, distance, inclination, v_flat
- Baryonic Tully-Fisher relation parameters from SPARC

## Procedure

1. Load SIM119 outputs. Confirm 65 PASS / 96 FAIL / 0 missing.
2. For each of the 161 galaxies, construct a feature vector:
   - Morphological type (Sa–Sm, dIrr, etc.)
   - log₁₀(M_HI/M_☉), log₁₀(M_*/M_☉), log₁₀(M_bary/M_☉)
   - Mean effective surface brightness μ_eff
   - v_flat (observed)
   - R_eff (effective radius)
   - Inclination
   - Distance
3. For each feature, compute:
   - Pass/fail distribution (histogram, KS test between pass and fail subsamples)
   - Correlation with SIM119 χ² residual
4. Multivariate: logistic regression of pass/fail against all features. Report coefficients with p-values.
5. If a structural pattern emerges, identify the subpopulation where χ-DM works. If not, report feature distributions overlap.

## Success criteria

- **Structural pattern found:** one or more features show KS p < 0.01 between pass and fail, or logistic regression finds a variable with |β| > 2σ significance. → χ-DM is a genuine mechanism for a physical subclass. Tier 2 work proceeds with that subclass in mind for any matter-sector coupling tests.
- **No pattern found:** all features overlap (KS p > 0.05), logistic regression uninformative. → χ-DM mechanism is incomplete. Flag for Phase 5 reconsideration; does not block Tier 2.

Both outcomes are acceptable conclusions. What would not be acceptable is a run that does not produce a clear verdict.

## Failure modes to watch

- Inclination bias: poorly inclined (i < 30° or i > 80°) galaxies have unreliable rotation curves. Check that the SIM119 sample was already cut for this; if not, redo the analysis excluding them.
- Distance uncertainty: galaxies with >20% distance error should be flagged separately
- m₂₂ tension with Lyα: SIM121 required f_FDM < 0.14. Check if the 65 PASS galaxies cluster at m₂₂ values consistent with this bound.

## Deliverables

- `output.json` per schema in RESEARCH_RULES.md §3.2
- `RESULT.md` with:
  - Verdict: STRUCTURAL_PATTERN | NO_PATTERN | INCONCLUSIVE
  - Tables of feature distributions (pass vs fail)
  - Logistic regression coefficients and p-values
  - Recommendation for Tier 2 scope
- Figures: pass/fail feature histograms (one per feature), correlation matrix
- Update `PHASE4_ROADMAP.md` SIM137 status

## Estimated time

Small. Pure data analysis, no new solvers.
