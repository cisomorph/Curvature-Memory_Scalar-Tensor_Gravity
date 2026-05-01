# SIM138 — DESI Y1 Per-Bin Sensitivity Decomposition

**Tier:** 1 (Diagnostic)
**Depends on:** SIM121C joint MCMC outputs (χ²_DESI = 18.26, tension 2.77σ)
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VII §2

---

## Question

Which of the 6 DESI Y1 redshift bins drives the 2.77σ tension floor? Is the tension localized (→ favours a sharp transition mechanism) or distributed (→ favours smooth evolution)?

## Motivation

Paper VII established that Phase 3 curvature-sourced scalars cannot reduce the DESI floor. The no-go theorem is global — it does not tell us *where in z* the tension lives. Phase 4 mechanism design depends on this:

- Localized tension at z ~ 0.5–1.0 → SIM140 symmetron-type transition is the right tool
- Smooth tension across all bins → SIM142 Galileon tracker is the right tool
- Tension concentrated at a single bin → possible DESI systematic, not a theoretical problem

## Action spec

No new physics. Decompose existing χ² against Phase 1 canonical prediction.

## Inputs

- Phase 1 canonical parameters: H₀ = 67.59, Ω_m = 0.312, Λ₀ = 0.003, F_eff = 0.521, r_d = 147.56 Mpc
- DESI Y1 BAO measurements at z ∈ {0.295, 0.510, 0.706, 0.930, 1.317, 2.330}
- DESI full covariance matrix (Adame et al. 2024, arXiv:2404.03002)
- CMSTG Friedmann solver from Phase 1 (used in SIM90, SIM121C)

## Procedure

1. Reproduce SIM121C χ²_DESI = 18.26 as a sanity check (must match within 0.01).
2. Compute model prediction D_M/r_d, D_H/r_d (and D_V/r_d where applicable) at each of the 6 DESI z-bins for:
   - (a) Phase 1 canonical CMSTG
   - (b) ΛCDM best-fit to same data (reference)
3. Compute per-bin residuals: Δ_i = model_i − observation_i for each measurement
4. Compute per-bin χ² contribution:
   - Diagonal: Δ_i² / σ_i²
   - Full: diagonal sum to get a "naive" contribution, plus reconstruct full quadratic form so the sum equals the total χ²
5. Rank bins by Δχ² contribution.
6. Identify the dominant bins (those accounting for >50% of the total χ² above what ΛCDM contributes).

## Success criteria

One of the following must be produced:

- **LOCALIZED:** one or two z-bins account for ≥70% of the CMSTG-excess χ². → SIM140 prioritized, transition z set to bracket the offending bin(s).
- **DISTRIBUTED:** all 6 bins contribute within a factor of 2 of each other. → SIM142 prioritized.
- **SINGLE_BIN_DOMINANT:** one bin accounts for >80% of the excess. → Flag as possible DESI systematic; cross-check with BOSS+eBOSS (used in SIM87) to see if the same z-range is discrepant.

## Failure modes to watch

- Covariance handling: the bins are correlated. The off-diagonal structure matters. Do NOT sum diagonal-only contributions and treat as independent.
- r_d degeneracy: CMSTG and ΛCDM are r_d-degenerate to 21 ppm (SIM101). Make sure the comparison holds r_d fixed from the SIM121C joint fit, not recomputed.
- Sign convention: Δ_i = model − obs. Document sign clearly — confusion here inverts the diagnosis.

## Deliverables

- `output.json` with:
  ```json
  {
    "per_bin_residuals": { "z_0.295": ..., ... },
    "per_bin_chi2_contribution": { ... },
    "total_chi2": 18.26,
    "total_chi2_reconstruction_check": ...,
    "verdict": "LOCALIZED" | "DISTRIBUTED" | "SINGLE_BIN_DOMINANT",
    "dominant_bins": [ ... ],
    "recommended_tier2_priority": "SIM140" | "SIM142" | "defer"
  }
  ```
- `RESULT.md` with per-bin table, figure of residuals vs z, and Tier 2 recommendation
- Figure: model vs DESI observations with error bars, both CMSTG and ΛCDM overlaid

## Estimated time

Small–medium. The solver already exists in Phase 1 code; the work is covariance-correct χ² decomposition and plotting.
