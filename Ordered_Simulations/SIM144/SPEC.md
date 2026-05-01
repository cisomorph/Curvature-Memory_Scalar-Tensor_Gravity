# SIM144 — Joint CMB+BAO+RSD MCMC on Winning Mechanism

**Tier:** 3 (Validation)
**Depends on:** At least one of SIM140–143 produced a PASS candidate
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, the winner's SPEC and RESULT, SIM121C infrastructure

---

## Question

When the Phase 4 winning mechanism is varied simultaneously with the standard cosmological parameters, does it still reduce the DESI tension below 2σ, and what are the joint parameter posteriors?

## Motivation

Tier 2 sims fix most parameters at Phase 1 canonical values and only scan the mechanism-specific parameters. SIM144 opens the full parameter space. Often, parameter shifts can reveal (a) the mechanism fitting data via degeneracies rather than genuine new physics, or (b) tensions that didn't appear in the restricted scan.

This is the analog of SIM121C for Phase 4.

## Action spec

Whichever of SIM140–143 produced PASS. Use the full Lagrangian from that sim's SPEC.

## Parameters (varied in MCMC)

Standard:
- H₀
- Ω_m h²
- Ω_b h²
- n_s
- A_s
- τ

CMSTG-specific (from winning mechanism):
- If SIM140: V₀, Ψ_c, σ
- If SIM141: β (or γ, a_trans, σ depending on form)
- If SIM142: c₃
- If SIM143: U₀, λ (or n)

Plus any derived quantities to track: w_0, w_a, F_eff(z=0), σ₈, S₈.

## Likelihoods

- Planck 2018 plikHM TTTEEE (via `clipy`, as in SIM97/98)
- Planck 2018 lowE
- DESI Y1 BAO with full covariance (arXiv:2404.03002)
- BOSS+eBOSS+Lyα BAO (for cross-check)
- 9-point fσ₈ RSD compilation
- Optional: Pantheon+ SNe (check if it sharpens constraints)

## Procedure

1. Port SIM121C MCMC infrastructure to accept the new Lagrangian.
2. Verify that at ΛCDM-equivalent point (mechanism parameters off), the MCMC reproduces Planck + DESI ΛCDM best-fits within known tolerances.
3. Run full joint MCMC. Use emcee or Cobaya. Target ~50k samples after burn-in, R < 1.02.
4. Diagnostic convergence: Gelman-Rubin statistic, trace plots, autocorrelation time.
5. Extract:
   - Best-fit parameters
   - 1D and 2D marginalized posteriors
   - ΔlnZ via nested sampling (MultiNest or dynesty) for Bayesian evidence vs Phase 1 canonical and vs ΛCDM
   - Tension with DESI: compute (χ²_joint − χ²_DESI_alone) / √(2 × N_DESI)
6. Compare to reference:
   - Phase 1 canonical: ln B = −0.71 (inconclusive), tension 2.77σ
   - Target: tension < 2σ AND ln B > 0 (CMSTG mildly preferred)

## Success criteria

- (a) MCMC converged: R < 1.02, ESS > 10× parameter count
- (b) DESI tension < 2σ at best-fit
- (c) Bayesian evidence ln B(winner vs Phase 1 canonical) > 1 (positive preference)
- (d) No parameter at boundary of prior (would indicate tuning pressure)
- (e) Standard parameters (H₀, Ω_m) consistent with Phase 1 canonical within 1σ — large shifts would indicate the mechanism is fitting via reparametrization
- (f) H₀ tension with SH0ES does NOT worsen beyond Phase 1 levels

All (a)-(f): SIM144 PASS → Paper VIII winning-mechanism paper can be drafted.
Partial pass: document which criteria failed and interpret accordingly.

## Failure modes to watch

- **Fitting via H₀ shift:** mechanisms that reduce DESI tension by shifting H₀ to higher values often worsen SH0ES tension in compensation. Check both.
- **Fitting via reparametrization:** if the "winner" just reproduces ΛCDM at a shifted Ω_m, it's not new physics. Check that the mechanism parameters are non-zero at >3σ.
- **Prior dependence:** Bayesian evidence is prior-sensitive. Run with 2–3 prior choices on the mechanism parameters. Report range.
- **Nested sampling failures:** high-dim MultiNest runs can miss modes. Cross-check lnZ with harmonic mean estimator or thermodynamic integration.
- **Nuisance parameters:** Planck plikHM has many nuisance parameters (foregrounds, calibration). Either marginalize or use compressed likelihood. SIM98 used full nuisance set.

## Deliverables

- `output.json` with best-fit, posterior summaries, lnZ, tension metrics
- `RESULT.md` with:
  - Convergence diagnostics
  - Best-fit table (all parameters)
  - Key 2D posterior figures (especially mechanism param vs H₀, vs Ω_m, vs σ₈)
  - Bayesian evidence comparison
  - Tension reduction statement
- Triangle plot (corner.py) of full posterior
- Draft Paper VIII results table

## Estimated time

Large. Joint MCMC with new Lagrangian is computationally expensive. Plan for 2–4 weeks of walltime depending on compute.
