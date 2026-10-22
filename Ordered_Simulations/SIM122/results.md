# SIM122 — Results

**Phase:** Phase 2  
**Title:** Phase 3 Direction A: Ψ² coupling  
**Verdict:** FAIL

## What was tested

Non-minimal Ψ² curvature coupling; FAIL — exceeds GR bound.

## Key numerical results

- DESI tension = 2.462

## Structural diagnosis

Posterior runs to Λ₀ upper boundary (0.010). CMB θ_* forces F₀→0.569, H₀→76, systematically overshooting DESI H(z) ~10%. w₀/wₐ identical to SIM121C (freezing). Unlocking Λ₀ does not resolve F₀ normalisation problem — the tension is deeper than the Λ₀ lock.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** PARTIAL — Λ0 unconstrained (published: FAIL)

- The χ² profile over Λ0 ∈ [1e-4, 0.01] is flat to 1.5e-8 (DR1) and 5e-8 (DR2). The Λ0 posterior width equals the prior width.
- The fit gain is entirely (w0, wa): Δχ² = +10.53 (DR1), +23.78 (DR2) vs ΛCDM.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
