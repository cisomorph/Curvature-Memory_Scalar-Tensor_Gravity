# SIM144 — Results

**Phase:** Phase 4  
**Title:** Tier 3: Joint MCMC on winner (not run)  
**Verdict:** Not run (spec only)

## What was planned

Would run if Tier 2 produced a PASS candidate; deferred — all Tier 2 failed.

## Status

This sim was not executed. See SPEC.md for the planned mechanism and success criteria.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** faithful FAIL / θ*-matched PASS (published: FAIL/STRUCTURAL) — Cassini 66–135×

- β_m sign error fixed (energy −2β_mφρ_m).
- ξ_R = β_m = 0.01: Δχ² = +9.07 (DR1) / +19.94 (DR2). Every cosmology pass fails Cassini.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
