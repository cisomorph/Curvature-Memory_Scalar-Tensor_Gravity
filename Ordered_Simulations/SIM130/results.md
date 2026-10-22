# SIM130 — Results

**Phase:** Phase 3  
**Title:** Curvature-Sourced Ψ Evolution  
**Verdict:** FAIL

## What was tested

Ψ driven by curvature source; FAIL — M-field programme closed.

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** faithful FAIL / θ*-matched PASS (published: FAIL) — but Cassini ~50×

- Ψ_ini = 2.62: Ψ0 = 2.88, H0 = 66.64, Δχ² = +5.43 (DR1) / +13.94 (DR2). Robust across three integrators and z_ini = 1e5–1e7.
- Ψ_ini = 3.5: Δχ² = +7.89 / +19.68.
- Cassini 49–85×; Ġ/G ≈ −1.2e-13 /yr (about 2.5σ vs LLR).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
