# SIM133 — Results

**Phase:** Phase 3  
**Title:** Phase 3 Direction D: Gauss-Bonnet coupling  
**Verdict:** FAIL

## What was tested

αΨ G coupling to Gauss-Bonnet invariant G; FAIL — structural.

## Key numerical results

- Λ₀ = 0.003

## Action tested

```
S=∫d⁴x√g[½(1+2Λ₀Ψ²)R - ½(∂Ψ)² + α̃·Ψ·𝒢] + S_SM
```

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** FAIL (published: FAIL)

- With the derived Gauss–Bonnet Friedmann term (8α_GBΨ'H⁴), the constraint has no physical root for any α̃ ≥ 1e-8.
- 𝒢 = −24H⁴ in the radiation era drives a runaway.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
