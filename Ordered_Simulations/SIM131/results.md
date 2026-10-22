# SIM131 — Results

**Phase:** Phase 3  
**Title:** Phase 3 Direction C: ξΨR coupling (ξ=1/6)  
**Verdict:** FAIL

## What was tested

Conformal coupling F_eff=(1+2Λ₀Ψ²)/2+ξΨ; FAIL — θ* violation 5.97σ.

## Action tested

```
S = ∫d⁴x√g [F_eff(Ψ)·R − ½(∂Ψ)²] + S_SM, F_eff = (1+2Λ₀Ψ²)/2 + ξΨ
```

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** FAIL / FAIL (published: FAIL)

- For ξ ≥ 1/6 there is no θ* solution. Best case ξ = 1/20: Δχ² = −49 (DR1).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
