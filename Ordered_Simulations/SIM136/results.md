# SIM136 — Results

**Phase:** Phase 3  
**Title:** Horndeski G^{μν}∂Ψ² coupling  
**Verdict:** FAIL

## What was tested

G^{μν}∂_μΨ∂_νΨ kinetic Horndeski term; FAIL — G₃_X=0 equivalent to rescaled kinetic.

## Key numerical results

- Λ₀ = 0.003

## Action tested

```
S=∫d⁴x√g[(1+2Λ₀Ψ²)/2·R − ½(∂Ψ)² + α·G^{μν}∂_μΨ∂_νΨ]+S_SM
```

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** FAIL — c_T ≠ c (published: FAIL)

- Corrected Horndeski equations: energy ½Ψ̇²(1 + 18αH²); the field equation carries the (1 + 6αH²) factors.
- c_T² − 1 = 2.4e-11 … 2.3e-8 today, far beyond the GW170817 bound (~6e-15).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
