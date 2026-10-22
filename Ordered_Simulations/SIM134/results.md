# SIM134 — Results

**Phase:** Phase 3  
**Title:** Gauss-Bonnet parameter scan  
**Verdict:** FAIL

## What was tested

Scan over GB coupling strength α; FAIL — CMB violation systematic.

## Action tested

```
S=∫d⁴x√g[F_eff(Ψ)R - ½(∂Ψ)²]+S_SM, F_eff=½exp(-2Λ₀Ψ²)+ξΨ
```

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** faithful FAIL / θ*-matched PASS (published: FAIL) — Solar-System safe

- Λ0 = 0.2: Δχ² = +8.88 (DR1) / +19.20 (DR2), robust. Ψ0 = 0.4527 at the attractor; F' ≈ −1.7e-4.
- Passes Cassini and LLR. Same open issues as SIM132 (BBN, full CMB, bare Λ).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
