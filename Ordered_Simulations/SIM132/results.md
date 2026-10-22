# SIM132 — Results

**Phase:** Phase 3  
**Title:** ξΨR scan over ξ  
**Verdict:** FAIL

## What was tested

Scan over ξ ∈ [0.01, 1]; FAIL — all ξ violate CMB θ*.

## Action tested

```
S=∫d⁴x√g[F_eff(Ψ)R - ½(∂Ψ)²]+S_SM, F_eff=½(1-2Λ₀Ψ²)+ξΨ
```

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** faithful PARTIAL / θ*-matched PASS (published: FAIL) — Solar-System safe

- Λ0 = 0.2: Δχ² = +9.27 (DR1) / +21.24 (DR2), robust to ±0.3. Ψ0 = 0.4163 sits at the attractor Ψ* = 0.4167, where F' ≈ 1.5e-4.
- |γ−1| = 7.9e-8 (passes Cassini); Ġ/G ≈ −4e-18 /yr (passes LLR).
- Open: G is about 7% higher at BBN; the full CMB spectrum (hi_class) is untested; Λ0 = 0.2 is at the grid edge; the DE is a bare Λ.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
