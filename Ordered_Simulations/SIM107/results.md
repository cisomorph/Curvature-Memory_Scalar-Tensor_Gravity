# SIM107 — Results

**Phase:** Phase 1  
**Title:** Perturbation Theory: δΨ  
**Verdict:** N/A

## What was tested

Linear perturbations δΨ; no tachyon (ω²>0 for all modes).

## Key numerical results

- Λ₀ = 0.003

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** FAIL (published: PASS)

- At the canonical Ψ̄ = 2.62: |γ−1| = 9.48e-4, which is 41× the Cassini bound 2.3e-5. The published run only evaluated Ψ0 = 1.0.
- Correct Einstein-frame coupling: α² = F̃'²/(3F̃'² + 2F̃) with F̃ = 1 + 2Λ0Ψ². Cassini requires Ψ0 < 0.40, i.e. F0 − ½ < 4.8e-4.
- No screening: m_in·R_sun = 2.6e-4 and Khoury–Weltman ΔR/R = 1.3e7. The derived matter term in m_eff² is negative.
- Unit bugs fixed: G·M_sun (1.476e-12 → 4.78e-20 Mpc) and 1 AU (4.848e-9 → 4.848e-12 Mpc).
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), χ²/p-value/GOF σ, and Δχ² vs ΛCDM in the same pipeline.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/`.
