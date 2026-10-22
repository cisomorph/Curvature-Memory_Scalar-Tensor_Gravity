# SIM124 — Results

**Phase:** Phase 3  
**Title:** Kinetic Ψ Modification  
**Verdict:** FAIL

## What was tested

Non-canonical kinetic term f(X)∂Ψ²; FAIL — no viable solutions.

## Structural diagnosis

K(Ψ)=exp(-αΨ²) does NOT resolve the frozen-field problem. 'Viable' solutions require Ψ_ini ∈ [1.84, 2.57] M_Pl — the field starts within 0.78 M_Pl of the target Ψ₀=2.62 and barely rolls. This is fine-tuning, not dynamics. The hilltop bottleneck (Ψ≈0, K≈1, V'≈0, H large) is unchanged by the kinetic m

## Context

Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc.


## Audit rerun (2026-10-22)

**Audited verdict:** FAIL (published: FAIL)

- Published bugs: Friedmann without a GR limit; the F'R curvature source was missing from the field equation; r_s integrated over 0 → z*; non-DESI H(z) table.
- No roll from the hilltop. Reaching Ψ0 = 2.62 needs Ψ_ini = 2.403 (α = 0) or 2.246 (α = 0.1), i.e. fine-tuning. For α ≥ 0.2 the field runs away near the Ψ0 = 2.62 locus (unresolved; 20 s watchdog per case).
- The solutions that reach Ψ0 = 2.62 (θ*-matched): Δχ² vs ΛCDM = +4.24 / +10.99 (α = 0) and +6.09 / +15.60 (α = 0.1) for DR1 / DR2, but they fail Cassini by 41× and 82×.
- Data/statistics: official DESI DR1/DR2 likelihood, r_s to z = ∞, θ* = r_s(z*)/D_M(z*), proper statistics.
- Details: `AUDIT_NOTES.md` (repo root); log: `audit_logs/SIM124.log`.
