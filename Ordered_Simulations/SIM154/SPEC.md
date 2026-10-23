# SIM154 — SPEC (pre-registered 2026-10-22, before run.py was written)

## Question
Do the Solar-System-safe attractor models found in the audit (SIM132 subtractive F, SIM134 exponential F)
have a region of parameter space that also improves on ΛCDM for DESI and survives BBN?

## Models (Jordan frame, 8πG = 1, GR is F = ½; canonical Ψ; bare Λ)
- Family A (subtractive, SIM132): F = ½(1 − 2Λ0Ψ²) + ξΨ. Attractor Ψ* = ξ/(2Λ0), F* = ½ + ξ²/(4Λ0).
- Family B (exponential, SIM134): F = ½e^{−2Λ0Ψ²} + ξΨ. Attractor where F' = 0, if ξ < √Λ0·e^{−1/2}.
- Initial condition: Ψ = 0, Ψ' = 0 at z = 1e5 (radiation era, field frozen). So F = ½ at BBN.

## Grid (both families)
Λ0 ∈ {0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5} × ξ ∈ {1/24, 1/12, 1/8, 1/6, 1/4}: 35 points per family, 70 in total.

## Method
- Background: corrected equations (`common/st_solver.py`, validated by `tests/`). H0 is solved from θ*.
- Data: DESI DR1 and DR2 (official likelihoods); Planck 100θ* = 1.04110 ± 0.00031; Planck ω_b, ω_m.
- Local gravity at z = 0: Einstein-frame α² from F(Ψ0), F'(Ψ0). Ġ/G = −Ḟ/F.
- BBN: speed-up S = √(G_BBN/G_N) = √(F0/F_BBN) = √(2F0). Kneller & Steigman (2004) fits:
  - Y_p = 0.2485 + 0.0016[(η10 − 6) + 100(S − 1)]
  - 10⁵ D/H = 2.64 [6/(η10 − 6(S − 1))]^1.6
  - with η10 = 273.9 ω_b. At S = 1 these give Y_p ≈ 0.249 and 10⁵ D/H ≈ 2.55.

## Survival criteria (copied from ROADMAP.md; all must hold)
- C1 cosmology: DESI DR2 GOF < 2σ (dof = 13 − 2) **and** Δχ²(DR2) vs ΛCDM ≥ +4 (θ*-matched)
- C2 Solar System: |γ − 1| < 2.3e-5 **and** |Ġ/G − 7.1e-14| ≤ 3 × 7.6e-14 /yr
- C3 BBN: Y_p within 2σ of 0.245 ± 0.003 (⊕ 0.0006 fit error) **and** D/H within 2σ of (2.547 ± 0.025)e-5 (⊕ 3% fit error)
- C4 stability: F > 0 along the whole history (no gravity flip)

## Outputs
`output.json` (every point, every criterion), `RESULT.md`, `results.md`, `figures/`.
Survivors (if any) go to ROADMAP Phase 3 (full CMB with hi_class + Cobaya).
