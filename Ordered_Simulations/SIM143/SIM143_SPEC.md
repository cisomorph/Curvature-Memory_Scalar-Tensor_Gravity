# SIM143 — Bi-Scalar: Phase 1 Ψ + Potential-Driven φ (Decoupled from R)

**Tier:** 2 (Mechanism)
**Depends on:** Tier 1 complete; SIM140–142 results reviewed
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VI SIM135 (distinction below)

---

## Question

Can an independent scalar field φ, governed by a self-potential V(φ) and NOT coupled to the Ricci scalar R, produce dynamical dark energy with w_0 ≈ −0.76, w_a ≈ −0.79 (DESI-preferred) while the Phase 1 canonical Ψ sector remains intact?

## Motivation and distinction from SIM135

SIM135 tried a bi-scalar approach where φ was **sourced by curvature**. The no-go theorem killed it: curvature sourcing + φ_init = 0 + R > 0 → φ grows → F_eff grows → DESI worsens.

SIM143 tries the explicit loophole: φ is NOT sourced by R. It has its own potential and its own dynamics, uncoupled from the gravitational sector except through the metric. This is standard quintessence layered on top of Phase 1 canonical CMSTG.

The Phase 2 programme tried and failed to make Ψ itself play the role of quintessence (SIM112–113 SSB hilltop). SIM143 accepts that Ψ is structurally a memory/curvature-coupling field and introduces a **separate** DE scalar that is a dedicated quintessence.

## Action spec

```
L = (1 + 2Λ₀Ψ²)/2 · R                         (Phase 1 Ψ sector: unchanged)
    − ½(∂Ψ)² − ½m₀²Ψ²                         (Phase 1 Ψ kinetics + mass)
    − ½(∂φ)² − U(φ)                           (new DE scalar φ)
    + L_matter
```

with Λ₀ = 0.003 and Ψ frozen at Ψ̄ = 2.62 M_Pl throughout (Phase 1 attractor). The new field φ is minimally coupled to gravity; ∴ no direct R-sourcing term.

Potential forms to scan:

- **Exponential:** U(φ) = U₀ exp(−λφ/M_Pl), λ ∈ {0.5, 1.0, 2.0}
- **Double exponential:** U(φ) = U₀ [exp(−λ₁φ/M_Pl) + exp(−λ₂φ/M_Pl)], known tracker
- **Power-law:** U(φ) = U₀ (M_Pl/φ)^n, n ∈ {0.5, 1, 2}
- **Hilltop:** U(φ) = U₀ [1 − (φ/M_Pl)²]² (freezing; for comparison with DESI preference which is thawing — should fail)

## Theoretical checks (required before PASS)

1. **GR recovery:** U(0) = 0 and U(φ → 0) does not source gravity. Only issue is the late-time late-time non-zero U acting as effective Λ. OK by design.
2. **No tachyon:** U''(φ) ≥ 0 at all points on the trajectory.
3. **c_T = c:** φ is minimally coupled — no kinetic mixing with gravitons. c_T = c exactly. OK.
4. **Ward identity:** preserved (no new graviton-scalar vertices beyond minimal).
5. **UV:** U(φ) introduces new vertices. Bounded potentials have no quartic divergences. Power-law and exponential are bounded in the relevant field range. Document.
6. **Phase 1 sector independence:** the critical check. Ψ must remain at Ψ̄ = 2.62 M_Pl throughout the φ evolution. Since φ does not couple to Ψ and does not source Ψ, this is automatic — but verify numerically that Ψ doesn't drift due to background H(z) changes from φ.

## Inputs

- Phase 1 Friedmann solver + Phase 1 Ψ EOM (both as in SIM82–86)
- New φ EOM: φ̈ + 3Hφ̇ + U'(φ) = 0
- Modified Friedmann: 3H²·F_eff = ρ_m + ρ_r + ½φ̇² + U(φ) + (½m₀²Ψ̄² from Ψ — frozen, so constant)
- Planck plikHM, DESI Y1, RSD, GW datasets

## Procedure

1. Implement φ EOM alongside Phase 1 solver. Verify U(φ) = 0 reduces to Phase 1 SIM90 exactly.
2. For each potential form, scan (U₀, λ or n) to find combinations producing:
   - Sufficient DE density at z = 0: Ω_φ(0) ≈ 0.7 − Ω_Λ_bare(Phase 1 contribution)
   - w(z) evolution close to DESI-preferred w_0 ≈ −0.76, w_a ≈ −0.79
3. For promising candidates, compute CMB, BAO, RSD χ².
4. Verify Ψ stays at Ψ̄ throughout (sanity check for sector independence).
5. Cross-check: the DE-DM link from SIM120 (κ ≈ 3.84 × 10⁻¹⁰³) depends on Ψ̄. Does adding φ change this? Should NOT since Ψ is frozen — verify.

## Success criteria (all must hold)

- (a) CMB preserved: standard quintessence CMB impact well-understood; should be minor for trackers. Check 100θ*, Δ(−2 ln L_plik) < +5.
- (b) DESI tension < 2σ (joint χ² < 12)
- (c) RSD fσ₈ χ²/N < 1.5
- (d) Phase 1 Ψ sector unaffected: Ψ(z) − Ψ̄ < 0.01 M_Pl at all z, SIM120 DE-DM link unchanged
- (e) No tachyon, c_T = c, Ward preserved
- (f) Theoretical honesty: explicitly identify this as introducing a NEW degree of freedom not derivable from the Phase 1 action. SIM143 is a "what-if" rather than a "from-first-principles" result.

If (a)-(e) pass: PASS_PHENOMENOLOGICAL. Not a PASS in the same sense as a first-principles result — this is a working-mechanism-demonstrated result pending derivation.

## Failure modes to watch

- **This is standard quintessence.** The literature on exponential and power-law quintessence is vast. Draw from it (Copeland, Sami, Tsujikawa 2006 review; Tsujikawa 2013). Do not reinvent. The CMSTG-specific part is that φ coexists with the Phase 1 Ψ sector, not that the φ physics is novel.
- **Thawing vs freezing:** DESI prefers thawing (w_a < 0). Hilltop potentials give freezing. If only freezing potentials pass, the sim has not resolved the DESI tension — it has just shifted it. Ensure the scan includes genuinely thawing cases.
- **Two Ω contributions:** ρ_φ + ρ_Λ_bare = total DE density ≈ 0.7 × 3H₀²M_Pl². The Phase 1 Λ_bare is already set. Either reduce it by hand when adding φ (reparametrization) or let φ be sub-dominant DE. Document choice clearly.
- **"Decoupled" honesty:** φ not coupling to Ψ or R is a tuning. From an effective-field-theory perspective, any scalar at this energy scale should have couplings to everything unless forbidden by a symmetry. A truly first-principles Phase 4 mechanism needs such a symmetry argument. Flag this as a Phase 5 question if SIM143 passes.

## Deliverables

- `output.json` with best potential form, best parameters, all χ², Ψ-independence check
- `RESULT.md` with verdict including PHENOMENOLOGICAL caveat
- Figures: w(z), Ω_φ(z), Ω_m(z), Ψ(z) (should be flat)

## Estimated time

Small–medium. Standard quintessence — well-worn infrastructure. Mostly new code for coupling to Phase 1 solver.

## Deferral rule

Same as SIM142: if SIM140 or SIM141 produces a clean PASS, defer SIM143 to a future cycle.
