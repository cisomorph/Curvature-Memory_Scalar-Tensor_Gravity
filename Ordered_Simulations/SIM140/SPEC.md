# SIM140 — P4-C: Symmetron-Type Step Potential

**Tier:** 2 (Mechanism)
**Depends on:** Tier 1 complete; TIER1_SUMMARY.md reviewed
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VII §7 direction P4-C, Sim 8 history

---

## Question

Can a step potential V(Ψ) that holds Ψ frozen at the Phase 1 canonical value until a transition redshift z_trans, then releases it to roll to smaller Ψ, produce F_eff(z < z_trans) < F_eff(z_CMB) without breaking the CMB acoustic scale?

## Motivation

The Phase 3 no-go theorem established that any scalar sourced monotonically by R > 0 and initialized at zero grows F_eff and fails. The symmetron-type loophole is a **potential-driven** mechanism: Ψ does not grow from zero — it sits at the Phase 1 attractor Ψ̄ = 2.62 M_Pl and later rolls **downward**. This inverts the direction of F_eff evolution without re-sourcing from curvature.

This is the most tractable of the three Paper VII directions because it preserves the Phase 1 backbone and introduces a single transition.

## Action spec

```
L = (1 + 2Λ₀Ψ²)/2 · R  −  ½(∂Ψ)²  −  V(Ψ)  +  L_matter
```

with Λ₀ = 0.003 (Phase 1 locked) and

```
V(Ψ) = V₀ · [1 − tanh((Ψ² − Ψ_c²) / σ²)]
```

or equivalent smooth-step form. Parameters to scan:
- V₀: height of the plateau
- Ψ_c: transition field value (related to z_trans through background evolution)
- σ: width of the transition

The transition redshift z_trans is an output, not an input — determined by when V'(Ψ) becomes large enough to release Ψ from the Phase 1 attractor.

## Theoretical checks (required before PASS)

1. **GR recovery:** V(Ψ → 0) = 0 ✗ (fails by construction — this is a late-time potential, not a GR-recovery potential). Must verify instead that the scalar sector at z ≫ z_trans reproduces the Phase 1 canonical F_eff = 0.521 to < 1% precision.
2. **No tachyon:** V''(Ψ̄) > 0 at the initial (frozen) attractor. If V''(Ψ̄) < 0, Ψ is tachyonic at the Phase 1 canonical value.
3. **c_T = c:** the action is minimally coupled in the tensor sector (no G⁴ or G⁵ Horndeski term added), so c_T = c automatically. Document this.
4. **Ward identity:** not affected by V(Ψ) addition; inherited from Phase 1 (SIM104, SIM106).
5. **UV behaviour:** V(Ψ) introduces no new quartic divergences as long as it is polynomially bounded. tanh form is bounded. OK.
6. **Causality:** retarded Green's function structure preserved; V is local.

## Inputs

- Phase 1 FLRW solver (SIM82–86)
- Planck 2018 plikHM TTTEEE likelihood via `clipy` (SIM97/98 infrastructure)
- DESI Y1 BAO full covariance (SIM138 machinery, updated with Tier 1 per-bin weighting if SIM138 found localization)
- Phase 1 canonical parameters as initial conditions at z_init ≫ z_CMB

## Procedure

1. Implement V(Ψ) as above; wire into Phase 1 Friedmann + scalar ODE system.
2. Integrate the background equations from z_init to z = 0 for a grid of (V₀, Ψ_c, σ).
3. For each (V₀, Ψ_c, σ):
   - Verify Ψ ≈ Ψ̄ at z_CMB (tolerance < 1%)
   - Extract z_trans (z where dΨ/dt exceeds 10% of Hubble rate)
   - Compute F_eff(z) trajectory
   - Compute H(z), D_M(z), D_H(z)
4. For configurations passing the CMB-frozen check, run:
   - Full CLASS Boltzmann (TT, EE, TE) → χ²_CMB
   - DESI Y1 BAO χ² with full covariance
5. Find the (V₀, Ψ_c, σ) minimizing joint χ²_CMB + χ²_DESI.
6. Bonus if time: verify against RSD fσ₈ (9-point dataset) to ensure structure growth not broken.

## Success criteria (all must hold)

- (a) **CMB preserved:** 100θ* ∈ [1.0408, 1.0414] and Planck plikHM Δ(−2 ln L) < +5 relative to Phase 1 canonical
- (b) **DESI tension reduced:** joint DESI χ² below 18.26 (Phase 1 canonical benchmark); target < 12 (≈ 2σ)
- (c) **No tachyon, no ghost:** explicitly verified
- (d) **Phase 1 passes preserved:** BBN G_eff, Solar System γ, GW speed unchanged (since late-time transition doesn't affect these)

If all four hold: PASS → proceed to SIM144.
If (a) and (c) and (d) hold but (b) partial: PARTIAL → document and continue to SIM141.
If (a) or (c) fails: FAIL → document the obstruction.

## Failure modes to watch (the Sim 8 lesson)

- **The step potential must actually produce a transition.** Just writing down V(Ψ) with a step shape does not guarantee Ψ finds the descent path — the EOM solution may stick to the metastable minimum forever. Verify numerically that dΨ/dt becomes non-zero at the expected z_trans.
- **Transition too fast = oscillations.** If σ is too small, Ψ overshoots and oscillates around the new minimum. Oscillations produce dark-energy-like EOS spikes incompatible with DESI. Enforce σ large enough to ensure overdamped descent.
- **Transition too slow = no effect.** If σ is too large, the transition smears over most of cosmic history and F_eff barely changes. The parameter space between these extremes is the target.
- **Phase 1 attractor destabilization.** If V''(Ψ̄) < 0 at the supposed frozen value, Ψ is already rolling at z_CMB and the CMB check will fail before any late-time physics engages. The tanh form with Ψ_c > Ψ̄ protects against this — verify.

## Deliverables

- `output.json` with best-fit (V₀, Ψ_c, σ), F_eff(z) trajectory, all χ² contributions, theoretical check flags
- `RESULT.md` with verdict and diagnosis
- Figures: F_eff(z), H(z) vs ΛCDM and Phase 1 canonical, CMB residuals, BAO residuals
- If PASS: draft of Paper VIII mechanism section

## Estimated time

Medium–large. New potential, full Boltzmann + BAO chain, parameter scan.
