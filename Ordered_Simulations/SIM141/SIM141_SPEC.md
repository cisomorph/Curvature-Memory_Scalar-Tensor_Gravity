# SIM141 — P4-B: Running Λ₀(a)

**Tier:** 2 (Mechanism)
**Depends on:** Tier 1 complete; SIM140 optionally complete (either result informs this)
**Prerequisite reads:** RESEARCH_RULES.md, PHASE4_ROADMAP.md, Paper VII direction P4-B, SIM105 RG flow results

---

## Question

Can a time-varying coupling Λ₀(a), decreasing after recombination, produce F_eff(z < 1) < F_eff(z_CMB) while preserving the CMB acoustic scale and respecting the SIM105 UV structure (negative beta function, Λ₀ fixed point)?

## Motivation

Phase 1 established (SIM105) that Λ₀ runs under RG flow with Λ₀ = 0 as a UV fixed point and Λ₀ ≈ 0.003 as the IR attractor. The running was computed at the level of quantum corrections (high-energy → low-energy). Cosmologically, Λ₀ has been treated as a constant.

The P4-B direction asks: if Λ₀ carries cosmological-scale a-dependence on top of the RG running, does this unlock a window in which F_eff weakens at late times? This is a Brans-Dicke-like mechanism phrased in CMSTG-native language.

**Critical caveat:** a-dependent Λ₀ in the cosmological background is NOT automatically compatible with the Lagrangian — it has to descend from some underlying dynamics (e.g., a secondary field that couples to Λ). This sim is exploratory; any PASS result becomes the starting point for a Phase 5 first-principles derivation, not a claim that Λ₀(a) is derived.

## Action spec

Effective Lagrangian (**phenomenological — label as such**):

```
L_eff = (1/(16πG) + Λ₀(a)·Ψ²) · R  −  ½(∂Ψ)²  −  ½m₀²Ψ²  +  L_matter
```

with Λ₀(a) scanned over three functional forms:

- **Linear:** Λ₀(a) = Λ₀_CMB · (1 − β(a − a_CMB))
- **Exponential:** Λ₀(a) = Λ₀_CMB · exp(−β(a − a_CMB))
- **tanh transition:** Λ₀(a) = Λ₀_CMB · [1 − γ · ½(1 + tanh((a − a_trans)/σ))]

All with Λ₀_CMB = 0.003 (Phase 1 locked).

## Theoretical checks (required before PASS)

1. **GR recovery:** Λ₀(a → 0) → 0.003 (Phase 1 value at CMB epoch) — high-z behaviour must match Phase 1.
2. **UV consistency with SIM105:** the phenomenological a-dependence must not contradict the RG flow. At high energy (early universe), Λ₀ should tend toward 0 (UV fixed point). Check that Λ₀(a → 0) is consistent with RG extrapolation.
3. **No tachyon:** m₀² > 0 preserved.
4. **c_T = c:** non-minimal coupling (1/16πG + Λ₀Ψ²) R gives c_T = c exactly when G⁵ = 0. Document.
5. **Ward identity:** if Λ₀ is time-varying, Π_hh(0) may pick up corrections. The Ward identity derivation in SIM104 assumed constant Λ₀. **Flag as open:** full UV recheck (SIM145) required if this sim passes.
6. **Causality:** Λ₀(a) depends on cosmic time but not on spatial position → no issue.

## Inputs

- Phase 1 Friedmann solver, modified to accept Λ₀ as a function of a
- CLASS Boltzmann solver (parameters passed per-z if feasible, or effective Λ₀ at the relevant scale)
- Planck plikHM, DESI Y1, 9-point fσ₈ datasets
- SIM105 RG flow result: analytic formula 1/Λ₀² = 1/Λ₀_obs² + (k_m² − k_m,nat²)/(16π² m₀²)

## Procedure

1. Modify Phase 1 Friedmann solver to accept Λ₀(a) input. Verify that with Λ₀ = constant it reproduces SIM90 exactly.
2. For each of the three functional forms, scan one or two parameters:
   - Linear: β ∈ {0.01, 0.1, 1.0}
   - Exponential: β ∈ {0.01, 0.1, 1.0}
   - tanh: (a_trans, σ, γ) on a small grid; γ controls amplitude, a_trans controls when, σ controls width
3. For each point, compute F_eff(z), H(z), D_M(z), D_H(z).
4. Check CMB acoustic scale (100θ*) preservation.
5. Run full Planck likelihood + DESI likelihood for a few best candidates.
6. Consistency check with SIM105: at Λ₀ values reached during the evolution, verify β-function sign has not been inverted.

## Success criteria (all must hold)

- (a) CMB preserved: 100θ* ∈ [1.0408, 1.0414], Δ(−2 ln L_plik) < +5
- (b) DESI tension < 2σ (joint χ² < 12)
- (c) SIM105 RG structure not contradicted (Λ₀ flow remains monotone with negative beta function)
- (d) No tachyon; BBN G_eff bound satisfied
- (e) Phase 1 passes preserved (BBN, Solar System, GW speed)

If (a)-(e) all hold: PASS (conditional on SIM145 UV recheck).
If (a), (c), (d), (e) hold but (b) partial: PARTIAL.
If (c) fails: FAIL. The parametrization is phenomenologically viable but theoretically inconsistent with Phase 1 foundations — not acceptable.

## Failure modes to watch

- **Phenomenological vs derived:** Λ₀(a) as written is NOT derived from the Lagrangian. Any PASS result must be explicitly flagged as phenomenological and linked to a Phase 5 derivation task. Do NOT let this parametrization leak into a paper as if it were derived.
- **Coupling with Ward identity:** as noted in theoretical check 5. If this sim passes, SIM145 becomes essential, not optional.
- **CMB at z ~ 1090 vs initial conditions at z_init:** Λ₀(a_CMB) is not the same as Λ₀(a → 0) strictly. Integrate from z_init and check boundary behaviour carefully.
- **Degeneracy with H₀:** any mechanism that alters H(z) can fit data by shifting H₀. Do not mistake a H₀ shift for a genuine DESI tension reduction. Fix H₀ at SIM121C best-fit and report ΔH₀ at the end.

## Deliverables

- `output.json` with best functional form, best parameters, all checks
- `RESULT.md` with verdict, per-form comparison, SIM105 consistency statement, flag on phenomenological status
- Figures: Λ₀(a), F_eff(z), H(z) vs data

## Estimated time

Medium. Mostly solver modifications + parameter scans. Less physics complexity than SIM140.
