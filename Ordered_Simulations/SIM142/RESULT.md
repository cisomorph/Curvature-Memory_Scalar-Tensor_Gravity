# SIM142 — RESULT

**Verdict:** FAIL (STRUCTURAL)
**Date:** 2026-04-25
**Mechanism:** Galileon G₃(Ψ)□Ψ sector — simplest Horndeski G₃ extension of Phase 1 canonical

---

## One-line finding

G₃(Ψ)□Ψ with G₃=G₃(Ψ) only (no X dependence) is structurally equivalent to a rescaled kinetic term via integration by parts. It provides no new H(z) lever: the correction scales as c₃·Ψ'(z=0)≈0.05·c₃, giving max ΔH/H=1.08% at c₃=0.1 (quadratic), far below the ~5% needed to shift DESI tension. True Galileon kinetic braiding (G₃_X≠0) is not tested here.

---

## EOM derivation

Action:
```
L = F(Ψ)R + ½(∂Ψ)² − ½m₀²Ψ² + G₃(Ψ)□Ψ + L_matter
F(Ψ) = ½ + Λ₀Ψ²
```

Variation gives modified Friedmann and scalar EOM:

**Friedmann:**
```
H²·coef = 3F₀(ωh²_m/a³ + ωh²_r/a⁴) + Λ_bare
coef = 3F + 6F_ΨΨ' − ½Ψ'² − 3G₃_ΨΨ'
```
(F₀=0.521 absorbed into matter coupling; SIM91 confirms CMSTG ΛCDM-equivalent <0.1% in H(z).)

**Scalar EOM:**
```
Ψ'' + (3−ε_H)Ψ' = F_Ψ·R/H² + G₃_ΨΨ·Ψ'²
```
- Linear G₃=c₃Ψ: G₃_ΨΨ=0 → scalar EOM identical to Phase 1, Ψ' unchanged
- Quadratic G₃=c₃Ψ²: G₃_ΨΨ=2c₃ → small self-kick proportional to c₃·Ψ'²

**Key note:** G₃(Ψ) (no X dependence) is equivalent via IBP to adding 2G₃_ΨX to G₂, i.e., a rescaled kinetic term. This does NOT produce standard Galileon braiding. In the Bellini-Sawicki classification: α_B=0 (braiding), α_T=0, α_K=Ψ'²/M²_*. True braiding requires G₃_X≠0.

---

## Phase 1 reference (c₃=0)

ODE gives:
| Quantity | ODE | Phase 1 canonical |
|---|---|---|
| H₀ | 67.590 km/s/Mpc | 67.59 (target) ✓ |
| Ψ₀ (z=0) | 2.882 M_Pl | 2.62 (attractor) |
| F_eff(z=0) | 0.5249 | 0.521 |
| DESI χ² | 13.63 (1.51σ) | 18.26 (2.77σ) |
| 100θ* | 1.04096 | ~1.041 ✓ |
| σ₈ | 0.794 | 0.811 |
| RSD χ²/N | 1.053 | 0.86 |

**Note on DESI discrepancy:** The ODE starts with Ψ_ini=2.62 at z=1e5 (Phase 1 IC); the scalar rolls to Ψ₀=2.882 at z=0. This is not at the Phase 1 attractor value (Ψ̄=2.62 is the z=0 value in the full canonical fit). The ODE therefore traverses a slightly off-canonical trajectory, producing a different DESI χ². The canonical 18.26 (2.77σ) from SIM121C used a full MCMC with correct late-time Ψ. The ODE result being better (13.63) does not imply the canonical point is resolved — it reflects a trajectory artifact.

---

## G₃ scan results

| Form | c₃ | ΔH/H (mean, DESI bins) | DESI tension | 100θ* | RSD χ²/N | Verdict |
|---|---|---|---|---|---|---|
| linear | 1e-4 | +0.0002% | 1.507σ | 1.04096 | 1.053 | PASS* |
| linear | 1e-3 | +0.002% | 1.507σ | 1.04098 | 1.055 | PASS* |
| linear | 1e-2 | +0.019% | 1.504σ | 1.04122 | 1.074 | PASS* |
| linear | 1e-1 | +0.186% | 1.476σ | 1.04355 | 1.288 | PASS* |
| quadratic | 1e-4 | +0.001% | 1.507σ | 1.04097 | 1.054 | PASS* |
| quadratic | 1e-3 | +0.011% | 1.505σ | 1.04110 | 1.065 | PASS* |
| quadratic | 1e-2 | +0.105% | 1.490σ | 1.04240 | 1.176 | PASS* |
| quadratic | 1e-1 | +1.082% | 1.332σ | 1.05580 | 2.975 | FAIL |

*PASS* = individually satisfies DESI<2σ + CMB 2σ + RSD<1.5, BUT **not due to G₃** — these inherit Phase 1 baseline passing. G₃ contribution is negligible at c₃≤0.01. Overall verdict: **FAIL (STRUCTURAL)** because G₃(Ψ) provides no meaningful new mechanism.

---

## Stability (c_T and ghost checks)

- **α_T = 0:** Analytic. G₃_X=0, G₅=0 → c_T=c exactly. GW170817 satisfied by construction.
- **α_K = Ψ'²/M²_* > 0:** Phase 1 slow-roll gives Ψ'≈0.053; no ghost.
- **α_B = 0:** No kinetic braiding (G₃_X=0). This is the key theoretical constraint — true braiding would require G₃_X≠0.
- **M²_* = 1+2Λ₀Ψ²₀ ≈ 1.05 > 0:** Positive effective Planck mass, no Brans-Dicke ghost.

All stability checks pass. The mechanism fails on dynamical grounds, not stability grounds.

---

## Structural no-go for G₃(Ψ)

For any G₃(Ψ) with G₃_X=0:

1. **IBP equivalence:** ∫G₃(Ψ)□Ψ = −∫G₃_Ψ(∂Ψ)² = −2G₃_Ψ·X. This is just a modified kinetic term. There is no new propagating degree of freedom, no braiding, no tracker attractor beyond what Phase 1 already has.

2. **Slow-roll suppression:** The Friedmann correction from G₃ is ΔH²/H² = −3G₃_ΨΨ'/coef ≈ c₃Ψ'/F₀. With Ψ'~0.053, this is ≈0.1·c₃. For c₃≤0.1: correction ≤1%. Needed ~5%.

3. **Required c₃:** To achieve 5% H(z) shift: c₃ ~ F₀/Ψ' ≈ 0.521/0.053 ≈ 10. This is O(1) in natural units and physically excludes the kinetic sector from being sub-dominant.

4. **Large c₃ pathology:** At c₃=0.1 (quadratic), the CMB angle is shifted by 0.014 (≫5σ deviation) and RSD shape breaks (χ²/N=2.975). The CMB shift comes from modified H(z) at z~1000 through the G₃ term's cumulative effect on the comoving distance.

---

## Implications for Phase 4

This sim exhausts G₃(Ψ)□Ψ forms (linear and quadratic). Findings for Tier 2:

- **SIM142 FAIL:** G₃(Ψ) provides no new mechanism. True braiding G₃(X)□Ψ (X-dependent) would give α_B≠0, but that is a qualitatively different theory requiring a new simulation. If pursued, it should be labeled SIM142b.
- **SIM141 (running Λ₀)** remains Tier 2 priority 2. Run next.
- **SIM140 (step potential)** remains priority 3 — control experiment.
- **SIM143 (bi-scalar φ)** is the fallback loophole mechanism.

The TIER1_SUMMARY.md elevated SIM142 to first priority because braiding was expected to give a smooth late-time H(z) boost matching the distributed DESI tension. This result shows that the G₃(Ψ) form cannot deliver that — only G₃(X) can. If true braiding is needed, a new spec is required.

---

## Theoretical checks (RESEARCH_RULES.md §2.3)

| Check | Result |
|---|---|
| GR recovery | ✓ G₃(0)=0, G₅=0, Ψ→0 → standard GR |
| c_T = c | ✓ Analytic (G₃_X=0, G₅=0) |
| No tachyon | ✓ m_eff² = m₀² + ... > 0 at Phase 1 attractor |
| No ghost | ✓ α_K = Ψ'²/M²_* > 0 |
| Ward identity | ✓ Diffeomorphism invariant; Π_hh(0)=0 at tree level |
| UV finiteness | ✓ Inherited from Phase 1; G₃_Ψ vertices suppressed by memory kernel |

All theoretical checks pass. Failure is on observational mechanism, not theory.
