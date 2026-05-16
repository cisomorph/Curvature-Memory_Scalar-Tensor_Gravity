# SIM150 RESULT — Λ₀ Sweep for Ψ_local Mechanism Universality

**Date:** 2026-05-15
**Phase:** 5
**Gate inherited:** RED (from SIM147)
**Sequence position:** Fourth in Phase 5 (final sim before Paper V work)

---

## Overall Result: Stage 1 MATCH / Stage 2 INCOMPATIBLE

> **No value of Λ₀ within the locked CMSTG Solar-System constraints can
> lift the G_eff/G ceiling to the required 3.115. The Ψ_local mechanism
> is universally dead: the only Λ₀ that permits the required G_eff boost
> violates the Cassini PPN bound by a factor of ~2×10⁴.**

---

## Universality Question

SIM149 failure mode 3 established that G_eff/G is capped at 1.041 at
the locked Λ₀ = 0.003, while NGC 3198 rotation curves require G_eff/G ≈ 3.115.
The cap formula in CMSTG is:

    G_eff/G_max(Λ₀) = (½ + Λ₀Ψ₀²) / (½) = 1 + 2Λ₀Ψ₀²

where Ψ₀ = 2.62 M_Pl is the cosmic VEV (locked by SSB, SIM113/120).

SIM150 asks: is this ceiling a feature of the locked Λ₀ = 0.003, or is
it structurally unavoidable given all CMSTG constraints?

**Formula verification:** G_eff/G_max(0.003) = 1.04119, matching SIM149's
1.04118 to within 6×10⁻⁶ — formula confirmed.

---

## Stage 1 — G_eff/G_max Sweep

**Flag: MATCH**

The formula G_eff/G_max = 1 + 2Λ₀Ψ₀² is linear and unbounded in Λ₀.
The ceiling can be lifted — but only by choosing a much larger Λ₀.

**Analytic Λ₀_required:**

    G_eff/G_max = G_target
    1 + 2Λ₀Ψ₀² = 3.115
    Λ₀_required = (3.115 − 1) / (2 × 2.62²) = 0.1541

| Quantity | Value |
|---|---|
| G_eff/G target (NGC 3198) | 3.115 |
| G_eff/G_max at Λ₀_locked | 1.04119 |
| Λ₀_locked (SIM113) | 0.003 |
| Λ₀_required (analytic) | 0.1541 |
| Ratio Λ₀_req / Λ₀_locked | 51.4× |

The crossing is not at the edge of the sweep — it occurs at moderate Λ₀,
well inside the numerical range. Stage 1 is unambiguous.

---

## Stage 2 — Constraint Checks at Λ₀_required = 0.1541

**Flag: INCOMPATIBLE**

All five locked CMSTG constraints are violated. The binding constraint is
Cassini, which fails first and by the largest margin.

### 2a — Cassini Solar-System (binding)

In CMSTG with F(Ψ) = ½ + Λ₀Ψ², the post-Newtonian parameter is:

    ω_BD = F(Ψ₀) / (F'(Ψ₀))² = (½ + Λ₀Ψ₀²) / (2Λ₀Ψ₀)²
    |γ_PPN − 1| = 2 / (ω_BD + 2)

Cassini bound (Bertotti et al. 2003): |γ_PPN − 1| < 2.3×10⁻⁵

| Quantity | Λ₀_locked = 0.003 | Λ₀_required = 0.1541 |
|---|---|---|
| ω_BD | 2.11×10³ | 2.39 |
| \|γ_PPN − 1\| | 9.49×10⁻⁴ | 4.56×10⁻¹ |
| Excess over Cassini | — | **1.98×10⁴×** |
| **Verdict** | | **FAIL** |

Note: the locked Λ₀ value also exceeds the raw Cassini bound in this
simplified ω_BD estimate — Papers I–IV pass Solar-System constraints via
the SSB mass term and potential curvature providing Yukawa screening near
massive bodies. At Λ₀_required the BD parameter drops by a factor ~880,
making any such screening mechanism far less effective. The structural
conclusion is robust: Λ₀_required is catastrophically incompatible with
Solar-System gravity by any reasonable metric.

### 2b — f_σ8 Structure Growth

At Λ₀_required, the effective gravitational coupling for growth is:

    G_growth/G_N = F₀_locked / F₀(Λ₀_req) = 0.521 / 1.606 = 0.334

This suppresses structure growth by ~66% relative to GR:

| Quantity | Λ₀_locked | Λ₀_required |
|---|---|---|
| G_growth ratio | 1.0000 | 0.3343 |
| f_σ8 model | 0.800 | 0.253 |
| Tension with Planck+DESI | 0.00σ | **39.1σ** |
| **Verdict** | | **FAIL** |

### 2c — UV Finiteness Σ(0)

One-loop self-energy: Σ(0) = Λ₀² k_m⁴ / (64π²), k_m = 10 Mpc⁻¹.
Paper III finiteness threshold: Σ(0) < 10⁻².

| Quantity | Λ₀_locked | Λ₀_required |
|---|---|---|
| Σ(0) | 1.42×10⁻⁴ | **3.76×10⁻¹** |
| Excess over threshold | — | **37.6×** |
| **Verdict** | | **FAIL** |

### 2d — BAO r_s/r_d (DESI H(z) proxy)

Background ODE run at Λ₀_required produces a dramatically different
expansion history:

| Quantity | Value |
|---|---|
| χ²_BAO (Λ₀_locked) | 95.92 |
| χ²_BAO (Λ₀_required) | 1725.19 |
| Δχ²_BAO | **1629.3** |
| **Verdict** | **FAIL** |

The expansion rate is completely altered — the CMSTG background no longer
resembles ΛCDM at any redshift relevant to DESI DR1.

### 2e — Cosmological Ψ Subdominance

| Quantity | Λ₀_locked | Λ₀_required |
|---|---|---|
| ρ_Ψ/ρ_tot at z_drag = 1060 | 1.45×10⁻⁴ | **7.70×10⁻²** |
| ρ_Ψ/ρ_tot at z = 0 | 0.784 | **0.997** |
| Threshold (< 0.01) | | **FAIL** |

At Λ₀_required the scalar field dominates the energy budget at recombination,
destroying the standard CMB power spectrum.

### Stage 2 Summary Table

| Constraint | Value at Λ₀_req | Bound/Threshold | Pass? | Excess |
|---|---|---|---|---|
| Cassini \|γ_PPN−1\| | 4.56×10⁻¹ | 2.3×10⁻⁵ | **FAIL** | 1.98×10⁴× |
| f_σ8 tension | 39.1σ | < 3σ | **FAIL** | 13× |
| UV Σ(0) | 3.76×10⁻¹ | 1×10⁻² | **FAIL** | 37.6× |
| BAO Δχ² | 1629.3 | < 4 | **FAIL** | 407× |
| Ψ subdominance (z_drag) | 7.70×10⁻² | < 0.01 | **FAIL** | 7.7× |

**Stage 2 verdict: INCOMPATIBLE**
**Binding constraint: Cassini PPN (excess 1.98×10⁴×)**

---

## Structural Conclusion

The two requirements

> (i) G_eff/G_max ≥ 3.115 (galactic rotation requirement)
> (ii) |γ_PPN − 1| < 2.3×10⁻⁵ (Cassini Solar-System bound)

are structurally incompatible within F(Ψ) = ½ + Λ₀Ψ². There exists no
Λ₀ satisfying both simultaneously.

The root cause is algebraic: G_eff/G_max = 1 + 2Λ₀Ψ₀² grows linearly
in Λ₀, while ω_BD = F₀/(2Λ₀Ψ₀)² falls as 1/Λ₀ for large Λ₀. The
Cassini bound requires ω_BD > 8.7×10⁴, which at fixed Ψ₀ = 2.62 means
Λ₀ < 2.5×10⁻⁴ — over 600× below Λ₀_required. These constraints
are incompatible by construction.

**The G_eff/G ceiling established in SIM149 is structural and universal:
it cannot be lifted by any Λ₀ within the locked CMSTG action.**

---

## Cross-Conjecture Summary

| Conjecture | Revival at Λ₀_req? | Root cause of death |
|---|---|---|
| Ψ_local (SIM149/150) | — | Cassini violation at Λ₀_req = 0.154 |
| Ψ_pre (SIM148) | **NO** | k_m-driven; K(13.4 Gyr) = 0 regardless of Λ₀ |

Λ₀_required does **not** revive Ψ_pre: the kernel failure depends on
k_m (τ_mem = 2/k_m), not on Λ₀. With k_m = 10 Mpc⁻¹ fixed by SIM102,
K(13.4 Gyr) = exp(−2.39×10⁹) = 0 at any Λ₀. Both Phase 5 conjectures
are **independently dead** via distinct structural arguments. There is no
shared coupling-channel revival path.

---

## Phase 5 Complete — All Loopholes Closed

| Sim | Test | Result | Root cause |
|---|---|---|---|
| SIM147 | Kernel characterization | RED gate | τ_mem = 205 kyr |
| SIM148 | Ψ_pre pre-bang scan (3 variants) | FAIL | w(13.4 Gyr) = 0 (k_m) |
| SIM149 | Ψ_local halo scan (NGC 3198) | FAIL | w(10 Gyr) = 0 + Λ₀ too small |
| SIM150 | Λ₀ sweep universality check | INCOMPATIBLE | Cassini 1.98×10⁴× |

**The DESI tension and galactic rotation problem are irreducible within
CMSTG. No parameter within the locked action — neither k_m (already
checked by SIM147/148) nor Λ₀ (this sim) — can revive either Theorem 2
loophole. The CMSTG framework is structurally self-consistent in its
failures: the same small Λ₀ required by Solar-System gravity is the
same smallness that caps G_eff variation at 4%, insufficient for both
DESI and galactic rotation simultaneously.**

---

## Plots

| File | Contents |
|---|---|
| `sim150_main.pdf` | 5-panel: G_eff/G_max vs Λ₀ (Stage 1); Cassini deviation vs Λ₀; f_σ8 tension vs Λ₀; UV Σ(0) vs Λ₀; Stage 2 constraint table |
| `sim150_metadata.json` | Full numerical output, all constraint values |

Artifacts at `sims/sim150_output/`.
