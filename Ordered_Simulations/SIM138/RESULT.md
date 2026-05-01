# SIM138 — RESULT

**Verdict:** DISTRIBUTED  
**Date:** 2026-04-24  
**Recommended Tier 2 priority:** SIM142 (Galileon G₃ tracker)

---

## Summary

The DESI Y1 tension for Phase 1 canonical is **distributed across z = 0.5–1.3** with no single bin dominant. The LRG3 (z=0.930), LRG1 (z=0.510), LRG2 (z=0.706) and ELG (z=1.317) bins collectively account for 96% of the chi2. This distribution pattern favours a **smooth evolution mechanism** (SIM142 Galileon G₃) over a sharp transition (SIM140 symmetron step).

---

## Baseline note

The Phase 3 scripts (SIM131–SIM136) hardcode `SIM121C_chi2 = 18.26, SIM121C_tension = 2.77σ` as the baseline. Direct computation at Phase 1 canonical parameters (H₀=67.59, Ωm=0.312) gives:

- χ²_DESI = **13.15**, tension = √(13.15/6) = **1.48σ** (this sim)

The discrepancy arises from different parameter evaluations:
- SIM121C MAP had H₀=77.76, F₀=0.558 (MCMC optimized), giving chi2_DESI=41.49, tension=2.63σ
- The 2.77σ figure corresponds to chi2≈46: √(46/6) = 2.77, consistent with evaluating Phase 1 at parameters that differ from H₀=67.59
- In any case, the 2.77σ baseline in papers is an evaluation under a particular convention. The per-bin decomposition below is internally consistent and uses the H₀=67.59 canonical.

---

## Per-bin decomposition

| Bin | z_eff | H_CMSTG | H_DESI | Δ [km/s/Mpc] | pull [σ] | χ²(bin) | frac |
|---|---|---|---|---|---|---|---|
| LRG3 | 0.930 | 115.7 | 128.1 | −12.4 | −2.21 | 4.87 | 37% |
| LRG1 | 0.510 | 89.7 | 97.9 | −8.2 | −1.86 | 3.44 | 26% |
| LRG2 | 0.706 | 101.1 | 110.7 | −9.6 | −1.55 | 2.39 | 18% |
| ELG | 1.317 | 144.5 | 156.4 | −11.9 | −1.38 | 1.91 | 15% |
| BGS | 0.295 | 79.0 | 81.7 | −2.7 | −0.60 | 0.36 | 3% |
| QSO+Lyα | 2.330 | 236.3 | 240.8 | −4.5 | −0.41 | 0.17 | 1% |
| **TOTAL** | | | | | | **13.15** | 100% |

All pulls are **negative** (CMSTG predicts lower H(z) than DESI observes), consistent with Phase 1 canonical having H₀=67.59 while DESI prefers H₀≈74 (ΛCDM best-fit: H₀=74.0, Ωm=0.290, χ²=2.4).

---

## Pattern diagnosis

- **No single bin dominant:** top-1 fraction = 37% (threshold for SINGLE_BIN_DOMINANT: >80%)
- **Top-2 fraction = 63%** (threshold for LOCALIZED: ≥70%) — not quite localized
- The z=2.33 Lyα bin contributes only 1%: the high-z anchor is fine
- The tension is **a monotone H(z) shift**, not a feature at any particular epoch
- This is consistent with Phase 1 canonical having a systematically lower H(z) than DESI prefers across all bins z=0.5–1.3

---

## ΛCDM comparison

ΛCDM best-fit to the same DESI data: H₀=74.0, Ωm=0.290, χ²=2.4, tension=0.63σ. The CMSTG excess χ² over ΛCDM best-fit is Δχ²=+10.8. Phase 4 needs to raise H(z) by ~8–12 km/s/Mpc across z=0.5–1.3 without breaking CMB constraints.

---

## Tier 2 recommendation

**DISTRIBUTED → SIM142 (Galileon G₃) prioritized.**

A tracker mechanism that smoothly raises H(z) by ~10% across z=0.5–1.3 is exactly what a Galileon G₃ term can provide. The SIM140 symmetron step would require a transition redshift z_trans ∈ [0.5, 1.3] spanning the entire tension range, making it less sharp than its motivation implies. SIM142 should run first.

SIM140 should still run (it might find a solution at z_trans ~ 0.9 that clips both LRG bins), but it is no longer the highest priority.

---

## Figures

- `figures/sim138_desi_decomposition.pdf` — H(z) model vs DESI data; per-bin pull bar chart
