# SIM139 — RESULT

**Verdict:** SHAPE_FIXABLE (phenomenological)  
**Date:** 2026-04-24  
**Note:** The fix requires a two-parameter coupling that is phenomenological — it must be derived from a Lagrangian before counting as a theoretical result.

---

## Summary

The SIM128 RSD shape failure is **structurally interpretable**: the model predicts too little fσ₈ at low-z (z < 0.5, deficit) and too much at high-z (z ≥ 0.5, excess). This is a **systematic slope** in the residuals, not a single outlier. The shape is fixable with a two-parameter phenomenological modification, but the one-parameter family (varying p or adding a cutoff) is insufficient. Phase 4 work is still needed to derive the shape fix from a Lagrangian.

---

## Sanity check

| Quantity | SIM128 | Reproduced |
|---|---|---|
| S₈ | 0.7532 | 0.7532 ✓ |
| χ²_RSD | 31.42 | 31.421 ✓ |
| χ²/N | 2.244 | 2.244 ✓ |

---

## Per-z residuals (ranked)

| z | Survey | fσ₈ model | fσ₈ obs | pull | χ²(z) | frac |
|---|---|---|---|---|---|---|
| 0.570 | BOSS | 0.505 | 0.427 | +3.38 | 11.43 | 36% |
| 0.020 | 2MTF | 0.309 | 0.428 | −2.54 | 6.44 | 20% |
| 0.170 | SDSS | 0.370 | 0.510 | −2.33 | 5.41 | 17% |
| 0.067 | 6dFGRS | 0.329 | 0.423 | −1.72 | 2.95 | 9% |
| 1.400 | eBOSS QSO | 0.645 | 0.482 | +1.40 | 1.97 | 6% |
| 0.800 | eBOSS LRG | 0.559 | 0.470 | +1.12 | 1.25 | 4% |
| … | … | … | … | … | … | … |

Top contributor (z=0.57, BOSS): 36% of total chi2. Excluding it: chi2/N = (31.42−11.43)/13 = 1.54 — still >1.2, so not a LOCALIZED_OUTLIER.

---

## Shape pattern

- **Low-z (z<0.5) mean pull:** −0.71σ (model systematically too low)
- **High-z (z≥0.5) mean pull:** +1.40σ (model systematically too high)

The M̂(a) function rises from 0 at early times to 1 today, meaning Ĝ_eff/G = 1−νM̂ falls monotonically as a→1. This **suppresses growth too much at late times** (the memory field has accumulated). At early times (z>0.5), M̂ is still small, so growth is barely screened — but the data wants stronger screening there. At low-z, M̂≈1 gives maximum screening — but the data wants less. This is the slope inversion.

---

## Modification scan

| Modification | ν_opt | S₈ | dS₈(σ) | χ²/N | Verdict |
|---|---|---|---|---|---|
| p=0.5 (phenom) | 0.375 | 0.722 | −1.54 | 1.66 | PARTIAL |
| p=1.5 (phenom) | 0.825 | 0.757 | −0.09 | 2.79 | PARTIAL |
| p=2.0 (phenom) | 1.050 | 0.763 | +0.17 | 3.38 | PARTIAL |
| Cutoff z_c=0.3 (phenom) | 0.600 | 0.747 | −0.51 | 2.14 | PARTIAL |
| Cutoff z_c=0.5 (phenom) | 0.625 | 0.745 | −0.57 | 2.00 | PARTIAL |
| Cutoff z_c=0.7 (phenom) | 0.675 | 0.743 | −0.68 | 1.83 | PARTIAL |
| **Two-param (phenom)** | **1.85** | **0.735** | **−1.02** | **0.82** | **PASS** |

**All "phenomenological" labels are mandatory per RESEARCH_RULES.md §2.2.**

The two-parameter form Ĝ_eff/G = 1 − νM̂ + ηM̂² (ν=1.85, η=2.56) achieves χ²/N=0.82 and S₈=0.735 (within KiDS-1000 range). The η>0 term adds a positive correction that partially restores growth at late times (where M̂≈1, M̂²≈1 too, so the net effect is 1−1.85+2.56=1.71... wait, that's >1 which means enhanced gravity). 

Recheck: at a=1 (z=0), M̂=1: Ĝ_eff/G = 1−1.85+2.56 = 1.71 > 1. At early times M̂→0: Ĝ_eff/G→1. So the coupling **enhances** gravity at low-z and is neutral at high-z. This inverts the original M̂ structure. A Lagrangian realization of this form would need to be derived.

---

## Verdict justification

**SHAPE_FIXABLE** because:
1. A mathematical solution exists (two-parameter phenomenological coupling)
2. The solution is physically interpretable (need enhanced, not suppressed, late-time gravity to simultaneously match low-z fσ₈ data and high-z BOSS measurement)

**But:** SHAPE_FIXABLE here means the data is compatible with a modified coupling — it does not mean Phase 4 is avoided. The fix requires a new Lagrangian term that naturally produces Ĝ_eff/G > 1 at late times. This is a positive Geff enhancement, opposite to the original memory-screened gravity picture.

---

## Recommendation

1. **Note in Paper VIII:** the SIM128 shape failure points toward a **late-time gravity enhancement**, not additional screening. This is potentially derivable from a Galileon G₃ term or a running Λ₀(a) — both of which can modify the effective gravitational coupling in the required direction.

2. **SIM142 (Galileon G₃)** is the correct Tier 2 candidate to realize the shape fix from a Lagrangian.

3. SIM128 remains PARTIAL, not upgraded to PASS, until a derived Lagrangian term reproduces the phenomenological fix.

---

## Figures

- `figures/sim139_rsd_diagnostic.pdf` — fσ₈ model vs data; per-z pull bar chart
