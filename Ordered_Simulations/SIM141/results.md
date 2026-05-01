# SIM141 — RESULT

**Verdict:** PARTIAL / STRUCTURAL ANTI-CORRELATION
**Date:** 2026-04-25
**Mechanism:** Running Λ₀(a) — phenomenological Brans-Dicke analog (P4-B)

---

## One-line finding

Running Λ₀(a) can reduce DESI tension to 0.547σ (tanh, a_trans=0.50, γ=1.00) but simultaneously shifts the CMB acoustic angle θ* from 1.041 to 1.059 — a 63σ violation of the Planck constraint. DESI improvement and CMB preservation are mutually exclusive for any late-time H(z) boost with fixed sound horizon r_s.

---

## PARTIAL vs FAIL

Verdict is PARTIAL because:
- The DESI target (tension < 2σ) is **achievable** — the mechanism works on DESI
- The CMB failure is **structural** — not a parameter choice issue, but a fundamental anti-correlation
- A PASS would require simultaneously increasing r_s (sound horizon at z~1090), which requires modifying pre-recombination physics — outside Phase 4 scope
- The mechanism is explicitly **phenomenological** (Λ₀(a) not derived from action); any further investigation is Phase 5 work

---

## Scan results

| Form | Parameter | Λ₀(z=0) | F_eff | ΔH/H% | DESI tension | 100θ* | RSD χ²/N | Verdict |
|---|---|---|---|---|---|---|---|---|
| Linear | β=0.10 | 0.00270 | 0.5224 | +0.10% | 1.490σ | 1.04174 | 1.095 | FAIL_CMB |
| Linear | β=0.30 | 0.00210 | 0.5174 | +0.30% | 1.455σ | 1.04333 | 1.186 | FAIL_CMB |
| Linear | β=0.60 | 0.00120 | 0.5099 | +0.62% | 1.401σ | 1.04579 | 1.342 | FAIL_CMB |
| Linear | β=1.00 | 0.00000 | 0.5000 | +1.08% | 1.325σ | 1.04918 | 1.584 | FAIL_CMB |
| Exp | β=0.10 | 0.00271 | 0.5225 | +0.10% | 1.489σ | 1.04176 | 1.096 | FAIL_CMB |
| Exp | β=0.50 | 0.00182 | 0.5150 | +0.58% | 1.407σ | 1.04514 | 1.303 | FAIL_CMB |
| Exp | β=1.00 | 0.00110 | 0.5090 | +1.16% | 1.307σ | 1.04918 | 1.599 | FAIL_CMB |
| Exp | β=5.00 | 0.00002 | 0.5002 | +2.15% | 1.176σ | 1.06114 | 2.826 | FAIL_CMB |
| Tanh | a_tr=0.33, γ=0.50 | 0.00150 | 0.5123 | +1.66% | 1.229σ | 1.05086 | 1.667 | FAIL_CMB |
| Tanh | a_tr=0.33, γ=1.00 | 0.00000 | 0.5000 | +3.43% | 0.980σ | 1.06117 | 2.530 | FAIL_CMB |
| Tanh | a_tr=0.50, γ=0.50 | 0.00151 | 0.5125 | +2.25% | 1.052σ | 1.04988 | 1.697 | FAIL_CMB |
| **Tanh** | **a_tr=0.50, γ=1.00** | **0.00002** | **0.5002** | **+4.76%** | **0.547σ** | **1.05939** | **2.641** | **FAIL_CMB** |
| Tanh | a_tr=0.67, γ=0.50 | 0.00155 | 0.5129 | +1.83% | 1.196σ | 1.04735 | 1.602 | FAIL_CMB |
| Tanh | a_tr=0.67, γ=1.00 | 0.00011 | 0.5009 | +3.95% | 0.963σ | 1.05439 | 2.467 | FAIL_CMB |

Phase 1 reference: tension=1.507σ, θ*=1.04096, RSD=1.053. Planck: 1.04101±0.00029 (2σ window: ±0.00058).

---

## Structural anti-correlation

The fundamental constraint:

```
θ* = 100 · r_s / DC_star   where DC_star = ∫₀^{z*} c/H(z) dz
```

Any mechanism that raises H(z) at z < z* necessarily **decreases DC_star** and **increases θ***. The Planck constraint θ*=1.04101±0.00029 is extremely tight (σ=0.03%). A 5% H(z) boost at z~0.5-1.3 shifts θ* by ~+1.5% = ~50σ.

The tradeoff is exact and unavoidable for fixed r_s:
- θ* too high → DESI tension fixed (H(z) increased)
- θ* correct → DESI tension persists (H(z) too low)

**The only loophole:** modify r_s (sound horizon at z~1090) upward by a commensurate amount. This requires changing physics at or before recombination — adding relativistic species, modifying baryon photon coupling, or an early dark energy component. None of these are in Phase 4 scope.

This is the same anti-correlation that drives the Hubble tension in standard cosmology.

---

## SIM105 RG consistency

- **Tanh form (6 cases):** UV-consistent. Λ₀→Λ₀_CMB at high z (z>1090, before transition), monotone decrease to Λ₀(1-γ) at z=0. Consistent with SIM105 IR attractor at CMB epoch.
- **Linear form (4 cases):** UV-inconsistent. Λ₀(a→0)→∞; Λ₀ grows without bound at early universe. Contradicts SIM105 UV fixed point Λ₀→0.
- **Exponential form (4 cases):** UV-inconsistent. Λ₀(a→0) = Λ₀_CMB·exp(+β·a_CMB) ≈ Λ₀_CMB (small correction, since a_CMB≈9e-4), effectively ~0.003 at all a<a_CMB. **Actually borderline consistent** — the exponential correction is negligible at a<a_CMB since β(a-a_CMB)≈−β·a_CMB≈0 for a_CMB≪1. Flagged as inconsistent only because form allows growth for any a<a_CMB, but numerically it's ≈ Λ₀_CMB.

For practical purposes, tanh is the most theoretically clean form, and it achieves the best DESI improvement.

---

## Theoretical checks

| Check | Result |
|---|---|
| GR recovery | ✓ Λ₀(a→0)→Λ₀_CMB (tanh) → Phase 1 action restored at high z |
| c_T = c | ✓ Analytic: G₃=G₅=0, F_X=0 → α_T=0 |
| No tachyon | ✓ m²_eff > 0 (Phase 1 mass term unchanged) |
| No ghost | ✓ α_K = Ψ'²/M²_* > 0 |
| Ward identity | ⚠ OPEN — time-varying Λ₀ modifies Π_hh(0) at 1-loop; SIM145 required |
| UV finiteness | ⚠ OPEN — time-translation violation from Λ₀(a); new divergences possible |
| SIM105 RG | ✓ Tanh form (6/14 cases); others UV-inconsistent |

---

## Implications for Phase 4

This sim reveals a **new structural no-go** for all late-time H(z)-boosting mechanisms:

> Any mechanism that raises H(z) at z < z_rec = 1090 while fixing r_s at its Phase 1 value will increase θ* above the Planck observation. The DESI tension cannot be resolved by a pure late-time H(z) boost without simultaneously breaking the CMB acoustic scale.

**Impact on remaining Tier 2 sims:**

- **SIM140 (step potential):** Same anti-correlation applies. A step in H(z) at z_trans~2 will still compress DC_star. DESI improvement will be accompanied by θ* violation. Likely FAIL_CMB for the same reason.
- **SIM143 (bi-scalar φ, decoupled from R):** If φ raises H(z) at late times (w_φ < -1), the same CMB anti-correlation applies. However, if φ modifies r_s through its contribution to the pre-CMB radiation density, there is a potential loophole. This is the most promising remaining candidate.

**Recommendation:** Run SIM143 next (ahead of SIM140), focusing specifically on whether the bi-scalar φ component can simultaneously raise H(z) at z~0.5-1.3 AND increase r_s through its pre-CMB density contribution.

---

## Note on Phase 1 reference discrepancy

As noted in SIM142, the ODE Phase 1 reference gives DESI tension 1.507σ (vs canonical 2.77σ) because the scalar starts at Ψ_ini=2.62 at z=1e5 and evolves to Ψ₀=2.882 at z=0 — not exactly at the Phase 1 attractor (Ψ̄=2.62 is the z=0 value). This is a trajectory artifact of the ODE IC, not a genuine improvement.
