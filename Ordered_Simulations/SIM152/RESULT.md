# SIM152 — Paper V Option B: SPARC Rotation Curve Test

**Phase:** Paper V (Option B branch)  
**Date:** 2026-05-15  
**Mode:** Two-stage pass/fail (Stage 1 ≤ 2.0 χ²/dof)  
**Status:** FAIL — structural shape mismatch, both stages

---

## Stage 1 — NGC 3198 (43 SPARC data points)

### Results

| Test | χ²/dof | Verdict |
|------|--------|---------|
| Uniform β₀ = 1.018, E = 3.115 everywhere | 91.9 | FAIL (threshold ≤ 2.0) |
| Radial β(ρ_b), best ρ_sc* sweep | 91.9 | FAIL |
| Free β₀ at ρ_sc* | 91.4 (β₀ = 0.999) | FAIL |
| Diagnostic ρ_sc = 3.5×10⁻²⁴ g/cm³ (mid-range) | 572.0 | FAIL (chameleon variation makes it worse) |

**ρ_sc* pinned at lower boundary** (10⁻²⁶ g/cm³, below all galactic densities) — the optimal is the uniform limit (ρ_sc → 0). Any positive ρ_sc introduces chameleon density-variation that reduces enhancement in the outer disk where V_obs is flat and V_bary is falling. This confirms that the chameleon shape variation is structurally unhelpful.

### Shape Diagnostic — Key Finding

**Pearson correlation between E_model(r) and E_required(r): r = −0.811**

| Radius | E_required = V_obs²/V_bary² | E_model (chameleon ρ_sc = 3.5×10⁻²⁴) |
|--------|------------------------------|---------------------------------------|
| 2 kpc | ~1.8 | ~3.1 |
| 10 kpc | ~2.6 | ~1.07 |
| 20 kpc | ~3.8 | ~1.04 |
| 40 kpc | ~5.2 | ~1.04 |

E_required is **monotonically increasing** with radius (outer disk needs more gravity).  
E_model is **monotonically decreasing** with radius (chameleon gives less gravity where less baryons).  
The two are anti-correlated (r = −0.81). No parameter choice can reconcile this.

### Structural Diagnosis

The β(ρ_b) = β₀ tanh²(ρ_b/ρ_sc) profile gives:
- **Inner disk** (high ρ_b → ρ_b >> ρ_sc): β → β₀, E → 3.115 (LARGE enhancement)
- **Outer disk** (low ρ_b → ρ_b << ρ_sc): β → 0, E → 1.041 (no enhancement)

Flat rotation curves require:
- **Inner disk**: modest enhancement (V_bary ≈ 0.7 V_obs → E_req ≈ 2.0)
- **Outer disk**: LARGE enhancement (V_bary ≈ 0.45 V_obs → E_req ≈ 5.2)

The coupling provides large enhancement exactly where it is not needed, and no enhancement where it is most needed. This is a **structural inversion** — not a parameter mismatch.

The minimum χ²/dof achievable by any choice of (β₀, ρ_sc, Υ_disk) is bounded from below by the anti-correlated shape requirement. Even with free Υ_disk marginalization, the shape mismatch persists (V_bary shape is set by the SPARC-measured Vdisk and Vgas, independent of Υ_disk normalization).

---

## Stage 2 — NOT REACHED

Per spec: Stage 1 FAIL stops Stage 2. "Option B has the wrong shape and Stage 2 will not rescue it."

The shape mismatch is galaxy-morphology-independent (flat rotation curves universally require outer-disk enhancement) so no ensemble of galaxies will reverse this verdict.

---

## Structural No-Go Theorem

**Theorem (SIM152):** A monotonically increasing coupling β(ρ_b) cannot produce flat galactic rotation curves within CMSTG Option B.

**Proof sketch:** For a flat rotation curve V_obs ≈ V_flat = const,  the required enhancement factor is:
```
E_required(r) = V_flat² / V_bary²(r)
```
Since V_bary(r) falls for r > R_disk (past the disk peak), E_required is monotonically increasing. The model produces E(r) = 1.041 + 2β²(ρ_b(r)). Since ρ_b(r) is also monotonically decreasing (exponential disk), and tanh² is increasing, β(r) is monotonically decreasing, and hence E(r) is monotonically decreasing. E_required ↗ and E_model ↘ are incompatible. No choice of ρ_sc (which controls where in r the transition from β₀ to 0 occurs) can reverse this because the ordering of the radii is fixed. ∎

**Corollary:** A coupling β(ρ_b) that DECREASES with density (anti-chameleon) would give E_model ↗, potentially matching E_required ↗. Example: β ∝ sech²(ρ_b/ρ_sc) or β ∝ exp(−ρ_b/ρ_sc). However, such a profile would also make β LARGE at cosmological densities (ρ_b ~ ρ_crit ≪ ρ_sc), threatening the SIM151 β_∞ ≲ 2.4×10⁻¹⁰ bound. A non-monotone profile would be required (large at intermediate galactic densities, small at both high inner-disk and low cosmological densities) — a qualitatively different architecture not derivable from the SIM151 ansatz.

---

## Conclusion

**Option B (monotone chameleon coupling β ∝ tanh²(ρ)) is not viable for galactic rotation curves.** The architecture is structurally inverted. χ²/dof = 91.9 at best (vs threshold 2.0), anti-correlation r = −0.81 between model and required enhancement. Stage 2 not reached.

The finding qualifies as a structural no-go for Option B as specified. It does not rule out all β(ρ) architectures — it rules out the specific monotonically-increasing profile. Options:
1. **Anti-chameleon profile** (β decreasing with ρ): potentially viable for rotation curves but requires re-examination of cosmological constraints (β_∞ at ρ_cosmo).
2. **Option A** (m(Ψ) mass modification): already FAIL (Phase 2, SIM112-type).
3. **Option C** (sextic condensate/Bose-star, Paper III §Option C): not yet tested.

---

## Outputs

- **PDF:** `sims/sim152_output/sim152_main.pdf` — 3-panel figure
- **JSON:** `sims/sim152_output/sim152_metadata.json`
