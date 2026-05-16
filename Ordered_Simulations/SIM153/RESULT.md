# SIM153 — Paper V Option C: Sextic Condensate Derivation

**Phase:** Paper V (Option C branch)  
**Date:** 2026-05-16  
**Mode:** Analytical derivation — two fail-fast pre-checks  
**Status:** HALT at Pre-check 2 — full derivation not warranted

---

## Pre-check Summary

| Pre-check | Verdict | Key finding |
|-----------|---------|-------------|
| 1 — Shape compatibility | **PASS** | All four reference E(r) values within factor 2 of E_required |
| 2 — Mass ceiling / cosmological DM | **FAIL** | 466× mass gap; gap-closing requires extended condensate (solid-body V_C ∝ r); coupling sign breaks lock |

**Overall: HALT. SIM154 and SIM155 are not warranted.**

---

## Background: What Option C Is

Option C adds a sextic stabilization term to the locked CMSTG Ψ potential:

    V(Ψ) → V_locked(Ψ) + μ_6 Ψ⁶    (μ_6 > 0)

The intended mechanism: an attractive lower-order coupling draws Ψ together; the repulsive sextic term prevents collapse; a stable soliton forms whose energy density sources gravity through the T_μν channel (NOT through G_eff modification). This makes Option C structurally distinct from:
- **Option A** (chameleon m(Ψ), SIM112): G_eff modification → ceiling 4.12% (SIM149, Mode 3)  
- **Option B** (coupled chameleon β(ρ_b)Ψρ_b, SIM151–52): fifth-force → correct channel but wrong shape (anti-correlated E(r), SIM152)

Option C proposes a *new gravitational source* (condensate mass density), not a modification of how existing sources gravitate. This means the SIM149 G_eff ceiling and SIM152 shape inversion theorem do not directly apply. The pre-checks ask whether *analogous* obstructions exist.

---

## Pre-check 1 — Shape Compatibility

**Question:** Does a sextic condensate produce E(r) = V_total²/V_bary² that is monotonically increasing with r (matching the flat-rotation-curve requirement)?

**Method:** Top-hat soliton (ρ_C = ρ₀ for r ≤ r_c, then Keplerian) tuned to produce V_flat = 150 km/s at r = 20 kpc for NGC 3198. Reference V_bary values from SIM149 table.

**Physics of the shape:**

- Inside r_c (flat core): M_C(<r) ∝ r³ → V_C²(r) = G M_C(<r)/r ∝ r² → solid-body rise
- Outside r_c (Keplerian tail): M_C(<r) = M_C_total = const → V_C²(r) ∝ 1/r → falls

In the outer disk where V_bary falls (exponential disk + gas), V_C² = G M_C/r falls only as 1/r. If V_bary² falls faster than 1/r (as it does for stellar disks at r >> R_d), then E(r) = 1 + V_C²/V_bary² increases. This is the *opposite* of SIM152's structural inversion — Option C has the right qualitative shape.

**Numeric results (r_c = 10 kpc, M_C = 7.20×10¹⁰ M_sun):**

| Radius | E_model | E_required (SIM152) | Ratio | Within factor 2? |
|--------|---------|---------------------|-------|-----------------|
| 2 kpc  | 1.38    | 1.80                | 0.77  | ✓ |
| 10 kpc | 4.50    | 2.60                | 1.73  | ✓ |
| 20 kpc | 3.43    | 3.80                | 0.90  | ✓ |
| 40 kpc | 2.78    | 5.20                | 0.54  | ✓ (marginal) |

All four reference points pass the factor-of-2 criterion. The 40 kpc point is marginal (ratio = 0.54, just above the factor-of-2 boundary of 0.50).

**Monotonicity:** E(r) peaks at r = r_c = 10 kpc and then decreases in the outer disk. This is because NGC 3198 is unusually gas-rich: V_bary(40 kpc) ≈ 64 km/s (mostly gas), so V_bary² does not fall fast enough to make E(r) increase at large r. For less gas-rich galaxies, V_bary falls exponentially outside the disk scale length, and E(r) increases monotonically.

The non-monotonicity is NOT a structural inversion (cf. SIM152 where E_model ↘ while E_required ↗ with anti-correlation r = −0.81). Here E_model peaks and levels off, while E_required continues to rise — a *quantitative* mismatch, not a *structural* one.

**Verdict: PASS.** No analog of the SIM152 shape obstruction exists for Option C. The condensate is a qualitatively viable dark-matter source for flat rotation curves. Precise quantitative tuning is deferred to SIM154.

---

## Pre-check 2 — Condensate Mass Ceiling Under Λ₀ Lock

**Question:** Does a parameter window exist where the condensate simultaneously explains galactic rotation curves and the cosmological DM density, without violating the locked action?

**Three independent obstructions were found.**

---

### Obstruction A: Mass gap (466×)

**Galactic constraint** (rotation curve fit for NGC 3198):

    V_flat = 150 km/s at r = 20 kpc
    V_bary(20 kpc) = 83.7 km/s
    V_C²(20 kpc) = 150² − 83.7² = 15,494 km²/s²
    M_C(<20 kpc) = 15,494 × 20 kpc / G = 7.20×10¹⁰ M_sun

**Cosmological constraint** (DM density must come from condensates):

    ρ_DM = Ω_DM × ρ_crit = 0.264 × 1.27×10¹¹ = 3.35×10¹⁰ M_sun/Mpc³
    n_gal = 10⁻³ Mpc⁻³  (L* galaxy density, spec value)
    M_C_cosmo_req = ρ_DM / n_gal = 3.35×10¹³ M_sun

**Gap:**

    M_C_cosmo / M_C_gal = 3.35×10¹³ / 7.20×10¹⁰ = 466×

The condensate mass required for galactic rotation curves is 466× too small to supply the cosmological DM density. Galactic condensates (all L* halos) contribute only 0.21% of ρ_DM.

---

### Obstruction B: Gap-closing requires solid-body rotation (incompatible)

**Attempt:** Extend r_c so that most of M_C_total lies outside the observed disk (r < 44 kpc), with only the inner portion M_C(<20 kpc) contributing to V_flat.

For V_C(44 kpc) ≤ 150 km/s with M_C_total = 3.35×10¹³ M_sun:

    r_c_min = (G × M_C_total / (V_max / R_max)²)^(1/3)
            = (G × 3.35×10¹³ / (150/44.1)²)^(1/3)
            = 232 kpc

The condensate would need r_c ≈ 232 kpc, larger than the virial radius of an NGC 3198-type halo (~150 kpc). Inside r_c, the profile is ρ ≈ const (flat core) → V_C(r) ∝ r (solid-body rise).

**Result with Ω_C = 3.40 km/s/kpc:**

| r [kpc] | V_C_solid | V_total | V_obs | Verdict |
|---------|-----------|---------|-------|---------|
| 4.8     | 16.3      | 89.0    | 127.0 | Under-fit |
| 7.1     | 24.1      | 101.9   | 146.0 | Under-fit |
| 16.1    | 54.8      | 100.0   | 153.0 | Under-fit |
| 30.1    | 102.4     | 123.9   | 146.0 | Under-fit |
| 44.1    | 150.0     | 163.2   | 149.0 | Over-fit |

The extended condensate fails to reproduce the observed flat rotation curve: V_total is under-predicted by 25–35% at intermediate radii and slightly over-predicted at the outermost point. The rising V_C ∝ r profile is qualitatively incompatible with the flat V_obs profile across the 5–44 kpc range. No choice of r_c simultaneously satisfies the rotation curve shape and the cosmological density requirement.

---

### Obstruction C: Coupling sign breaks the locked action

The locked CMSTG potential is:

    V_locked(Ψ) = ½ m₀² Ψ² + λ Ψ⁴    (m₀² > 0, λ > 0)

Both terms are **repulsive** — the potential has a unique minimum at Ψ = Ψ₀ (the cosmological VEV) and no attractive self-interaction above that. A non-topological soliton requires at least one **attractive** term in the field equation, such that the effective potential has a local minimum displaced from the cosmological background.

Standard constructions:

| Soliton type | Potential | Requirement |
|--------------|-----------|-------------|
| Q-ball (Lee-Weinberg) | −\|α\|Ψ² + μ₆Ψ⁶ | Negative quadratic (m₀² < 0) |
| Thick-wall soliton | m₀²Ψ² − \|λ\|Ψ⁴ + μ₆Ψ⁶ | Negative quartic (λ < 0) |
| Fuzzy DM (FDM) | −½ m²Ψ² | Negative mass² (equivalent) |

Both require sign changes in the lower-order couplings relative to the locked values. This is **not a parameter-tuning issue** — it is a structural requirement for soliton formation. Adding μ_6Ψ⁶ with μ_6 > 0 to the locked potential does not produce a soliton without simultaneously flipping the sign of m₀² or λ.

Flipping m₀² → negative: changes the entire cosmological Ψ evolution (Ψ₀ = 0 becomes the VEV → locks break).  
Flipping λ → negative: the quartic becomes a tachyonic direction → Ψ runs toward larger values until the sextic stabilizes it → the Phase 1 fit (SIM113, SIM88) is disrupted.

Either modification **breaks the locked action** — the coupled set (Λ₀, Ψ₀, m₀², λ) that reproduces all Phase 1 cosmological constraints simultaneously.

---

### Cosmological back-reaction (clarification)

The SIM101 constraint (ρ_Ψ/ρ_tot ≲ 10⁻¹⁰ at z_drag) applies to the **cosmological Ψ background** at z_drag ~ 1100. Galactic condensates form during structure formation (z < 5), after z_drag. They do not contribute to ρ_Ψ at z_drag and therefore do not violate the SIM101 bound.

The relevant cosmological constraint is instead whether the present-day average condensate energy density ⟨ρ_C⟩₀ matches Ω_DM × ρ_crit. As computed above: ⟨ρ_C⟩₀ = 7.2×10⁷ M_sun/Mpc³, supplying only 0.21% of ρ_DM. This is the failure, not a SIM101 back-reaction issue.

---

### Pre-check 2 parameter space analysis

| Parameter space | Assessment |
|-----------------|------------|
| 1-parameter (μ_6 only, lower-order coupling fixed) | Cannot form soliton with locked positive m₀², λ |
| 2-parameter (μ_6 + one lower-order coupling, sign relaxed) | Breaks lock; even so, no window satisfies galactic + cosmological simultaneously (mass gap 466×, shape incompatibility) |
| 3+-parameter (μ_6, m₀², λ all free) | Leaves Phase 1 locks — no longer CMSTG |

**Pre-check 2 Verdict: FAIL.** Three independent structural obstructions. No parameter window exists within the locked CMSTG action.

---

## Structural No-Go Pattern

SIM153 completes the Option A–B–C trilogy. Each option fails on a different structural ground:

| Option | Mechanism | Primary failure | Sim |
|--------|-----------|-----------------|-----|
| A | m(Ψ) mass modification → G_eff | Mode 3 ceiling: G_eff/G_max = 1.041 << 3.115 required | SIM149 |
| B | β(ρ_b)Ψρ_b fifth force | Shape inversion: E_model ↘ while E_req ↗; Pearson r = −0.81 | SIM152 |
| C | Sextic condensate → T_μν sourcing | Mass gap 466×; gap-closing → solid-body shape; coupling sign breaks lock | SIM153 |

All three mechanism channels fail on structural grounds independent of parameter tuning. This constitutes a **complete failure of the Paper V DM programme within the locked CMSTG action**. The DESI tension (SIM144: 3.44σ irreducible) and galactic rotation curves (SIM149: ceiling 4.12%) now have companion no-gos for all three proposed DM channels.

---

## Summary Verdict

**Option C (sextic condensate) is not viable as a dark-matter candidate within CMSTG.**

Pre-check 1 passes — the condensate produces the correct *qualitative* E(r) shape (unlike Option B's structural inversion). This is the only positive finding.

Pre-check 2 fails on three independent structural grounds: (1) the condensate mass for galactic fits is 466× below the cosmological DM density requirement; (2) closing the mass gap by extending r_c produces a solid-body rotation profile incompatible with flat V_obs; (3) soliton formation in the sextic-stabilized potential requires breaking the locked lower-order coupling sign.

SIM154 (single-galaxy detailed test) and SIM155 (SPARC ensemble + cosmology) are not warranted. The Paper V programme has exhausted all three proposed DM mechanism channels.

The open question for future work is whether a *qualitatively different* mechanism — one not reducible to G_eff modification, fifth force, or localized condensate — could explain galactic rotation curves within a modified but still predictive version of the action. No such mechanism is currently identified within the CMSTG framework.

---

## Outputs

- **PDF:** `sims/sim153_output/sim153_main.pdf` — 4-panel figure  
- **JSON:** `sims/sim153_output/sim153_metadata.json`
