# No-Go Theorems

Two structural no-go theorems (Paper I) together prove that the 2.77σ DESI Y1 BAO residual cannot be reduced within the CMSTG action class by any late-time modification, no matter how it is engineered.

---

## Theorem 1 — Curvature-Sourced Scalars Grow the Wrong Way

> **Any scalar field φ sourced by curvature R > 0 and initialized at zero grows monotonically, increasing F_eff and thereby suppressing H(z) at DESI redshifts.**

### Plain language

The DESI tension requires H(z) to be *higher* than ΛCDM predicts at z ≈ 0.3–2.3. To raise H(z) you need F_eff to *decrease* (less effective gravitational coupling → faster expansion).

But in CMSTG — and in any scalar-tensor theory of this class — a scalar coupled to curvature via R > 0 will grow as it is sourced. Growing Ψ means growing Ψ², which means F_eff = (1 + 2Λ₀Ψ²)/2 *increases*. This is the opposite of what is needed.

Every curvature-sourced mechanism tested in Phase 3 (SIM131–136) worsened the tension — some dramatically (SIM133: Gauss-Bonnet, 9.68σ).

### Technical statement

For any mechanism of the form `□φ ⊃ f(Ψ, R)` with f > 0 when R > 0 and φ(t₀) = 0:

```
dφ/dt > 0  for all t > t₀
  ⟹  d/dt[(1 + 2Λ₀(Ψ+φ)²)/2] > 0
  ⟹  H²(z) decreases  (wrong direction for DESI)
```

**Established by:** SIM131–136 (Phase 3) · Proved in Paper I §3

---

## Theorem 2 — DESI Fix and CMB Preservation Are Mutually Exclusive

> **Any mechanism that raises H(z) on z ∈ [0, z_rec] with fixed sound horizon r_s necessarily increases θ* = r_s / D_C* above the Planck bound 100θ* = 1.04101 ± 0.00029. Pure late-time modifications cannot satisfy both constraints simultaneously.**

### Plain language

The CMB acoustic peak position fixes the angular scale θ* = r_s / D_C* with extraordinary precision. The sound horizon r_s is set at recombination — it is essentially the same in CMSTG as in ΛCDM (given the tiny Λ₀ = 0.003).

To reduce the DESI BAO tension you need H(z) to be higher at late times. But raising H(z) while keeping r_s fixed *shrinks* D_C* (the comoving distance to last scattering). A smaller D_C* means a larger θ*, which violates the Planck bound. This is a geometric constraint independent of any specific mechanism.

The most extreme instance: SIM141 (running Λ₀(a) / BD analog) violated θ* by **63σ**.

### Technical statement

With r_s fixed at its CMSTG value:

```
δH(z) > 0  on  z ∈ [0, z_rec]
  ⟹  δD_C* = −∫ δH/H² dz/(1+z) < 0
  ⟹  δθ* = −r_s δD_C* / D_C*² > 0
  ⟹  100θ* exceeds Planck bound
```

**Established by:** SIM137–144 (Phase 4), especially SIM141 (63σ) · Proved in Paper I §4

---

## Together: The Irreducibility Proof

| Route | Theorem | Outcome |
|-------|---------|---------|
| Raise H(z) via curvature-sourced scalar | Theorem 1 | Scalar grows wrong way — H(z) decreases |
| Raise H(z) via any other late-time mechanism | Theorem 2 | CMB θ* violated |
| Both constraints simultaneously | 1 + 2 | Impossible within CMSTG action class |

**Conclusion:** The 2.77σ DESI Y1 BAO residual is a structural floor of the CMSTG Phase 1 action. It is irreducible by any late-time modification.

---

## The Surviving Loophole

Neither theorem rules out a *pre-recombination* solution. If new physics operating before recombination increases the sound horizon r_s, then both constraints can in principle be satisfied simultaneously — a larger r_s at fixed θ* allows larger D_C*, which permits higher H(z) at z < z_rec.

**Candidate mechanisms:**

| Mechanism | Description | Key constraint |
|-----------|-------------|----------------|
| Early dark energy (EDE) | Sub-dominant energy component peaking near z_rec | Planck CMB power spectrum |
| ΔNeff > 0 | Extra relativistic species at BBN/recombination | BBN, Planck Neff |
| Modified recombination | Changes to hydrogen recombination history | CMB peak structure |

Each constitutes a distinct theoretical programme, separately constrained by Planck CMB power spectra.

**DESI Y3 trigger:** If tension ≥ 3σ in ≥ 2 tracer samples with no systematic explanation → Phase 5 activates and Paper VII is written.

---

## Simulation Evidence

| SIM | Test | σ | Note |
|-----|------|---|------|
| SIM111 | m₀ scan — minimum tension | 3.44σ | Structural floor confirmed |
| SIM131 | Additive ξΨR | 5.97σ | Theorem 1 instance |
| SIM133 | Gauss-Bonnet sourcing | 9.68σ | Theorem 1 — worst case |
| SIM136 | Horndeski G_μν∂Ψ² kinetic | 3.53–3.75σ | Theorem 1 instance |
| SIM141 | Running Λ₀(a) / BD analog | 63σ (θ*) | Theorem 2 — most extreme |
| SIM144 | Mixed-source completeness probe | FAIL | Closes Tier 2 mechanism space |

See [Simulation Programme](Simulation-Programme) for full tables.
