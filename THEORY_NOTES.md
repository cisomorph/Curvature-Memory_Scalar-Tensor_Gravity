# CMSTG theory notes (ROADMAP Phase 4, started 2026-10-23)

## 4.1 Why the attractor models work, and whether the attractor form can be derived

The attractor families (SIM132, SIM134, SIM154):

    S = ∫√-g [ F(Ψ)R − ½(∂Ψ)² ] + S_m + Λ_bare
    family A:  F = ½(1 − 2Λ0Ψ²) + ξΨ
    family B:  F = ½e^{−2Λ0Ψ²} + ξΨ

In the matter era the scalar equation is Ψ'' + (3 − ε)Ψ' = F'(Ψ)·R/H², with R/H² = 6(2 − ε) ≈ 3 > 0. So Ψ moves
up the slope of F. It behaves as if it had an effective potential V_eff ∝ −F(Ψ)·R, which is minimised where F is
**largest**, i.e. where F'(Ψ*) = 0:

- Family A: Ψ* = ξ/(2Λ0), F* = ½ + ξ²/(4Λ0).
- In the radiation era R ≈ 0, so Ψ stays frozen. BBN therefore sees F(Ψ_ini) = ½.

Ingredients, mapped to the original CMSTG language:
- The linear term ξΨR is the curvature-memory source: it lets curvature switch Ψ on from Ψ = 0 (SIM131's idea).
- The quadratic term must have the **opposite sign to Phase 1** (−Λ0Ψ² instead of +Λ0Ψ²). Otherwise F has no maximum and Ψ runs away (SIM131 FAIL). SIM132's docstring calls this a "local Deser–Woodard analog"; that link has not been worked out.
- Open question: is there a symmetry or principle that fixes ξ²/Λ0, the total change of G? BBN requires ξ²/(4Λ0) ≲ 0.012 (SIM154). The best DESI fit sits right at that edge. A principled reason to be there would be the real theoretical result; right now it is a fit.

## 4.2 Relation to Damour–Nordtvedt (GR as a cosmological attractor)
Damour & Nordtvedt (PRL 70, 2217, 1993; PRD 48, 3436, 1993) showed that scalar-tensor theories whose Einstein-frame
coupling function A(φ) has a minimum are driven toward GR during cosmic expansion. The present-day coupling
α = d ln A/dφ ends up small but nonzero, typically |γ − 1| ~ 1e-7 … 1e-5.

The CMSTG attractor models are a member of this class:
- In the Einstein frame A² ∝ 1/(2F), and α ∝ F'(Ψ). The fixed point F' = 0 is their α = 0.
- The SIM154 survivors have |γ − 1| ≈ 1e-7, squarely in the DN range.

What CMSTG adds: the source term (ξΨR) that moves Ψ *away* from Ψ = 0, so G changes by about 2% between BBN and today.
That change, not the attractor itself, is what improves the DESI fit.
**Refined by SIM155 (full CMB, see 4.5):** at the full-CMB best fits the change happens between BBN and recombination,
not at late times, and it improves the CMB damping tail and the BAO scale rather than mimicking dynamical dark energy. Any paper must cite DN and position the work
relative to it. There is also later literature on nonminimally coupled scalars and DESI that should be surveyed
before claiming novelty.

## 4.3 Dark energy is still a bare Λ
- Every model that survives (SIM130–154) needs Λ_bare.
- SIM123 replaces it with a hilltop potential V = λ(Ψ² − v²)², but with v = 13.16 ≫ Ψ that is ≈ λv⁴, a cosmological constant in disguise.
- So the README's "dark energy emerges from the Ψ–geometry coupling without a bare cosmological constant" is not supported by any audited simulation. That remains open.

## 4.4 Paper III (UV finiteness, Λ0 fixed point)
Not audited.
- Two-loop finiteness of a scalar-tensor theory and a negative β-function for Λ0 are strong claims; they need an independent check.
- If the fixed point is Λ0 → 0 (GR), the attractor models' Λ0 ≈ 0.15 needs an explanation.
- The sign of the quadratic term is flipped relative to Phase 1, so Paper III's results may not carry over at all.

## 4.5 What the full-CMB fits show (SIM155, 2026-10-24): an early-universe effect, and a structural tension
Full hi_class + Planck 2018 + DESI DR2 + SN fits (Ordered_Simulations/SIM155/RESULT.md):
- **When G changes.** In the best-fit region the field reaches its attractor before recombination.
  - At the D1 best fit, G at BBN is 7–8% stronger than today, 0.2% stronger at recombination (z ≈ 1090), and equal to today from z ≈ 10 on.
  - The reason: in the matter era the effective mass of Ψ around Ψ* is of order √(Λ0 R/F) ~ H for Λ0 ~ 0.3, so the field settles within about an e-fold of matter–radiation equality.
- **What it changes.** CMSTG is therefore an *early* modification of gravity. After recombination its expansion is ΛCDM-like (with bare Λ) at shifted parameters.
  - It improves the fit through the CMB damping tail and the sound horizon r_d: plik-lite −3.4, low-ℓ EE −1.1, DESI −1.9 in χ².
  - It does not help the supernovae (+0.6). w0waCDM gains −2.6 there.
  - So the model cannot reproduce the late-time w(z) signal that drives DESI's preference for dynamical dark energy.
- **BBN.** The early change of G is exactly what BBN measures (S = √(G_BBN/G0) = √(2F0)). The improvement grows with ΔF = ξ²/(4Λ0) up to ΔF ≈ 0.04 (Δχ² ≈ +4.7), but BBN allows only ΔF ≲ 0.015 (Δχ² ≈ +3).
- **The structural tension (any Damour–Nordtvedt-class model of this type).**
  - Local gravity needs the field at the attractor today: |γ − 1| ∝ α0² ∝ F'(Ψ0)².
  - A late-time dark-energy-like signal needs the field still evolving at z ≲ 2, i.e. F' ≠ 0 recently.
  - Fast convergence (large Λ0) satisfies Cassini but makes the effect early (BBN-limited, no late-time signal). Slow convergence (small Λ0) gives late-time evolution but fails Cassini (SIM154: Λ0 ≤ 0.05 fails).
  - Escaping this needs a mechanism that separates local tests from cosmological evolution, i.e. screening (chameleon, symmetron or Vainshtein-type terms such as G3 or G4(X)). Alternatively the late-time signal must come from a separate dark-energy sector, which brings back the bare-Λ issue of 4.3.
  - This is the most useful theoretical question for CMSTG now: can "curvature memory" be given a screening form?
