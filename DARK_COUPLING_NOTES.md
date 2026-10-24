# Dark curvature memory: CMSTG coupled to the dark sector only (ROADMAP Phase 9, started 2026-10-24)

Option 1 of SCREENING_NOTES §4, chosen by the user together with option 2.

**Idea.** Keep CMSTG's coupling function, but let only dark matter feel it. Baryons, photons and the Solar System see
plain GR with a constant G, so Cassini and LLR are satisfied by construction. This also solves option 2's
problem: ordinary matter never feels a drift in G.

## 1. The model
Einstein frame, 8πG = 1:

    S = ∫√-g [ R/2 − ½(∂φ)² − V(φ) ] + S_b,γ,ν[g] + S_c[A²(φ) g]

- Dark matter lives in the CMSTG "Jordan frame": A(φ) = (2F(φ))^{-1/2}, so its mass is m_c(φ) = A(φ) m_c0.
- For family A around its attractor, F = ½ − Λ0 δφ² (F* = ½), so A = (1 − 2Λ0δφ²)^{-1/2} ≈ 1 + Λ0δφ².
- The coupling is β(φ) ≡ d ln A/dφ = 2Λ0δφ/(1 − 2Λ0δφ²).

**"Dark pinned thawing" (DPT)**, the dark-sector analogue of SIM156's PT:
- V = V0 + v1 δφ, with V0 the bare Λ fixed by closure.
- Initial conditions δφ = δφ' = 0 deep in the radiation era (z = 3000).
- **Early on** the dark-matter density pins the field at the minimum of A, where β = 0: V_eff = V + ρ_c A, and the restoring force is 2Λ0ρ_c δφ. So nothing happens before recombination, and the CMB and BBN are standard.
- **Late**, ρ_c dilutes and the tilt pushes the field away. A grows, so **dark matter gets heavier**. That is the direction that mimics DESI's phantom crossing (Khoury, Lin & Trodden, PRL 135, 181001, 2025; arXiv:2604.08449; arXiv:2601.05235).

## 2. Equations (N = ln a; H in units of 100 km/s/Mpc; ρ in units with 3H² = ρ)
Background:

    ρ_c = 3 ω_c a⁻³ A(φ)                (ω_c = early-time value; today ω_c A(φ0))
    H² (1 − φ'²/6) = ω_b a⁻³ + ω_ν a⁻³ + ω_c A a⁻³ + ω_r a⁻⁴ + Λ_b + V/3
    φ'' + (3 − ε) φ' = −(V_φ + β ρ_c) / H²
    ε = −H'/H = [3(ω_b + ω_ν + ω_c A) a⁻³ + 4 ω_r a⁻⁴ + H² φ'²] / (2H²)

Linear growth (sub-horizon, quasi-static, light field). Derivation:
- Dark-matter momentum obeys v̇ + (H + βφ̇) v = −∇(Φ + βδφ)/a.
- The heavier particles slow down; the βφ̇ term is the *mass-change drag*.
- The field responds as ∇²δφ = a² β ρ_c δ_c, so the scalar adds a dark-matter-only force 2β² × gravity.
- Baryons feel only gravity.

    δ_c'' + (2 − ε + β φ') δ_c' = (3/2) [ (1 + 2β²) Ω_c δ_c + Ω_b δ_b ]
    δ_b'' + (2 − ε) δ_b'         = (3/2) [ Ω_c δ_c + Ω_b δ_b ]          (Ω_i = ρ_i / 3H²)

- **Two competing effects.** When dark matter gains mass (βφ' > 0) the drag *suppresses* growth, while the fifth force (1 + 2β²) *enhances* it. Which wins sets the sign of the S8 shift.
- **Validity of the light-field limit.** It assumes the Compton wavelength of δφ exceeds the scales of interest. Near the attractor m² ≈ 2Λ0ρ_c, i.e. m/H ~ √(6Λ0Ω_c). That is ≲ 2 for Λ0 ≲ 1, so the fifth force acts on all sub-horizon scales. For larger Λ0 this limit overstates the fifth force.

**Apparent dark energy** (diagnostic). An observer who assumes uncoupled dark matter with today's density infers

    ρ_DE^eff(a) = 3H²(a) − 3(ω_b + ω_ν + ω_c A0) a⁻³ − 3 ω_r a⁻⁴,   w_eff = −1 − (1/3) d ln ρ_DE^eff / d ln a

Dark matter that was lighter in the past makes ρ_DE^eff smaller in the past, which reads as w_eff < −1.

## 3. What is CMSTG about it
The coupling function A = (2F)^{-1/2} is CMSTG's conformal factor; only the set of matter that feels it changes.
- The attractor (maximum of F, minimum of A) is what pins the field early. That is "curvature memory" of the dark-matter curvature: the field is held by the dark part of R and released when it fades.
- Compared with generic coupled quintessence (constant β), this coupling **vanishes in the early universe automatically**. The literature shows that is required to evade CMB constraints (arXiv:2604.08449: the field must start "frozen deep in the radiation era").

## 4. Pre-registered test: SIM157
See Ordered_Simulations/SIM157/SPEC.md and ROADMAP Phase 9.

## References (verified 2026-10-24)
- J. Khoury, M.-X. Lin, M. Trodden, *Apparent w < −1 and a lower S8 from dark axion and dark baryons interactions*, PRL 135, 181001 (2025), arXiv:2503.16415.
- arXiv:2604.08449, *Coupled dark energy and dark matter for DESI: an effective guide to the phantom divide* (2026).
- arXiv:2604.12032, *Constraints on coupled dark energy in the DESI era* (2026): |β| ~ 0.03 peak with Planck CamSpec + DESI DR2 + SNe; no coupling excluded at ~95%; not strongly preferred over ΛCDM.
- arXiv:2601.05235, *Mimicking phantom dark energy with evolving dark matter mass* (2026).
- L. Amendola, *Coupled quintessence*, PRD 62, 043511 (2000); *Linear and non-linear perturbations in dark energy models*, PRD 69, 103524 (2004).
- S. Alam et al. (eBOSS), PRD 103, 083533 (2021), Table III. RSD-only fσ8: 0.53 ± 0.16 (z = 0.15), 0.500 ± 0.047 (0.38), 0.455 ± 0.039 (0.51), 0.448 ± 0.043 (0.70), 0.315 ± 0.095 (0.85), 0.462 ± 0.045 (1.48). Read from the PDF.
