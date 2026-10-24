# Can "curvature memory" be screened? (ROADMAP Phase 8, started 2026-10-24)

**Question** (THEORY_NOTES 4.5). Local gravity needs the CMSTG field at its attractor today. A late-time,
DESI-like signal needs the field still evolving. Screening can in principle separate the two by hiding the field
in dense places (the Solar System) while letting it evolve cosmologically. Can CMSTG be given a screened form that
produces a DESI-sized late-time effect?

**Short answer** (derived below for CMSTG's own F(Ψ)R form, consistent with the literature):
- **Thin-shell screening** (symmetron/chameleon) is the natural screened form of the CMSTG attractor. It can screen the Solar System only if the cosmological effect is ≲ 1e-5. No-go.
- **Vainshtein screening** (adding a genuine cubic Galileon term) screens Cassini easily, but it cannot screen the *time* variation of G. Lunar Laser Ranging then caps the running of the Planck mass today at |α_M| ≲ 0.002 per e-fold. A late-time effect driven by the F(Ψ)R coupling is then limited to ~0.1% of H(z). Larger effects require either tuning (F freezing exactly now) or that the coupling plays no role, in which case the model is just quintessence.
- **Ways out need a different theory**: non-universal (dark-sector-only) coupling, theories in which the time variation of G is also screened, or accepting fine-tuned narrow bands.

The quantitative check is SIM156 (pre-registered in ROADMAP Phase 8).

Units: 8πG = 1 (reduced Planck mass M_pl = 1), GR is F = ½, Ψ canonical. Field equation from
S = ∫√-g [F(Ψ)R − ½(∂Ψ)² − V(Ψ)] + S_m:

    □Ψ = ∂V_eff/∂Ψ,   V_eff(Ψ; R) ≡ V(Ψ) − F(Ψ) R

For non-relativistic matter R ≈ ρ/(2F). So R is the local matter density in curvature units: ≈ 1e5 × cosmic in the
Milky Way, ≈ 1e30 × cosmic in the Sun.

---

## 1. Thin-shell screening: the CMSTG attractor *is* a symmetron in the curvature
Family A is exactly quadratic around its attractor:

    F(Ψ) = ½ − Λ0Ψ² + ξΨ = F* − Λ0 δΨ²,   δΨ = Ψ − Ψ*,  Ψ* = ξ/(2Λ0),  F* = ½ + ξ²/(4Λ0).

Add the simplest potential that destabilises Ψ* at low curvature, V = V0 − ½μ²δΨ² + ¼λδΨ⁴:

    V_eff = V0 − F*R + (Λ0 R − ½μ²) δΨ² + ¼λ δΨ⁴

This is the symmetron (Hinterbichler & Khoury 2010) with the matter density replaced by R:
- For R > R_crit = μ²/(2Λ0) the minimum is δΨ = 0. There F' = 0: GR, no fifth force (screened).
- For R < R_crit the symmetry breaks: δΨ_min² = (μ² − 2Λ0R)/λ. The field couples, and G changes by ΔF = Λ0 δΨ².

"Curvature memory" in this form: the field is held at the GR point while the universe (or a region) is dense, and
released when the curvature drops below R_crit. That is the late-time trigger we want, and it is screened in dense places.

**Einstein-frame dictionary.** A² = 1/(2F), so ln A ≈ const + Λ0δΨ²/(2F*).
- This is the symmetron coupling A ≈ 1 + δΨ²/(2M²) with **M² = F*/Λ0 ≈ 1/(2Λ0)** (in M_pl² units).
- The vacuum coupling strength is β = δΨ_vac/M² ≈ 2Λ0 δΨ_vac. The cosmological change of G is ΔF ≈ 2F*·Δln A = β²F*²/Λ0 ≈ β²/(4Λ0).

**Screening condition (thin shell).**
- A body of Newtonian potential Φ_N is screened if Δφ/(6βΦ_N) ≪ 1, with Δφ = δΨ_vac (Khoury & Weltman 2004).
- For the symmetron coupling β ∝ δΨ_vac, so δΨ_vac cancels: **screened iff Φ_N ≳ M²/6 = F*/(6Λ0) ≈ 1/(12Λ0).**
- The Sun has Φ_N = 2.1e-6 and the Milky Way ~1e-6. Screening them needs **Λ0 ≳ 8e4**.

**Consequence (the WHK no-go, made explicit for CMSTG).**
- With Λ0 ≳ 8e4 and a coupling no stronger than gravity in unscreened regions (β ≲ 1), the cosmological change of G is ΔF ≈ β²/(4Λ0) ≲ **3e-6**.
- A DESI-relevant ΔF ~ 1e-2 would need β ≈ √(4Λ0ΔF) ≈ 57. That is a fifth force 2β² ≈ 6500× gravity in every unscreened environment within the Compton wavelength.
- That Compton wavelength is m ≈ √(4Λ0R_crit) ~ 2000 H0, i.e. ≈ 2 Mpc, as WHK require. Voids, dwarf galaxies and the dynamics of unscreened bodies exclude such a force.
- This reproduces theorems (i) and (ii) of Wang, Hui & Khoury (PRL 109, 241301, 2012): Compton wavelength ≲ 1 Mpc, and the conformal factor constant over the last Hubble time.
- **Thin-shell-screened curvature memory cannot change the cosmological background.** The same holds for chameleon (runaway-potential) versions.

## 2. Vainshtein screening: CMSTG + a genuine cubic Galileon
Note that the project's earlier "Galileon" SIM142 used G3(Ψ)□Ψ. That equals −G3'(Ψ)(∂Ψ)² up to a total derivative, a
rescaled kinetic term with **no** screening (AUDIT_NOTES; its Cassini failure stands). A genuine cubic Galileon depends on X:

    L ⊃ −(c3/Λ3³) X □Ψ,   X = −½(∂Ψ)²,   Λ3³ = M_pl H0²  (= H0² in our units)

- **c_T = 1.** G4 = F(Ψ) has no X dependence and G3 does not change c_T. GW170817 is satisfied.
- **Cassini is screened.** The Vainshtein radius of the Sun is r_V ~ (c3 β r_s/H0²)^{1/3} ≈ 120 pc × (c3β)^{1/3}.
  - Inside it the fifth force is suppressed by ~(r/r_V)^{3/2}. At 1 AU that is ~1e-11/(c3β)^{1/2}.
  - So γ − 1 is negligible for any c3 ≳ 1e-12. The attractor condition F'(Ψ0) = 0 is *no longer needed* for Cassini.
- **Braiding is small.** α_B ~ c3 (dΨ/dN)³ (H/H0)². For c3 ≲ 1 and dΨ/dN ≲ 0.2 this gives α_B ≲ 1e-2. That is the bound above which gravitational waves from binaries trigger dark-energy instabilities (Creminelli, Tambalo, Vernizzi & Yingcharoenrat, JCAP 05 (2020) 002).
- **Ġ is not screened.** Vainshtein hides static, short-distance effects, but not the cosmological time evolution of the field (Babichev, Deffayet & Esposito-Farèse, PRL 107, 251102, 2011).
  - Inside the Solar System G_local(t) ∝ 1/F(Ψ_cosmo(t)).
  - Lunar Laser Ranging then gives **|α_M(z=0)| = |d ln F/d ln a|₀ ≲ 0.002** (Burrage & Dombrowski, JCAP 07 (2020) 060; this applies to Vainshtein and kinetic, K-mouflage-type, screening with shift symmetry).

So with Vainshtein screening the binding local constraint moves from Cassini (F' = 0 today) to LLR (Ḟ ≈ 0 today).

## 3. What |α_M| ≤ 0.002 leaves for a late-time curvature-memory effect
Three ways to make F(Ψ)R produce late-time dynamics:

**(a) Adiabatic late release.** A heavy field tracks its moving minimum. With a tilt V = V0 + V1 δΨ, the minimum is δΨ = −V1/(2Λ0R), so ΔF = Λ0δΨ² ∝ R⁻².
- Today d ln R/dN = −3Ω_m/(Ω_m + 4Ω_Λ) = −0.30, so α_M ≈ −1.2 ΔF.
- **LLR gives ΔF₀ ≲ 0.0017.** Since R(z=1)/R₀ = 1.71, F changes by only 0.66ΔF₀ ≲ 0.0011 between z = 1 and today.
- That is a ~0.1% effect on H(z): invisible to DESI.

**(b) Light (thawing) field.** The field rolls with y = dΨ/dN.
- A DESI-sized w₀ + 1 ≈ 0.15 needs y₀² ≈ 3Ω_DE(1 + w₀) ≈ 0.3.
- Then α_M = F'y/F ≤ 0.002 forces |F'| ≲ 0.002 along the trajectory: the field must roll where F is flat to 0.2%.
- The non-minimal coupling then plays no role. The model is minimally coupled thawing quintessence, which **cannot cross the phantom divide**, the part of the DESI w0wa preference that needs α_M or α_B ≠ 0.

**(c) Transient.** F evolves at z ~ 0.5–2 and freezes by today. That needs α_M to fall from O(0.01–0.1) to < 0.002 within about half an e-fold, i.e. F reaching a turning point *now*.
- The natural realisation is an oscillation about the attractor, which leaves an oscillating α_M ∝ δΨ·dδΨ/dN.
- After a transition at z_t ~ 1–2 its amplitude keeps α_M ≫ 0.002 today, except at a phase-tuned zero crossing. A coincidence.

**Conclusion.** For a universally coupled F(Ψ)R field, combining SIM155 with the above:
- *early* curvature memory is BBN-limited (Δχ² ≈ +3);
- *late* curvature memory is limited by WHK (thin-shell screening) or by LLR (Vainshtein/kinetic screening) to sub-percent effects;
- unless it is either tuned or reduced to ordinary quintessence.

This agrees with the recent literature, which finds that non-minimally coupled dark energy fitting DESI and passing
Solar System and Ġ bounds lives in **very narrow** parameter bands (arXiv:2606.14014). Alternatively it survives
without screening by keeping G_eff within 0.3% of G throughout (arXiv:2609.10133). See also Wolf, Ferreira &
García-García, PRD 111, L041303 (2025), on the improved DESI fit of non-minimally coupled dark energy.

## 4. Ways out (each a different theory, not a tweak)
1. **Dark-sector-only coupling.** If Ψ couples to dark matter but not to baryons, Cassini and LLR do not apply ("coupled dark energy"). Constraints then come from the CMB and large-scale structure. This drops CMSTG's universal F(Ψ)R form.
2. **Screening of the time variation.** Theories (beyond Horndeski / DHOST) in which the local G decouples from the cosmological evolution. A research frontier, with GW and stability constraints to check.
3. **Tuned bands.** Accept a coincidence (F stationary today) and quantify the tuning, as the 2026 papers above do.

## 5. Quantitative check: SIM156 (pre-registered in ROADMAP Phase 8)
The cleanest test of the open branch (Vainshtein-screened, LLR-limited) is a model that *is* curvature memory with a
late trigger: **"pinned thawing" CMSTG**.
- The field starts at the attractor Ψ* (GR early on, so BBN and the CMB are untouched).
- A tilt V1 pushes it off once the curvature restoring force 2Λ0δΨR weakens at late times.
- Fit it to background data: compressed Planck 2018 distance priors (Chen, Huang & Wang, JCAP 02 (2019) 028) + DESI DR2 + Pantheon+.
- Compare with ΛCDM, w0waCDM and minimally coupled thawing quintessence, under the LLR bound.
- Section 3 predicts that within |α_M| ≤ 0.002 it does no better than thawing quintessence.

## References (verified 2026-10-24)
- J. Wang, L. Hui, J. Khoury, *No-go theorems for generalized chameleon field theories*, PRL 109, 241301 (2012), arXiv:1208.4612.
- E. Babichev, C. Deffayet, G. Esposito-Farèse, *Constraints on shift-symmetric scalar-tensor theories with a Vainshtein mechanism from bounds on the time variation of G*, PRL 107, 251102 (2011), arXiv:1107.1569.
- C. Burrage, J. Dombrowski, *Constraining the cosmological evolution of scalar-tensor theories with local measurements of the time variation of G*, JCAP 07 (2020) 060, arXiv:2004.14260.
- P. Creminelli, G. Tambalo, F. Vernizzi, V. Yingcharoenrat, *Dark-energy instabilities induced by gravitational waves*, JCAP 05 (2020) 002, arXiv:1910.14035.
- K. Hinterbichler, J. Khoury, *Screening long-range forces through local symmetry restoration*, PRL 104, 231301 (2010).
- J. Khoury, A. Weltman, *Chameleon fields*, PRL 93, 171104 (2004).
- L. Chen, Q.-G. Huang, K. Wang, *Distance priors from Planck final release*, JCAP 02 (2019) 028, arXiv:1808.05724: R = 1.7502 ± 0.0046, l_A = 301.471 ± 0.090, ω_b = 0.02236 ± 0.00015; correlations R–l_A +0.46, R–ω_b −0.66, l_A–ω_b −0.33 (read from Table I).
- W. J. Wolf, P. G. Ferreira, C. García-García, PRD 111, L041303 (2025), arXiv:2409.17019.
- arXiv:2606.14014, *Constraints on stable scalar–tensor dark energy from DESI data and Solar System tests* (2026).
- arXiv:2609.10133, *Safe phantom divide crossing from unscreened non-minimal coupling to gravity* (2026).
