// CMSTG repository data, extracted faithfully from README.md
// Source: github.com/cisomorph/Curvature-Memory_Scalar-Tensor_Gravity

window.CMSTG_DATA = {
  meta: {
    title: "Curvature-Memory Scalar-Tensor Gravity",
    short: "CMSTG",
    tagline: "A Lagrangian-based framework for curvature memory and modified gravity",
    author: "Christopher Robert Barrick Wilson",
    contact: "C.Isomorph@gmail.com",
    repo: "https://github.com/cisomorph/Curvature-Memory_Scalar-Tensor_Gravity",
    version: "1.2.1",
    released: "2026-05-14",
    license: "CC BY-NC 4.0",
  },

  // Canonical parameters from README
  parameters: [
    { sym: "\\Lambda_0",        val: "0.003",  unit: "M_\\mathrm{Pl}^{-2}", status: "locked",     note: "Non-minimal coupling" },
    { sym: "\\bar\\Psi",         val: "2.62",   unit: "M_\\mathrm{Pl}",      status: "attractor",  note: "Cosmological VEV" },
    { sym: "F_\\mathrm{eff}",    val: "0.521",  unit: "",                    status: "locked",     note: "(1+2\\Lambda_0\\bar\\Psi^2)/2" },
    { sym: "H_0",                val: "67.59",  unit: "km/s/Mpc",            status: "best-fit",   note: "Joint CMB+BAO" },
    { sym: "\\sigma_\\mathrm{DESI}", val: "2.77", unit: "\\sigma",          status: "structural floor", note: "DESI Y1 BAO residual" },
    { sym: "\\ln B",             val: "-0.71",  unit: "",                    status: "inconclusive", note: "Bayes factor vs ΛCDM" },
  ],

  // 4 papers + master
  papers: [
    {
      id: "I",
      slug: "paper1_nogo",
      repoPath: "Papers/paper1_nogo",
      title: "Two Structural No-Go Theorems for Late-Time Modifications",
      sims: "SIM131–144",
      verdict: "FAIL / NO-GO",
      verdictKind: "fail",
      explainer: "explainers/paper1-nogo.html",
      summary: "Establishes two structural no-go theorems that together prove the 2.77σ DESI Y1 BAO residual cannot be reduced within the CMSTG action class without pre-recombination new physics.",
    },
    {
      id: "II",
      slug: "paper2_framework",
      repoPath: "Papers/paper2_framework",
      title: "Framework, Field Equations, and Observational Tests",
      sims: "SIM80–111",
      verdict: "PASS (2.77σ floor)",
      verdictKind: "pass",
      explainer: "explainers/framework.html",
      summary: "The Phase 1 canonical action: variation, field equations, and 22 observational tests covering BAO, CMB, RSD, GW speed, and Solar System bounds.",
    },
    {
      id: "III",
      slug: "paper3_uv",
      repoPath: "Papers/paper3_uv",
      title: "UV Finiteness and the Λ₀ Fixed Point",
      sims: "SIM102–106",
      verdict: "PASS",
      verdictKind: "pass",
      explainer: "explainers/paper3-uv.html",
      summary: "Resummed memory propagator, two-loop UV finiteness, negative β-function driving Λ₀→0 in the UV, and the kernel-frame scope addressed in §2.2(d).",
    },
    {
      id: "IV",
      slug: "paper4_galactic",
      repoPath: "Papers/paper4_galactic",
      title: "Galactic-Scale Constraints",
      sims: "SIM99–100, 103",
      verdict: "FAIL",
      verdictKind: "fail",
      explainer: "explainers/paper4-galactic.html",
      summary: "Ψ-sector rotation curves: the Ψ field cannot replace dark matter at galactic scales. Phase 2 χ-DM is the surviving route.",
    },
  ],

  theorems: [
    {
      n: "1",
      phase: "Paper I — Phase 3",
      title: "Curvature-sourced scalars grow the wrong way",
      body: "Any scalar field φ sourced by curvature R > 0 and initialized at zero grows monotonically, increasing F_eff and thereby suppressing H(z) at DESI redshifts. This is the wrong direction.",
    },
    {
      n: "2",
      phase: "Paper I — Phase 4",
      title: "DESI fix and CMB preservation are mutually exclusive",
      body: "Any mechanism that raises H(z) on z ∈ [0, z_rec] with fixed sound horizon r_s necessarily increases θ* = r_s/D_C* above the Planck bound 100θ* = 1.04101 ± 0.00029. Pure late-time modifications cannot satisfy both.",
    },
  ],

  loophole: "Simultaneously increase r_s via pre-recombination physics (early dark energy, extra relativistic species). Separately constrained by Planck CMB power spectra; constitutes a distinct theoretical programme (Phase 5).",

  // Phases of simulation work, faithful to README structure
  phases: [
    {
      key: "phase1",
      label: "Phase 1",
      title: "Locked Action",
      range: "SIM80–SIM111",
      summary: "Phase 1 canonical action established. 22 observational tests; 2.77σ DESI Y1 residual identified as structural floor.",
      sims: [
        { id: "SIM80",      topic: "Resummed memory propagator / UV falloff",          verdict: "PASS" },
        { id: "SIM82",      topic: "Retarded kernel causality and Yukawa profile",     verdict: "PARTIAL" },
        { id: "SIM83–86",   topic: "Memory cutoff, coupling survey, Friedmann stability", verdict: "PASS" },
        { id: "SIM87",      topic: "BAO full-covariance refit (BOSS + eBOSS + Lyα)",   verdict: "PASS" },
        { id: "SIM88",      topic: "GW speed cᴛ = c; Solar System G_eff",              verdict: "PASS" },
        { id: "SIM89",      topic: "Planck acoustic-peak shift diagnostic",            verdict: "PASS" },
        { id: "SIM90",      topic: "Joint CMB+BAO parameter fit (Δχ² ≈ 0)",            verdict: "PASS" },
        { id: "SIM91–96",   topic: "RSD, σ₈, lensing, BBN, neutrinos, Ly-α",           verdict: "PASS" },
        { id: "SIM97",      topic: "Full Planck plikHM TTTEEE (BAO-only params)",     verdict: "FAIL", detail: "Δχ² = +7.0" },
        { id: "SIM98",      topic: "Joint plikHM + BAO refit",                         verdict: "PASS", detail: "Δ(−2lnL) = −0.013" },
        { id: "SIM99–100",  topic: "Ψ-sector rotation curves (two approaches)",        verdict: "FAIL" },
        { id: "SIM101–106", topic: "BAO rᴅ, two-loop UV, Ward identity, RG flow",      verdict: "PASS" },
        { id: "SIM107–108", topic: "Solar System, GW speed",                           verdict: "PASS" },
        { id: "SIM109–110", topic: "DESI w(z), inflation",                             verdict: "FAIL" },
        { id: "SIM111",     topic: "m₀ scan — DESI tension floor 3.44σ",               verdict: "STRUCTURAL" },
      ],
    },
    {
      key: "phase2",
      label: "Phase 2",
      title: "Matter Sector",
      range: "SIM112–SIM130",
      summary: "Matter-sector extensions exhausted. χ-DM oscillating field passes SPARC 161-galaxy fit; joint DE+DM converges on 2.77σ floor.",
      sims: [
        { id: "SIM112–113",  topic: "λΨ⁴ and hilltop quintessence",                 verdict: "PARTIAL" },
        { id: "SIM114–115",  topic: "χ-DM condensate and gradient soliton",          verdict: "FAIL" },
        { id: "SIM116–119",  topic: "χ-DM oscillating field, halos, SPARC 161-gal",  verdict: "PASS" },
        { id: "SIM120–121C", topic: "Joint DE+DM, Ly-α, H₀, DESI+Planck MCMC",       verdict: "PARTIAL", detail: "2.77σ floor" },
        { id: "SIM122–130",  topic: "Phase 3 probes (matter sector)",                verdict: "FAIL" },
      ],
    },
    {
      key: "phase3",
      label: "Phase 3",
      title: "Curvature-Memory Scalars",
      range: "SIM131–SIM136",
      summary: "All six curvature-sourced mechanisms fail. Structural no-go (Theorem 1, Paper I) established.",
      sims: [
        { id: "SIM131", topic: "Additive ξΨR (conformal)",                  verdict: "FAIL", detail: "5.97σ" },
        { id: "SIM132", topic: "Deser–Woodard analog",                      verdict: "FAIL", detail: "3.59–3.83σ" },
        { id: "SIM133", topic: "Gauss–Bonnet sourcing",                     verdict: "FAIL", detail: "9.68σ" },
        { id: "SIM134", topic: "Dilaton / exponential coupling",            verdict: "FAIL", detail: "3.59–3.84σ" },
        { id: "SIM135", topic: "Bi-scalar (Ψ frozen, φ sourced by R)",      verdict: "FAIL", detail: "3.73–5.60σ" },
        { id: "SIM136", topic: "Horndeski Gᵐᵘ∂Ψ² kinetic",                  verdict: "FAIL", detail: "3.53–3.75σ" },
      ],
    },
    {
      key: "phase4",
      label: "Phase 4",
      title: "Late-Time Extensions",
      range: "SIM137–SIM144",
      summary: "Second structural no-go (Theorem 2) established. Tier 2 mechanism space exhaustively covered by SIM144 completeness probe.",
      sims: [
        { id: "SIM137", tier: "T1", topic: "SPARC failure-mode analysis",                                 verdict: "STRUCTURAL_PATTERN" },
        { id: "SIM138", tier: "T1", topic: "DESI per-bin sensitivity",                                    verdict: "DISTRIBUTED" },
        { id: "SIM139", tier: "T1", topic: "RSD shape diagnostic",                                        verdict: "SHAPE_FIXABLE" },
        { id: "SIM140", tier: "T2", topic: "Step potential (P4-A)",                                       verdict: "PREDICTED FAIL", detail: "not run" },
        { id: "SIM141", tier: "T2", topic: "Running Λ₀(a) / BD analog (P4-B)",                            verdict: "FAIL", detail: "θ* 63σ" },
        { id: "SIM142", tier: "T2", topic: "Galileon G₃(Ψ)□Ψ (P4-C)",                                     verdict: "FAIL" },
        { id: "SIM143", tier: "T2", topic: "Bi-scalar Ψ+φ quintessence (P4-D)",                           verdict: "FAIL" },
        { id: "SIM144", tier: "T2", topic: "Mixed-source φ (ξ_R + β_m) — completeness probe (P4-E)",      verdict: "FAIL" },
        { id: "SIM145–146", tier: "T3", topic: "UV recheck / predictions",                                verdict: "DEFERRED", detail: "no Tier 2 PASS" },
      ],
    },
  ],

  predictions: [
    {
      n: "01",
      title: "DESI Y3 BAO tension",
      milestone: "≈ 2027",
      body: "If the 2.77σ floor persists at ≥3σ with more data, CMSTG requires pre-recombination new physics (Phase 5). If tension falls below 2σ, Phase 1 is fully validated.",
    },
    {
      n: "02",
      title: "χ-DM mass",
      milestone: "next-gen stellar kinematics",
      body: "m₂₂ ≈ 0.28–0.34 for low-mass gas-dominated dwarfs and LSBs (v_flat ≲ 120 km/s). Testable by next-generation stellar kinematics.",
    },
    {
      n: "03",
      title: "G_eff / G",
      milestone: "future GW standard sirens",
      body: "G_eff / G = 1 − 15 ppm at z = 0 — negligible at current precision but in principle testable by future GW standard-siren surveys.",
    },
  ],
};
