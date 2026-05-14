#!/usr/bin/env python3
"""
Generate results.md for every simulation in Ordered_Simulations.
Reads existing JSON/README data from each sim's Outputs/ folder.
"""
import json, os, glob, re

BASE = os.path.dirname(os.path.abspath(__file__))

# ──────────────────────────────────────────────────────────────────────────────
# Sim catalogue: id → (phase, short_title, description)
# ──────────────────────────────────────────────────────────────────────────────
SIM_META = {
    "SIM80":  ("Phase 1", "Recursive UV Finiteness",
               "UV falloff and bubble-integral cutoff suppression with memory self-energy"),
    "SIM82":  ("Phase 1", "Retarded Green's Function Kernel",
               "Causality test: retarded kernel propagation, Yukawa falloff, memory accumulation"),
    "SIM83":  ("Phase 1", "Ψ Field Halo Solver",
               "Static NFW/halo Ψ profile with memory source; gradient stability check"),
    "SIM84":  ("Phase 1", "Ward Identity / Graviton Sector",
               "Π_hh(0)=0 Ward identity; massless graviton; c_T verification"),
    "SIM85":  ("Phase 1", "CMB TT Power Spectrum",
               "CMB TT Cℓ from CLASS with CMSTG Ψ coupling vs Planck 2018"),
    "SIM86":  ("Phase 1", "EE/TE Polarisation Spectrum",
               "CMB EE+TE spectra; F_eff(z_CMB)=0.521 consistency"),
    "SIM87":  ("Phase 1", "BAO Likelihood (12-point full-covariance)",
               "Full-covariance BAO chi2 (BOSS+eBOSS+DESI); rd and DV/rd model"),
    "SIM88":  ("Phase 1", "G_eff Solar System / GW Speed",
               "G_eff/G correction at lab scale; GW170817 c_T=c bound satisfied"),
    "SIM89":  ("Phase 1", "Planck Acoustic-Peak Shift Test",
               "Δ_l shift from CMSTG F_eff; 30.1σ tension at BAO-only best-fit"),
    "SIM90":  ("Phase 1", "Joint CMB+BAO Parameter Fit",
               "Joint fit resolves SIM89 tension: CMSTG and ΛCDM indistinguishable (Δχ²≈0)"),
    "SIM91":  ("Phase 1", "RSD / Growth Rate",
               "fσ₈(z) predictions vs BOSS/WiggleZ/2dF; growth-rate consistency"),
    "SIM92":  ("Phase 1", "σ₈ / S₈ Tension",
               "S₈ = σ₈√(Ω_m/0.3) at CMSTG canonical parameters"),
    "SIM93":  ("Phase 1", "Lensing Power Spectrum",
               "CMB lensing Cℓ^κκ and weak-lensing convergence power spectrum"),
    "SIM94":  ("Phase 1", "BBN Helium Abundance",
               "Y_p(CMSTG) consistency with BBN; H(z_BBN) constraint"),
    "SIM95":  ("Phase 1", "Neutrino Mass Sum",
               "Σm_ν marginalisation at CMSTG canonical; compatibility with Planck"),
    "SIM96":  ("Phase 1", "Ly-α Forest Power Spectrum",
               "P(k) from Ly-α flux at CMSTG Ψ; UV damping consistency"),
    "SIM97":  ("Phase 1", "Full plikHM TTTEEE Likelihood",
               "Planck plikHM TTTEEE likelihood at CMSTG best-fit; Δ(-2lnL) vs ΛCDM"),
    "SIM98":  ("Phase 1", "Joint plikHM + BAO Refit",
               "Optimised joint plikHM+BAO: CMSTG Δ(-2lnL)=−0.013 vs ΛCDM → PASS"),
    "SIM99":  ("Phase 1", "Rotation Curves — NFW DM + Ψ halo",
               "SPARC rotation curves with NFW+Ψ sub-component; FAIL for most galaxies"),
    "SIM100": ("Phase 1", "Rotation Curves — Ψ-only halo",
               "Pure Ψ soliton rotation curves; FAIL — insufficient DM density"),
    "SIM101": ("Phase 1", "Memory Kernel Convergence",
               "Convergence test on retarded kernel integrator; step-halving check"),
    "SIM102": ("Phase 1", "Two-loop UV Finiteness",
               "Two-loop self-energy diagram; memory kernel suppresses quartic divergences"),
    "SIM103": ("Phase 1", "Λ₀ RG Flow / UV Fixed Point",
               "β-function for Λ₀; UV fixed point at Λ₀=0 confirmed"),
    "SIM104": ("Phase 1", "Four-diagram UV Check",
               "Four 1PI diagrams; Ward identity Π_hh(0)=0 verified at 1-loop"),
    "SIM105": ("Phase 1", "Λ₀ Running / RG Flow and Fixed Point",
               "Λ₀(μ) RG flow; negative beta function; Λ₀ fixed point at Λ₀=0 (GR limit)"),
    "SIM106": ("Phase 1", "Phase Space / Attractor",
               "Ψ phase-space portrait; late-time attractor at Ψ̄=2.62 M_Pl"),
    "SIM107": ("Phase 1", "Perturbation Theory: δΨ",
               "Linear perturbations δΨ; no tachyon (ω²>0 for all modes)"),
    "SIM108": ("Phase 1", "GR Weak-Field Limit",
               "Post-Newtonian expansion; GR recovered as Λ₀→0, V→0"),
    "SIM109": ("Phase 1", "DESI Y1 BAO Full Fit",
               "DESI Y1 6-bin BAO fit at Phase 1 canonical; 2.77σ residual identified"),
    "SIM110": ("Phase 1", "DESI BAO w₀wₐ Contours",
               "w₀–wₐ contours from DESI+CMB; CMSTG vs ΛCDM model comparison"),
    "SIM111": ("Phase 1", "Dynamical m₀ Scan / DESI Floor",
               "Scan over m₀ to test DESI w₀wₐ tension; minimum tension=3.44σ — structural floor"),
    "SIM112": ("Phase 2", "Quintessence DE (λΨ⁴)",
               "Phase 2 opening: quintessence from λΨ⁴ potential; PARTIAL — DESI 7.74σ"),
    "SIM113": ("Phase 2", "Quintessence + Ψ Tilt",
               "λ(Ψ²−v²)² potential; best-fit w₀=−0.973, wₐ=−0.41; PARTIAL"),
    "SIM114": ("Phase 2", "χ-DM Condensate (structural trilemma)",
               "Ultra-light χ condensate DM; FAIL — structural trilemma"),
    "SIM115": ("Phase 2", "Gradient Soliton DM",
               "Gradient-dominated soliton χ; FAIL — H₀² suppression structural"),
    "SIM116": ("Phase 2", "χ-DM Oscillating Field",
               "Oscillating χ field DM with κ coupling to Ψ; PASS density match"),
    "SIM117": ("Phase 2", "χ-DM Galactic Soliton Profile",
               "NFW+soliton χ halo; mass-radius relation vs observational constraint"),
    "SIM118": ("Phase 2", "χ-DM Halo Mass Function",
               "Halo mass function from χ-DM; consistency with large-scale structure"),
    "SIM119": ("Phase 2", "SPARC 161-galaxy χ-DM Fit",
               "χ-DM rotation curve fit on full SPARC dataset; 65 PASS / 95 MARGINAL"),
    "SIM120": ("Phase 2", "Joint DE+DM Background Consistency",
               "SIM113 Ψ quintessence + SIM119 χ-DM background; DE-DM decoupling confirmed"),
    "SIM121": ("Phase 2", "Ly-α Power Spectrum (χ-DM)",
               "Ly-α P(k) with χ-DM; CONDITIONAL — f_FDM=1 ruled out"),
    "SIM121B": ("Phase 2", "H₀ Tension via G_eff",
               "G_eff(z) variation to address H₀ tension; FAIL — CMB worsens"),
    "SIM121C": ("Phase 2", "Joint DESI+Planck MCMC",
               "Full joint DESI+Planck MCMC at Phase 2 canonical; PARTIAL — 2.77σ floor"),
    "SIM122": ("Phase 2", "Phase 3 Direction A: Ψ² coupling",
               "Non-minimal Ψ² curvature coupling; FAIL — exceeds GR bound"),
    "SIM123": ("Phase 2", "Phase 3 Direction B: Inverse Ψ",
               "1/Ψ coupling extension; FAIL — GR recovery broken"),
    "SIM124": ("Phase 3", "Kinetic Ψ Modification",
               "Non-canonical kinetic term f(X)∂Ψ²; FAIL — no viable solutions"),
    "SIM125": ("Phase 3", "Memory Field M(a)",
               "Dynamical memory field M(a); FAIL — DESI unchanged"),
    "SIM126": ("Phase 3", "Exponential Memory Decay",
               "e^{−βM} coupling; FAIL — insufficient H(z) modification"),
    "SIM127": ("Phase 3", "G_eff/G = 1−νM̂(a) Coupling",
               "Linear memory-gravity coupling; PARTIAL — RSD improvement only"),
    "SIM128": ("Phase 3", "Power-law G_eff/G = 1−νM̂^p",
               "Power-law growth profile; PARTIAL — slope failure in RSD"),
    "SIM129": ("Phase 3", "Memory-modulated Λ_eff",
               "Λ_eff=Λ₀(1−γM̂); FAIL — M-field DOF exhausted"),
    "SIM130": ("Phase 3", "Curvature-Sourced Ψ Evolution",
               "Ψ driven by curvature source; FAIL — M-field programme closed"),
    "SIM131": ("Phase 3", "Phase 3 Direction C: ξΨR coupling (ξ=1/6)",
               "Conformal coupling F_eff=(1+2Λ₀Ψ²)/2+ξΨ; FAIL — θ* violation 5.97σ"),
    "SIM132": ("Phase 3", "ξΨR scan over ξ",
               "Scan over ξ ∈ [0.01, 1]; FAIL — all ξ violate CMB θ*"),
    "SIM133": ("Phase 3", "Phase 3 Direction D: Gauss-Bonnet coupling",
               "αΨ G coupling to Gauss-Bonnet invariant G; FAIL — structural"),
    "SIM134": ("Phase 3", "Gauss-Bonnet parameter scan",
               "Scan over GB coupling strength α; FAIL — CMB violation systematic"),
    "SIM135": ("Phase 3", "Phase 3 Direction E: k-essence Ψ",
               "k-essence P(X) kinetic modification; FAIL — no new mechanism"),
    "SIM136": ("Phase 3", "Horndeski G^{μν}∂Ψ² coupling",
               "G^{μν}∂_μΨ∂_νΨ kinetic Horndeski term; FAIL — G₃_X=0 equivalent to rescaled kinetic"),
    "SIM137": ("Phase 4", "SPARC Failure-Mode Analysis",
               "Tier 1 diagnostic: χ-DM failure morphology; structural pattern in m₂₂"),
    "SIM138": ("Phase 4", "DESI Y1 Per-Bin Sensitivity",
               "Tier 1 diagnostic: DESI tension distributed across z=0.5–1.3 (LRG1/2/3+ELG)"),
    "SIM139": ("Phase 4", "SIM128 RSD Shape Diagnostic",
               "Tier 1 diagnostic: slope failure in fσ₈(z); phenomenological 2-param fix"),
    "SIM140": ("Phase 4", "Symmetron Step Potential (not run)",
               "P4-C: symmetron V(Ψ) step; predicted FAIL by SIM141/143 structural argument"),
    "SIM141": ("Phase 4", "Running Λ₀(a) — Brans-Dicke Analog",
               "P4-B: phenomenological Λ₀(a) decay; PARTIAL — θ* 63σ CMB violation"),
    "SIM142": ("Phase 4", "Galileon G₃(Ψ)□Ψ Sector",
               "P4-A: Horndeski G₃ term; FAIL — G₃_X=0 suppressed by slow-roll"),
    "SIM143": ("Phase 4", "Bi-scalar Ψ + decoupled φ",
               "P4-D: independent quintessence φ not sourced by R; FAIL — same anti-correlation"),
    "SIM144": ("Phase 4", "Tier 3: Joint MCMC on winner (not run)",
               "Would run if Tier 2 produced a PASS candidate; deferred — all Tier 2 failed"),
    "SIM145": ("Phase 4", "Tier 3: UV consistency recheck (not run)",
               "Would recheck UV finiteness of winning mechanism; deferred"),
    "SIM146": ("Phase 4", "Tier 3: Distinctive prediction extraction (not run)",
               "Would extract 2–3 testable predictions for DESI Y3/Euclid; deferred"),
}

# ──────────────────────────────────────────────────────────────────────────────
# Phase 1 canonical (for reference in results.md)
# ──────────────────────────────────────────────────────────────────────────────
CANONICAL_NOTE = "Phase 1 canonical: Λ₀=0.003, Ψ̄=2.62 M_Pl, F_eff(z_CMB)=0.521, H₀≈67.59 km/s/Mpc."


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def find_json(sim_dir):
    """Return the first JSON found in Outputs/, preferring *_results.json."""
    out = os.path.join(sim_dir, "Outputs")
    if not os.path.isdir(out):
        # Phase 4 sims have flat structure
        j = os.path.join(sim_dir, "output.json")
        if os.path.isfile(j):
            return j
        return None
    # Prefer *_results.json over *_diagnostics.json
    for pat in ["*_results.json", "*_diagnostics.json", "*.json"]:
        matches = glob.glob(os.path.join(out, pat))
        if matches:
            return matches[0]
    return None


def extract_verdict(d):
    if isinstance(d, dict):
        v = d.get("verdict")
        if isinstance(v, dict):
            v = v.get("overall", v.get("pass", "N/A"))
        if v:
            return str(v)
        # Try nested
        if "passed" in d:
            return "PASS" if d["passed"] else "FAIL"
        if "all_pass" in d:
            return "PASS" if d["all_pass"] else "FAIL"
        if "viable" in d:
            return "VIABLE" if d["viable"] else "FAIL"
        ac = d.get("acceptance_checks", {})
        if ac:
            all_p = all(v.get("passed", False) for v in ac.values() if isinstance(v, dict))
            return "PASS" if all_p else "PARTIAL"
    return "N/A"


def extract_key_numbers(d, sim_id):
    """Extract the most important numbers from the JSON."""
    nums = []
    if not isinstance(d, dict):
        return nums
    # Common keys
    for key, label in [
        ("tension", "DESI tension"), ("DESI_tension_sigma", "DESI tension"),
        ("chi2_DESI", "DESI χ²"), ("chi2_dof", "χ²/N"),
        ("chi2_total", "Total χ²"), ("delta_chi2", "Δχ²"),
        ("delta_m2lnL_plikHM_cmstg_minus_lcdm", "Δ(−2lnL) plikHM"),
        ("H0", "H₀"), ("Lambda0", "Λ₀"),
        ("Feff0", "F_eff(z=0)"), ("Feff_CMB", "F_eff(z_CMB)"),
        ("theta_100", "100θ*"), ("theta_star_100", "100θ*"),
        ("s8", "S₈"), ("sigma8", "σ₈"),
    ]:
        if key in d:
            val = d[key]
            if isinstance(val, (int, float)):
                nums.append(f"{label} = {val:.4g}")
    # Comparison sub-dict
    comp = d.get("comparison", {})
    if isinstance(comp, dict):
        for k, lbl in [("tension_sigma", "DESI tension"), ("delta_chi2", "Δχ²")]:
            if k in comp:
                nums.append(f"{lbl} = {comp[k]:.4g}")
    # CMSTG sub-dict
    cmstg = d.get("cmstg", d.get("CMSTG_joint", {}))
    if isinstance(cmstg, dict):
        for k, lbl in [("H0", "H₀(CMSTG)"), ("Lambda0", "Λ₀(best-fit)"),
                       ("chi2_total", "Total χ²(CMSTG)")]:
            if k in cmstg:
                nums.append(f"{lbl} = {cmstg[k]:.4g}")
    # Verdict sub-dict
    vd = d.get("verdict", {})
    if isinstance(vd, dict) and "summary" in vd:
        nums.append(f"Summary: {vd['summary'][:120]}")
    return nums[:8]


def read_readme(sim_dir):
    """Read README from Outputs or root of sim dir."""
    for pat in [
        os.path.join(sim_dir, "Outputs", "*.md"),
        os.path.join(sim_dir, "*.md"),
        os.path.join(sim_dir, "Outputs", "README*"),
    ]:
        matches = glob.glob(pat)
        for m in matches:
            if "results.md" not in m.lower():
                try:
                    with open(m) as f:
                        return f.read()
                except Exception:
                    pass
    return None


def make_results_md_from_json(sim_id, meta, json_data):
    phase, title, desc = meta
    verdict = extract_verdict(json_data)
    nums = extract_key_numbers(json_data, sim_id)

    # failure_mode / structural finding
    fm = json_data.get("failure_mode", json_data.get("structural_finding", ""))
    note = json_data.get("note", json_data.get("interpretation", ""))

    lines = [
        f"# {sim_id} — Results",
        "",
        f"**Phase:** {phase}  ",
        f"**Title:** {title}  ",
        f"**Verdict:** {verdict}",
        "",
        "## What was tested",
        "",
        desc + ".",
        "",
    ]

    if nums:
        lines += ["## Key numerical results", ""]
        for n in nums:
            lines.append(f"- {n}")
        lines.append("")

    if fm:
        lines += ["## Structural diagnosis", "", str(fm)[:300], ""]

    if note and len(note) > 10:
        lines += ["## Notes", "", str(note)[:400], ""]

    # Action spec if present
    action = json_data.get("action", json_data.get("action_spec", ""))
    if action and action != "Phase 1 canonical":
        lines += ["## Action tested", "", f"```\n{action}\n```", ""]

    lines += [
        "## Context",
        "",
        CANONICAL_NOTE,
        "",
    ]
    return "\n".join(lines)


def make_results_md_from_readme(sim_id, meta, readme_text):
    phase, title, desc = meta
    lines = [
        f"# {sim_id} — Results",
        "",
        f"**Phase:** {phase}  ",
        f"**Title:** {title}  ",
        "",
        "## What was tested",
        "",
        desc + ".",
        "",
        "## Summary (from README)",
        "",
        readme_text[:1200].strip(),
        "",
        "## Context",
        "",
        CANONICAL_NOTE,
        "",
    ]
    return "\n".join(lines)


def make_results_md_stub(sim_id, meta):
    phase, title, desc = meta
    lines = [
        f"# {sim_id} — Results",
        "",
        f"**Phase:** {phase}  ",
        f"**Title:** {title}  ",
        f"**Verdict:** Not run (spec only)",
        "",
        "## What was planned",
        "",
        desc + ".",
        "",
        "## Status",
        "",
        "This sim was not executed. See SPEC.md for the planned mechanism and success criteria.",
        "",
        "## Context",
        "",
        CANONICAL_NOTE,
        "",
    ]
    return "\n".join(lines)


def write_results_md(sim_id, content):
    sim_dir = os.path.join(BASE, sim_id)
    if not os.path.isdir(sim_dir):
        print(f"  SKIP {sim_id}: directory not found")
        return
    out_path = os.path.join(sim_dir, "results.md")
    with open(out_path, "w") as f:
        f.write(content)
    print(f"  WROTE {out_path}")


def copy_existing_result(sim_id):
    """For Phase 4 sims that already have RESULT.md — copy to results.md."""
    sim_dir = os.path.join(BASE, sim_id)
    src = os.path.join(sim_dir, "RESULT.md")
    dst = os.path.join(sim_dir, "results.md")
    if os.path.isfile(src):
        with open(src) as f:
            content = f.read()
        # Add header if not present
        if not content.startswith("# SIM"):
            content = f"# {sim_id} — Results\n\n" + content
        with open(dst, "w") as f:
            f.write(content)
        print(f"  COPIED RESULT.md → results.md for {sim_id}")
        return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# Old simulation_XX series catalogue
# ──────────────────────────────────────────────────────────────────────────────
OLD_SIM_META = {
    "simulation_01_finalize_lagrangian_recursive_potential":
        ("Pre-SIM80", "Lagrangian Finalisation", "Finalise recursive potential V(Ψ); establish Lagrangian form for CMSTG action"),
    "simulation_02_define_equation_of_state_psi":
        ("Pre-SIM80", "Ψ Equation of State", "Define w(Ψ) for the recursive scalar field; early DE parametrisation"),
    "simulation_03_finalize_mass_quantization_rule":
        ("Pre-SIM80", "Mass Quantisation", "Establish mass quantisation rule from CMSTG field spectrum"),
    "simulation_04_run_gradient_stability_sweep":
        ("Pre-SIM80", "Gradient Stability Sweep", "Sweep parameter space for gradient stability of Ψ field"),
    "simulation_05_validate_GR_weak_field_limit":
        ("Pre-SIM80", "GR Weak-Field Limit", "Validate Newtonian limit and post-Newtonian corrections"),
    "simulation_06_verify_divergence_suppression_mechanism":
        ("Pre-SIM80", "Divergence Suppression", "Verify memory kernel suppresses UV divergences at 1-loop"),
    "simulation_07_cmstg_evolve_FRW":
        ("Pre-SIM80", "FRW Evolution", "Evolve Ψ in flat FRW background; attractor solution"),
    "simulation_08_cmstg_evolve_statichalo":
        ("Pre-SIM80", "Static Halo", "Static halo Ψ profile; NFW comparison"),
    "simulation_09_cmstg_evolve_voids":
        ("Pre-SIM80", "Voids Evolution", "Ψ field in void regions; underdense environment"),
    "simulation_10_cmstg_phase_basinmap":
        ("Pre-SIM80", "Phase Basin Map", "Phase-space basin of attraction for Ψ"),
    "simulation_11_cmstg_gradient_stability":
        ("Pre-SIM80", "Gradient Stability", "Gradient stability check with perturbations"),
    "simulation_12_cmstg_Psi_spectrum":
        ("Pre-SIM80", "Ψ Spectrum", "Power spectrum of Ψ perturbations"),
    "simulation_13_cmstg_CLASS_bridge":
        ("Pre-SIM80", "CLASS Bridge", "Interface CMSTG modified equations to CLASS Boltzmann solver"),
    "simulation_13_2_cmstg_SDE1_kernel":
        ("Pre-SIM80", "SDE Kernel 1", "Stochastic DE kernel — first version"),
    "simulation_13_3_cmstg_SDE2_kernel":
        ("Pre-SIM80", "SDE Kernel 2", "Stochastic DE kernel — second version"),
    "simulation_13_4_cmstg_dark_sector_field":
        ("Pre-SIM80", "Dark Sector Field", "Dark sector field coupling in CMSTG action"),
    "simulation_14_cmstg_C_ell_generator":
        ("Pre-SIM80", "Cℓ Generator", "CMB angular power spectrum generator via CLASS"),
    "simulation_15_cmstg_BAO_shellspace":
        ("Pre-SIM80", "BAO Shell Space", "BAO in 3D shell space; r_d measurement"),
    "simulation_16_cmstg_BAO_projected":
        ("Pre-SIM80", "BAO Projected", "Projected 2D BAO; angular diameter distance"),
    "simulation_17_cmstg_low_ell_CMB":
        ("Pre-SIM80", "Low-ℓ CMB", "Low-multipole CMB; quadrupole and octupole anomalies"),
    "simulation_18_cmstg_graviton_emit":
        ("Pre-SIM80", "Graviton Emission", "Graviton emission rate from Ψ field coupling"),
    "simulation_19_cmstg_graviton_energy":
        ("Pre-SIM80", "Graviton Energy", "Graviton energy density from Ψ → hh process"),
    "simulation_20_cmstg_tensor_mode_propagation":
        ("Pre-SIM80", "Tensor Mode Propagation", "Tensor mode propagation; GW speed verification"),
    "simulation_21_cmstg_damping_vs_emission":
        ("Pre-SIM80", "Damping vs Emission", "Memory damping vs graviton emission rate"),
    "simulation_22_cmstg_voidlens_HEALPix":
        ("Pre-SIM80", "Void Lensing (HEALPix)", "Void lensing signal on HEALPix grid"),
    "simulation_23_cmstg_echo_delaymap":
        ("Pre-SIM80", "Echo Delay Map", "Memory echo delay map for photon propagation"),
    "simulation_24_cmstg_void_shelllens":
        ("Pre-SIM80", "Void Shell Lensing", "Shell-lensing in CMSTG voids"),
    "simulation_25_cmstg_collapse_solver_1D":
        ("Pre-SIM80", "Collapse Solver 1D", "1D gravitational collapse with Ψ field"),
    "simulation_26_cmstg_attractor_saturation":
        ("Pre-SIM80", "Attractor Saturation", "Ψ attractor saturation at Ψ̄=2.62 M_Pl"),
    "simulation_27_cmstg_soft_horizon":
        ("Pre-SIM80", "Soft Horizon", "Soft-horizon structure from Ψ memory feedback"),
    "simulation_28_cmstg_energy_density_profiles":
        ("Pre-SIM80", "Energy Density Profiles", "Ψ energy density ρ(r) profiles"),
    "simulation_29_SM_fit_step1_step2":
        ("Pre-SIM80", "SM Fit Steps 1–2", "Standard Model mass fitting via CMSTG eigenmodes"),
    "simulation_30_SM_fit_step3_to_step5_bridge_particle":
        ("Pre-SIM80", "SM Fit Steps 3–5", "SM particle bridge: charge, spin, mass hierarchy"),
    "simulation_31_cmstg_SM_eigenmodes":
        ("Pre-SIM80", "SM Eigenmodes", "CMSTG field eigenmodes mapped to SM particle content"),
    "simulation_32_cmstg_SM_massfit":
        ("Pre-SIM80", "SM Mass Fit", "CMSTG prediction of SM particle masses"),
    "simulation_33_cmstg_charge_winding":
        ("Pre-SIM80", "Charge Winding", "Topological winding numbers as charge analogues"),
    "simulation_34_cmstg_spin_parity":
        ("Pre-SIM80", "Spin Parity", "Spin and parity from CMSTG field topology"),
    "simulation_35_cmstg_plot_shell_spacing":
        ("Pre-SIM80", "Shell Spacing Plot", "Visualisation of shell spacing in CMSTG structure"),
    "simulation_36_cmstg_plot_CMB_vs_Planck":
        ("Pre-SIM80", "CMB vs Planck Plot", "CMB Cℓ comparison against Planck 2018 data"),
    "simulation_37_cmstg_plot_voidlens":
        ("Pre-SIM80", "Void Lens Plot", "Void lensing signal visualisation"),
    "simulation_38_cmstg_plot_emission_spectra":
        ("Pre-SIM80", "Emission Spectra Plot", "Graviton emission spectra from CMSTG coupling"),
    "simulation_39_cmstg_plot_entropy_growth":
        ("Pre-SIM80", "Entropy Growth Plot", "Entropy growth in CMSTG field evolution"),
    "simulation_40_document_all_modules_GitHub":
        ("Pre-SIM80", "Documentation", "Document all simulation modules for GitHub release"),
    "simulation_41_finalize_overleaf_sections":
        ("Pre-SIM80", "Paper Finalisation", "Finalise paper sections for Overleaf submission"),
    "simulation_42_archive_to_zenodo":
        ("Pre-SIM80", "Zenodo Archive", "Archive simulation outputs to Zenodo"),
    "simulation_43_submit_to_arxiv":
        ("Pre-SIM80", "arXiv Submission", "Submit CMSTG paper I to arXiv"),
    "simulation_44_journal_submissions":
        ("Pre-SIM80", "Journal Submission", "Submit to Phys. Rev. D / JCAP"),
    "simulation_45_cmstg_field_snapshot_exporter":
        ("Pre-SIM80", "Field Snapshot Exporter", "Export field snapshots for visualisation"),
    "simulation_46_cmstg_energy_tracking_diagnostics":
        ("Pre-SIM80", "Energy Tracking", "Track total energy conservation in CMSTG evolution"),
    "simulation_47_cmstg_tensor_eigensystem_solver":
        ("Pre-SIM80", "Tensor Eigensystem", "Solve tensor eigensystem for field decomposition"),
    "simulation_48_cmstg_nonlinear_damping_profile_tuner":
        ("Pre-SIM80", "Nonlinear Damping Tuner", "Tune nonlinear damping profile in memory kernel"),
    "simulation_49_cmstg_mode_transition_mapping":
        ("Pre-SIM80", "Mode Transition", "Map mode transitions in Ψ field evolution"),
    "simulation_50_cmstg_field_memory_decay_checker":
        ("Pre-SIM80", "Memory Decay Checker", "Verify memory field decay behaviour"),
}


def process_old_sims():
    for sim_id, meta in OLD_SIM_META.items():
        sim_dir = os.path.join(BASE, sim_id)
        if not os.path.isdir(sim_dir):
            continue
        out_path = os.path.join(sim_dir, "results.md")
        # Try README first
        readme = read_readme(sim_dir)
        if readme:
            content = make_results_md_from_readme(sim_id, meta, readme)
        else:
            phase, title, desc = meta
            content = f"# {sim_id} — Results\n\n**Phase:** {phase}  \n**Title:** {title}  \n\n## What was tested\n\n{desc}.\n\n## Status\n\nPre-framework simulation. See scripts/ and outputs/ in this directory for full details.\n\n## Context\n\n{CANONICAL_NOTE}\n"
        with open(out_path, "w") as f:
            f.write(content)
        print(f"  WROTE {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main():
    print("=== Generating results.md for all SIM### sims ===")
    phase4_with_result = {"SIM137", "SIM138", "SIM139", "SIM141", "SIM142", "SIM143"}
    unrun = {"SIM140", "SIM144", "SIM145", "SIM146"}

    for sim_id, meta in SIM_META.items():
        sim_dir = os.path.join(BASE, sim_id)
        if not os.path.isdir(sim_dir):
            print(f"  SKIP {sim_id}: not found")
            continue

        if sim_id in phase4_with_result:
            if not copy_existing_result(sim_id):
                content = make_results_md_stub(sim_id, meta)
                write_results_md(sim_id, content)
            continue

        if sim_id in unrun:
            content = make_results_md_stub(sim_id, meta)
            write_results_md(sim_id, content)
            continue

        json_path = find_json(sim_dir)
        if json_path:
            d = load_json(json_path)
            if d:
                content = make_results_md_from_json(sim_id, meta, d)
                write_results_md(sim_id, content)
                continue

        # Fallback: README
        readme = read_readme(sim_dir)
        if readme:
            content = make_results_md_from_readme(sim_id, meta, readme)
            write_results_md(sim_id, content)
        else:
            content = make_results_md_stub(sim_id, meta)
            write_results_md(sim_id, content)

    print("\n=== Processing old simulation_XX series ===")
    process_old_sims()

    print("\nDone.")


if __name__ == "__main__":
    main()
