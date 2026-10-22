"""
st_solver.py — corrected single-field scalar-tensor background solver (audit 2026-10-22).

Action (signature -+++, 8πG = 1, GR is F = 1/2):
    S = ∫√-g [ F(Ψ, N) R − ½(∂Ψ)² − V(Ψ) + f(Ψ) 𝒢 ] + S_m + bare Λ

Friedmann (lapse variation; ρ in GR units 3H² = ρ; densities written as ω, ρ = 3ω):
    H²(6F + 6(F_Ψ Ψ' + F_N) − ½Ψ'²) + 24 f_Ψ Ψ' H⁴ = ρ + V
  → H²(2F + 2(F_Ψ Ψ' + F_N) − Ψ'²/6) + 8 f_Ψ Ψ' H⁴ = rho_scale·(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b + V/3
Scalar EOM:
    Ψ'' + (3 − ε)Ψ' = [F_Ψ R + f_Ψ 𝒢 − V_Ψ]/H²,  R/H² = 6(2 − ε),  𝒢/H² = 24H²(1 − ε)
ε is taken from the constraint along the trajectory (cmstg_bg.eps_from_constraint). For
covariant F(Ψ) this equals the Raychaudhuri closed form (validated to 2e-10).
With F depending explicitly on N the model is phenomenological (not from a covariant
action); the constraint-consistent ε is then the natural closure.

rho_scale: 1 for GR-normalised densities (Planck ω enters as in ΛCDM when F = ½);
           2F_ref for the Phase-4 "Planck-calibrated" convention (early-time H equal to
           ΛCDM when F = F_ref).
"""

import numpy as np
from scipy.integrate import solve_ivp

import cmstg_obs as obs
import cmstg_bg as bgm

H100 = 100.0


class _Breakdown(Exception):
    """Raised when the Friedmann constraint has no physical (H² > 0) solution."""


class ScalarTensorBG:
    def __init__(self, F, F_psi, F_N=None, V=None, V_psi=None, f_gb_psi=None,
                 rho_scale=1.0, omega_m=obs.OMEGA_M, omega_r=obs.OMEGA_R,
                 z_ini=1e5, n_points=4000, max_step=0.05, coef_floor=1e-3):
        self.F, self.F_psi = F, F_psi
        self.F_N = F_N or (lambda p, N: 0.0)
        self.V = V or (lambda p: 0.0)
        self.V_psi = V_psi or (lambda p: 0.0)
        self.f_gb_psi = f_gb_psi          # None or callable f'(Ψ) [H100⁻² units]
        self.rho_scale = rho_scale
        self.om, self.orr = omega_m, omega_r
        self.z_ini, self.n_points, self.max_step = z_ini, n_points, max_step
        self.coef_floor = coef_floor

    # ── constraint ────────────────────────────────────────────────────────────
    def coef(self, p, y, N):
        return 2.0 * self.F(p, N) + 2.0 * (self.F_psi(p, N) * y + self.F_N(p, N)) - y * y / 6.0

    def rho_hat(self, N, Lb, p):
        a = np.exp(N)
        return self.rho_scale * (self.om / a**3 + self.orr / a**4) + Lb + self.V(p) / 3.0

    def H2(self, p, y, N, Lb):
        c = self.coef(p, y, N)
        r = self.rho_hat(N, Lb, p)
        if self.f_gb_psi is None:
            return r / c if c > self.coef_floor else np.nan
        b = 8.0 * self.f_gb_psi(p) * y            # b X² + c X − r = 0, X = H²
        if abs(b) < 1e-300:
            return r / c if c > self.coef_floor else np.nan
        disc = c * c + 4.0 * b * r
        if disc < 0:
            return np.nan
        X = 2.0 * r / (c + np.sqrt(disc))          # root continuous with r/c as b → 0
        return X if X > 0 else np.nan

    # ── dynamics ──────────────────────────────────────────────────────────────
    def AB(self, p, y, N, H2):
        """Scalar EOM written as Ψ'' = A + B ε."""
        Fp = self.F_psi(p, N)
        A = 12.0 * Fp - 3.0 * y - self.V_psi(p) / H2           # V, V' in ρ-units (3H² = ρ in GR)
        B = y - 6.0 * Fp
        if self.f_gb_psi is not None:
            g = 24.0 * self.f_gb_psi(p) * H2
            A += g
            B -= g
        return A, B

    def eps(self, p, y, N, Lb):
        H2 = self.H2(p, y, N, Lb)
        if not np.isfinite(H2):
            return np.nan, np.nan
        A, B = self.AB(p, y, N, H2)
        e = bgm.eps_from_constraint(lambda pp, yy, nn: self.H2(pp, yy, nn, Lb), p, y, N, A, B)
        return e, H2

    def rhs(self, N, s, Lb):
        p, y = s
        e, H2 = self.eps(p, y, N, Lb)
        if not np.isfinite(e):
            raise _Breakdown(N)              # constraint has no physical solution
        A, B = self.AB(p, y, N, H2)
        return [y, A + B * e]

    def integrate(self, psi_ini, y_ini, Lb, psi_max=1e3):
        Ni = np.log(1.0 / (1.0 + self.z_ini))
        bad = lambda N, s, Lb_: (self.coef(s[0], s[1], N) - 2 * self.coef_floor)
        bad.terminal = True
        blow = lambda N, s, Lb_: psi_max - abs(s[0])
        blow.terminal = True
        try:
            sol = solve_ivp(self.rhs, (Ni, 0.0), [psi_ini, y_ini], args=(Lb,),
                            t_eval=np.linspace(Ni, 0.0, self.n_points), method='RK45',
                            rtol=1e-8, atol=1e-11, max_step=self.max_step, events=[bad, blow])
        except _Breakdown:
            return None
        if not sol.success or sol.t[-1] < -1e-9:
            return None
        N, p, y = sol.t, sol.y[0], sol.y[1]
        H2 = np.array([self.H2(p[i], y[i], N[i], Lb) for i in range(len(N))])
        if not np.all(np.isfinite(H2)) or np.any(H2 <= 0):
            return None
        F = np.array([self.F(p[i], N[i]) for i in range(len(N))])
        z = np.exp(-N) - 1.0
        o = np.argsort(z)
        return dict(z=z[o], N=N[o], psi=p[o], y=y[o], H=H100 * np.sqrt(H2[o]), F=F[o],
                    psi0=float(p[-1]), y0=float(y[-1]), F0=float(F[-1]), H0=float(H100 * np.sqrt(H2[-1])),
                    Lb=float(Lb))

    # ── H0 calibration (as the published sims) ────────────────────────────────
    def calibrated(self, psi_ini, y_ini, h, n_iter=12, tol=1e-9):
        """Iterate Λ_b so that H(z=0) = 100 h exactly."""
        Lb = h * h * 2.0 * self.F(psi_ini, 0.0) - self.rho_scale * (self.om + self.orr)
        bg = None
        for _ in range(n_iter):
            bg = self.integrate(psi_ini, y_ini, Lb)
            if bg is None:
                return None
            p0, y0 = bg['psi0'], bg['y0']
            # solve H2(p0, y0, 0, Lb_new) = h² for Lb_new (linear in Lb when f_GB = 0)
            X = h * h
            c = self.coef(p0, y0, 0.0)
            gb = 8.0 * self.f_gb_psi(p0) * y0 * X * X if self.f_gb_psi is not None else 0.0
            Lb_new = c * X + gb - self.rho_scale * (self.om + self.orr) - self.V(p0) / 3.0
            if abs(Lb_new - Lb) < tol:
                Lb = Lb_new
                break
            Lb = Lb_new
        bg = self.integrate(psi_ini, y_ini, Lb)
        if bg is not None:
            bg['h_target'] = h
        return bg

    # ── observables ───────────────────────────────────────────────────────────
    def H_func(self, bg):
        return obs.H_from_table(bg['z'], bg['H'], self.rho_scale * self.om, self.rho_scale * self.orr)

    def observe(self, bg, n_params=0):
        Hf = self.H_func(bg)
        rd = obs.r_drag(Hf)
        th = obs.theta_star_100(Hf)
        out = dict(theta_star_100=th, chi2_theta=obs.chi2_theta(th), r_drag=rd,
                   H0=bg['H0'], psi0=bg['psi0'], F0=bg['F0'], Lb=bg['Lb'])
        for rel in ('DR1', 'DR2'):
            lik = _lik(rel)
            c2 = lik.chi2(Hf, rd)
            ref = obs.lcdm_reference(rel)
            s = obs.summarize_fit(c2, lik.n, n_params, out['chi2_theta'])
            d = (ref['chi2_DESI'] + ref['chi2_theta']) - (c2 + out['chi2_theta'])
            s.update(delta_chi2_tot_vs_LCDM=d,
                     improvement_sigma=obs.delta_chi2_to_sigma(d, max(n_params, 1)))
            out[rel] = s
        return out

    def theta_matched(self, psi_ini, y_ini, h_guess=0.674):
        """Profile over H0: solve h so that 100θ* = Planck (fair comparison)."""
        cache = {}
        def f(h):
            if not (0.40 < h < 1.0):
                return np.nan
            bg = self.calibrated(psi_ini, y_ini, h)
            if bg is None:
                return np.nan
            Hf = self.H_func(bg)
            cache[h] = bg
            return obs.theta_star_100(Hf) - obs.THETA_STAR_100
        h = bgm.secant_match(f, h_guess, h_guess * 1.01, tol=1e-6)
        if not np.isfinite(h):
            return None
        return cache.get(h) or self.calibrated(psi_ini, y_ini, h)


_LIKS = {}
def _lik(rel):
    if rel not in _LIKS:
        _LIKS[rel] = obs.DesiBAO(rel)
    return _LIKS[rel]


def verdict_from(o_faithful, o_matched, n_params, release='DR1'):
    """
    Uniform verdict used by the audited Phase 3/4 sims (same thresholds as published):
      PASS    : GOF < 2σ and |Δθ*| < 2σ
      PARTIAL : beats ΛCDM (Δχ²_tot > 0) in the same pipeline
      FAIL    : otherwise
    Evaluated in both modes: 'faithful' (H0 fixed as published) and 'θ*-matched' (H0 solved).
    """
    def one(o):
        if o is None:
            return 'FAIL (no solution)'
        s = o[release]
        th_ok = o['chi2_theta'] < 4.0
        if s['gof_sigma'] < 2.0 and th_ok:
            return 'PASS'
        if s['delta_chi2_tot_vs_LCDM'] > 0:
            return 'PARTIAL'
        return 'FAIL'
    return dict(faithful=one(o_faithful), theta_matched=one(o_matched))
