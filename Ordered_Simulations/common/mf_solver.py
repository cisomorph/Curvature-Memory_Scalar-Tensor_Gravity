"""
mf_solver.py — corrected multi-field scalar-tensor background + growth (audit 2026-10-22).

Action (signature -+++, 8πG = 1, GR is F = 1/2), fields φ_i:
    S = ∫√-g [ F(φ, N) R − Σ_i ½K_i(φ_i)(∂φ_i)² − V(φ) ] + S_m + bare Λ  (+ V_int(φ, N))
V_int is an optional phenomenological field–matter term (energy density V_int, force −∂V_int/∂φ).

Friedmann (lapse variation; ρ in GR units 3H² = ρ; written in ω units, ρ = 3ω):
    H² [2F + 2(Σ F_i φ_i' + F_N) − Σ K_i φ_i'²/6] = rho_scale·(ω_m a⁻³ + ω_r a⁻⁴) + Λ_b + (V + V_int)/3
Field equations (' = d/dN):
    φ_i'' = [F_i·6(2−ε) − (V_i + V_int,i)/H² − ½K_i' φ_i'²]/K_i − (3−ε)φ_i'   ≡ A_i + B_i ε
ε from the Friedmann constraint along the trajectory (generalises cmstg_bg.eps_from_constraint).

Growth (sub-horizon, quasi-static, massless scalars):
    δ'' + (2 − ε)δ' = (3/2)·[rho_scale ω_m a⁻³ / H²]·g_eff·δ,
    g_eff = (1/(2F))·(2F + 4S)/(2F + 3S),  S = Σ F_i²/K_i     (GR: F = ½ → 1)
    plus an optional extra fifth-force factor for direct matter couplings.
"""

import numpy as np
from scipy.integrate import solve_ivp

import cmstg_obs as obs

H100 = 100.0


class _Breakdown(Exception):
    pass


class MultiFieldBG:
    def __init__(self, n, F, F_i, F_N=None, K=None, K_i=None, V=None, V_i=None,
                 V_int=None, V_int_i=None, rho_scale=1.0, omega_m=obs.OMEGA_M, omega_r=obs.OMEGA_R,
                 z_ini=1e5, n_points=4000, max_step=0.05, coef_floor=1e-3, fifth_force=None):
        self.n = n
        self.F, self.F_i = F, F_i                         # F(f, N), F_i(f, N) -> array(n)
        self.F_N = F_N or (lambda f, N: 0.0)
        self.K = K or (lambda f: np.ones(n))               # K_i(f_i) as array(n)
        self.K_i = K_i or (lambda f: np.zeros(n))          # dK_i/df_i
        self.V = V or (lambda f: 0.0)
        self.V_i = V_i or (lambda f: np.zeros(n))
        self.V_int = V_int or (lambda f, N: 0.0)
        self.V_int_i = V_int_i or (lambda f, N: np.zeros(n))
        self.rho_scale = rho_scale
        self.om, self.orr = omega_m, omega_r
        self.z_ini, self.n_points, self.max_step = z_ini, n_points, max_step
        self.coef_floor = coef_floor
        self.fifth_force = fifth_force                     # optional g_extra(f, N)

    def coef(self, f, y, N):
        return (2.0 * self.F(f, N) + 2.0 * (np.dot(self.F_i(f, N), y) + self.F_N(f, N))
                - np.dot(self.K(f), y * y) / 6.0)

    def H2(self, f, y, N, Lb):
        c = self.coef(f, y, N)
        if not (c > self.coef_floor):
            return np.nan
        a = np.exp(N)
        r = self.rho_scale * (self.om / a**3 + self.orr / a**4) + Lb + (self.V(f) + self.V_int(f, N)) / 3.0
        X = r / c
        return X if X > 0 else np.nan

    def AB(self, f, y, N, H2):
        Fi, K, Ki = self.F_i(f, N), self.K(f), self.K_i(f)
        Vi = self.V_i(f) + self.V_int_i(f, N)
        A = (12.0 * Fi - Vi / H2 - 0.5 * Ki * y * y) / K - 3.0 * y
        B = y - 6.0 * Fi / K
        return A, B

    def eps(self, f, y, N, Lb):
        H2 = self.H2(f, y, N, Lb)
        if not np.isfinite(H2):
            return np.nan, np.nan
        A, B = self.AB(f, y, N, H2)
        ln = lambda ff, yy, nn: np.log(self.H2(ff, yy, nn, Lb))
        dN = 1e-5
        aN = (ln(f, y, N + dN) - ln(f, y, N - dN)) / (2 * dN)
        num, den = aN, 1.0
        for i in range(self.n):
            e = np.zeros(self.n); e[i] = 1.0
            hf = 1e-6 * max(1.0, abs(f[i])); hy = 1e-6 * max(1.0, abs(y[i]))
            af = (ln(f + hf * e, y, N) - ln(f - hf * e, y, N)) / (2 * hf)
            ay = (ln(f, y + hy * e, N) - ln(f, y - hy * e, N)) / (2 * hy)
            num += af * y[i] + ay * A[i]
            den += 0.5 * ay * B[i]
        return -0.5 * num / den, H2

    def rhs(self, N, s, Lb):
        f, y = s[:self.n], s[self.n:]
        e, H2 = self.eps(f, y, N, Lb)
        if not np.isfinite(e):
            raise _Breakdown(N)
        A, B = self.AB(f, y, N, H2)
        return np.concatenate([y, A + B * e])

    def integrate(self, f_ini, y_ini, Lb, f_max=1e3):
        Ni = np.log(1.0 / (1.0 + self.z_ini))
        s0 = np.concatenate([np.asarray(f_ini, float), np.asarray(y_ini, float)])
        bad = lambda N, s, Lb_: self.coef(s[:self.n], s[self.n:], N) - 2 * self.coef_floor
        bad.terminal = True
        blow = lambda N, s, Lb_: f_max - np.max(np.abs(s[:self.n]))
        blow.terminal = True
        try:
            sol = solve_ivp(self.rhs, (Ni, 0.0), s0, args=(Lb,), t_eval=np.linspace(Ni, 0.0, self.n_points),
                            method='RK45', rtol=1e-8, atol=1e-11, max_step=self.max_step, events=[bad, blow])
        except _Breakdown:
            return None
        if not sol.success or sol.t[-1] < -1e-9:
            return None
        N = sol.t
        f, y = sol.y[:self.n], sol.y[self.n:]
        H2 = np.array([self.H2(f[:, k], y[:, k], N[k], Lb) for k in range(len(N))])
        if not np.all(np.isfinite(H2)):
            return None
        z = np.exp(-N) - 1.0
        o = np.argsort(z)
        Fv = np.array([self.F(f[:, k], N[k]) for k in range(len(N))])
        return dict(z=z[o], N=N[o], f=f[:, o], y=y[:, o], H=H100 * np.sqrt(H2[o]), F=Fv[o],
                    f0=f[:, -1].copy(), y0=y[:, -1].copy(), F0=float(Fv[-1]),
                    H0=float(H100 * np.sqrt(H2[-1])), Lb=float(Lb))

    def calibrated(self, f_ini, y_ini, h, n_iter=12, tol=1e-9):
        f_ini = np.asarray(f_ini, float); y_ini = np.asarray(y_ini, float)
        Lb = h * h * 2.0 * self.F(f_ini, 0.0) - self.rho_scale * (self.om + self.orr)
        bg = None
        for _ in range(n_iter):
            bg = self.integrate(f_ini, y_ini, Lb)
            if bg is None:
                return None
            f0, y0 = bg['f0'], bg['y0']
            Lb_new = (self.coef(f0, y0, 0.0) * h * h - self.rho_scale * (self.om + self.orr)
                      - (self.V(f0) + self.V_int(f0, 0.0)) / 3.0)
            if abs(Lb_new - Lb) < tol:
                Lb = Lb_new
                break
            Lb = Lb_new
        bg = self.integrate(f_ini, y_ini, Lb)
        if bg is not None:
            bg['h_target'] = h
        return bg

    def H_func(self, bg):
        return obs.H_from_table(bg['z'], bg['H'], self.rho_scale * self.om, self.rho_scale * self.orr)

    def theta_matched(self, f_ini, y_ini, h_guess=0.674):
        import cmstg_bg as bgm
        cache = {}
        def g(h):
            if not (0.40 < h < 1.0):
                return np.nan
            bg = self.calibrated(f_ini, y_ini, h)
            if bg is None:
                return np.nan
            cache[h] = bg
            return obs.theta_star_100(self.H_func(bg)) - obs.THETA_STAR_100
        h = bgm.secant_match(g, h_guess, h_guess * 1.01, tol=1e-6)
        if not np.isfinite(h):
            return None
        return cache.get(h) or self.calibrated(f_ini, y_ini, h)

    # ── growth ────────────────────────────────────────────────────────────────
    def g_eff(self, f, N):
        F = self.F(f, N)
        S = float(np.sum(self.F_i(f, N)**2 / self.K(f)))
        g = (1.0 / (2.0 * F)) * (2.0 * F + 4.0 * S) / (2.0 * F + 3.0 * S)
        if self.fifth_force is not None:
            g *= self.fifth_force(f, N)
        return g

    def fsigma8(self, bg, z_out, sigma8_ref=0.811):
        """fσ8(z) with σ8 normalised to sigma8_ref × D(0)/D_ΛCDM(0) at equal early amplitude
        (same convention as the published sims, now with G_eff(z) in the growth source)."""
        Hf = self.H_func(bg)
        Nz = np.log(1.0 / (1.0 + bg['z'][::-1]))          # ascending N
        fz = bg['f'][:, ::-1]
        def src(N, use_geff):
            a = np.exp(N); z = 1.0 / a - 1.0
            H2 = (Hf(np.array([z]))[0] / H100)**2
            g = 1.0
            if use_geff:
                k = np.searchsorted(Nz, N).clip(1, len(Nz) - 1)
                fv = fz[:, k - 1] + (fz[:, k] - fz[:, k - 1]) * (N - Nz[k - 1]) / (Nz[k] - Nz[k - 1])
                g = self.g_eff(fv, N)
            return 1.5 * self.rho_scale * self.om / a**3 / H2 * g
        def ode(N, s, use_geff):
            a = np.exp(N); z = 1.0 / a - 1.0; d = 1e-4
            lnH = lambda zz: np.log(Hf(np.array([zz]))[0])
            eps = -(lnH(1.0 / np.exp(N + d) - 1.0) - lnH(1.0 / np.exp(N - d) - 1.0)) / (2 * d)
            return [s[1], src(N, use_geff) * s[0] - (2.0 - eps) * s[1]]
        Ns = np.linspace(np.log(1e-3), 0.0, 600)
        out = solve_ivp(ode, (Ns[0], 0.0), [1e-3, 1e-3], args=(True,), t_eval=Ns, rtol=1e-8, atol=1e-11)
        D, Dp = out.y
        # ΛCDM reference with the same ω_m, H0 from θ* (pipeline reference)
        ref = obs.lcdm_reference('DR1')
        Hl = obs.H_lcdm(ref['h'])
        def ode_l(N, s):
            a = np.exp(N); z = 1.0 / a - 1.0; d = 1e-4
            lnH = lambda zz: np.log(Hl(np.array([zz]))[0])
            eps = -(lnH(1.0 / np.exp(N + d) - 1.0) - lnH(1.0 / np.exp(N - d) - 1.0)) / (2 * d)
            H2 = (Hl(np.array([z]))[0] / H100)**2
            return [s[1], 1.5 * obs.OMEGA_M / a**3 / H2 * s[0] - (2.0 - eps) * s[1]]
        ol = solve_ivp(ode_l, (Ns[0], 0.0), [1e-3, 1e-3], t_eval=Ns, rtol=1e-8, atol=1e-11)
        s8 = sigma8_ref * D[-1] / ol.y[0][-1]
        zN = 1.0 / np.exp(Ns) - 1.0
        f_growth = Dp / D
        fs8 = np.interp(z_out, zN[::-1], (f_growth * s8 * D / D[-1])[::-1])
        return fs8, s8
