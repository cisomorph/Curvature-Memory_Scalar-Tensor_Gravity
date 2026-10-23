"""
Robust maximum-posterior optimiser for noisy, badly conditioned cosmological fits (SIM155).

    x, fx, info = minimize(neg_logpost, x0, lo, hi, prop)

1. coarse : Py-BOBYQA inside the box [lo, hi], each axis scaled to [0, 1]
2. curvature: central-difference Hessian H at the current point, steps h = HSTEP·prop, in proposal units
   Hs = D H D (D = diag(prop)); eigenvalues are folded to |w| and floored at W_FLOOR, so no direction is
   stretched more than 1/sqrt(W_FLOOR) proposal widths (noise or non-convexity cannot blow it up)
3. fine   : Py-BOBYQA in whitened coordinates u, x = x_c + D V diag(w)^(-1/2) u, so −log(post) ≈ ½|u|² + const
Steps 2–3 repeat until a round improves −log(post) by less than `tol` (at most `rounds` times).

Noise: with hi_class at default precision, −log(post) jitters by ≈ 0.035 RMS between points 0.004 proposal
widths apart (discrete k-sampling; SIM155 diagnostics.py). The defaults are set for that: Hessian steps of one
proposal width (≈ 1σ, Δf ≈ 0.5 ≫ noise), stopping radii ≈ 0.1σ, and tol = 0.1. A best-fit −log(post) is then
good to ≈ ±0.05 (±0.1 in χ²).

Outside the box the objective is f(clip(x)) + WALL·|(x − clip(x))/prop|², which is continuous, so the
rotated stage can approach a prior edge (e.g. the GR limit ξ = 0 of CMSTG) without a jump.

A failed evaluation (non-finite, e.g. a hi_class 0/0 stability-test hit at a field turning point) is retried
up to RETRIES times at x + RETRY_STEP·prop·N(0,1); the value at the jittered point stands in for f(x)
(error ≈ |∇f|·1e-3 proposal widths). Only if every retry fails does it return best-so-far + FAIL_STEP (a wall).
Hessian stencils that hit a wall are not used: a failed diagonal falls back to unit curvature in proposal
units and a failed off-diagonal to zero. Without these two rules 5% random failures left fits 2–12 in
−log(post) above the minimum in synthetic tests; with them the result is unaffected (tests/test_fitcore.py).
"""
import numpy as np
import pybobyqa

WALL = 50.0         # −log(post) cost of stepping one proposal width outside the prior box
FAIL_STEP = 50.0    # failed evaluation = best value so far + this
W_FLOOR = 1e-3      # eigenvalue floor in proposal units (max stretch ≈ 32 proposal widths)
HSTEP = 1.0         # Hessian finite-difference step, in proposal widths
RETRIES = 3         # re-tries of a failed evaluation at a jittered point
RETRY_STEP = 1e-3   # jitter, in proposal widths


class Objective:
    """Wraps neg_logpost with the box wall, failure handling and bookkeeping."""

    def __init__(self, neg_logpost, lo, hi, prop, seed=0):
        self.g, self.lo, self.hi, self.prop = neg_logpost, np.asarray(lo, float), np.asarray(hi, float), np.asarray(prop, float)
        self.n = self.fail = self.retries = self.rescued = 0
        self.best_f, self.best_x = np.inf, None
        self.last_failed = False
        self.rng = np.random.default_rng(10_000 + seed)

    def __call__(self, x):
        x = np.asarray(x, float)
        self.n += 1
        xc = np.clip(x, self.lo, self.hi)
        wall = WALL * float(np.sum(((x - xc) / self.prop) ** 2))
        xe, v = xc, self.g(xc)
        k = 0
        while not np.isfinite(v) and k < RETRIES:
            k += 1
            self.retries += 1
            xe = np.clip(xc + RETRY_STEP * self.prop * self.rng.normal(size=len(xc)), self.lo, self.hi)
            v = self.g(xe)
        self.last_failed = not np.isfinite(v)
        if self.last_failed:
            self.fail += 1
            v = (self.best_f if np.isfinite(self.best_f) else 1e4) + FAIL_STEP
        else:
            self.rescued += k > 0
            if v < self.best_f:
                self.best_f, self.best_x = float(v), xe.copy()
        return float(v) + wall


def hessian(f, x, h, f0):
    """Central-difference Hessian and gradient of f at x with steps h (2n² evaluations).
    Stencils containing a failed evaluation (f.last_failed, if f has it) are replaced: diagonal → 1/h_i²·¼
    (unit curvature in units of 2h), off-diagonal → 0, gradient → 0. Returns (H, g, n_bad)."""
    n = len(x)
    E = np.eye(n)
    failed = lambda: bool(getattr(f, 'last_failed', False))
    fp, fm, okp, okm = np.zeros(n), np.zeros(n), np.ones(n, bool), np.ones(n, bool)
    for i in range(n):
        fp[i] = f(x + h[i] * E[i]); okp[i] = not failed()
        fm[i] = f(x - h[i] * E[i]); okm[i] = not failed()
    ok = okp & okm
    g = np.where(ok, (fp - fm) / (2 * h), 0.0)
    H = np.diag(np.where(ok, (fp - 2 * f0 + fm) / h**2, 0.25 / h**2))
    bad = int(np.sum(~ok))
    for i in range(n):
        for j in range(i + 1, n):
            vals, good = [], True
            for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                vals.append(f(x + si * h[i] * E[i] + sj * h[j] * E[j]))
                good &= not failed()
            if good:
                H[i, j] = H[j, i] = (vals[0] - vals[1] - vals[2] + vals[3]) / (4 * h[i] * h[j])
            else:
                bad += 1
    return H, g, bad


def minimize(neg_logpost, x0, lo, hi, prop, save=None, maxfun1=400, rhoend1=5e-3, maxfun3=300, rhoend3=0.1,
             rounds=3, tol=0.1, seed=0):
    """Returns (x_best, f_best, info). `save(stage, x, f, info)` is called after every stage."""
    f = Objective(neg_logpost, lo, hi, prop, seed=seed)
    lo, hi, prop = f.lo, f.hi, f.prop
    n = len(lo)
    D = np.diag(prop)
    info = dict(rounds=[])

    s1 = pybobyqa.solve(f, np.clip(np.asarray(x0, float), lo, hi), bounds=(lo, hi), scaling_within_bounds=True,
                        rhobeg=0.05, rhoend=rhoend1, maxfun=maxfun1, do_logging=False)
    x, fx = np.clip(s1.x, lo, hi), float(s1.f)
    info['coarse'] = dict(flag=int(s1.flag), nf=int(s1.nf), f=fx)
    if f.best_x is None:
        raise RuntimeError('no successful evaluation in the coarse stage: check the model set-up')
    if save:
        save('coarse', x, fx, info)

    cov = None
    for r in range(rounds):
        H, g, nbad = hessian(f, x, HSTEP * prop, fx)
        w, V = np.linalg.eigh(D @ H @ D)
        w_reg = np.clip(np.abs(w), W_FLOOR, None)
        gs = prop * g
        pred = 0.5 * float(gs @ V @ ((V.T @ gs) / w_reg))
        T = D @ V @ np.diag(1.0 / np.sqrt(w_reg))
        cov = T @ T.T                                   # ≈ inverse Hessian (Laplace covariance) in raw units
        s3 = pybobyqa.solve(lambda u: f(x + T @ u), np.zeros(n), bounds=(-8.0 * np.ones(n), 8.0 * np.ones(n)),
                            rhobeg=0.5, rhoend=rhoend3, maxfun=maxfun3, do_logging=False)
        xr = x + T @ s3.x
        xn = np.clip(xr, lo, hi)
        fn = float(s3.f) if np.array_equal(xn, xr) else f(xn)
        gain = fx - fn
        if gain > 0:
            x, fx = xn, fn
        info['rounds'].append(dict(round=r + 1, eig_scaled=[float(v) for v in w], n_floored=int(np.sum(np.abs(w) < W_FLOOR)),
                                   bad_stencils=nbad,
                                   n_negative=int(np.sum(w < 0)), predicted_gain=pred, gain=float(gain),
                                   flag=int(s3.flag), nf=int(s3.nf), u_norm=float(np.linalg.norm(s3.x))))
        if save:
            save(f'round{r + 1}', x, fx, info)
        if gain < tol:
            break

    if f.best_x is not None and f.best_f < fx:         # an evaluation (e.g. in a Hessian stencil) beat the optimiser
        info['stencil_best_used'] = float(fx - f.best_f)
        x, fx = f.best_x, f.best_f
    info.update(evaluations=f.n, failed_evaluations=f.fail, retries=f.retries, rescued_by_retry=f.rescued,
                sigma=[float(s) for s in np.sqrt(np.diag(cov))] if cov is not None else None)
    return x, fx, info
