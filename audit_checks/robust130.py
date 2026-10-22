import sys, numpy as np
sys.path.insert(0, '/home/aion/CMSTG-fixes/Ordered_Simulations/common')
import cmstg_obs as obs, cmstg_bg as bgm
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
L0, om, orr = 0.003, obs.OMEGA_M, obs.OMEGA_R
F = lambda p: 0.5+L0*p*p; Fp = lambda p: 2*L0*p; Fpp = 2*L0
def solve(psi_ini, Lb, z_ini, method):
    def H2(p,y,N):
        a=np.exp(N); return (om/a**3+orr/a**4+Lb)/(2*F(p)+2*Fp(p)*y-y*y/6)
    def rhs(N,s):
        p,y=s; a=np.exp(N); h2=H2(p,y,N)
        eps=bgm.eps_closed_form(F(p),Fp(p),Fpp,y,h2,3*(om/a**3+4/3*orr/a**4))   # Raychaudhuri closed form
        return [y, Fp(p)*6*(2-eps)-(3-eps)*y]
    Ni=np.log(1/(1+z_ini))
    s=solve_ivp(rhs,(Ni,0),[psi_ini,0.0],method=method,rtol=1e-10,atol=1e-13,dense_output=True)
    N=np.linspace(Ni,0,6000); p,y=s.sol(N)
    H=100*np.sqrt([H2(p[i],y[i],N[i]) for i in range(len(N))])
    return np.exp(-N)-1, H, p[-1], y[-1]
def run(psi_ini, h, z_ini, method):
    Lb=2*F(psi_ini)*h*h-om-orr
    for _ in range(12):
        z,H,p0,y0=solve(psi_ini,Lb,z_ini,method)
        new=h*h*(2*F(p0)+2*Fp(p0)*y0-y0*y0/6)-om-orr
        if abs(new-Lb)<1e-12: break
        Lb=new
    return z,H,p0,y0
for z_ini,method in ((1e5,'DOP853'),(1e6,'LSODA'),(1e7,'Radau')):
    def th(h):
        z,H,_,_=run(2.62,h,z_ini,method); return obs.theta_star_100(obs.H_from_table(z,H))-obs.THETA_STAR_100
    h=brentq(th,0.64,0.69,xtol=1e-7)
    z,H,p0,y0=run(2.62,h,z_ini,method); Hf=obs.H_from_table(z,H)
    for rel in ('DR1','DR2'):
        lik=obs.DesiBAO(rel); c2=lik.chi2(Hf,obs.r_drag(Hf)); ref=obs.lcdm_reference(rel)['chi2_DESI']
        print(f"z_ini={z_ini:.0e} {method:6s} {rel}: H0={100*h:.3f} psi0={p0:.4f} y0={y0:.4f} chi2={c2:.3f} dchi2 vs LCDM={ref-c2:+.3f}")
