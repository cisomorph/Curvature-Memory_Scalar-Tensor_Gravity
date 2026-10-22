import sys, numpy as np
sys.path.insert(0, '/home/aion/CMSTG-fixes/Ordered_Simulations/common')
import cmstg_obs as obs, cmstg_bg as bgm
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
om, orr, XI = obs.OMEGA_M, obs.OMEGA_R, 1/6
def make(L0, kind):
    if kind == '132':
        F = lambda p: 0.5*(1-2*L0*p*p)+XI*p; Fp = lambda p: -2*L0*p+XI; Fpp = lambda p: -2*L0
    else:
        F = lambda p: 0.5*np.exp(-2*L0*p*p)+XI*p; Fp = lambda p: -2*L0*p*np.exp(-2*L0*p*p)+XI
        Fpp = lambda p: (-2*L0+8*L0*L0*p*p)*np.exp(-2*L0*p*p)
    return F, Fp, Fpp
def solve(F,Fp,Fpp,Lb,z_ini,method):
    H2 = lambda p,y,N: (om*np.exp(-3*N)+orr*np.exp(-4*N)+Lb)/(2*F(p)+2*Fp(p)*y-y*y/6)
    def rhs(N,s):
        p,y=s; a=np.exp(N); h2=H2(p,y,N)
        e=bgm.eps_closed_form(F(p),Fp(p),Fpp(p),y,h2,3*(om/a**3+4/3*orr/a**4))
        return [y, Fp(p)*6*(2-e)-(3-e)*y]
    Ni=np.log(1/(1+z_ini)); s=solve_ivp(rhs,(Ni,0),[0.0,0.0],method=method,rtol=1e-10,atol=1e-13,dense_output=True)
    N=np.linspace(Ni,0,6000); p,y=s.sol(N)
    return np.exp(-N)-1, 100*np.sqrt([H2(p[i],y[i],N[i]) for i in range(len(N))]), p[-1], y[-1]
def run(F,Fp,Fpp,h,z_ini,method):
    Lb=h*h-om-orr
    for _ in range(15):
        z,H,p0,y0=solve(F,Fp,Fpp,Lb,z_ini,method)
        new=h*h*(2*F(p0)+2*Fp(p0)*y0-y0*y0/6)-om-orr
        if abs(new-Lb)<1e-12: break
        Lb=new
    return z,H,p0
for kind in ('132','134'):
    F,Fp,Fpp = make(0.2, kind)
    for z_ini,method in ((1e5,'DOP853'),(1e6,'LSODA'),(1e7,'Radau')):
        th=lambda h: obs.theta_star_100(obs.H_from_table(*run(F,Fp,Fpp,h,z_ini,method)[:2]))-obs.THETA_STAR_100
        h=brentq(th,0.62,0.72,xtol=1e-7); z,H,p0=run(F,Fp,Fpp,h,z_ini,method); Hf=obs.H_from_table(z,H)
        d=[obs.lcdm_reference(r)['chi2_DESI']-obs.DesiBAO(r).chi2(Hf,obs.r_drag(Hf)) for r in ('DR1','DR2')]
        print(f"SIM{kind} L0=0.2 z_ini={z_ini:.0e} {method:6s}: H0={100*h:.3f} psi0={p0:.4f} F'(psi0)={Fp(p0):+.2e}  dchi2 DR1={d[0]:+.3f} DR2={d[1]:+.3f}")
