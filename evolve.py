#!/usr/bin/env python
"""1+1D time-domain evolution of the massive-scalar QBS on Kerr (Dolan 2012),
single-mode l=m=1, on the BL tortoise grid, with an ingoing-horizon PML.

The horizon absorber is Dolan's near-horizon PML (his Sec. on boundary
conditions): in the deep tortoise region (r ~ r_+, where V->0 and the full
equation reduces to  psi_tt = psi_xx + 2im Om(psi_x - psi_t) ) we evolve

  psi_t = pi
  pi_t  = psi_xx + 2im Om (psi_x - pi) - chi_x - im Om chi - g(pi + im Om psi)
  chi_t = g(psi_x - chi + im Om psi) - im Om chi

with attenuation g(x) ramping up toward x_min. With g=0 this is IDENTICAL to
the full equation near the horizon, so the match to the interior (full,
r-dependent-coefficient) equation at x_match (r = r_+ + eps) is smooth.

Extraction: probe time series -> Prony/matrix-pencil (handles tiny Gamma from a
moderate window) or a log-amplitude fit. Benchmarks against leaver.solve.
"""
import numpy as np, sys, argparse, os
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
from radial1d import build_grid, coeffs, d_dx, d2_dx2, horizons

def run(alpha=0.50, a=0.90, xmin=-70., xmax=220., N=4000, T=3000., cfl=0.4,
        x0=20., sigx=6., eps_match=0.02, gmax=1.0, probe_r=None, nrec=6000,
        nsnap=0):
    w_ref = complex(leaver.solve(alpha, a))
    rp,_ = horizons(a); OmH = a/(2*rp); m = 1
    x, r = build_grid(a, xmin, xmax, N); dx = x[1]-x[0]
    D,A2,A1,At,Vc,Delta,P = coeffs(r, a, alpha); invD = 1.0/D
    # PML region: x below the match radius r_+ + eps_match
    from radial1d import tortoise_x
    x_match = tortoise_x(rp+eps_match, a)
    pml = x < x_match
    # attenuation ramps 0 -> gmax from x_match down to xmin (quadratic)
    g = np.zeros_like(x)
    g[pml] = gmax*((x_match - x[pml])/(x_match - xmin))**2
    # mild outer sponge for numerical safety (bound state ~ 0 out there anyway)
    xo = xmax-40; out = x>xo; gout = np.zeros_like(x); gout[out]=0.8*((x[out]-xo)/40)**2
    imO = 1j*m*OmH
    if probe_r is None: probe_r = min(2.0/alpha**2, 0.6*r[-1])
    ip = int(np.argmin(np.abs(r - probe_r)))
    dt = cfl*dx; nsteps = int(T/dt)

    def rhs(psi, pi, chi):
        px = d_dx(psi,dx); pxx = d2_dx2(psi,dx); cx = d_dx(chi,dx)
        # interior (full, r-dependent) momentum RHS
        dpi_int = invD*(A2*pxx + A1*px + At*pi - Vc*psi)
        # PML momentum RHS (near-horizon reduced form + auxiliary chi + attenuation)
        dpi_pml = pxx + 2*imO*(px - pi) - cx - imO*chi - g*(pi + imO*psi)
        dpi = np.where(pml, dpi_pml, dpi_int) - gout*pi
        dpsi = pi - gout*psi
        dchi = np.where(pml, g*(px - chi + imO*psi) - imO*chi, 0.0)
        return dpsi, dpi, dchi

    psi = np.exp(-((x-x0)**2)/(2*sigx**2)).astype(complex)
    pi  = np.zeros_like(psi); chi = np.zeros_like(psi)
    trec = np.linspace(0, T, nrec); rec_i=(trec/dt).astype(int); rp_=0
    tsig=[]; psig=[]; energy=[]
    tsnap=np.linspace(0, T, nsnap) if nsnap>0 else np.array([])
    snap_i=(tsnap/dt).astype(int); sp_=0
    field_t=[]; field_psi=[]
    for n in range(nsteps):
        if rp_<nrec and n>=rec_i[rp_]:
            tsig.append(n*dt); psig.append(psi[ip]); rp_+=1
        if sp_<nsnap and n>=snap_i[sp_]:
            field_t.append(n*dt); field_psi.append(psi.astype(np.complex64).copy()); sp_+=1
        k1=rhs(psi,pi,chi)
        k2=rhs(psi+0.5*dt*k1[0], pi+0.5*dt*k1[1], chi+0.5*dt*k1[2])
        k3=rhs(psi+0.5*dt*k2[0], pi+0.5*dt*k2[1], chi+0.5*dt*k2[2])
        k4=rhs(psi+dt*k3[0], pi+dt*k3[1], chi+dt*k3[2])
        psi=psi+(dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        pi =pi +(dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        chi=chi+(dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        if not np.all(np.isfinite(psi)): print(f"  NaN at t={n*dt:.1f}"); break
    out = dict(w_ref=w_ref, OmH=OmH, t=np.array(tsig), psi=np.array(psig),
               probe_r=r[ip], dx=dx, dt=dt, x_match=x_match)
    if nsnap>0:
        out.update(field_t=np.array(field_t), field_psi=np.array(field_psi),
                    x=x, r=r, rp=rp)
    return out

def prony(t, psi, p=6, tmin_frac=0.35):
    """Matrix-pencil: fit psi(t) ~ sum_k c_k exp(s_k t); return the dominant
    (largest |c|) complex frequency omega = i*s = E + i*Gamma (e^{-i w t})."""
    m = t > tmin_frac*t[-1]; tt=t[m]-t[m][0]; y=psi[m]
    dt=np.median(np.diff(tt)); Nn=len(y); L=min(p+ p, Nn//2)
    # Hankel matrices for matrix pencil
    Y0=np.array([y[i:i+L] for i in range(Nn-L)])
    Y1=np.array([y[i+1:i+1+L] for i in range(Nn-L)])
    U,S,Vh=np.linalg.svd(Y0,full_matrices=False)
    r_=min(p,np.sum(S>1e-10*S[0]))
    U,S,V=U[:,:r_],S[:r_],Vh.conj().T[:,:r_]
    A=U.conj().T@Y1@V@np.diag(1/S); lam,_=np.linalg.eig(A)
    s=np.log(lam)/dt                    # psi ~ e^{s t}
    om=1j*s                             # e^{-i w t} => w = i s ; E=Re, Gamma=Im
    # amplitudes by least squares
    Van=np.exp(np.outer(tt,s)); c,_,_,_=np.linalg.lstsq(Van,y,rcond=None)
    k=np.argmax(np.abs(c)); return om[k].real, om[k].imag

def logfit(t,psi,tmin_frac=0.4):
    m=t>tmin_frac*t[-1]; tt=t[m]
    G=np.polyfit(tt,np.log(np.abs(psi[m])),1)[0]
    E=-np.polyfit(tt,np.unwrap(np.angle(psi[m])),1)[0]; return E,G

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--alpha",type=float,default=0.50); ap.add_argument("--a",type=float,default=0.90)
    ap.add_argument("--T",type=float,default=3000.); ap.add_argument("--N",type=int,default=4000)
    ap.add_argument("--xmin",type=float,default=-70.); ap.add_argument("--xmax",type=float,default=220.)
    ap.add_argument("--gmax",type=float,default=1.0); ap.add_argument("--save",type=str,default=None)
    ap.add_argument("--nsnap",type=int,default=0,
                    help="save this many full-field psi(x) snapshots (for a space-time heatmap)")
    q=ap.parse_args()
    o=run(q.alpha,q.a,q.xmin,q.xmax,q.N,q.T,gmax=q.gmax,nsnap=q.nsnap)
    w=o["w_ref"]; E1,G1=logfit(o["t"],o["psi"])
    try: E2,G2=prony(o["t"],o["psi"])
    except Exception as e: E2,G2=(float("nan"),float("nan")); print("prony failed:",e)
    print(f"alpha={q.alpha} a={q.a} N={q.N} T={q.T} gmax={q.gmax} (probe r={o['probe_r']:.1f}, dx={o['dx']:.3f}, x_match={o['x_match']:.1f})")
    print(f"  Leaver : E={w.real:.6f}  Gamma={w.imag:+.4e}")
    print(f"  logfit : E={E1:.6f}  Gamma={G1:+.4e}   (Gamma {G1/w.imag:+.3f}x)")
    print(f"  prony  : E={E2:.6f}  Gamma={G2:+.4e}   (Gamma {G2/w.imag:+.3f}x)")
    if q.save is not None:
        save_kw = dict(t=o["t"], psi=o["psi"], w_ref=w)
        if q.nsnap>0:
            save_kw.update(field_t=o["field_t"], field_psi=o["field_psi"],
                            x=o["x"], r=o["r"], rp=o["rp"], x_match=o["x_match"])
        np.savez(q.save, **save_kw); print("saved",q.save)

if __name__=="__main__": main()
