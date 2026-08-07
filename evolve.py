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
    # flux probe at the PML interface (Dolan's x_+, eq. 34/36): the boundary
    # of the "exterior" region his flux/energy diagnostics are defined on.
    ifl = int(np.argmin(np.abs(x - x_match))); rfl = r[ifl]
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
    tsig=[]; psig=[]; phi_fl=[]; phidot_fl=[]; phix_fl=[]
    tsnap=np.linspace(0, T, nsnap) if nsnap>0 else np.array([])
    snap_i=(tsnap/dt).astype(int); sp_=0
    field_t=[]; field_psi=[]
    for n in range(nsteps):
        if rp_<nrec and n>=rec_i[rp_]:
            tsig.append(n*dt); psig.append(psi[ip])
            phi_fl.append(psi[ifl]); phidot_fl.append(pi[ifl])
            phix_fl.append(d_dx(psi,dx)[ifl])
            rp_+=1
        if sp_<nsnap and n>=snap_i[sp_]:
            field_t.append(n*dt); field_psi.append(psi.astype(np.complex64).copy()); sp_+=1
        k1=rhs(psi,pi,chi)
        k2=rhs(psi+0.5*dt*k1[0], pi+0.5*dt*k1[1], chi+0.5*dt*k1[2])
        k3=rhs(psi+0.5*dt*k2[0], pi+0.5*dt*k2[1], chi+0.5*dt*k2[2])
        k4=rhs(psi+dt*k3[0], pi+dt*k3[1], chi+dt*k3[2])
        psi=psi+(dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        pi =pi +(dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        chi=chi+(dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        # Dolan-style Dirichlet pin: psi=0=pi on the final grid point at both
        # ends (d_dx/d2_dx2 already leave the derivative at these two points
        # undefined/zero, so this makes the edge treatment explicit instead
        # of an uncontrolled third thing).
        psi[0]=0; pi[0]=0; psi[-1]=0; pi[-1]=0
        if not np.all(np.isfinite(psi)): print(f"  NaN at t={n*dt:.1f}"); break
    out = dict(w_ref=w_ref, OmH=OmH, t=np.array(tsig), psi=np.array(psig),
               probe_r=r[ip], dx=dx, dt=dt, x_match=x_match,
               flux_r=rfl, phi_flux=np.array(phi_fl),
               phidot_flux=np.array(phidot_fl), phix_flux=np.array(phix_fl))
    if nsnap>0:
        out.update(field_t=np.array(field_t), field_psi=np.array(field_psi),
                    x=x, r=r, rp=rp)
    return out

def prony(t, psi, p=10, tmin_frac=0.35, anchor=None, full=False):
    """Matrix-pencil (Hua & Sarkar): fit psi(t) ~ sum_k c_k exp(s_k t) with up
    to p simultaneous complex exponentials -- built to separate genuinely
    beating modes (e.g. fundamental + slow overtone), unlike logfit's single
    straight-line slope.

    Mode selection: if `anchor` (a reference complex omega, e.g. Leaver's or
    a prior logfit estimate) is given, return the recovered mode CLOSEST to
    it in the complex plane -- this is the "physical root" selection the
    original version lacked (it picked argmax|c|, which is not reliably the
    physical mode when multiple comparable-amplitude modes are present, e.g.
    numerical/PML artifacts or the beating documented in floor_check.png).
    Falls back to argmax|c| if no anchor is given.

    If `full=True`, also return the full candidate list (om, c) for
    inspection -- useful for seeing what else the pencil found.
    """
    m = t > tmin_frac*t[-1]; tt=t[m]-t[m][0]; y=psi[m]
    dt=np.median(np.diff(tt)); Nn=len(y); L=min(4*p, Nn//2)
    # Hankel matrices for matrix pencil
    Y0=np.array([y[i:i+L] for i in range(Nn-L)])
    Y1=np.array([y[i+1:i+1+L] for i in range(Nn-L)])
    U,S,Vh=np.linalg.svd(Y0,full_matrices=False)
    r_=max(1, min(p,int(np.sum(S>1e-6*S[0]))))
    U,S,V=U[:,:r_],S[:r_],Vh.conj().T[:,:r_]
    A=U.conj().T@Y1@V@np.diag(1/S); lam,_=np.linalg.eig(A)
    keep=np.abs(lam)>1e-12; lam=lam[keep]        # drop exact-zero (log undefined)
    s=np.log(lam)/dt                    # psi ~ e^{s t}
    om=1j*s                             # e^{-i w t} => w = i s ; E=Re, Gamma=Im
    # amplitudes by least squares
    Van=np.exp(np.outer(tt,s)); c,_,_,_=np.linalg.lstsq(Van,y,rcond=None)
    if anchor is not None:
        k=np.argmin(np.abs(om-anchor))
    else:
        k=np.argmax(np.abs(c))
    if full: return om[k].real, om[k].imag, om, c
    return om[k].real, om[k].imag

def flux_gamma(t, phi, phidot, phix, r, a, m=1, tmin_frac=0.4):
    """Dolan eq. 34 (single l=1 term, no mode sum): energy flux through a
    fixed near-horizon surface, F_E = (r^2+a^2) Re(Phidot* . Phitilde'),
    Phitilde' = phix + i*a*m*phi/(r^2+a^2). For a quasi-stationary QBS this
    is quadratic in the field, F_E ~ exp(2 Gamma t), giving an independent
    Gamma estimate at the horizon interface (Dolan's x_+) instead of at the
    cloud-peak probe. Returns (Gamma_flux, F_E(t))."""
    P = r*r + a*a
    phit_tilde = phix + 1j*a*m*phi/P
    F_E = P*np.real(np.conj(phidot)*phit_tilde)
    msk = t > tmin_frac*t[-1]; tt = t[msk]
    slope = np.polyfit(tt, np.log(np.abs(F_E[msk])), 1)[0]
    return slope/2.0, F_E

def logfit(t,psi,tmin_frac=0.4):
    m=t>tmin_frac*t[-1]; tt=t[m]
    G=np.polyfit(tt,np.log(np.abs(psi[m])),1)[0]
    E=-np.polyfit(tt,np.unwrap(np.angle(psi[m])),1)[0]; return E,G

def freqfilter(t, psi, E_peak, kbins=5.0, keep_frac=0.5):
    """Dolan's frequency-filter (his Sec. IV.A.2, eq. 49): isolate ONE mode by
    a super-Gaussian bandpass around its spectral peak, then log-regress the
    reconstructed single-mode series. This is the tool that gives a CLEAN
    per-mode Gamma when modes are beating -- unlike logfit (fits through the
    beat) and unlike prony (nails E but is Gamma-resolution-limited for the
    ultra-slow superradiant growth). Requires a UNIFORMLY sampled series.

    E_peak: real frequency to centre the filter on (e.g. Leaver E, or a
      power-spectrum peak). kbins: filter half-width in FFT bins (Dolan ~5).
    keep_frac: fraction of the (central) time window used for the final fit,
      since the reconstruction is only faithful near t_max/2 (filter edges).
    Returns (E_fit, Gamma_fit).
    """
    N = len(t); dt = (t[-1]-t[0])/(N-1)
    w = np.fft.fft(psi)
    omega = 2*np.pi*np.fft.fftfreq(N, d=dt)      # angular
    dom = 2*np.pi/(t[-1]-t[0])                    # FFT angular resolution
    # psi ~ e^{-iEt} = e^{+i(-E)t}; in the np.fft basis a_j=sum_k A_k e^{+i w_k t_j}
    # that mode lives in the bin at w_k = -E_peak, NOT +E_peak.
    F = np.exp(-((omega + E_peak)/(kbins*dom))**4)
    psi_n = np.fft.ifft(w*F)                      # reconstructed single mode
    lo = int((0.5-keep_frac/2)*N); hi = int((0.5+keep_frac/2)*N)
    tt = t[lo:hi]; yy = psi_n[lo:hi]
    G = np.polyfit(tt, np.log(np.abs(yy)), 1)[0]
    E = -np.polyfit(tt, np.unwrap(np.angle(yy)), 1)[0]
    return E, G

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--alpha",type=float,default=0.50); ap.add_argument("--a",type=float,default=0.90)
    ap.add_argument("--T",type=float,default=3000.); ap.add_argument("--N",type=int,default=4000)
    ap.add_argument("--xmin",type=float,default=-70.); ap.add_argument("--xmax",type=float,default=220.)
    ap.add_argument("--gmax",type=float,default=1.0); ap.add_argument("--save",type=str,default=None)
    ap.add_argument("--eps_match",type=float,default=0.02)
    ap.add_argument("--nsnap",type=int,default=0,
                    help="save this many full-field psi(x) snapshots (for a space-time heatmap)")
    q=ap.parse_args()
    o=run(q.alpha,q.a,q.xmin,q.xmax,q.N,q.T,gmax=q.gmax,eps_match=q.eps_match,nsnap=q.nsnap)
    w=o["w_ref"]; E1,G1=logfit(o["t"],o["psi"])
    try: E2,G2=prony(o["t"],o["psi"],anchor=w)
    except Exception as e: E2,G2=(float("nan"),float("nan")); print("prony failed:",e)
    G3,_=flux_gamma(o["t"],o["phi_flux"],o["phidot_flux"],o["phix_flux"],o["flux_r"],q.a)
    print(f"alpha={q.alpha} a={q.a} N={q.N} T={q.T} gmax={q.gmax} (probe r={o['probe_r']:.1f}, dx={o['dx']:.3f}, x_match={o['x_match']:.1f})")
    print(f"  Leaver : E={w.real:.6f}  Gamma={w.imag:+.4e}")
    print(f"  logfit : E={E1:.6f}  Gamma={G1:+.4e}   (Gamma {G1/w.imag:+.3f}x)")
    print(f"  prony  : E={E2:.6f}  Gamma={G2:+.4e}   (Gamma {G2/w.imag:+.3f}x)")
    print(f"  flux   :            Gamma={G3:+.4e}   (Gamma {G3/w.imag:+.3f}x)  [at x_match, Dolan eq.34]")
    if q.save is not None:
        save_kw = dict(t=o["t"], psi=o["psi"], w_ref=w)
        if q.nsnap>0:
            save_kw.update(field_t=o["field_t"], field_psi=o["field_psi"],
                            x=o["x"], r=o["r"], rp=o["rp"], x_match=o["x_match"])
        np.savez(q.save, **save_kw); print("saved",q.save)

if __name__=="__main__": main()
