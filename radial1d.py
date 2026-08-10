#!/usr/bin/env python
"""1+1D time-domain massive-scalar QBS on Kerr, following Dolan 2012 (1212.1477).

Field: Phi = (1/r) psi(t,x) Y_lm e^{im phi}, x = BL tortoise (dx/dr=(r^2+a^2)/Delta).
Single-mode (l=m=1) truncation of Dolan Eq.(25) [c_ll = <11|cos^2 t|11> = 1/5;
the l->l+-2 coupling is O(c^2) small for a QBS and dropped in v1]:

  D psi_tt = (r^2+a^2)^2 psi_xx + [2iam(r^2+a^2) - 2a^2 Delta/r] psi_x
             - 4iamr psi_t - Vc psi
  D  = (r^2+a^2)^2 - (1-c_ll) a^2 Delta
  Vc = Delta[ l(l+1) + (2/r)(1 - a^2/r) + 2iam/r + mu^2(r^2 + c_ll a^2) ]   (M=1)

STAGE A here: validate the spatial operator by plugging in the Leaver (omega, R):
psi = r*R(r); the frequency-domain residual (psi_t->-i w psi, psi_tt->-w^2 psi)
must vanish. This checks the transcription before any time stepping.
"""
import numpy as np, mpmath as mp, sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver, leaver_eigenfunction as lef

def horizons(a): b=np.sqrt(1-a*a); return 1+b,1-b
def tortoise_x(r,a):
    rp,rm=horizons(a)
    return r + (2*rp/(rp-rm))*np.log(r-rp) - (2*rm/(rp-rm))*np.log(r-rm)

def build_grid(a, xmin, xmax, N):
    """Uniform x grid; invert to r(x) via a monotone spline of x(r)."""
    rp,_=horizons(a)
    # dense r-samples covering the needed range, mapped to x, then interpolate r(x)
    r_lo=rp+1e-6; r_hi=2.0
    # expand r_hi until x(r_hi) exceeds xmax
    while tortoise_x(r_hi,a)<xmax: r_hi*=1.5
    rr=np.linspace(r_lo,r_hi,400000); xx=tortoise_x(rr,a)
    x=np.linspace(xmin,xmax,N)
    r=np.interp(x,xx,rr)
    return x,r

def coeffs(r,a,mu,l=1,m=1,c_ll=0.2):
    Delta=r*r-2*r+a*a; P=r*r+a*a
    D  = P*P - (1-c_ll)*a*a*Delta
    A2 = P*P                                    # psi_xx coefficient
    A1 = 2j*a*m*P - 2*a*a*Delta/r               # psi_x  coefficient
    At = -4j*a*m*r                              # psi_t  coefficient (frame drag)
    Vc = Delta*( l*(l+1) + (2/r)*(1-a*a/r) + 2j*a*m/r + mu*mu*(r*r + c_ll*a*a) )
    return D,A2,A1,At,Vc,Delta,P

def d_dx(f,dx):    # 4th-order centered 1st derivative
    g=np.zeros_like(f); g[2:-2]=(-f[4:]+8*f[3:-1]-8*f[1:-3]+f[:-4])/(12*dx)
    g[1]=(f[2]-f[0])/(2*dx); g[-2]=(f[-1]-f[-3])/(2*dx); return g
def d2_dx2(f,dx):  # 4th-order centered 2nd derivative
    g=np.zeros_like(f); g[2:-2]=(-f[4:]+16*f[3:-1]-30*f[2:-2]+16*f[1:-3]-f[:-4])/(12*dx*dx)
    g[1]=(f[2]-2*f[1]+f[0])/dx**2; g[-2]=(f[-1]-2*f[-2]+f[-3])/dx**2; return g

def eigenmode_seed(x, r, alpha, a, w, l=1, m=1, Nseries=400, amplitude=1.0):
    """Clean single-mode initial data: psi = amplitude * normalized(r*R(r)), with
    R(r) the Leaver bound-state radial profile at the eigenfrequency w, AND the
    matching momentum pi = -i*w*psi so the state is a pure e^{-i w t} eigenmode.

    This is the transient-free alternative to the Gaussian seed (which starts with
    pi=0 and the wrong radial shape, so it rings as it sheds onto the true mode).
    On a linear (lambda=0) run a correct eigenmode seed reproduces the Leaver
    growth rate from t=0 with no shedding. Returns (psi, pi) as complex arrays on
    the grid r(x); normalized so max|psi| = amplitude (same peak convention as the
    Gaussian seed, so lambda_quartic*|psi|^2 strengths are comparable)."""
    coef = lef.series_coeffs(mp.mpc(w), alpha, a, l, m, Nseries)
    R = np.array([complex(lef.R_of_r(mp.mpf(ri), mp.mpc(w), alpha, a, coef))
                  for ri in r])
    psi = r * R
    psi = amplitude * psi / np.max(np.abs(psi))
    pi = -1j * w * psi
    return psi.astype(complex), pi.astype(complex)


def stageA(alpha=0.30, a=0.90):
    mp.mp.dps=30
    w=complex(leaver.solve(alpha,a)); print(f"Leaver: E={w.real:.8f} Gamma={w.imag:.4e}")
    x,r=build_grid(a,-40.0,120.0,6000); dx=x[1]-x[0]
    # Leaver eigenfunction on the grid: psi = r*R(r)
    coef=lef.series_coeffs(mp.mpc(w),alpha,a,1,1,400)
    R=np.array([complex(lef.R_of_r(mp.mpf(ri),mp.mpc(w),alpha,a,coef)) for ri in r])
    psi=r*R
    # normalize + zero the tiny-but-noisy far tail beyond where series is trusted
    psi=psi/np.max(np.abs(psi))
    D,A2,A1,At,Vc,Delta,P=coeffs(r,a,alpha)
    px=d_dx(psi,dx); pxx=d2_dx2(psi,dx)
    # frequency domain: psi_t->-i w psi, psi_tt->-w^2 psi
    res = A2*pxx + A1*px + At*(-1j*w)*psi - Vc*psi + D*(w*w)*psi
    # compare residual to the size of the individual terms (relative)
    scale=np.abs(A2*pxx)+np.abs(Vc*psi)+np.abs(D*w*w*psi)+1e-300
    m=(r>1.6)&(r<40)   # bulk region where R is well-resolved & series-converged
    relres=np.abs(res[m])/scale[m]
    print(f"grid: N={len(x)} dx={dx:.4f}, r in [{r[0]:.3f},{r[-1]:.1f}]")
    print(f"freq-domain residual (bulk r=1.6..40): median={np.median(relres):.2e} max={np.max(relres):.2e}")
    print("  -> if <~1e-3 the equation transcription is correct (limited by FD + series tail)")

if __name__=="__main__": stageA()
