import sys, glob, numpy as np, pandas as pd
files = sorted(glob.glob('results/spindown_*_diagnostics.csv') + glob.glob('results/fixed_background_*_diagnostics.csv')
               + glob.glob('results/omiya_fig6_*_diagnostics.csv'))
files = [f for f in files if 'spindown_free' not in f]
ALPHA0 = lambda f: 0.42 if 'alpha0.42' in f else 0.3
def analyse(f):
    d = pd.read_csv(f, comment='#')
    if 'b' not in d: d['b'] = 1.0
    if 'a' not in d: d['a'] = float(__import__('re').search(r'_a(0\.\d+)_', f).group(1))
    d = d[np.isfinite(d.A1) & np.isfinite(d.A2)].sort_values('time').drop_duplicates('time').reset_index(drop=True)
    tau = d.time.values
    lnA1, lnA2 = np.log(d.A1.values), np.log(d.A2.values)
    b = d.b.values; alpha = ALPHA0(f) * b
    n = len(d)
    r1 = np.full(n, np.nan); r2 = np.full(n, np.nan)
    for i in range(n):
        delta = 0.01 * max(tau[i], 1.0)
        j = np.searchsorted(tau, tau[i] - delta, side='right') - 1
        k = np.searchsorted(tau, tau[i] + delta, side='left')
        j = max(j, 0); k = min(k, n - 1)   # one-sided at the ends
        if k <= j: continue
        r1[i] = (lnA1[k] - lnA1[j]) / (tau[k] - tau[j]); r2[i] = (lnA2[k] - lnA2[j]) / (tau[k] - tau[j])
    # convert to current-mass time: d/dt_hat = b d/dtau
    r1 *= b; r2 *= b
    w1, w2 = d.Re_omega1.values, d.Re_omega2.values
    g1, g2 = d.Im_omega1.values, d.Im_omega2.values
    wq = 2*w2 - w1
    gap = wq - alpha
    with np.errstate(invalid='ignore', divide='ignore'):
        k_q = np.sqrt(np.where(gap > 0, wq**2 - alpha**2, np.nan)); vg = k_q / wq
        Gamma_q = np.abs(r1) + 2*np.abs(r2)
        R = 100.0 * (0.3 / alpha)**2
        eps = Gamma_q * R / vg
        assigned_conj = g1 + 2*g2           # |q1| g1 + |q2| g2 with q=(-1,2)
        assigned_signed = -g1 + 2*g2
        ratio_conj = assigned_conj / gap
        ratio_signed = assigned_signed / gap
    out = pd.DataFrame(dict(tau=tau, A1=d.A1, A2=d.A2, ratioA=d.A2/d.A1, b=b, a=d.a, gap=gap, k=k_q, vg=vg, Gq=Gamma_q, eps=eps,
                            rc=ratio_conj, rs=ratio_signed, r1=r1, r2=r2, g1=g1, g2=g2))
    return out
pd.set_option('display.width', 250); pd.set_option('display.float_format', lambda x: f'{x:.4g}')
for f in files:
    o = analyse(f)
    print('=====', f.split('/')[-1], f'rows={len(o)}  b:[{o.b.min():.4f},{o.b.max():.4f}] a:[{o.a.min():.4f},{o.a.max():.4f}]')
    lin = o[(o.A1<1e-4)&(o.tau>0)].dropna(subset=['r1']); print(f'  FD check (A1<1e-4): median|r1/g1-1|={np.median(np.abs(lin.r1/lin.g1-1)) if len(lin) else float("nan"):.2e}')
    print(f'  A1 range [{o.A1.min():.3e},{o.A1.max():.3e}]  A2 range [{o.A2.min():.3e},{o.A2.max():.3e}]  A2/A1 [{o.ratioA.min():.3e},{o.ratioA.max():.3e}]')
    print(f'  gap=w_q-alpha: min {o.gap.min():.4e} max {o.gap.max():.4e}  sign change: {bool((o.gap>0).any() and (o.gap<=0).any())}  frac rows gap<=0: {(o.gap<=0).mean():.3f}')
    for c in ['k','vg','Gq','eps','rc','rs']:
        v = o[c].replace([np.inf,-np.inf],np.nan)
        print(f'  {c:3s}: min {v.min():.4g}  max {v.max():.4g}  (nan {v.isna().sum()})')
    ev = o.dropna(subset=['eps'])
    # also summary of eps and rc extremes by absolute value
    for label, idx in [('argmax A2/A1', o.ratioA.idxmax()), ('argmax A1', o.A1.idxmax()), ('argmax A2', o.A2.idxmax())]:
        r = o.loc[idx]
        print(f'  @{label}: tau={r.tau:.4e} A=({r.A1:.4e},{r.A2:.4e}) b={r.b:.5f} a={r.a:.5f} gap={r.gap:.4e} k={r.k:.4e} vg={r.vg:.4e} Gq={r.Gq:.3e} eps={r.eps:.3e} ratio_conj={r.rc:.4g} ratio_signed={r.rs:.4g}  [row {idx}]')
