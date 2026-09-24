# Jacobian audit

Frozen cubic prediction at A0=10.33565530025236.

AD condition number: 6.626766710230434e9; after illustrative variable scaling: 5.050885085248954e10.

| Increment | Relative matrix change | Relative correction change | Worst column | Column error | AD linear-model residual ratio |
|---|---:|---:|---:|---:|---:|
| 1.0e-7 | 0.05055113188471363 | 2.3455409433387886 | 16 | 0.05745097291399937 | 0.002063758467792077 |
| 1.0e-8 | 0.009035382108778325 | 2.1141254208908653 | 16 | 0.010201620741910056 | 0.0020878761747743454 |
| 1.0e-9 | 0.00010845901017028924 | 0.01005694923425022 | 16 | 0.00012252034188683548 | 9.878197648523452e-6 |
| 1.0e-10 | 1.0866720888979e-6 | 0.00010165341787906827 | 16 | 1.2275708816538392e-6 | 2.775592354584236e-7 |
| 1.0e-11 | 1.1343317385353076e-8 | 6.19058187890018e-6 | 1 | 9.856634879976045e-6 | 4.344257531716627e-6 |

AD differentiates the discrete RK4 map. This audit tests derivative consistency, not radial discretization accuracy.
