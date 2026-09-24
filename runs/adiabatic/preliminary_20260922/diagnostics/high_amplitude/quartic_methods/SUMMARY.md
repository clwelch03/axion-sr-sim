# Quartic continuation method comparison

All cases start at A0=5.5837601798 and target 5.7 with the same grids,
1e-9 residual tolerance, 12-correction attempt limit, and 600-second
budget checked between attempts. SLURM allows 15 minutes per case.

| Method | Last A0 | Gain/min | Seconds | Accepted | Rejected | Status |
|---|---:|---:|---:|---:|---:|---|
| cubic_h0025_1e10 | 5.70000000 | 0.030688 | 227.3 | 4 | 0 | target_reached |
| cubic_h001_1e10 | 5.66280514 | 0.0075866 | 625.1 | 5 | 2 | budget_reached |
| cubic_h005_1e10 | 5.65139767 | 0.0062885 | 645.3 | 4 | 1 | budget_reached |
| cubic_h001_1e8 | 5.64249985 | 0.0057055 | 617.7 | 4 | 2 | budget_reached |
| secant_small_step | 5.62343666 | 0.003935 | 605.0 | 5 | 0 | budget_reached |
| secant_1e10 | 5.63017979 | 0.0036885 | 755.1 | 4 | 2 | budget_reached |
| secant_1e8 | 5.61877232 | 0.0030696 | 684.4 | 3 | 2 | budget_reached |

Results available for 7 of 7 cases.

Highest measured progress rate: **cubic_h0025_1e10**. Fewer than two cases reached the common endpoint; the ranking is provisional.

No production settings or checkpoints were changed by these experiments.
