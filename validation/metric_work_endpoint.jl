# Dependence of the (Boyer–Lindquist-fixed) metric work on the horizon endpoint: it grows like 1/s.
# Superseded by the action-based metric work (validation/action_balance.jl). Run from the package root:  julia --project=. validation/metric_work_endpoint.jl
using AxionCloud, Printf
const AC = AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
state = BackgroundState(1.0, background.a)
for rstar_min in (-60.0, -80.0, -100.0, -120.0)
    mesh = RadialMesh(background; rstar_min, rstar_max=180.0, step=0.4)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    inner_flux = AC.outward_flux(solution, 1)
    state_rate = BackgroundState(-inner_flux.energy, -state.mass * inner_flux.angular_momentum)
    work = boyer_lindquist_metric_work(solution, state, state_rate)
    @printf("r*_min = %6.1f  s_min = %.3e  w = %.6e  w/(E b') = %.6f\n", rstar_min, mesh.offsets[1], work,
            work / (cloud_charges(solution).energy * state_rate.mass))
    flush(stdout)
end
