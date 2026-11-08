using AxionCloud, LinearAlgebra, SparseArrays
const AC = AxionCloud
background = KerrBackground(0.99)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
mesh = RadialMesh(background; rstar_min=-100.0, rstar_max=200.0, step=0.4)
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
layout = AC.UnknownLayout(problem)
amplitudes = [1e-4, 1e-4]; scale = 1e-4
unknowns = AC.initial_unknowns(problem, layout, amplitudes, scale)
J = AC.jacobian(problem, layout, unknowns, amplitudes, scale)
println("jacobian total:"); @time AC.jacobian(problem, layout, unknowns, amplitudes, scale)
rows, cols, vals = findnz(J)
println("sparse() from triplets:"); @time sparse(rows, cols, vals, size(J)...)
println("lu:"); @time lu(J)
println("residual:"); @time AC.residual(problem, layout, unknowns, amplitudes, scale)
frequency_reals = unknowns[AC.parameter_range(layout)][AC.frequency_range(layout)]
println("node_derivs (751):"); @time [AC.rhs_derivatives(problem, mesh.offsets[n], unknowns[AC.node_range(layout, n)], frequency_reals, scale) for n in 1:layout.n_nodes]
println("node_rhs (751):"); @time AC.node_rhs_blocks(problem, layout, unknowns, frequency_reals, scale)
