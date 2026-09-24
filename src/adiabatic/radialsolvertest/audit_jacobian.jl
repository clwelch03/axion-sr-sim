# Diagnose derivatives at exactly the same frozen prediction as the solver cases.
include("model.jl")
include("correction.jl")
root = normpath(joinpath(@__DIR__, "../../.."))
directory = joinpath(root, "runs/adiabatic/radialsolvertest")
saved = deserialize(joinpath(directory, "seed_cosine.jls"))
context = benchmark_context(saved)
amplitude = last(saved.continuation).amplitude + 0.1
prediction = predict_nonlinear_mode_cubic(saved.continuation, amplitude; history_spacing=0.025)
problem = shooting_problem(context, amplitude, prediction)
state, residual = problem.initial, problem.residual
automatic = ForwardDiff.jacobian(residual, state)
mismatch = residual(state)
automatic_correction = -(automatic \ mismatch)
columns = norm.(eachcol(automatic))
rows = NamedTuple[]
for increment in (1e-7, 1e-8, 1e-9, 1e-10, 1e-11)
    matrix = difference_jacobian(residual, state, increment)
    correction = -(matrix \ mismatch)
    errors = [norm(matrix[:,index]-automatic[:,index])/max(columns[index], eps()) for index in eachindex(state)]
    record = (; increment, matrix_relative=norm(matrix-automatic)/norm(automatic),
        correction_relative=norm(correction-automatic_correction)/norm(automatic_correction),
        worst_column=argmax(errors), worst_column_relative=maximum(errors),
        predicted_by_ad=norm(mismatch+automatic*correction)/norm(mismatch),
        full_step_residual=norm(residual(state+correction)))
    push!(rows, record)
    println(record)
    flush(stdout)
end
scales = max.(abs.(state), 1e-12)
conditioning = (; unscaled=cond(automatic), variable_scaled=cond(automatic*Diagonal(scales)))
serialize(joinpath(directory,"jacobian_audit.jls"), (; amplitude, state, automatic, rows, conditioning))
open(joinpath(directory,"JACOBIAN_AUDIT.md"),"w") do stream
    println(stream,"# Jacobian audit\n\nFrozen cubic prediction at A0=$amplitude.\n")
    println(stream,"AD condition number: $(conditioning.unscaled); after illustrative variable scaling: $(conditioning.variable_scaled).\n")
    println(stream,"| Increment | Relative matrix change | Relative correction change | Worst column | Column error | AD linear-model residual ratio |")
    println(stream,"|---|---:|---:|---:|---:|---:|")
    for row in rows
        println(stream,"| $(row.increment) | $(row.matrix_relative) | $(row.correction_relative) | $(row.worst_column) | $(row.worst_column_relative) | $(row.predicted_by_ad) |")
    end
    println(stream,"\nAD differentiates the discrete RK4 map. This audit tests derivative consistency, not radial discretization accuracy.")
end
