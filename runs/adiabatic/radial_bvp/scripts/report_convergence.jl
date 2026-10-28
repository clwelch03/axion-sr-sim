using TOML
include("inputs.jl")

function convergence_report(paths)
    records = deserialize.(paths)
    reference = last(records)
    rows = map(zip(paths,records)) do (path,r)
        a = r.assessment
        row = Dict{String,Any}(
            "file"=>path,"sha256"=>file_hash(path),"case"=>r.case,
            "return_code"=>r.retcode,"successful"=>r.successful,
            "amplitudes"=>r.amplitudes,"mesh_points"=>length(r.mesh),
            "solve_seconds"=>r.elapsed,"requested_tolerance"=>r.settings.tolerance,
            "boundary_residual"=>a.boundary_residual,
            "off_mesh_defect"=>maximum(a.defect_by_state),
            "passed_residual_checks"=>r.successful && all(isfinite,a.states) &&
                a.boundary_residual <= r.settings.tolerance/100 &&
                maximum(a.defect_by_state) <= r.settings.tolerance,
            "parameter_drift"=>a.parameter_drift,
            "frequencies_real"=>real.(a.profile.frequencies),
            "frequencies_imaginary"=>imag.(a.profile.frequencies))
        hasproperty(r,:compilation_seconds) &&
            (row["compilation_seconds"]=r.compilation_seconds)
        if a.positions == reference.assessment.positions
            fields = reference.assessment.states
            row["relative_profile_difference_by_mode"] =
                vec(maximum(abs.(a.states-fields);dims=2)) ./
                max.(vec(maximum(abs.(fields);dims=2)),1e-30)
            row["maximum_frequency_difference"] = maximum(abs.(
                a.profile.frequencies-reference.assessment.profile.frequencies))
        end
        row
    end
    Dict("reference"=>last(paths),"runs"=>rows)
end

length(ARGS)>=2 || error("Usage: report_convergence.jl OUTPUT.toml RESULT.jls [...]")
open(ARGS[1],"w") do io
    TOML.print(io,convergence_report(ARGS[2:end]))
end
