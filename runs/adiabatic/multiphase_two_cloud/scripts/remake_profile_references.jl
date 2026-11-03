# Faithful cosine primary/secondary references and free modes on the new
# common domain, normalized to the original benchmark's per-cloud energies.
include("omiya_background_references.jl")

function remake_reference_main()
    length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY")
    inputs, output = abspath.(ARGS)
    ispath(output) && error("Preserve existing reference output: $output")
    config = TOML.parsefile(joinpath(inputs, "request.toml"))
    primary_path, secondary_path = joinpath(inputs,"primary.jls"), joinpath(inputs,"secondary.jls")
    saved_secondary = deserialize(secondary_path)
    old_setup = load_primary_secondary_setup(primary_path; secondary_modes=saved_secondary.secondary_modes)
    p = old_setup.primary
    @assert p.alpha == 0.3 && p.spin == 0.99 && Symbol(p.potential) === :cosine
    options = merge(p.options,(;inner=config["inner"],outer=config["outer"]))
    setup = merge(old_setup,(;primary=merge(p,(;options))))
    target = config["reference_energy"]
    mkpath(output)
    provenance = Dict("created_utc"=>string(now(UTC)),
        "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
        "slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "input_sha256"=>Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in
            (primary_path,secondary_path,joinpath(inputs,"request.toml"))))
    # Start at the original reference on its original domain, then solve again
    # after moving the endpoints. Never just crop and relabel the old profiles.
    old = primary_target(old_setup,target,joinpath(output,"original"))
    current = correct_reference_primary(setup,old,old.amplitude,joinpath(output,"domain"))
    current = match_reference_energy(setup,current,target,joinpath(output,"energy"))
    secondary = secondary_target(setup,saved_secondary,current,joinpath(output,"secondary"))
    case = Dict("name"=>"omiya_cosine","target_energies"=>[target,target])
    export_reference(setup,current,secondary,case,joinpath(output,"omiya_cosine"),provenance)
    export_linear_reference(setup,target,joinpath(output,"linear"),provenance)
end

abspath(PROGRAM_FILE) == @__FILE__ && remake_reference_main()
