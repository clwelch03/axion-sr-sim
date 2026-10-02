# julia --project=. src/adiabatic/drivers/secondary_energies.jl primary_branch.jls secondary_branch.jls output.csv
# Energy E2 (per unit (2,2) outer amplitude), boundary fluxes, and the energy
# balance -F/(2 Γ2 E2) of the averaged secondary at every saved branch point.
# Set AXION_MAX_POINTS to stop after that many points.
using Serialization, Printf
include(joinpath(@__DIR__, "continue_secondary.jl"))

function secondary_energies(primary_path, secondary_path, output_path; maximum_points=nothing)
    setup = load_secondary_setup(primary_path)
    (; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
    secondary_branch = deserialize(secondary_path)
    get(secondary_branch, :reading, "fixed_time") == "averaged" ||
        error("secondary energies are defined for the averaged reading only")

    rows = secondary_branch.rows
    last_index = maximum_points === nothing ? length(rows) : min(length(rows), maximum_points)

    temporary = output_path * ".tmp"
    open(temporary, "w") do io
        println(io, "primary_index,amplitude,omega2_real,omega2_imag,energy,",
                    "horizon_flux,outer_flux,balance_ratio")
        for row in rows[1:last_index]
            (; profiles) = primary_profiles_at(setup, row.primary_index)
            (; energy, horizon_flux, outer_flux) = secondary_energy_and_fluxes(
                row.solution, profiles, primary.modes, primary.spin, primary.alpha,
                secondary_modes, couplings, polar, azimuthal, grid,
            )
            frequency = row.solution.frequency
            balance_ratio = -(horizon_flux + outer_flux) / (2imag(frequency) * energy)
            @printf(io, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n",
                    row.primary_index, row.amplitude, real(frequency), imag(frequency),
                    energy, horizon_flux, outer_flux, balance_ratio)
            flush(io)
            @printf("%3d A1=%-9.5g E2=%.6e balance=%.8f\n",
                    row.primary_index, row.amplitude, energy, balance_ratio)
            flush(stdout)
        end
    end
    mv(temporary, output_path; force=true)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 3 || error("require primary_branch.jls secondary_branch.jls output.csv")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    secondary_energies(ARGS...; maximum_points)
end
