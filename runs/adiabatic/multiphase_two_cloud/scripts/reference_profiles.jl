using Serialization
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

#=
joint_input
├── problem                         Fixed physical and numerical settings
│   ├── spin                        0.99
│   ├── alpha                       0.3
│   ├── components                  Retained Fourier/angular components
│   ├── potential                   :cosine
│   ├── integration                 Radial settings
│   │   ├── step                    0.1
│   │   ├── inner                   -100.0
│   │   ├── outer                   300.0
│   │   ├── match                   5.0
│   │   ├── horizon_order           3
│   │   └── outer_order             7
│   └── …                           Angular tables and other settings
│
├── profile_guess                   Initial frequencies and boundary data
│   ├── frequencies                 Two complex frequencies
│   ├── horizon_coefficients        14 complex coefficients
│   └── outer_coefficients          14 complex coefficients
│
├── cloud_amplitudes                [A₁, A₂], the fixed outer amplitudes
├── packed_guess                    Same guess in 56 real Newton coordinates
│
├── reference                       Information about the Omiya reference
│   ├── primary_index               Position in primary.continuation; currently 72
│   ├── secondary_index             Position in secondary checkpoint.rows; currently 72
│   ├── reference_energy            650.5419114886128
│   ├── secondary_energy_per_amplitude_squared
│   └── normalization_description   Text explaining the amplitude choice
│
├── matching_tolerance              1e-5
├── maximum_newton_iterations        12
└── provenance                      Saved source/input provenance text
=#

campaign_directory = normpath(joinpath(@__DIR__, ".."))

# Read joint run input record to reuse its physical settings, A1, A2, reference point
joint_run_input = deserialize(joinpath(campaign_directory,
    "results/direct_fig7_20261026_v5-input.jls"))

joint_problem_settings = joint_run_input.problem
primary_outer_amplitude, secondary_outer_amplitude = joint_run_input.cloud_amplitudes
@show joint_problem_settings.spin joint_problem_settings.alpha primary_outer_amplitude secondary_outer_amplitude

# Load the primary and secondary branches
reference_input_directory = joinpath(campaign_directory, "work/direct_fig7_20261026_v5/inputs")
reference_setup = load_primary_secondary_setup(joinpath(reference_input_directory,
    "alpha030_cosine_branch.jls"))
secondary_branch_checkpoint = deserialize(joinpath(reference_input_directory,
    "alpha030_cosine_secondary.jls"))

# Grab the primary/secondary index and the corresponding ??? from joint_input
#=
primary_continuation_point
├── amplitude = 0.036390480353373555       Fixed primary outer coefficient A₁
└── solution
    ├── frequency                        Complex fundamental frequency ω₁
    ├── horizon_amplitudes               6 complex horizon coefficients
    ├── outer_amplitudes                 6 complex outer coefficients
    ├── residual                         12 complex matching mismatches
    ├── residual_norm                    2.0333e-12
    ├── scaled_residual_norm             5.5874e-11
    └── iterations                       6 Newton corrections
=#
primary_continuation_index = joint_run_input.reference.primary_index
secondary_continuation_index = joint_run_input.reference.secondary_index
primary_continuation_point = reference_setup.primary.continuation[primary_continuation_index]
secondary_continuation_point = secondary_branch_checkpoint.rows[secondary_continuation_index]

@assert secondary_continuation_point.primary_index == primary_continuation_index
@assert secondary_branch_checkpoint.reading == "averaged"

omiya_primary_solution = primary_continuation_point.solution
omiya_secondary_solution = secondary_continuation_point.solution
@show omiya_primary_solution.outer_amplitudes[1] omiya_secondary_solution.outer_amplitudes[1]

primary_reconstruction = primary_profiles_at(reference_setup, primary_continuation_index)
omiya_primary_profiles = primary_reconstruction.profiles

primary_matching_norm = norm(omiya_primary_profiles.residual) / primary_outer_amplitude
@show primary_matching_norm primary_reconstruction.residual_change

unit_secondary_profiles = secondary_matching_residual(
    omiya_secondary_solution.frequency,
    omiya_secondary_solution.horizon_amplitudes,
    omiya_secondary_solution.outer_amplitudes,
    omiya_primary_profiles, reference_setup.primary.modes,
    joint_problem_settings.spin, joint_problem_settings.alpha, reference_setup.secondary_modes,
    reference_setup.couplings, reference_setup.polar, reference_setup.azimuthal,
    reference_setup.grid;
    save_profile=true, reading=:averaged,
)

omiya_secondary_profiles = (
    from_horizon=(
        positions=unit_secondary_profiles.from_horizon.positions,
        states=secondary_outer_amplitude .* unit_secondary_profiles.from_horizon.states,
    ),
    from_outer=(
        positions=unit_secondary_profiles.from_outer.positions,
        states=secondary_outer_amplitude .* unit_secondary_profiles.from_outer.states,
    ),
    residual=secondary_outer_amplitude .* unit_secondary_profiles.residual,
)

secondary_matching_norm =
    norm(omiya_secondary_profiles.residual) / secondary_outer_amplitude
@show secondary_matching_norm
@assert secondary_matching_norm <= joint_run_input.matching_tolerance