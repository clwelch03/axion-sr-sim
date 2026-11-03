using ForwardDiff

@testset "Chiral potential: vacuum, mass, and quartic coupling" begin
    potential(field) = potential_value(field, :chiral)
    remainder(field) = potential_remainder(field, :chiral)
    first_derivative(field) = ForwardDiff.derivative(potential, field)
    second_derivative(field) = ForwardDiff.derivative(first_derivative, field)
    third_derivative(field) = ForwardDiff.derivative(second_derivative, field)

    # Independent small-field prediction from the quark-mass ratio:
    # U = phi^2/2 - c4*phi^4/24 + O(phi^6).
    ratio = UP_QUARK_MASS / DOWN_QUARK_MASS
    quartic_coefficient = (1 - ratio + ratio^2) / (1 + ratio)^2
    @test potential(0.0) == 0.0
    @test first_derivative(0.0) == 0.0
    @test second_derivative(0.0) ≈ 1.0
    @test third_derivative(0.0) == 0.0
    @test isapprox(ForwardDiff.derivative(third_derivative, 0.0),
        -quartic_coefficient; rtol=1e-13)
    @test remainder(0.0) == 0.0
    @test ForwardDiff.derivative(remainder, 0.0) == 0.0

    # Relative checks preserve sensitivity to a nonzero cubic interaction even
    # where subtracting phi - U'(phi) in Float64 would give zero.
    for field in (1e-60, -1e-20, 1e-7)
        @test isapprox(remainder(field) / field^3,
            quartic_coefficient / 6; rtol=1e-13)
        @test isapprox(ForwardDiff.derivative(remainder, field) / field^2,
            quartic_coefficient / 2; rtol=1e-13)
    end
    println("Chiral mass normalization U''(0) = ", second_derivative(0.0),
        "; z = ", ratio, "; c4 = ", quartic_coefficient)
end

@testset "Chiral energy and source against the original high-precision potential" begin
    potential(field) = potential_value(field, :chiral)
    remainder(field) = potential_remainder(field, :chiral)
    maximum_relative_errors = zeros(4)

    # Use the original subtractive expression at high precision, independently
    # of the rationalized Float64 implementation. 768 bits resolve its severe
    # cancellation for fields as small as 1e-60, including in the source slope.
    setprecision(BigFloat, 768) do
        up, down = BigFloat(UP_QUARK_MASS), BigFloat(DOWN_QUARK_MASS)
        kappa = 4up * down / (up + down)^2
        for magnitude in (1e-60, 1e-20, 1e-10, 1e-4, 0.03,
                          prevfloat(0.1), 0.1, nextfloat(0.1), 0.7, Float64(pi), 3.8, 8.0),
                sign in (-1, 1)
            field = sign * magnitude
            reference_field = BigFloat(field)
            root = sqrt(1 - kappa * sin(reference_field / 2)^2)
            force = sin(reference_field) / root
            curvature = cos(reference_field) / root +
                kappa * sin(reference_field)^2 / (4root^3)
            references = (
                4 / kappa * (1 - root),
                reference_field - force,
                force,
                1 - curvature,
            )
            calculated = (
                potential(field),
                remainder(field),
                ForwardDiff.derivative(potential, field),
                ForwardDiff.derivative(remainder, field),
            )
            for index in eachindex(references)
                reference = Float64(references[index])
                @test isapprox(calculated[index], reference; rtol=5e-13, atol=0)
                maximum_relative_errors[index] = max(maximum_relative_errors[index],
                    abs(calculated[index] / reference - 1))
            end
        end
    end
    println("Chiral maximum relative errors (U, remainder, U', remainder'): ",
        maximum_relative_errors)
end
