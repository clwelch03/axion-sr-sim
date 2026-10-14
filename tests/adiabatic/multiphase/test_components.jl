@testset "Component frequencies and azimuthal numbers" begin
    cloud_frequencies = [0.405, 0.415]
    cloud_m = [1, 2]

    for (q, expected_frequency, expected_m) in (
        ((1, 0), 0.405, 1),
        ((0, 1), 0.415, 2),
        ((2, -1), 0.395, 0),
        ((-1, 2), 0.425, 3),
    )
        result = component_frequency_and_m(q, cloud_frequencies, cloud_m)
        @test result.frequency ≈ expected_frequency
        @test result.m == expected_m
    end

    # Equal m does not imply equal frequency for independent phases.
    third_primary = component_frequency_and_m((3, 0), cloud_frequencies, cloud_m)
    radiative = component_frequency_and_m((-1, 2), cloud_frequencies, cloud_m)
    @test third_primary.m == radiative.m
    @test third_primary.frequency ≉ radiative.frequency

    @test_throws DimensionMismatch component_frequency_and_m((1,), cloud_frequencies, cloud_m)
    @test_throws DimensionMismatch component_frequency_and_m((1, 0), cloud_frequencies, [1])
end
