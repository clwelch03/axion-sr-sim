@testset "Adaptive amplitude continuation" begin
    point(amplitude, count=3) = (; amplitude, solution=(; iterations=count))
    seed = point(1.0)
    attempts, saved = [], []
    correct = function(previous, target, attempt)
        push!(attempts, (previous.amplitude, target, attempt))
        attempt == 1 && error("Newton step could not reduce the residual (1.0)")
        return point(target, attempt == 3 ? 12 : 3)
    end
    result = continue_amplitude(seed; amplitude=p -> p.amplitude, correct,
        save=(p, attempt) -> push!(saved, (p.amplitude, attempt)),
        target_amplitude=1.2, fractional_step=0.1, log_io=devnull)
    @test result.amplitude == 1.2
    @test seed.amplitude == 1.0
    @test attempts[1] == (1.0, 1.1, 1)
    @test attempts[2][1] == 1.0  # Retry from the accepted point, not the failed target.
    @test attempts[2][2] ≈ 1.05
    @test attempts[3][2] ≈ 1.05 * 1.075  # Easy correction grows the fractional step.
    @test attempts[4][2] ≈ attempts[3][2] * 1.0375  # Costly correction halves it.
    @test first(saved)[2] == 2
    @test length(saved) == length(attempts) - 1
    @test all(p -> p[1] <= 1.2, saved)

    # A failed shortened final step halves the increment actually attempted.
    empty!(attempts)
    empty!(saved)
    result = continue_amplitude(seed; amplitude=p -> p.amplitude, correct,
        save=(p, attempt) -> push!(saved, (p.amplitude, attempt)),
        target_amplitude=1.03, fractional_step=0.1, log_io=devnull)
    @test attempts[1][2] == 1.03
    @test attempts[2][2] ≈ 1.015
    @test result.amplitude == 1.03

    calls = Ref(0)
    solve = (previous, target, attempt) -> (calls[] += 1; point(target))
    options = (; amplitude=p -> p.amplitude, correct=solve,
        save=(p, attempt) -> nothing, target_amplitude=1.1, log_io=devnull)
    @test continue_amplitude(seed; options..., target_amplitude=1.0) === seed
    @test calls[] == 0
    @test_throws ArgumentError continue_amplitude(seed; options..., fractional_step=NaN)
    @test_throws ArgumentError continue_amplitude(seed; options..., target_amplitude=0.9)
    @test_throws ArgumentError continue_amplitude(seed; options..., fractional_step=0.5)
    @test_throws ErrorException continue_amplitude(seed; options...,
        correct=(args...) -> error("programming error"))
    @test_throws DomainError continue_amplitude(seed; options...,
        validate=p -> throw(DomainError(p.amplitude)))
    @test calls[] == 1  # Even a numerical error during validation is not retried.
    @test_throws ErrorException continue_amplitude(seed; options...,
        save=(args...) -> error("disk write failed"))
    @test calls[] == 2
    @test_throws ErrorException continue_amplitude(seed; options...,
        correct=(args...) -> seed)  # Wrong target must not be saved.

    for exception in (SingularException(1), DomainError(1.0),
            ErrorException("Newton iteration limit reached; residual=1.0"),
            ErrorException("nonfinite Newton residual"))
        calls[] = 0
        failed = (args...) -> (calls[] += 1; throw(exception))
        @test_throws typeof(exception) continue_amplitude(seed; options...,
            correct=failed, fractional_step=0.01, minimum_fractional_step=0.004)
        @test calls[] == 2
    end
end

@testset "Continuation checkpoints and branch records" begin
    profile = TwoPhaseProfile(ComplexF64[0.4], ComplexF64[1], ComplexF64[0.1])
    point = ContinuationPoint(7, (; amplitude=0.1), profile,
        (; tolerance=1e-5, iterations=3), (), nothing)
    branch = ContinuationBranch((; model=:example), [point], (; source="test"))
    mktempdir() do directory
        path = joinpath(directory, "branch.jls")
        write_continuation_checkpoint(path, branch; overwrite=false)
        @test !ispath(path * ".tmp")
        @test only(read_branch(path).points).point_id == 7
        @test only(read_branch(path).points).profile.frequencies == profile.frequencies
        @test_throws ErrorException write_continuation_checkpoint(path, nothing; overwrite=false)
        @test only(read_branch(path).points).point_id == 7
        write_continuation_checkpoint(path, (; branch, restart=(; attempt=4)))
        @test read_branch(path).settings == branch.settings

        write_continuation_checkpoint(path, (; legacy=true))
        @test_throws ArgumentError read_branch(path)
        adapter = (saved; checkpoint_path) -> begin
            @test saved.legacy
            @test checkpoint_path == abspath(path)
            branch
        end
        @test read_branch(path; checkpoint_adapter=adapter) === branch
        @test_throws ArgumentError read_branch(path;
            checkpoint_adapter=(saved; checkpoint_path) -> nothing)
        write_continuation_checkpoint(path, (; branch=nothing))
        @test_throws ArgumentError read_branch(path)

        callback = continuation_newton_checkpoint(path,
            (; cloud_amplitudes=[0.1, 0.2], residual_weights=ones(2)); log_io=devnull)
        callback(0, [1.0, 2.0], [0.2, 0.1])
        callback(1, [1.1, 2.1], [0.01, 0.02])
        saved = deserialize(path)
        @test saved.iteration == 1
        @test saved.state == [1.1, 2.1]
        @test saved.mismatch == [0.01, 0.02]
        @test saved.cloud_amplitudes == [0.1, 0.2]
        @test !ispath(path * ".tmp")
    end
end
