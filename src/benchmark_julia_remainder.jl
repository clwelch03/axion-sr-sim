# Called by benchmark_julia_remainder.py with a temporary data directory.
# File I/O, startup, and compilation are outside the timed regions.
include(joinpath(@__DIR__, "sine_remainder.jl"))

function main(directory)
    @assert Threads.nthreads() == 1
    @assert !get_zero_subnormals()
    for name in ("small", "mixed", "tiny", "cloud", "edges")
        field = reinterpret(Float64, read(joinpath(directory, name * ".input")))
        saved = copy(field)
        output = similar(field)
        sine_remainder(field)  # Compile and warm both entry points.
        sine_remainder!(output, field)
        allocating_times = Float64[]
        inplace_times = Float64[]
        for repeat in 1:11
            if isodd(repeat)
                push!(allocating_times, @elapsed output = sine_remainder(field))
                push!(inplace_times, @elapsed sine_remainder!(output, field))
            else
                push!(inplace_times, @elapsed sine_remainder!(output, field))
                push!(allocating_times, @elapsed output = sine_remainder(field))
            end
        end
        @assert isequal(field, saved)
        @assert isequal(output, sine_remainder(field))
        write(joinpath(directory, name * ".output"), output)
        allocation = @allocated sine_remainder(field)
        inplace_allocation = @allocated sine_remainder!(output, field)
        println(name, " ", sort(allocating_times)[6], " ", sort(inplace_times)[6],
                " ", allocation, " ", inplace_allocation)
    end
    @assert isempty(sine_remainder(Float64[]))
    scalar = fill(0.01)
    @assert size(sine_remainder(scalar)) == ()
    field = collect(range(-1.0, 1.0, length=31))
    @assert sine_remainder(view(field, 1:2:31)) == sine_remainder(field[1:2:31])
end

main(ARGS[1])
