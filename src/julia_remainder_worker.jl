# Experimental persistent worker: array data stays in shared mapped files;
# stdin/stdout carry only the element count and completion acknowledgement.
using Mmap
include(joinpath(@__DIR__, "sine_remainder.jl"))

function main(input_path, output_path, capacity)
    field = Mmap.mmap(input_path, Vector{Float64}, capacity; grow=false)
    output_file = open(output_path, "r+")
    output = Mmap.mmap(output_file, Vector{Float64}, capacity; grow=false)
    sine_remainder!(view(output, 1:1), view(field, 1:1)) # Compile before ready.
    println("ready")
    flush(stdout)
    for line in eachline(stdin)
        line == "quit" && break
        count = parse(Int, line)
        0 <= count <= capacity || error("invalid element count")
        sine_remainder!(view(output, 1:count), view(field, 1:count))
        println("done")
        flush(stdout)
    end
    close(output_file)
end

main(ARGS[1], ARGS[2], parse(Int, ARGS[3]))
