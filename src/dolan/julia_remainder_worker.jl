# Experimental persistent worker: array data stays in shared mapped files;
# stdin/stdout carry only the element count and completion acknowledgement.
using Mmap
include(joinpath(@__DIR__, "sine_remainder.jl"))

"""
    main(input_path, output_path, capacity)

Serve sine-remainder requests using memory-mapped `Float64` buffers.

Map `capacity` elements from each existing file and overwrite requested prefixes
of `output_path`. Compile the kernel, print `ready`, then read element counts
from standard input. Each count must lie in `0:capacity`; print and flush `done`
after each request. Stop on `quit` or end of input and close the output file.
"""
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
