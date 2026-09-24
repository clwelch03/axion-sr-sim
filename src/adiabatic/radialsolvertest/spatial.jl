"""
    hermite_simpson(left, right, frequency, geometry, spacing, context)

Return the fourth-order Lobatto collocation defect divided by interval length.

Reconstruct the midpoint with cubic Hermite interpolation and enforce Simpson
integration of the ODE. This is radial collocation, independent of the angular
quadrature used to evaluate the nonlinear source.
"""
function hermite_simpson(left, right, frequency, geometry, spacing, context)
    first = radial_rhs(left, frequency, geometry[1], context)
    last = radial_rhs(right, frequency, geometry[3], context)
    middle = (left+right)/2 + spacing/8*(first-last)
    center = radial_rhs(middle, frequency, geometry[2], context)
    return (right-left)/spacing - (first+4center+last)/6
end

"""
    spatial_problem(shooting; method=:multiple, resolution=2, align_mesh=false)

Build a sparse multiple-shooting or radial-collocation system.

Keep all 24 real boundary parameters and introduce the 24 real state components
at every mesh node, scaled by A0. Match boundary series at both ends. Multiple
shooting integrates toward the original matching point from each side, using
`resolution` segments per side. Collocation uses `resolution` as maximum spacing.
Assemble local AD Jacobians into a sparse matrix; never differentiate a dense
global system. The initial profile comes from the same prediction as shooting.
Set `align_mesh=true` to place multiple-shooting endpoints on the existing RK4
grid and retain its step counts, rather than changing the discrete integrator.
"""
function spatial_problem(shooting; method=:multiple, resolution=2, align_mesh=false)
    context, amplitude = shooting.context, shooting.amplitude
    options = context.options
    left_count = method === :multiple ? Int(resolution) : ceil(Int, (options.match-options.inner)/resolution)
    right_count = method === :multiple ? Int(resolution) : ceil(Int, (options.outer-options.match)/resolution)
    positions = vcat(collect(range(options.inner, options.match; length=left_count+1)),
        collect(range(options.match, options.outer; length=right_count+1))[2:end])
    interval_counts = nothing
    if method === :multiple && align_mesh
        inner_indices = round.(Int, range(0, shooting.inner_mesh.count; length=left_count+1))
        outer_indices = round.(Int, range(0, shooting.outer_mesh.count; length=right_count+1))
        inner_positions = options.inner .+ inner_indices .* shooting.inner_mesh.spacing
        outer_positions = reverse(options.outer .+ outer_indices .* shooting.outer_mesh.spacing)
        positions = vcat(inner_positions, outer_positions[2:end])
        interval_counts = vcat(diff(inner_indices), reverse(diff(outer_indices)))
    end
    nodes = length(positions)
    width = length(shooting.initial)
    frequency_indices = [1, width÷2+1]
    frequency, horizon, outer = shooting.boundaries(shooting.initial)
    values = zeros(width, nodes)
    values[:,1], values[:,end] = horizon/amplitude, outer/amplitude
    meshes = [integration_mesh(index <= left_count ? positions[index] : positions[index+1],
        index <= left_count ? positions[index+1] : positions[index], options.step, context;
        count_override=isnothing(interval_counts) ? nothing : interval_counts[index])
        for index in 1:nodes-1]
    for index in 1:left_count
        values[:,index+1] = integrate_cached(amplitude*values[:,index], frequency,
            meshes[index], context)/amplitude
    end
    for index in nodes-1:-1:left_count+1
        candidate = integrate_cached(amplitude*values[:,index+1], frequency,
            meshes[index], context)/amplitude
        index == left_count+1 || (values[:,index] = candidate)
    end
    geometries = [[tortoise_geometry(position, context.spin) for position in
        (positions[index], (positions[index]+positions[index+1])/2, positions[index+1])]
        for index in 1:nodes-1]
    node_indices(index) = (width+(index-1)*width+1):(width+index*width)
    local_columns(index) = vcat(collect(node_indices(index)), collect(node_indices(index+1)), frequency_indices)
    interval_residual = function(local_state, index)
        left, right = local_state[1:width], local_state[width+1:2width]
        frequency = complex(local_state[end-1], local_state[end])
        if method === :multiple
            return index <= left_count ?
                integrate_cached(amplitude*left, frequency, meshes[index], context)/amplitude-right :
                integrate_cached(amplitude*right, frequency, meshes[index], context)/amplitude-left
        end
        return hermite_simpson(amplitude*left, amplitude*right, frequency,
            geometries[index], positions[index+1]-positions[index], context)/amplitude
    end
    boundary_columns = vcat(collect(1:width), collect(node_indices(1)), collect(node_indices(nodes)))
    boundary_residual = function(local_state)
        _, horizon, outer = shooting.boundaries(local_state[1:width])
        return vcat(local_state[width+1:2width]-horizon/amplitude,
            local_state[2width+1:3width]-outer/amplitude)
    end
    initial = vcat(shooting.initial, vec(values))
    residual = function(state)
        blocks = [interval_residual(state[local_columns(index)], index) for index in 1:nodes-1]
        return vcat(blocks..., boundary_residual(state[boundary_columns]))
    end
    jacobian = function(state)
        row_indices, column_indices, entries = Int[], Int[], Float64[]
        function append_block(matrix, columns, row_offset)
            for column in axes(matrix,2), row in axes(matrix,1)
                value = matrix[row,column]
                iszero(value) && continue
                push!(row_indices, row_offset+row)
                push!(column_indices, columns[column])
                push!(entries, value)
            end
        end
        for index in 1:nodes-1
            columns = local_columns(index)
            matrix = ForwardDiff.jacobian(local_state -> interval_residual(local_state, index), state[columns])
            append_block(matrix, columns, (index-1)*width)
        end
        matrix = ForwardDiff.jacobian(boundary_residual, state[boundary_columns])
        append_block(matrix, boundary_columns, (nodes-1)*width)
        return sparse(row_indices, column_indices, entries, length(state), length(state))
    end
    return (; initial, residual, jacobian, positions, width, method, resolution)
end
