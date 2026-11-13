"""
Direct solver for the Newton Jacobian that uses its block structure.

Let `n` be the node block size, `N` the number of nodes, and split the
parameter block into the mode frequencies (`n_f` reals), the horizon
coefficients (`n_h` reals) and the free outer coefficients (`n_o` reals). The
Jacobian rows are

* one interval block of `n` rows per mesh interval, touching the two node blocks
  of the interval and the frequency columns;
* the horizon block of `n` rows, touching the first node block, the frequency
  columns and the horizon-coefficient columns;
* the outer block of `n` rows, touching the last node block, the frequency
  columns and the outer-coefficient columns.

Elimination is Gaussian elimination with partial pivoting for the column order
horizon coefficients, node 1, ..., node `N`, frequencies, outer coefficients,
done in stages that touch only the rows that have entries in the columns being
eliminated. First the horizon coefficients are eliminated from the horizon block,
which leaves `c = n - n_h` rows. Stage `j` stacks the `c` rows left by the previous
stage on the `n` rows of interval `j` and eliminates the `n` columns of node `j`
with row exchanges among those `c + n` rows. It passes on `c` rows that involve only
node `j + 1` and the frequencies; the frequency columns are carried along as `n_f`
extra trailing columns. The last stage is a dense `(c + n)`-square solve for the
last node block, the frequencies and the outer coefficients. Cost is linear in `N`,
and nearly all of it is dense matrix-matrix work on blocks of about 50 x 70.

This is not the recursion `x_{j+1} = M_j x_j` from one end: that is unstable for
a boundary-value problem with growing and decaying solutions, and has no pivoting.
"""

using LinearAlgebra: BlasInt

"""
    BlockFactorization

LU factors of a Jacobian with the structure above: the dense panel of each
stage after elimination, with its pivot vector. The unknowns are in the
column order of the Jacobian (nodes, then frequencies, horizon coefficients, outer
coefficients). `factors \\ b` solves for a vector or a matrix of right-hand sides.
"""
struct BlockFactorization
    node_size::Int
    carried_rows::Int
    frequency_count::Int
    horizon_count::Int
    outer_count::Int
    first_panel::Matrix{Float64}        # n x (n_h + n + n_f): horizon rows; columns [horizon coefficients, node 1, frequencies]
    first_pivots::Vector{BlasInt}
    middle_panels::Array{Float64,3}     # (c + n) x (2n + n_f) x (N - 1): columns [node j, node j + 1, frequencies]
    middle_pivots::Matrix{BlasInt}      # n x (N - 1)
    last_panel::Matrix{Float64}         # (c + n) square: columns [node N, frequencies, outer coefficients]
    last_pivots::Vector{BlasInt}
end

"""
    block_dimensions(layout)

`(node_size, node_count, frequency_count, horizon_count, outer_count, carried_rows)` of the
Jacobian blocks described at the top of this file.
"""
function block_dimensions(layout::UnknownLayout)
    node_size = node_block_size(layout)
    frequency_count = length(frequency_range(layout))
    horizon_count = 2layout.n_channels
    outer_count = 2length(layout.free_outer_channels)
    carried_rows = node_size - horizon_count
    carried_rows == frequency_count + outer_count || error("parameter block does not match the node block size")
    return node_size, layout.n_nodes, frequency_count, horizon_count, outer_count, carried_rows
end

"""
    panel_location(row, column, dimensions)

Where the Jacobian entry `(row, column)` goes: `(panel, stage, i, j)` with `panel` 1 (first), 2 (middle, interval
`stage`), 3 (last), or `panel = 0` for a position outside the block structure.
"""
function panel_location(row::Int, column::Int, dimensions)
    node_size, node_count, frequency_count, horizon_count, _, carried_rows = dimensions
    row_block, row_in_block = fldmod1(row, node_size)
    if column ≤ node_count * node_size
        node, column_in_block = fldmod1(column, node_size)
        row_block < node_count && node == row_block && return (2, row_block, carried_rows + row_in_block, column_in_block)
        row_block < node_count && node == row_block + 1 && return (2, row_block, carried_rows + row_in_block, node_size + column_in_block)
        row_block == node_count && node == 1 && return (1, 1, row_in_block, horizon_count + column_in_block)
        row_block == node_count + 1 && node == node_count && return (3, 1, carried_rows + row_in_block, column_in_block)
        return (0, 0, 0, 0)
    end
    parameter = column - node_count * node_size
    if parameter ≤ frequency_count
        row_block < node_count && return (2, row_block, carried_rows + row_in_block, 2node_size + parameter)
        row_block == node_count && return (1, 1, row_in_block, horizon_count + node_size + parameter)
        return (3, 1, carried_rows + row_in_block, node_size + parameter)
    elseif parameter ≤ frequency_count + horizon_count
        return row_block == node_count ? (1, 1, row_in_block, parameter - frequency_count) : (0, 0, 0, 0)
    end
    return row_block == node_count + 1 ? (3, 1, carried_rows + row_in_block, node_size + parameter - horizon_count) : (0, 0, 0, 0)
end

"""Apply the row exchanges `pivots` (LAPACK convention) to the rows of `block`, in order."""
function apply_row_exchanges!(block::AbstractMatrix, pivots::AbstractVector)
    for (row, partner) in enumerate(pivots)
        row == partner && continue
        for column in 1:size(block, 2)
            block[row, column], block[partner, column] = block[partner, column], block[row, column]
        end
    end
    return block
end

"""
    eliminate_leading_columns!(panel, pivot_count)

Gaussian elimination with partial pivoting of the first `pivot_count` columns of `panel`. On return the
leading columns hold `L` below the diagonal and `U` on and above it, the other columns of the pivot
rows hold `U` rows, and the remaining rows hold the Schur complement of the other columns. Returns the pivots.
"""
function eliminate_leading_columns!(panel::AbstractMatrix, pivot_count::Int)
    pivots = LAPACK.getrf!(view(panel, :, 1:pivot_count))[2]
    trailing_columns = pivot_count + 1:size(panel, 2)
    trailing = view(panel, :, trailing_columns)
    apply_row_exchanges!(trailing, pivots)
    pivot_rows = view(trailing, 1:pivot_count, :)
    BLAS.trsm!('L', 'L', 'N', 'U', 1.0, view(panel, 1:pivot_count, 1:pivot_count), pivot_rows)
    mul!(view(trailing, pivot_count + 1:size(panel, 1), :), view(panel, pivot_count + 1:size(panel, 1), 1:pivot_count),
         pivot_rows, -1.0, 1.0)
    return pivots
end

"""
    stamp_panels(matrix, layout, row_weights, column_weights)

Dense `(first_panel, middle_panels, last_panel)` holding the entries of `D_r matrix D_c` (see `BlockFactorization`
for the panel layouts), before elimination. Throws if `matrix` has a nonzero entry outside the block structure.
"""
function stamp_panels(matrix::SparseMatrixCSC, layout::UnknownLayout, row_weights::AbstractVector,
                      column_weights::AbstractVector)
    dimensions = block_dimensions(layout)
    node_size, node_count, frequency_count, horizon_count, _, carried_rows = dimensions
    first_panel = zeros(node_size, horizon_count + node_size + frequency_count)
    middle_panels = zeros(carried_rows + node_size, 2node_size + frequency_count, node_count - 1)
    last_panel = zeros(carried_rows + node_size, carried_rows + node_size)
    panels = (first_panel, middle_panels, last_panel)
    row_indices, values = rowvals(matrix), nonzeros(matrix)
    node_columns = node_count * node_size
    for column in 1:size(matrix, 2)
        node, column_in_block = column ≤ node_columns ? fldmod1(column, node_size) : (0, 0)
        for position in nzrange(matrix, column)
            value = values[position]
            iszero(value) && continue
            row = row_indices[position]
            scaled = row_weights[row] * value * column_weights[column]
            # Fast path without integer division for the entries of node columns in interval rows (almost all entries):
            # a node column meets interval `node - 1` (rows offset 1:n) and interval `node` (offset n + 1:2n).
            offset = row - (node - 2) * node_size
            if node > 0 && 1 ≤ offset ≤ node_size
                middle_panels[carried_rows + offset, node_size + column_in_block, node - 1] = scaled
            elseif node > 0 && node < node_count && node_size < offset ≤ 2node_size
                middle_panels[carried_rows + offset - node_size, column_in_block, node] = scaled
            else
                panel, stage, i, j = panel_location(row, column, dimensions)
                panel == 0 && error("Jacobian entry ($row, $column) lies outside the block structure")
                panel == 2 ? (middle_panels[i, j, stage] = scaled) : (panels[panel][i, j] = scaled)
            end
        end
    end
    return panels
end

"""
    stamp_panels(blocks::JacobianBlocks, layout, row_weights, column_weights)

The panels of `stamp_panels` filled from the Jacobian blocks of `jacobian_blocks!` (no sparse matrix). The
interval panels are filled in parallel over the Julia threads. Throws if the horizon rows have an entry in an
outer-coefficient column or the outer rows in a horizon-coefficient column.
"""
function stamp_panels(blocks::JacobianBlocks, layout::UnknownLayout, row_weights::AbstractVector,
                      column_weights::AbstractVector)
    node_size, node_count, frequency_count, horizon_count, outer_count, carried_rows = block_dimensions(layout)
    any(!iszero, view(blocks.horizon_wrt_parameters, :, frequency_count + horizon_count + 1:size(blocks.horizon_wrt_parameters, 2))) &&
        error("Jacobian: horizon rows depend on outer coefficients, outside the block structure")
    any(!iszero, view(blocks.outer_wrt_parameters, :, frequency_count + 1:frequency_count + horizon_count)) &&
        error("Jacobian: outer rows depend on horizon coefficients, outside the block structure")
    first_panel = zeros(node_size, horizon_count + node_size + frequency_count)
    middle_panels = zeros(carried_rows + node_size, 2node_size + frequency_count, node_count - 1)
    last_panel = zeros(carried_rows + node_size, carried_rows + node_size)
    last_node_offset = (node_count - 1) * node_size   # first row of the horizon rows and first column of node N, minus 1
    parameter_offset = node_count * node_size         # likewise for the outer rows and the parameter columns
    # Weighted copy of `block[:, source_columns]` into `panel[panel_row_offset .+ (1:node_size), panel_columns]`.
    function stamp!(panel, panel_columns, block, source_columns, row_offset, column_offset, panel_row_offset)
        for (panel_column, source_column) in zip(panel_columns, source_columns), row in 1:node_size
            panel[panel_row_offset + row, panel_column] =
                row_weights[row_offset + row] * block[row, source_column] * column_weights[column_offset + source_column]
        end
    end
    _foreach_in_chunks(1:node_count - 1, () -> nothing) do _, stage
        panel, row_offset = view(middle_panels, :, :, stage), (stage - 1) * node_size
        stamp!(panel, 1:node_size, view(blocks.wrt_left, :, :, stage), 1:node_size, row_offset, row_offset, carried_rows)
        stamp!(panel, node_size + 1:2node_size, view(blocks.wrt_right, :, :, stage), 1:node_size, row_offset,
               row_offset + node_size, carried_rows)
        stamp!(panel, 2node_size + 1:2node_size + frequency_count, view(blocks.wrt_frequencies, :, :, stage),
               1:frequency_count, row_offset, parameter_offset, carried_rows)
    end
    # Horizon rows (row block N): columns [horizon coefficients, node 1, frequencies].
    stamp!(first_panel, 1:horizon_count, blocks.horizon_wrt_parameters, frequency_count + 1:frequency_count + horizon_count,
           last_node_offset, parameter_offset, 0)
    stamp!(first_panel, horizon_count + 1:horizon_count + node_size, blocks.horizon_wrt_node, 1:node_size, last_node_offset, 0, 0)
    stamp!(first_panel, horizon_count + node_size + 1:size(first_panel, 2), blocks.horizon_wrt_parameters, 1:frequency_count,
           last_node_offset, parameter_offset, 0)
    # Outer rows (row block N + 1): columns [node N, frequencies, outer coefficients].
    stamp!(last_panel, 1:node_size, blocks.outer_wrt_node, 1:node_size, parameter_offset, last_node_offset, carried_rows)
    stamp!(last_panel, node_size + 1:node_size + frequency_count, blocks.outer_wrt_parameters, 1:frequency_count,
           parameter_offset, parameter_offset, carried_rows)
    stamp!(last_panel, node_size + frequency_count + 1:node_size + frequency_count + outer_count, blocks.outer_wrt_parameters,
           frequency_count + horizon_count + 1:frequency_count + horizon_count + outer_count, parameter_offset, parameter_offset,
           carried_rows)
    return first_panel, middle_panels, last_panel
end

"""
    BlockFactorization(jacobian, layout, row_weights, column_weights)

Factors of `D_r J D_c` for a Jacobian given as a `SparseMatrixCSC` or as `JacobianBlocks`, with the weights
`D_r = Diagonal(row_weights)`, `D_c = Diagonal(column_weights)`.
"""
function BlockFactorization(jacobian::Union{SparseMatrixCSC, JacobianBlocks}, layout::UnknownLayout,
                            row_weights::AbstractVector, column_weights::AbstractVector)
    node_size, node_count, frequency_count, horizon_count, outer_count, carried_rows = block_dimensions(layout)
    first_panel, middle_panels, last_panel = stamp_panels(jacobian, layout, row_weights, column_weights)
    first_pivots = eliminate_leading_columns!(first_panel, horizon_count)
    middle_pivots = zeros(BlasInt, node_size, node_count - 1)
    # Columns of the next panel that the carried rows occupy: node j + 1 (its leading columns) and the frequencies
    carried_columns = [1:node_size; 2node_size + 1:2node_size + frequency_count]
    carried = view(first_panel, horizon_count + 1:node_size, horizon_count + 1:horizon_count + node_size + frequency_count)
    for stage in 1:node_count - 1
        panel = view(middle_panels, :, :, stage)
        panel[1:carried_rows, carried_columns] .= carried
        middle_pivots[:, stage] = eliminate_leading_columns!(panel, node_size)
        carried = view(panel, node_size + 1:size(panel, 1), node_size + 1:size(panel, 2))
    end
    last_panel[1:carried_rows, 1:node_size + frequency_count] .= carried
    last_pivots = LAPACK.getrf!(last_panel)[2]
    return BlockFactorization(node_size, carried_rows, frequency_count, horizon_count, outer_count,
                              first_panel, first_pivots, middle_panels, middle_pivots, last_panel, last_pivots)
end

function Base.:\(factors::BlockFactorization, right_hand_side::AbstractVector)
    return vec(factors \ reshape(right_hand_side, :, 1))
end

"""
    forward_eliminate(factors, right_hand_side)

Apply the row exchanges and `L` factors of every stage to `right_hand_side` (original row order). Returns the
pivot-row parts of the first stage and of each middle stage, and the right-hand side of the last dense stage.
"""
function forward_eliminate(factors::BlockFactorization, right_hand_side::AbstractMatrix)
    node_size, carried_rows, horizon_count = factors.node_size, factors.carried_rows, factors.horizon_count
    node_count = size(factors.middle_panels, 3) + 1
    block_rows(block) = (block - 1) * node_size + 1:block * node_size   # interval, horizon (block N) and outer (block N + 1) rows
    rhs_count = size(right_hand_side, 2)

    work = right_hand_side[block_rows(node_count), :]
    apply_row_exchanges!(work, factors.first_pivots)
    first_panel = factors.first_panel
    BLAS.trsm!('L', 'L', 'N', 'U', 1.0, view(first_panel, 1:horizon_count, 1:horizon_count), view(work, 1:horizon_count, :))
    mul!(view(work, horizon_count + 1:node_size, :), view(first_panel, horizon_count + 1:node_size, 1:horizon_count),
         view(work, 1:horizon_count, :), -1.0, 1.0)
    first_solved = work[1:horizon_count, :]
    carried = work[horizon_count + 1:node_size, :]

    middle_solved = Array{Float64}(undef, node_size, rhs_count, node_count - 1)
    stage_work = Matrix{Float64}(undef, carried_rows + node_size, rhs_count)
    for stage in 1:node_count - 1
        stage_work[1:carried_rows, :] = carried
        stage_work[carried_rows + 1:end, :] = right_hand_side[block_rows(stage), :]
        apply_row_exchanges!(stage_work, view(factors.middle_pivots, :, stage))
        panel = view(factors.middle_panels, :, :, stage)
        pivot_rows = view(stage_work, 1:node_size, :)
        BLAS.trsm!('L', 'L', 'N', 'U', 1.0, view(panel, 1:node_size, 1:node_size), pivot_rows)
        mul!(view(stage_work, node_size + 1:size(stage_work, 1), :), view(panel, node_size + 1:size(panel, 1), 1:node_size),
             pivot_rows, -1.0, 1.0)
        middle_solved[:, :, stage] = pivot_rows
        carried = stage_work[node_size + 1:end, :]
    end
    last_work = vcat(carried, right_hand_side[block_rows(node_count + 1), :])
    apply_row_exchanges!(last_work, factors.last_pivots)
    return first_solved, middle_solved, last_work
end

function Base.:\(factors::BlockFactorization, right_hand_side::AbstractMatrix)
    node_size, frequency_count, horizon_count = factors.node_size, factors.frequency_count, factors.horizon_count
    node_count = size(factors.middle_panels, 3) + 1
    first_solved, middle_solved, last_work = forward_eliminate(factors, right_hand_side)

    ldiv!(UnitLowerTriangular(factors.last_panel), last_work)
    ldiv!(UpperTriangular(factors.last_panel), last_work)
    solution = Matrix{Float64}(undef, size(right_hand_side))
    node_rows(node) = (node - 1) * node_size + 1:node * node_size
    parameter_rows(offset, count) = node_count * node_size + offset + 1:node_count * node_size + offset + count
    frequencies = last_work[node_size + 1:node_size + frequency_count, :]
    solution[node_rows(node_count), :] = last_work[1:node_size, :]
    solution[parameter_rows(0, frequency_count), :] = frequencies
    solution[parameter_rows(frequency_count + horizon_count, factors.outer_count), :] = last_work[node_size + frequency_count + 1:end, :]

    for stage in node_count - 1:-1:1
        panel = view(factors.middle_panels, :, :, stage)
        node = middle_solved[:, :, stage]
        mul!(node, view(panel, 1:node_size, node_size + 1:2node_size), view(solution, node_rows(stage + 1), :), -1.0, 1.0)
        mul!(node, view(panel, 1:node_size, 2node_size + 1:2node_size + frequency_count), frequencies, -1.0, 1.0)
        ldiv!(UpperTriangular(view(panel, 1:node_size, 1:node_size)), node)
        solution[node_rows(stage), :] = node
    end
    first_panel = factors.first_panel
    horizon = first_solved
    mul!(horizon, view(first_panel, 1:horizon_count, horizon_count + 1:horizon_count + node_size), view(solution, node_rows(1), :), -1.0, 1.0)
    mul!(horizon, view(first_panel, 1:horizon_count, horizon_count + node_size + 1:size(first_panel, 2)), frequencies, -1.0, 1.0)
    ldiv!(UpperTriangular(view(first_panel, 1:horizon_count, 1:horizon_count)), horizon)
    solution[parameter_rows(frequency_count, horizon_count), :] = horizon
    return solution
end
