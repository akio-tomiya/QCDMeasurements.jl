using LinearAlgebra
using Random
using StaticArrays

function dense_projected_bilinear_reference(
    block1,
    block2,
    left,
    right;
    axis,
    origin,
    momentum,
    parity_mask,
    coefficient,
)
    NS = length(block1)
    NC = size(block1[1], 1)
    lattice_size = size(block1[1])[3:end]
    result = zeros(ComplexF64, lattice_size[axis])
    for site in CartesianIndices(lattice_size)
        position = Tuple(site)
        relative = position .- origin
        separation = mod(relative[axis], lattice_size[axis]) + 1
        phase = cis(-2pi * sum(
            momentum[d] * relative[d] / lattice_size[d] for d in 1:4))
        staggered_sign = iseven(sum(
            parity_mask[d] * relative[d] for d in 1:4)) ? 1 : -1
        for sink_color in 1:NC
            P = [block1[source_spin][sink_color, sink_spin, position...]
                 for sink_spin in 1:NS, source_spin in 1:NS]
            Q = [block2[source_spin][sink_color, sink_spin, position...]
                 for sink_spin in 1:NS, source_spin in 1:NS]
            result[separation] += coefficient * staggered_sign * phase *
                tr(left * P * right * Q')
        end
    end
    return result
end

@testset "LatticeMatrices projected meson kernel" begin
    @test QCDMeasurements.projected_bilinear_slices ===
        LatticeMatrices.projected_bilinear_slices
    @test QCDMeasurements.set_global_component! ===
        LatticeMatrices.set_global_component!

    Random.seed!(0x4d45534f4e)
    NS = 4
    NC = 2
    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    lattice_size = (2 * number_of_processes, 3, 2, 4)
    process_grid = (number_of_processes, 1, 1, 1)
    arrays1 = ntuple(
        _ -> randn(ComplexF64, NC, NS, lattice_size...), NS)
    arrays2 = ntuple(
        _ -> randn(ComplexF64, NC, NS, lattice_size...), NS)
    fields1 = ntuple(
        i -> LatticeMatrix(arrays1[i], 4, process_grid; nw=0), NS)
    fields2 = ntuple(
        i -> LatticeMatrix(arrays2[i], 4, process_grid; nw=0), NS)
    left = randn(ComplexF64, NS, NS)
    right = randn(ComplexF64, NS, NS)

    for (axis, origin, momentum, parity_mask, coefficient) in (
        (4, (2, 2, 1, 3), (1, -1, 0, 0), (1, 0, 1, 0), -0.75 + 0.2im),
        (2, (1, 3, 2, 1), (1, 0, -1, 1), (0, 1, 1, 0), 1.25 - 0.1im),
    )
        expected = dense_projected_bilinear_reference(
            arrays1,
            arrays2,
            left,
            right;
            axis,
            origin,
            momentum,
            parity_mask,
            coefficient,
        )
        observed = QCDMeasurements.projected_bilinear_slices(
            fields1,
            fields2,
            left,
            right;
            axis,
            origin,
            momentum,
            parity_mask,
            coefficient,
        )
        @test observed ≈ expected rtol=2e-13 atol=2e-13
    end

    source_field = LatticeMatrix(
        zeros(ComplexF64, 3, 4, lattice_size...), 4, process_grid; nw=0)
    QCDMeasurements.set_global_component!(
        source_field, 2 - 3im, 2, 3, (2, 2, 1, 4))
    source_array = Array(source_field.A)
    global_sum = MPI.Allreduce(sum(source_array), MPI.SUM, MPI.COMM_WORLD)
    global_nonzeros = MPI.Allreduce(
        count(!iszero, source_array), MPI.SUM, MPI.COMM_WORLD)
    @test global_sum == 2 - 3im
    @test global_nonzeros == 1
end
