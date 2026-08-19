using Random

# Independent transcription of SIMULATeQCD's measureHadrons contraction and
# M1--M8 phase_factor switch:
# src/modules/measureHadrons/measureHadrons.cpp (contractPropagators) and
# src/modules/measureHadrons/measureHadrons.h (phase_factor).
const SIMULATEQCD_MESON_MASKS = (
    (1, 1, 1),
    (0, 0, 0),
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)

function simulateqcd_hadron_reference(first_mass, second_mass, axis, source_axis)
    lattice_size = size(first_mass[1])[3:end]
    transverse = Tuple(d for d in 1:4 if d != axis)
    output = [zeros(ComplexF64, lattice_size[axis]) for _ in 1:8]
    for site in CartesianIndices(lattice_size)
        position = Tuple(site)
        separation = mod(position[axis] - source_axis, lattice_size[axis]) + 1
        contraction = sum(
            conj(first_mass[source_color][sink_color, 1, position...]) *
            second_mass[source_color][sink_color, 1, position...]
            for source_color in 1:3, sink_color in 1:3)
        transverse_zero_based = ntuple(
            i -> position[transverse[i]] - 1, 3)
        for channel in 1:8
            parity = sum(SIMULATEQCD_MESON_MASKS[channel][i] *
                         transverse_zero_based[i] for i in 1:3)
            output[channel][separation] += (iseven(parity) ? 1 : -1) * contraction
        end
    end
    return output
end

@testset "SIMULATeQCD M1--M8 meson reference" begin
    Random.seed!(0x53494d514344)
    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    lattice_size = (2 * number_of_processes, 4, 2, 4)
    process_grid = (number_of_processes, 1, 1, 1)
    first_mass = [randn(ComplexF64, 3, 1, lattice_size...) for _ in 1:3]
    second_mass = [randn(ComplexF64, 3, 1, lattice_size...) for _ in 1:3]
    first_fields = [LatticeMatrix(array, 4, process_grid; nw=0)
                    for array in first_mass]
    second_fields = [LatticeMatrix(array, 4, process_grid; nw=0)
                     for array in second_mass]
    channels = simulateqcd_staggered_channels(4)

    expected = simulateqcd_hadron_reference(first_mass, second_mass, 4, 1)
    observed = [zeros(ComplexF64, lattice_size[4]) for _ in 1:8]
    for source_color in 1:3, channel in 1:8
        # SIMULATeQCD computes dot_prod(first, second) = conj(first)*second.
        # The QCDMeasurements primitive is P*Q†, hence the intentionally
        # reversed block order below.
        observed[channel] .+= QCDMeasurements.projected_bilinear_slices(
            (second_fields[source_color],),
            (first_fields[source_color],),
            ones(ComplexF64, 1, 1),
            ones(ComplexF64, 1, 1);
            axis=4,
            origin=(1, 1, 1, 1),
            momentum=(0, 0, 0, 0),
            parity_mask=channels[channel].parity_mask,
            coefficient=1,
        )
    end

    @test getfield.(channels, :parity_mask) == [
        (1, 1, 1, 0),
        (0, 0, 0, 0),
        (0, 1, 1, 0),
        (1, 0, 1, 0),
        (1, 1, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
    ]
    for channel in 1:8
        @test observed[channel] ≈ expected[channel] rtol=2e-13 atol=2e-13
    end
end
