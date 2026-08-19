using Gaugefields
using LatticeDiracOperators

function direct_wilson_pion_correlator(U; κ=0.1)
    source = Initialize_pseudofermion_fields(U[1], "Wilson")
    solution = similar(source)
    parameters = Dict(
        "Dirac_operator" => "Wilson",
        "κ" => κ,
        "faster version" => true,
        "boundarycondition" => [1, 1, 1, -1],
        "eps_CG" => 1e-14,
        "MaxCGstep" => 10_000,
        "method_CG" => "bicg",
        "verbose_level" => 0,
    )
    D = Dirac_operator(U, source, parameters)
    correlator = zeros(Float64, U[1].NT)

    for source_color in 1:3, source_spin in 1:4
        clear_fermion!(source)
        clear_fermion!(solution)
        setindex_global!(
            source,
            1.0,
            source_color,
            1,
            1,
            1,
            1,
            source_spin,
        )
        solve_DinvX!(solution, D, source)
        for t in 1:U[1].NT, z in 1:U[1].NZ, y in 1:U[1].NY,
                x in 1:U[1].NX, sink_color in 1:3, sink_spin in 1:4
            correlator[t] +=
                abs2(solution[sink_color, x, y, z, t, sink_spin])
        end
    end

    return correlator
end

@testset "Wilson pion correlator uses sink color-spin indices" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 4; condition="cold", verbose_level=0
    )
    expected = direct_wilson_pion_correlator(U)
    measurement = Pion_correlator_measurement(
        U;
        fermiontype="Wilson",
        κ=0.1,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        verbose_level=0,
    )
    observed = get_value(measure(measurement, U))

    generalized = Meson_correlator_measurement(
        U;
        channels=["pseudoscalar"],
        fermiontype="Wilson",
        κ=0.1,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        verbose_level=0,
    )
    generalized_result = get_value(measure(generalized, U))

    @test observed ≈ expected rtol=1e-12 atol=1e-12
    @test generalized_result[:pseudoscalar][:, 1] ≈ expected rtol=1e-12 atol=1e-12
    @test generalized_result.axis == 4
    @test generalized_result.source_position == (1, 1, 1, 1)

    dictionary_measurement = prepare_measurement_from_dict(
        U,
        Dict(
            "methodname" => "Meson_correlator",
            "fermiontype" => "Wilson",
            "hop" => 0.1,
            "channels" => ["pseudoscalar", "vector_1"],
            "momenta" => [[0, 0, 0]],
            "printvalues" => false,
        ),
    )
    @test dictionary_measurement isa Meson_correlator_measurement
    @test getfield.(dictionary_measurement.channels, :name) ==
        [:pseudoscalar, :vector_1]
end

@testset "Wilson pion MPILattice contraction uses the LM kernel" begin
    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    global_nx = 2 * number_of_processes
    process_grid = (number_of_processes, 1, 1, 1)
    U = Initialize_Gaugefields(
        3,
        1,
        global_nx,
        2,
        2,
        4;
        condition="cold",
        isMPILattice=true,
        PEs=process_grid,
        verbose_level=0,
    )
    U_serial = Initialize_4DGaugefields(
        3, 0, global_nx, 2, 2, 4; condition="cold", verbose_level=0,
    )
    measurement_keywords = (
        fermiontype="Wilson",
        κ=0.1,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
    measurement = Pion_correlator_measurement(
        U; measurement_keywords...,
    )
    serial_measurement = Pion_correlator_measurement(
        U_serial; measurement_keywords...,
    )

    observed = get_value(measure(measurement, U))
    expected = get_value(measure(serial_measurement, U_serial))
    @test observed ≈ expected rtol=5e-12 atol=5e-13
    if number_of_processes == 1
        @test observed ≈ [
            18.986971373596035,
            1.3697586050068575,
            0.5488341661984579,
            1.369758605006858,
        ] rtol=5e-12 atol=5e-13
    end
end

@testset "Staggered M2 agrees with the existing pion correlator" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 4; condition="cold", verbose_level=0
    )
    pion_measurement = Pion_correlator_measurement(
        U;
        fermiontype="Staggered",
        mass=0.5,
        eps_CG=1e-14,
        MaxCGstep=10_000,
        verbose_level=0,
    )
    generalized_measurement = Meson_correlator_measurement(
        U;
        fermiontype="Staggered",
        mass=0.5,
        channels=["M2"],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        verbose_level=0,
    )

    expected = get_value(measure(pion_measurement, U))
    observed = get_value(measure(generalized_measurement, U))[
        :M2_pseudoscalar][:, 1]
    @test observed ≈ expected rtol=1e-12 atol=1e-12
end
