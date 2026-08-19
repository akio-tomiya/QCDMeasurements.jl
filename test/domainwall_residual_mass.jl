using Gaugefields

@testset "Grid Shamir domain-wall residual-mass reference" begin
    # Independent reference generated with Grid develop commit
    # 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1 by
    # references/grid/qcdm_domainwall_mres_reference.cc.  The field is cold
    # except for U_1(0)=diag(exp(i*1.1), exp(-i*1.1), 1), so this also avoids
    # validating only the trivial free field.
    grid_PP = ComplexF64[
        0.14592222922302317 + 5.4286636189814794e-20im,
        0.097088345874940141 + 2.041030642822601e-20im,
        0.059072177389470064 + 6.7600199723056984e-21im,
        0.097088345874940141 + 1.0307496632516971e-21im,
    ]
    grid_J5qP = ComplexF64[
        0.049655759640000491 + 1.6361798280871522e-21im,
        0.0028730807215859382 + 1.9668614644521858e-22im,
        0.00020791404200740403 + 4.3292054124976462e-23im,
        0.0028730807215859378 + 1.1579717795075363e-21im,
    ]

    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    4 % number_of_processes == 0 || error(
        "the Grid reference needs an MPI process count dividing NX=4")
    U = gauge_configuration(
        (4, 4, 4, 4);
        colors=3,
        halo=1,
        start=:cold,
        process_grid=(number_of_processes, 1, 1, 1),
        verbose=0,
    )
    QCDMeasurements.set_global_component!(
        U[1].U, cis(1.1), 1, 1, (1, 1, 1, 1))
    QCDMeasurements.set_global_component!(
        U[1].U, cis(-1.1), 2, 2, (1, 1, 1, 1))
    Gaugefields.set_wing_U!(U[1])

    measurement = DomainWallResidualMassMeasurement(
        U;
        mass=0.1,
        L5=4,
        M=-1,
        BoundaryCondition=[1, 1, 1, -1],
        # Grid specifies a residual norm; LDO compares its square.
        eps_CG=1e-28,
        MaxCGstep=100_000,
        method_CG="bicg",
        verbose_level=0,
    )
    result = redirect_stdout(devnull) do
        get_value(measure(measurement, U))
    end

    @test result[:PP] === result.pseudoscalar_correlator
    @test result[:J5qP] === result.midpoint_pseudoscalar_correlator
    @test result[:mres] === result.residual_mass
    @test result[:PP] ≈ grid_PP rtol=2e-11 atol=2e-13
    @test result[:J5qP] ≈ grid_J5qP rtol=2e-11 atol=2e-13
    @test result[:mres] ≈ grid_J5qP ./ grid_PP rtol=3e-11 atol=3e-13
    diagnostics = get_solver_diagnostics(measurement)
    @test length(diagnostics) == 12
    @test all(diagnostic ->
        diagnostic.recursive_residual_squared <=
        diagnostic.target_residual_squared,
        diagnostics,
    )

    configured = prepare_measurement_from_dict(
        U,
        Dict(
            "methodname" => "Domainwall_residual_mass",
            "L5" => 4,
            "M" => -1.0,
            "mass" => 0.1,
            "eps" => 1e-28,
            "MaxCGstep" => 100_000,
            "printvalues" => false,
            "verbose_level" => 0,
        ),
    )
    @test configured isa DomainWallResidualMassMeasurement
    @test supported_fermions(configured) == (:Domainwall,)

    @test_throws ArgumentError DomainWallResidualMassMeasurement(U; L5=3)
    serial_U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0)
    @test_throws ArgumentError DomainWallResidualMassMeasurement(serial_U)
    two_color_U = gauge_configuration(
        (2 * number_of_processes, 2, 2, 2);
        colors=2,
        halo=1,
        start=:cold,
        process_grid=(number_of_processes, 1, 1, 1),
        verbose=0,
    )
    @test_throws ArgumentError DomainWallResidualMassMeasurement(two_color_U)
end
