using Gaugefields

function localized_grid_gaugefield(; theta=1.1)
    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0)
    U[1][1, 1, 1, 1, 1, 1] = cis(theta)
    U[1][2, 2, 1, 1, 1, 1] = cis(-theta)
    return U
end

function grid_pcac_measurement(U; improvement_coefficient=0.0)
    return PCAC_mass_measurement(
        U;
        κ=0.1,
        improvement_coefficient,
        BoundaryCondition=[1, 1, 1, -1],
        # Grid uses a residual norm while LDO compares eps_CG with its square.
        eps_CG=1e-28,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
end

@testset "Wilson PCAC mass and Grid PP/AP reference" begin
    # Grid develop commit 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1.
    # The two channels are the {Gamma5,Gamma5} and
    # {Gamma5,GammaTGamma5} entries of Example_wall_wall_spectrum.cc.
    # Grid's D=(4+m)-H/2 and QCDM's D=1-kappa*H imply the common
    # correlator conversion C_QCDM=C_Grid/(4*kappa^2).
    kappa = 0.1
    grid_pp_raw = [
        0.54647487404684458,
        0.019335553983958836,
        0.0042604591777126418,
        0.019335553983958833,
    ]
    grid_ap_raw = [
        2.3742840420560186e-18,
        -0.015991082544102354,
        -5.3644686855748924e-19,
        0.01599108254410235,
    ]
    expected_pp = grid_pp_raw ./ (4kappa^2)
    expected_ap = grid_ap_raw ./ (4kappa^2)
    grid_pcac = [
        -0.014631123317420421,
        -3.7634439036868414e-17,
        1.8766853380212005,
        3.7634439036868421e-17,
    ]

    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0)
    measurement = grid_pcac_measurement(U)
    result = get_value(measure(measurement, U))

    @test real.(result[:PP]) ≈ expected_pp rtol=5e-12 atol=5e-13
    @test real.(result[:AP]) ≈ expected_ap rtol=5e-12 atol=5e-13
    @test result[:derivative] ≈ ComplexF64[
        (expected_ap[2] - expected_ap[4]) / 2,
        (expected_ap[3] - expected_ap[1]) / 2,
        (expected_ap[4] - expected_ap[2]) / 2,
        (expected_ap[1] - expected_ap[3]) / 2,
    ] rtol=5e-12 atol=5e-13
    @test result[:mass] ≈ result[:derivative] ./ (2 .* result[:PP])
    @test real.(result[:mass]) ≈ grid_pcac rtol=5e-12 atol=5e-13
    @test length(get_solver_diagnostics(measurement)) == 12

    dictionary_measurement = prepare_measurement_from_dict(
        U,
        Dict(
            "methodname" => "PCAC_mass",
            "hop" => kappa,
            "improvement_coefficient" => -0.02,
            "printvalues" => false,
        ),
    )
    @test dictionary_measurement isa PCAC_mass_measurement
    @test dictionary_measurement.improvement_coefficient == -0.02
end

@testset "Wilson PCAC mass non-cold Grid reference" begin
    # Independent Grid run on cold links with only U_1(0) replaced by
    # diag(exp(i*1.1), exp(-i*1.1), 1).
    kappa = 0.1
    grid_pp_raw = [
        0.54571523373660025,
        0.018892421105410458,
        0.0040817974287809429,
        0.018892421105410454,
    ]
    grid_ap_raw = [
        1.3504024247264962e-18,
        -0.015629101933822297,
        -1.5837150473941831e-19,
        0.015629101933822297,
    ]
    expected_pp = grid_pp_raw ./ (4kappa^2)
    expected_ap = grid_ap_raw ./ (4kappa^2)
    grid_pcac = [
        -0.014319832916159693,
        -1.996533320223616e-17,
        1.9144876989265531,
        1.9965333202236163e-17,
    ]

    U = localized_grid_gaugefield()
    result = get_value(measure(grid_pcac_measurement(U), U))
    @test real.(result[:PP]) ≈ expected_pp rtol=5e-12 atol=5e-13
    @test real.(result[:AP]) ≈ expected_ap rtol=5e-12 atol=5e-13
    @test real.(result[:mass]) ≈ grid_pcac rtol=5e-12 atol=5e-13
end

@testset "PCAC improvement term and validation" begin
    pp = ComplexF64[2, 3, 5, 7]
    ap = ComplexF64[11, 13, 17, 19]
    derivative = QCDMeasurements._periodic_symmetric_derivative(ap)
    laplacian = QCDMeasurements._periodic_laplacian(pp)
    @test derivative == ComplexF64[-3, 3, 3, -3]
    @test laplacian == ComplexF64[6, 1, 0, -7]

    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0)
    @test_throws ArgumentError PCAC_mass_measurement(
        U; fermiontype="Staggered")
    @test_throws ArgumentError PCAC_mass_measurement(
        U; improvement_coefficient=Inf)

    short_U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 2; condition="cold", verbose_level=0)
    @test_throws ArgumentError PCAC_mass_measurement(short_U)
end
