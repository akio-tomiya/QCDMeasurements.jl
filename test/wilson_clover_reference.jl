using Gaugefields

function localized_wilson_clover_gaugefield(process_grid)
    U = Initialize_Gaugefields(
        3, 1, 4, 4, 4, 4;
        condition="cold",
        isMPILattice=true,
        PEs=process_grid,
        verbose_level=0,
    )
    QCDMeasurements.set_global_component!(
        U[1].U, cis(1.1), 1, 1, (1, 1, 1, 1))
    QCDMeasurements.set_global_component!(
        U[1].U, cis(-1.1), 2, 2, (1, 1, 1, 1))
    Gaugefields.set_wing_U!(U[1])
    return U
end

@testset "Wilson-clover Grid reference on MPILattice" begin
    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    if 4 % number_of_processes == 0
        # The pinned Grid values describe the global 4^4 field, so they must
        # be independent of whether that field is held on one rank or split
        # over MPI ranks in the x direction.
        U = localized_wilson_clover_gaugefield(
            (number_of_processes, 1, 1, 1))
        measurement = PCACMassMeasurement(
            U;
            fermiontype="WilsonClover",
            κ=0.1,
            cSW=1.2,
            BoundaryCondition=[1, 1, 1, -1],
            # Grid's tolerance is a residual norm, whereas LDO's eps_CG is
            # compared with the squared residual.
            eps_CG=1e-28,
            MaxCGstep=10_000,
            method_CG="bicg",
            verbose_level=0,
        )
        result = get_value(measure(measurement, U))

        # Grid develop commit 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1,
        # examples/qcdm_wilson_clover_reference.cc. Grid uses
        # D=(4+m)-H/2, so at m=1/(2κ)-4 its correlators are divided by
        # 4κ^2 to match QCDMeasurements' D=1-κH normalization.
        grid_pp_raw = [
            0.54585899959258521,
            0.018936150732242110,
            0.0040886011626595367,
            0.018936150732242106,
        ]
        grid_ap_raw = [
            1.3984569142974939e-18,
            -0.015634849086767995,
            -3.5795140852803518e-19,
            0.015634849086767995,
        ]
        grid_pcac = [
            -0.014321325743861909,
            -2.3188560701448903e-17,
            1.9120046765087135,
            2.3188560701448909e-17,
        ]
        normalization = 4 * 0.1^2
        @test isapprox(
            real.(result[:PP]), grid_pp_raw ./ normalization;
            rtol=2e-11, atol=2e-12,
        )
        @test isapprox(
            real.(result[:AP]), grid_ap_raw ./ normalization;
            rtol=2e-11, atol=2e-12,
        )
        @test isapprox(
            real.(result[:mass]), grid_pcac;
            rtol=2e-11, atol=2e-12,
        )
        @test measurement.meson_measurement.first_quark.D.D.cSW == 1.2

        dictionary_measurement = prepare_measurement(
            U,
            Dict(
                "methodname" => "Pion_correlator",
                "Dirac_operator" => "WilsonClover",
                "hop" => 0.1,
                "cSW" => 1.2,
                "method_CG" => "bicg",
                "printvalues" => false,
            ),
        )
        @test dictionary_measurement.D.D.cSW == 1.2
        @test_throws ArgumentError PionCorrelatorMeasurement(
            U;
            fermiontype="WilsonClover",
            method_CG="preconditiond_bicgstab",
        )
    else
        @test_skip 4 % number_of_processes == 0
    end

    serial_U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 4; condition="cold", verbose_level=0)
    @test_throws ArgumentError PionCorrelatorMeasurement(
        serial_U; fermiontype="WilsonClover")
end
