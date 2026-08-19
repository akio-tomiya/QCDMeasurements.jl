using CUDA
using JACC
JACC.@init_backend

using Gaugefields
using LatticeMatrices
using MPI
using QCDMeasurements
using Test

MPI.Initialized() || MPI.Init()

@testset "QCDMeasurements H100 CUDA parity" begin
    @test CUDA.functional()
    @test JACC.backend == "cuda"
    @info "QCDMeasurements GPU validation" device = CUDA.name(CUDA.device()) capability = CUDA.capability(CUDA.device())

    U = gauge_configuration(
        (2, 2, 2, 4);
        colors=3,
        halo=1,
        start=:cold,
        process_grid=(1, 1, 1, 1),
        verbose=0,
    )
    @test gauge_backend(U) isa LatticeMatricesBackend
    @test getproperty(U[1].U, :A) isa CUDA.CuArray

    plaquette = get_value(measure(PlaquetteMeasurement(U), U))
    polyakov = get_value(measure(PolyakovMeasurement(U), U))
    energy = get_value(measure(EnergyDensityMeasurement(U), U))
    topological_charge =
        get_value(measure(TopologicalChargeMeasurement(U), U))

    @test plaquette ≈ 1.0 atol=2e-14 rtol=0
    @test polyakov ≈ 3.0 + 0.0im atol=2e-14 rtol=0
    @test energy ≈ 1.0 atol=2e-14 rtol=0
    @test all(value -> isapprox(value, 0.0; atol=2e-14),
        values(topological_charge))

    measurement_keywords = (
        fermiontype="Wilson",
        κ=0.1,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
    pion = PionCorrelatorMeasurement(U; measurement_keywords...)
    observed_pion = get_value(measure(pion, U))
    expected_pion = [
        18.986971373596035,
        1.3697586050068575,
        0.5488341661984579,
        1.369758605006858,
    ]
    @test observed_pion ≈ expected_pion rtol=5e-12 atol=5e-13

    meson = MesonCorrelatorMeasurement(
        U;
        channels=["pseudoscalar"],
        momenta=[(0, 0, 0)],
        measurement_keywords...,
    )
    observed_meson =
        get_value(measure(meson, U))[:pseudoscalar][:, 1]
    @test observed_meson ≈ expected_pion rtol=5e-12 atol=5e-13

    # Reconstruct the same non-cold one-link SU(3) field used by the frozen
    # Grid PCAC reference without host-side scalar indexing.  The LM setter is
    # a JACC kernel, so this also exercises the intended CUDA construction path.
    U_nontrivial = gauge_configuration(
        (4, 4, 4, 4);
        colors=3,
        halo=1,
        start=:cold,
        process_grid=(1, 1, 1, 1),
        verbose=0,
    )
    QCDMeasurements.set_global_component!(
        U_nontrivial[1].U, cis(1.1), 1, 1, (1, 1, 1, 1))
    QCDMeasurements.set_global_component!(
        U_nontrivial[1].U, cis(-1.1), 2, 2, (1, 1, 1, 1))
    @test get_value(measure(PlaquetteMeasurement(U_nontrivial), U_nontrivial)) ≈
          0.998577073232879 rtol=2e-13 atol=2e-14

    # CPU/CUDA parity for a nontrivial Wilson-flow trajectory.  The flow and
    # energy definitions are independently checked against SIMULATeQCD by the
    # normal suite; these frozen values specifically guard the H100 backend.
    flow_history = redirect_stdout(devnull) do
        get_value(measure(
            GradientFlowScaleMeasurement(
                U_nontrivial;
                flow_step_size=0.01,
                number_of_flow_steps=2,
                energy_methods=["plaquette", "clover"],
                measure_topological_charge=false,
            ),
            U_nontrivial,
        ))
    end
    @test flow_history.plaquette ≈ [
        0.998577073232879,
        0.9987541739598155,
        0.9989133245620224,
    ] rtol=2e-12 atol=2e-13
    @test flow_history.energy_density["plaquette"] ≈ [
        0.05122536361635488,
        0.044849737446643,
        0.03912031576719288,
    ] rtol=3e-10 atol=5e-13
    @test flow_history.energy_density["clover"] ≈ [
        0.004653811866959022,
        0.004435109063318777,
        0.004177219274613179,
    ] rtol=3e-10 atol=5e-13

    # Grid develop commit 0ac72cb6: Gamma5/Gamma5 and
    # Gamma5/GammaTGamma5 on the same non-cold field.  Grid's operator
    # normalization gives C_QCDM=C_Grid/(4*kappa^2).
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
    grid_pcac = [
        -0.014319832916159693,
        -1.996533320223616e-17,
        1.9144876989265531,
        1.9965333202236163e-17,
    ]
    pcac_result = redirect_stdout(devnull) do
        get_value(measure(
            PCACMassMeasurement(
                U_nontrivial;
                κ=kappa,
                BoundaryCondition=[1, 1, 1, -1],
                eps_CG=1e-28,
                MaxCGstep=10_000,
                method_CG="bicg",
                verbose_level=0,
            ),
            U_nontrivial,
        ))
    end
    @test real.(pcac_result[:PP]) ≈ grid_pp_raw ./ (4kappa^2) rtol=5e-12 atol=5e-13
    @test real.(pcac_result[:AP]) ≈ grid_ap_raw ./ (4kappa^2) rtol=5e-12 atol=5e-13
    @test real.(pcac_result[:mass]) ≈ grid_pcac rtol=5e-12 atol=5e-13

    # Wilson--clover uses the same LatticeMatrices/JACC operator on CPU and
    # CUDA.  Compare the H100 result directly with the frozen Grid clover run
    # on this nonzero-field-strength configuration.
    clover_grid_pp_raw = [
        0.54585899959258521,
        0.018936150732242110,
        0.0040886011626595367,
        0.018936150732242106,
    ]
    clover_grid_ap_raw = [
        1.3984569142974939e-18,
        -0.015634849086767995,
        -3.5795140852803518e-19,
        0.015634849086767995,
    ]
    clover_grid_pcac = [
        -0.014321325743861909,
        -2.3188560701448903e-17,
        1.9120046765087135,
        2.3188560701448909e-17,
    ]
    clover_measurement = PCACMassMeasurement(
        U_nontrivial;
        fermiontype="WilsonClover",
        κ=kappa,
        cSW=1.2,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-28,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
    clover_result = redirect_stdout(devnull) do
        get_value(measure(clover_measurement, U_nontrivial))
    end
    @test real.(clover_result[:PP]) ≈
          clover_grid_pp_raw ./ (4kappa^2) rtol=2e-11 atol=2e-12
    @test real.(clover_result[:AP]) ≈
          clover_grid_ap_raw ./ (4kappa^2) rtol=2e-11 atol=2e-12
    @test real.(clover_result[:mass]) ≈
          clover_grid_pcac rtol=2e-11 atol=2e-12
    @test clover_measurement.meson_measurement.first_quark.D.D.cSW == 1.2

    # The same non-cold field, Shamir parameters, and physical/midpoint
    # contractions were evaluated independently with Grid's
    # ImportPhysicalFermionSource, ExportPhysicalFermionSolution, and
    # ContractJ5q APIs.
    grid_domainwall_PP = ComplexF64[
        0.14592222922302317 + 5.4286636189814794e-20im,
        0.097088345874940141 + 2.041030642822601e-20im,
        0.059072177389470064 + 6.7600199723056984e-21im,
        0.097088345874940141 + 1.0307496632516971e-21im,
    ]
    grid_domainwall_J5qP = ComplexF64[
        0.049655759640000491 + 1.6361798280871522e-21im,
        0.0028730807215859382 + 1.9668614644521858e-22im,
        0.00020791404200740403 + 4.3292054124976462e-23im,
        0.0028730807215859378 + 1.1579717795075363e-21im,
    ]
    domainwall_measurement = DomainWallResidualMassMeasurement(
        U_nontrivial;
        mass=0.1,
        L5=4,
        M=-1,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-28,
        MaxCGstep=100_000,
        method_CG="bicg",
        verbose_level=0,
    )
    domainwall_result = redirect_stdout(devnull) do
        get_value(measure(domainwall_measurement, U_nontrivial))
    end
    @test domainwall_result[:PP] ≈
          grid_domainwall_PP rtol=2e-11 atol=2e-13
    @test domainwall_result[:J5qP] ≈
          grid_domainwall_J5qP rtol=2e-11 atol=2e-13
    @test domainwall_result[:mres] ≈
          grid_domainwall_J5qP ./ grid_domainwall_PP rtol=3e-11 atol=3e-13
    @test all(diagnostic ->
        diagnostic.recursive_residual_squared <=
        diagnostic.target_residual_squared,
        get_solver_diagnostics(domainwall_measurement),
    )

    CUDA.synchronize()
end
