using Gaugefields

@testset "Grid Wilson pseudoscalar reference" begin
    # Grid develop commit 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1,
    # examples/qcdm_wilson_meson_reference.cc in the adjacent independent clone.
    # The contraction follows channel 0 of Grid's official
    # Example_wall_wall_spectrum.cc:
    # https://github.com/paboyle/Grid/blob/0ac72cb6a30ccdc41d664e7e0759f0c8833078f1/examples/Example_wall_wall_spectrum.cc
    #
    # Grid uses D = (4 + m) - H/2 and QCDMeasurements uses D = 1 - κH.
    # At m = 1/(2κ)-4, C_QCDM = C_Grid/(4κ²).
    κ = 0.1
    grid_raw = [
        0.54647487404684458,
        0.019335553983958836,
        0.0042604591777126418,
        0.019335553983958833,
    ]
    expected = grid_raw ./ (4κ^2)

    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0
    )
    measurement = Meson_correlator_measurement(
        U;
        channels=["pseudoscalar"],
        momenta=[(0, 0, 0)],
        fermiontype="Wilson",
        κ,
        BoundaryCondition=[1, 1, 1, -1],
        eps_CG=1e-14,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
    observed = real.(get_value(measure(measurement, U))[:pseudoscalar][:, 1])

    @test observed ≈ expected rtol=5e-12 atol=5e-13
end

@testset "Grid Wilson pseudoscalar non-cold reference" begin
    # This is a compact non-cold configuration that can be reconstructed
    # independently without storing a large gauge-field fixture.  Starting
    # from the cold field, Grid and QCDMeasurements both replace U_1(0) by
    # diag(exp(iθ), exp(-iθ), 1), θ=1.1.  The corresponding Grid program is
    # examples/qcdm_wilson_meson_hot_reference.cc in the adjacent clone, run
    # with QCDM_GRID_FIELD=localized at the same pinned commit as above.
    κ = 0.1
    θ = 1.1
    grid_plaquette = 0.99857707323287903
    grid_raw = [
        0.54571523373660025,
        0.018892421105410458,
        0.0040817974287809429,
        0.018892421105410454,
    ]
    expected = grid_raw ./ (4κ^2)

    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0
    )
    U[1][1, 1, 1, 1, 1, 1] = cis(θ)
    U[1][2, 2, 1, 1, 1, 1] = cis(-θ)

    observed_plaquette =
        get_value(measure(Plaquette_measurement(U; verbose_level=0), U))
    measurement = Meson_correlator_measurement(
        U;
        channels=["pseudoscalar"],
        momenta=[(0, 0, 0)],
        fermiontype="Wilson",
        κ,
        BoundaryCondition=[1, 1, 1, -1],
        # Grid's CG tolerance is a residual norm. LDO compares eps_CG with
        # the squared residual, hence (1e-14)^2 here.
        eps_CG=1e-28,
        MaxCGstep=10_000,
        method_CG="bicg",
        verbose_level=0,
    )
    observed = real.(get_value(measure(measurement, U))[:pseudoscalar][:, 1])

    @test observed_plaquette ≈ grid_plaquette rtol=5e-14 atol=5e-15
    @test observed ≈ expected rtol=5e-12 atol=5e-13
end
