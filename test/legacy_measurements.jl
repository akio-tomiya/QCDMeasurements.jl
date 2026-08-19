using Gaugefields
using LatticeDiracOperators
using LinearAlgebra

@testset "Legacy gluonic measurement regressions" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )

    wilson_loops = get_value(measure(
        Wilson_loop_measurement(
            U; Tmax=2, Rmax=2, verbose_level=0, printvalues=false,
        ),
        U,
    ))
    @test wilson_loops == ones(2, 2)

    plaquette_loop = [[(1, +1), (2, +1), (1, -1), (2, -1)]]
    origin_correlation = Correlation_measurement(
        U,
        plaquette_loop,
        plaquette_loop,
        [0, 0, 0, 0];
        originonly=true,
        verbose_level=0,
    )
    volume_correlation = Correlation_measurement(
        U,
        plaquette_loop,
        plaquette_loop,
        [0, 0, 0, 0];
        originonly=false,
        verbose_level=0,
    )
    gluonic_correlation = QCDMeasurements.Guluonic_correlators_measurement(
        U, plaquette_loop, plaquette_loop; verbose_level=0,
    )

    @test get_value(measure(origin_correlation, U)) == 9.0 + 0.0im
    @test get_value(measure(volume_correlation, U)) ==
          9 * prod(size(U[1])[3:end]) + 0.0im
    @test get_value(measure(
        gluonic_correlation,
        U,
        [1, 1, 1, 1],
        [0, 0, 0, 0],
    )) == 9.0 + 0.0im

    density_measurement =
        QCDMeasurements.Topological_charge_density_correlation_measurement(
            U;
            TC_methods=["plaquette", "clover"],
            improved_topological_charge_definition="both",
            verbose_level=0,
        )
    density = get_value(measure(
        density_measurement,
        U,
        [1, 1, 1, 1],
        [1, 0, 0, 0],
    ))
    @test Set(keys(density)) == Set((
        "plaquette",
        "clover",
        "clover improved",
        "clover improved bilson-thompson",
    ))
    @test all(iszero, values(density))
end

@testset "Nontrivial analytic Wilson loop" begin
    lattice_length = 4
    U = Initialize_4DGaugefields(
        3,
        0,
        lattice_length,
        lattice_length,
        lattice_length,
        lattice_length;
        condition="cold",
        verbose_level=0,
    )
    theta = 1.1
    U[1][1, 1, 1, 1, 1, 1] = cis(theta)
    U[1][2, 2, 1, 1, 1, 1] = cis(-theta)

    observed = get_value(measure(
        Wilson_loop_measurement(
            U; Tmax=1, Rmax=1, verbose_level=0, printvalues=false,
        ),
        U,
    ))[1, 1]
    volume = lattice_length^4
    # Only the two x-t plaquettes that contain the modified link differ from
    # identity.  The measurement averages over three spatial directions and
    # divides each trace by NC=3.
    expected = 1 + 4 * (cos(theta) - 1) / (9 * volume)
    @test observed ≈ expected atol=2.0e-15 rtol=0
end

@testset "Exact fermionic measurement limits" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )

    # At kappa=0 the Wilson operator is the identity.  Every Z4 vector then
    # gives tr(D^-n)/V = NC*Nspin = 12 exactly, so this regression does not
    # depend on a random seed or on stochastic convergence.
    chiral = Chiral_condensate_measurement(
        U;
        fermiontype="Wilson",
        κ=0.0,
        Nr=2,
        order=1,
        eps_CG=1.0e-20,
        MaxCGstep=100,
        verbose_level=0,
    )
    chiral_higher = Chiral_condensate_measurement(
        U;
        fermiontype="Wilson",
        κ=0.0,
        Nr=2,
        order=3,
        eps_CG=1.0e-20,
        MaxCGstep=100,
        verbose_level=0,
    )
    @test get_value(measure(chiral, U)) ≈ 12.0 atol=1.0e-14
    @test get_value(measure(chiral_higher, U)) ≈ fill(12.0, 3) atol=1.0e-14

    chiral_from_dictionary = prepare_measurement(U, Dict(
        :methodname => :Chiral_condensate,
        :fermiontype => :Wilson,
        :hop => 0.0,
        :Nr => 1,
        :eps => 1.0e-20,
        :MaxCGstep => 100,
        :verbose_level => 0,
        :printvalues => false,
    ))
    @test chiral_from_dictionary isa ChiralCondensateMeasurement
    @test supported_fermions(chiral_from_dictionary) == (:Wilson, :Staggered)
    @test get_value(measure(chiral_from_dictionary, U)) ≈ 12.0 atol=1.0e-14

    @test_throws ArgumentError PionCorrelatorMeasurement(
        U; fermiontype="Domainwall",
    )
    @test_throws ArgumentError ChiralCondensateMeasurement(
        U; fermiontype="Domainwall",
    )

    eigenmeasurement = QCDMeasurements.Eigenvalue_measurement(
        U;
        fermiontype="Wilson",
        κ=0.0,
        nev=4,
        solver="exact",
        verbose_level=0,
    )
    eigenvalues, eigenvectors = get_value(measure(eigenmeasurement, U))
    operator_matrix = construct_sparsematrix(eigenmeasurement.D(U))
    @test eigenvalues == ones(ComplexF64, 4)
    @test operator_matrix * eigenvectors ≈
          eigenvectors * Diagonal(eigenvalues) atol=1.0e-14
end

@testset "MdagM local spectral density" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )
    eta = 0.2
    spectrum_measurement = QCDMeasurements.MdagMspectrum_measurement(
        U;
        fermiontype="Wilson",
        κ=0.1,
        emin=0.0,
        emax=2.0,
        eta,
        numdatapoints=3,
        position=(1, 1),
        verbose_level=0,
    )
    energies, observed = get_value(measure(spectrum_measurement, U))

    D = construct_sparsematrix(spectrum_measurement.D(U))
    DdagD = Matrix(D' * D)
    source = zeros(ComplexF64, size(DdagD, 1))
    source[1] = 1
    expected = [
        -imag((((energy + eta * im) * I - DdagD) \ source)[1]) / pi
        for energy in energies
    ]
    @test observed ≈ expected atol=5.0e-11 rtol=5.0e-11
    @test all(>(0), observed)

    @test_throws ArgumentError QCDMeasurements.MdagMspectrum_measurement(
        U; numdatapoints=1,
    )
    @test_throws ArgumentError QCDMeasurements.MdagMspectrum_measurement(
        U; eta=0.0,
    )
end

@testset "Measurement parameter validation" begin
    @test_throws ArgumentError QCDMeasurements.initialize_fermion_parameters(
        "unsupported",
    )
    @test_throws ArgumentError QCDMeasurements.initialize_measurement_parameters(
        "unsupported",
    )
    @test QCDMeasurements.Guluonic_correlators_parameters().methodname ==
          "Guluonic_correlators"
end
