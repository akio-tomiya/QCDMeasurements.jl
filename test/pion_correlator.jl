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
    U = Initialize_4DGaugefields(3, 0, 2, 2, 2, 4; condition="cold")
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

    @test observed ≈ expected rtol=1e-12 atol=1e-12
end
