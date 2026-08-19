using Gaugefields

@testset "Basic gluonic measurements" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )

    plaquette = get_value(measure(Plaquette_measurement(U), U))
    polyakov = get_value(measure(Polyakov_measurement(U), U))
    energy = get_value(measure(Energy_density_measurement(U), U))
    topological_charge =
        get_value(measure(Topological_charge_measurement(U), U))

    @test plaquette ≈ 1.0
    @test polyakov ≈ 3.0 + 0.0im
    @test energy ≈ 1.0
    @test all(value -> isapprox(value, 0.0; atol=1.0e-14),
        values(topological_charge))
end
