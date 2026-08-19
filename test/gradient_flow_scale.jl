using Gaugefields
using Random

@testset "Gradient-flow history" begin
    Ucold = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )
    measurement = GradientFlowScale_measurement(
        Ucold;
        flow_step_size=0.01,
        number_of_flow_steps=3,
        flow_measure_every=2,
        energy_methods=["plaquette", "clover"],
        topological_charge_methods=["plaquette", "clover"],
    )
    history = get_value(measure(measurement, Ucold))

    @test history.flow_time == [0.0, 0.02, 0.03]
    @test history.plaquette ≈ ones(3)
    @test all(iszero, history.energy_density["plaquette"])
    @test all(value -> isapprox(value, 0.0; atol=1.0e-14),
        history.energy_density["clover"])
    @test all(values -> all(value -> isapprox(value, 0.0; atol=1.0e-14), values),
        values(history.topological_charge))
    @test all(iszero, history.t2_energy["plaquette"])

    cold_estimate = estimate_flow_scales([history])
    @test ismissing(cold_estimate.t0["plaquette"])
    @test ismissing(cold_estimate.w0["clover"])

    Random.seed!(12345)
    Uhot = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="hot", verbose_level=0
    )
    plaquette_measurement = Plaquette_measurement(Uhot)
    plaquette_before = get_value(measure(plaquette_measurement, Uhot))
    hot_history = get_value(measure(
        GradientFlowScale_measurement(
            Uhot;
            flow_step_size=0.01,
            number_of_flow_steps=5,
            measure_topological_charge=false,
        ),
        Uhot,
    ))
    plaquette_after = get_value(measure(plaquette_measurement, Uhot))

    @test plaquette_after ≈ plaquette_before atol=1.0e-14
    @test hot_history.plaquette[end] > hot_history.plaquette[1]
    @test all(isfinite, hot_history.energy_density["plaquette"])
    @test all(value -> value >= -1.0e-13, hot_history.energy_density["plaquette"])
    @test all(value -> value >= -1.0e-13, hot_history.energy_density["clover"])
    @test isempty(hot_history.topological_charge)
end

@testset "Ensemble gradient-flow scales" begin
    times = collect(0.1:0.1:0.6)
    # This gives F(t) = t^2 E(t) = t and W(t) = t dF/dt = t.
    energy = 1.0 ./ times
    history1 = GradientFlowHistory(times, Dict("clover" => energy); volume=16)
    history2 = GradientFlowHistory(times, Dict("clover" => energy); volume=16)
    estimate = estimate_flow_scales([history1, history2]; c=0.3)

    @test estimate.number_of_configurations == 2
    @test estimate.t0["clover"] ≈ 0.3
    @test estimate.w0["clover"] ≈ sqrt(0.3)
    @test estimate.t0_error["clover"] ≈ 0.0 atol=1.0e-14
    @test estimate.w0_error["clover"] ≈ 0.0 atol=1.0e-14
    @test estimate.t2_energy["clover"] ≈ times
    @test estimate.w_observable["clover"] ≈ times

    different_grid = GradientFlowHistory(
        times .+ 0.01, Dict("clover" => energy); volume=16
    )
    @test_throws ArgumentError estimate_flow_scales([history1, different_grid])
    @test_throws ArgumentError estimate_flow_scales(GradientFlowHistory[])
    @test_throws ArgumentError estimate_flow_scales([history1]; c=0.0)
end

@testset "Gradient-flow dictionary API" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )
    measurement = prepare_measurement_from_dict(U, Dict(
        "methodname" => "Gradient_flow_scale",
        "flow_step_size" => 0.02,
        "number_of_flow_steps" => 2,
        "flow_measure_every" => 2,
        "energy_methods" => ["clover"],
        "measure_topological_charge" => false,
        "printvalues" => false,
    ))
    @test measurement isa GradientFlowScale_measurement
    history = get_value(measure(measurement, U))
    @test history.flow_time == [0.0, 0.04]
    @test collect(keys(history.energy_density)) == ["clover"]
end

@testset "Gradient-flow input validation" begin
    U = Initialize_4DGaugefields(
        3, 0, 2, 2, 2, 2; condition="cold", verbose_level=0
    )
    @test_throws ArgumentError GradientFlowScale_measurement(
        U; flow_step_size=0.0
    )
    @test_throws ArgumentError GradientFlowScale_measurement(
        U; number_of_flow_steps=-1
    )
    @test_throws ArgumentError GradientFlowScale_measurement(
        U; energy_methods=["unknown"]
    )
end
