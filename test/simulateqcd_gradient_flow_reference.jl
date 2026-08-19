using Gaugefields
using LinearAlgebra

"""
Construct the non-Abelian SU(3) field used for the SIMULATeQCD v1.2.0
(git 767a1b1) cross-code comparison.  The formula is deterministic and avoids
shipping a binary configuration in the test suite.
"""
function simulateqcd_reference_field()
    U = Initialize_4DGaugefields(
        3, 0, 4, 4, 4, 4; condition="cold", verbose_level=0
    )
    generators = [
        ComplexF64[0 1 0; 1 0 0; 0 0 0],
        ComplexF64[0 -im 0; im 0 0; 0 0 0],
        ComplexF64[0 0 1; 0 0 0; 1 0 0],
        ComplexF64[0 0 -im; 0 0 0; im 0 0],
    ]
    for μ in 1:4, it in 1:4, iz in 1:4, iy in 1:4, ix in 1:4
        phase = 0.18 +
                0.07 * sinpi(((μ + 1) * ix + (2μ + 1) * iy + (μ + 3) * iz) / 2) +
                0.04 * cospi(((μ + 2) * it + ix + 2iz) / 2)
        link = exp(im * phase * generators[μ])
        for column in 1:3, row in 1:3
            U[μ][row, column, ix, iy, iz, it] = link[row, column]
        end
    end
    return U
end

@testset "SIMULATeQCD Wilson-flow reference" begin
    # Generated on an NVIDIA H100 with SIMULATeQCD v1.2.0, Wilson force,
    # fixed-step RK3, dt=0.01.  Its topology convention has the opposite sign
    # to QCDMeasurements, so the stored Q references below use our sign.
    reference_plaquette = [
        0.99640017,
        0.99677379,
        0.99709911,
        0.99738267,
        0.99763010,
        0.99784622,
    ]
    reference_clover_energy = [
        0.044165739,
        0.041663721,
        0.039431116,
        0.037436974,
        0.035653889,
        0.034057610,
    ]
    reference_clover_charge = -[
        1.86847372839524e-3,
        1.82574887807958e-3,
        1.78603462424029e-3,
        1.74907660463734e-3,
        1.71462703141252e-3,
        1.68244987512044e-3,
    ]
    reference_improved_charge = -[
        3.34620251122277e-4,
        3.08187799657249e-4,
        2.83614310731973e-4,
        2.60959772766252e-4,
        2.40212551377418e-4,
        2.21312332199286e-4,
    ]

    U = simulateqcd_reference_field()
    history = get_value(measure(
        GradientFlowScale_measurement(
            U;
            flow_step_size=0.01,
            number_of_flow_steps=5,
            energy_methods=["plaquette", "clover"],
            topological_charge_methods=["clover"],
            improved_topological_charge_definition="bilson_thompson",
        ),
        U,
    ))

    # The reference file prints plaquette and E to eight significant digits.
    @test history.plaquette ≈ reference_plaquette atol=5.0e-8 rtol=0
    @test history.energy_density["clover"] ≈ reference_clover_energy atol=5.0e-8 rtol=0
    @test history.topological_charge["clover"] ≈ reference_clover_charge atol=5.0e-10 rtol=0
    @test history.topological_charge["clover improved bilson-thompson"] ≈
          reference_improved_charge atol=5.0e-10 rtol=0

    # The site-local density used by the correlation measurement must integrate
    # to the same improved charge as the volume-wide implementation.
    local_measurement =
        QCDMeasurements.Topological_charge_density_correlation_measurement(U)
    local_improved_charge = 0.0
    for it in 1:4, iz in 1:4, iy in 1:4, ix in 1:4
        position = [ix, iy, iz, it]
        q_clover = QCDMeasurements.calculate_topological_charge_clover(
            U,
            position,
            local_measurement.temp_UμνTA,
            local_measurement._temporary_matrices,
        )
        local_improved_charge +=
            QCDMeasurements.calculate_topological_charge_bilson_thompson(
            U,
            position,
            local_measurement.temp_UμνTA,
            q_clover,
            local_measurement._temporary_matrices,
        )
    end
    @test local_improved_charge ≈
          history.topological_charge["clover improved bilson-thompson"][1] atol=1.0e-12
end

@testset "Improved topological-charge definitions" begin
    U = simulateqcd_reference_field()

    default_measurement = Topological_charge_measurement(U; TC_methods=["clover"])
    default_values = get_value(measure(default_measurement, U))
    @test default_measurement.improved_topological_charge_definition == "alexandrou"
    @test haskey(default_values, "clover improved")
    @test !haskey(default_values, "clover improved bilson-thompson")

    # Eq. (15)-(17) of Alexandrou, Athenodorou, and Jansen (2015):
    # q_imp = (5/3)q_clover - (1/12)q_rect, with the factor two in q_rect.
    q_clover = QCDMeasurements.calculate_topological_charge_clover(
        U,
        default_measurement.temp_UμνTA,
        default_measurement._temporary_gaugefields,
    )
    rectangle_loops = QCDMeasurements.calc_UμνTA!(
        default_measurement.temp_UμνTA,
        "rect",
        U,
        default_measurement._temporary_gaugefields,
    )
    q_rectangle = 2 * QCDMeasurements.calc_Q(
        default_measurement.temp_UμνTA,
        rectangle_loops,
        U,
    )
    @test default_values["clover improved"] ≈
          (5 / 3) * q_clover - (1 / 12) * q_rectangle atol=1.0e-14

    bilson_values = get_value(measure(
        Topological_charge_measurement(
            U;
            TC_methods=["clover"],
            improved_topological_charge_definition="bilson_thompson",
        ),
        U,
    ))
    @test !haskey(bilson_values, "clover improved")
    @test haskey(bilson_values, "clover improved bilson-thompson")
    @test bilson_values["clover improved bilson-thompson"] ≈
          -3.34620251122277e-4 atol=5.0e-10 rtol=0

    both_values = get_value(measure(
        Topological_charge_measurement(
            U;
            TC_methods=["clover"],
            improved_topological_charge_definition="both",
        ),
        U,
    ))
    @test both_values["clover improved"] ≈ default_values["clover improved"]
    @test both_values["clover improved bilson-thompson"] ≈
          bilson_values["clover improved bilson-thompson"]

    @test_throws ArgumentError Topological_charge_measurement(
        U;
        improved_topological_charge_definition="unknown",
    )

    dictionary_measurement = prepare_measurement_from_dict(U, Dict(
        "methodname" => "Topological_charge",
        "kinds_of_topological_charge" => ["clover"],
        "improved_topological_charge_definition" => "both",
        "printvalues" => false,
    ))
    @test dictionary_measurement.improved_topological_charge_definition == "both"

    density_measurement =
        QCDMeasurements.Topological_charge_density_correlation_measurement(
            U;
            TC_methods=["clover"],
            improved_topological_charge_definition="both",
        )
    density_values = get_value(measure(
        density_measurement,
        U,
        [1, 1, 1, 1],
        [0, 0, 0, 0],
    ))
    @test haskey(density_values, "clover improved")
    @test haskey(density_values, "clover improved bilson-thompson")
end
