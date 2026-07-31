using Gaugefields

U_solver_diagnostics = Initialize_Gaugefields(
    3,
    0,
    2,
    2,
    2,
    2;
    condition="cold",
)
measurement_solver_diagnostics = Pion_correlator_measurement(
    U_solver_diagnostics;
    fermiontype="Wilson",
    κ=0.10,
    r=1.0,
    eps_CG=1.0e-20,
    MaxCGstep=10_000,
    BoundaryCondition=[1, 1, 1, -1],
    method_CG="bicgstab",
    verbose_level=0,
    printvalues=false,
)
output_solver_diagnostics =
    get_value(measure(measurement_solver_diagnostics, U_solver_diagnostics))
diagnostics = get_solver_diagnostics(measurement_solver_diagnostics)
sorted_iterations =
    sort([diagnostic.iterations for diagnostic in diagnostics])
median_iterations =
    (sorted_iterations[6] + sorted_iterations[7]) / 2

@test length(output_solver_diagnostics) == 2
@test length(diagnostics) == 12
@test [diagnostic.source_number for diagnostic in diagnostics] == 1:12
@test all(diagnostic -> diagnostic.method === :bicgstab, diagnostics)
@test all(
    diagnostic -> 0 < diagnostic.iterations <
                  diagnostic.maximum_iterations,
    diagnostics,
)
@test all(
    diagnostic ->
        diagnostic.recursive_residual_squared <
        diagnostic.target_residual_squared,
    diagnostics,
)
@test maximum(
    diagnostic.true_relative_residual for diagnostic in diagnostics
) < 1.0e-10

@test_throws ArgumentError Pion_correlator_measurement(
    U_solver_diagnostics;
    fermiontype="Wilson",
    method_CG="unsupported",
)

parameter_measurement_solver_diagnostics = Pion_correlator_measurement(
    U_solver_diagnostics,
    QCDMeasurements.Pion_parameters(
        eps=1.0e-20,
        MaxCGstep=10_000,
        method_CG="bicgstab",
        fermion_parameters=QCDMeasurements.Wilson_parameters(hop=0.10),
        verbose_level=0,
        printvalues=false,
    ),
)
@test parameter_measurement_solver_diagnostics.D.method_CG == "bicgstab"

failed_measurement_solver_diagnostics = Pion_correlator_measurement(
    U_solver_diagnostics;
    fermiontype="Wilson",
    κ=0.10,
    eps_CG=0.0,
    MaxCGstep=0,
    method_CG="bicgstab",
    verbose_level=0,
    printvalues=false,
)
@test_throws ErrorException measure(
    failed_measurement_solver_diagnostics,
    U_solver_diagnostics,
)
@test isempty(
    get_solver_diagnostics(failed_measurement_solver_diagnostics),
)

@info "pion solver diagnostic metrics" minimum_iterations = minimum(
    diagnostic.iterations for diagnostic in diagnostics
) median_iterations maximum_iterations = maximum(
    diagnostic.iterations for diagnostic in diagnostics
) maximum_true_relative_residual = maximum(
    diagnostic.true_relative_residual for diagnostic in diagnostics
)
