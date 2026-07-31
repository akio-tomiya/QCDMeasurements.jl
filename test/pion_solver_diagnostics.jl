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
@test measurement_solver_diagnostics.S === nothing
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

preconditioned_measurement_solver_diagnostics =
    Pion_correlator_measurement(
        U_solver_diagnostics;
        fermiontype="Wilson",
        κ=0.10,
        method_CG="preconditiond_bicgstab",
        verbose_level=0,
        printvalues=false,
    )
@test preconditioned_measurement_solver_diagnostics.D.method_CG ==
      "preconditiond_bicgstab"
@test !occursin(
    "faster",
    string(typeof(preconditioned_measurement_solver_diagnostics.D)),
)

U_hot_solver_crosscheck = Initialize_Gaugefields(
    3,
    0,
    2,
    2,
    2,
    2;
    condition="hot",
    randomnumber="Reproducible",
)
normal_hot_measurement = Pion_correlator_measurement(
    U_hot_solver_crosscheck;
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
evenodd_hot_measurement = Pion_correlator_measurement(
    U_hot_solver_crosscheck;
    fermiontype="Wilson",
    κ=0.10,
    r=1.0,
    eps_CG=1.0e-20,
    MaxCGstep=10_000,
    BoundaryCondition=[1, 1, 1, -1],
    method_CG="preconditiond_bicgstab",
    verbose_level=0,
    printvalues=false,
)
normal_hot_correlator =
    get_value(measure(normal_hot_measurement, U_hot_solver_crosscheck))
evenodd_hot_correlator =
    get_value(measure(evenodd_hot_measurement, U_hot_solver_crosscheck))
normal_hot_diagnostics = get_solver_diagnostics(normal_hot_measurement)
evenodd_hot_diagnostics = get_solver_diagnostics(evenodd_hot_measurement)

@test length(normal_hot_diagnostics) == 12
@test length(evenodd_hot_diagnostics) == 12
@test all(
    diagnostic ->
        diagnostic.method === :bicgstab &&
        0 < diagnostic.iterations < diagnostic.maximum_iterations &&
        isfinite(diagnostic.recursive_residual_squared) &&
        diagnostic.recursive_residual_squared <
        diagnostic.target_residual_squared &&
        isfinite(diagnostic.true_relative_residual) &&
        diagnostic.true_relative_residual < 1.0e-10,
    normal_hot_diagnostics,
)
@test all(
    diagnostic ->
        diagnostic.method === :preconditiond_bicgstab &&
        0 < diagnostic.iterations < diagnostic.maximum_iterations &&
        isfinite(diagnostic.recursive_residual_squared) &&
        diagnostic.recursive_residual_squared <
        diagnostic.target_residual_squared &&
        isfinite(diagnostic.true_relative_residual) &&
        diagnostic.true_relative_residual < 1.0e-10,
    evenodd_hot_diagnostics,
)
hot_correlator_relative_difference = sqrt(
    sum(abs2, evenodd_hot_correlator - normal_hot_correlator) /
    sum(abs2, normal_hot_correlator),
)
@test isfinite(hot_correlator_relative_difference)
@test hot_correlator_relative_difference < 1.0e-8

normal_hot_iterations =
    [diagnostic.iterations for diagnostic in normal_hot_diagnostics]
evenodd_hot_iterations =
    [diagnostic.iterations for diagnostic in evenodd_hot_diagnostics]
normal_hot_maximum_true_relative_residual =
    maximum(
        diagnostic.true_relative_residual
        for diagnostic in normal_hot_diagnostics
    )
evenodd_hot_maximum_true_relative_residual =
    maximum(
        diagnostic.true_relative_residual
        for diagnostic in evenodd_hot_diagnostics
    )
@info "hot pion solver cross-check" hot_correlator_relative_difference normal_minimum_iterations =
    minimum(normal_hot_iterations) normal_maximum_iterations =
    maximum(normal_hot_iterations) evenodd_minimum_iterations =
    minimum(evenodd_hot_iterations) evenodd_maximum_iterations =
    maximum(evenodd_hot_iterations) normal_hot_maximum_true_relative_residual evenodd_hot_maximum_true_relative_residual

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
