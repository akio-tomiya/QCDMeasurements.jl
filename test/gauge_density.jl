include("gauge_density_helpers.jl")

@testset "Gauge density fields" begin
    @test QCDMeasurements.topological_charge_density !== Gaugefields.topological_charge_density
    number_of_processes = MPI.Comm_size(MPI.COMM_WORLD)
    if number_of_processes == 1
        @testset "Unit links" begin
            for backend in (LatticeMatricesBackend(), LegacyBackend())
                U = gauge_configuration((2, 3, 4, 5); colors=3, backend, start=:cold, verbose=0)
                for (method, _) in DENSITY_METHOD_KEYS
                    q = QCDMeasurements.topological_charge_density(U; method)
                    @test q isa Array{Float64,4}
                    @test size(q) == (2, 3, 4, 5)
                    @test all(iszero, q)
                end
                for method in ("plaquette", "clover")
                    @test all(iszero, action_density(U; method))
                end
            end
        end

        backend_reference = Dict{String,Array{Float64,4}}()
        @testset "Noncubic non-Abelian fields, $backend, Nc=$colors" for
            (backend, colors) in ((LatticeMatricesBackend(), 3), (LegacyBackend(), 3),
                                 (LatticeMatricesBackend(), 2))
            L = (3, 4, 5, 6)
            U = density_test_field(L; backend, colors)
            original = density_test_snapshot(U)
            workspace = GaugeDensityWorkspace(U)
            q = fill(NaN, L)
            measurement = TopologicalChargeMeasurement(U;
                TC_methods=["plaquette", "clover"],
                improved_topological_charge_definition="both")
            reference = get_value(measure(measurement, U))
            local_measurement = TopologicalChargeDensityCorrelationMeasurement(U)
            positions = ([1, 1, 1, 1], [3, 4, 5, 6], [2, 1, 3, 4],
                         [1, 2, 4, 3], [3, 3, 1, 5], [2, 4, 2, 6])

            for (method, key) in DENSITY_METHOD_KEYS
                @test topological_charge_density!(q, U, workspace; method) === q
                @test all(isfinite, q)
                @test maximum(abs, q) > 1e-7
                @test sum(q) ≈ reference[key] atol=2e-12 rtol=2e-12
                if colors == 3
                    if backend isa LatticeMatricesBackend
                        backend_reference[method] = copy(q)
                    else
                        @test q ≈ backend_reference[method] atol=2e-13 rtol=2e-12
                    end
                end
                allocated = QCDMeasurements.topological_charge_density(U; method=Symbol(method))
                @test q ≈ allocated atol=2e-14 rtol=2e-13
                for position in positions
                    @test q[position...] ≈ density_point_reference(
                        U, position, method, local_measurement) atol=2e-13 rtol=2e-12
                end
            end

            history = get_value(measure(GradientFlowScaleMeasurement(U;
                number_of_flow_steps=0, measure_topological_charge=false), U))
            for method in ("plaquette", "clover")
                @test action_density!(q, U, workspace; method) === q
                @test minimum(q) >= -2e-13
                @test sum(q) / length(q) ≈ history.energy_density[method][1] atol=2e-12 rtol=2e-12
                @test q ≈ action_density(U; method=Symbol(method)) atol=2e-14 rtol=2e-13
                if colors == 3
                    if backend isa LatticeMatricesBackend
                        backend_reference["action_" * method] = copy(q)
                    else
                        @test q ≈ backend_reference["action_" * method] atol=2e-12 rtol=2e-12
                    end
                end
                for position in positions
                    @test q[position...] ≈ action_point_reference(
                        U, position, method, local_measurement) atol=2e-12 rtol=2e-12
                end
            end
            @test density_test_snapshot(U) == original

            # Reuse after both raw-loop and improved-tensor calculations, on
            # changed links: a workspace must not cache a previous configuration.
            Ucold = gauge_configuration(L; colors, backend, start=:cold, verbose=0)
            @test all(iszero, topological_charge_density!(q, Ucold, workspace))
            @test all(iszero, action_density!(q, Ucold, workspace))
            @test_throws ArgumentError topological_charge_density!(q, U, workspace; method="unknown")
            @test_throws ArgumentError action_density!(q, U, workspace; method="alexandrou")
            @test_throws DimensionMismatch topological_charge_density!(zeros(2, 2, 2, 2), U, workspace)
            @test_throws DimensionMismatch action_density!(zeros(2, 2, 2, 2), U, workspace)
            @test_throws DimensionMismatch topological_charge_density!(q, U, GaugeDensityWorkspace(
                gauge_configuration((2, 2, 2, 2); colors, backend, verbose=0)))
            @test_throws ArgumentError GaugeDensityWorkspace(U[1:3])
            @test_throws ArgumentError QCDMeasurements.topological_charge_density(U; method="improved")
            @test_throws ArgumentError action_density(U; method="unknown")
        end

        @testset "Backend and geometry validation" begin
            U = gauge_configuration((3, 4, 5, 6); colors=3, verbose=0)
            workspace = GaugeDensityWorkspace(U)
            q = zeros(3, 4, 5, 6)
            legacy = gauge_configuration((3, 4, 5, 6); colors=3, backend=LegacyBackend(), verbose=0)
            @test_throws ArgumentError topological_charge_density!(q, legacy, workspace)
            U2d = gauge_configuration((4, 4); colors=3, backend=LegacyBackend(), verbose=0)
            @test_throws ArgumentError QCDMeasurements.topological_charge_density(U2d)
            @test_throws ArgumentError action_density(U2d)
            # Legacy types encode colors but not lattice dimensions in their type.
            legacy[2] = gauge_configuration((2, 2, 2, 2); colors=3, backend=LegacyBackend(), verbose=0)[2]
            @test_throws DimensionMismatch GaugeDensityWorkspace(legacy)
        end
    end
end
