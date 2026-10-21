using Gaugefields, LatticeMatrices, LinearAlgebra, Random, SHA, TOML

# Construct exactly the same SU(N) links used for the pre-change reference.
# Host construction is independent of process decomposition and device backend.
function density_reference_links(L, colors)
    links = Array{ComplexF64}(undef, colors, colors, L..., 4)
    rng = MersenneTwister(210426)
    for site in CartesianIndices(L), μ in 1:4
        A = randn(rng, ComplexF64, colors, colors)
        H = A + A'
        H -= tr(H) / colors * I
        links[:, :, Tuple(site)..., μ] = exp(0.25im * H)
    end
    return links
end

function density_backend_field(links, grid; halo=1, legacy=false)
    L = size(links)[3:6]
    if legacy
        U = Initialize_Gaugefields(size(links, 1), 0, L...;
            mpi=true, mpiinit=false, PEs=grid, condition="cold", verbose_level=0)
        for μ in 1:4
            ranges = ntuple(d ->
                (U[μ].myrank_xyzt[d] * U[μ].PN[d] + 1):
                ((U[μ].myrank_xyzt[d] + 1) * U[μ].PN[d]), 4)
            copyto!(U[μ].U, links[:, :, ranges..., μ])
        end
        return U
    end
    U = gauge_configuration(L; colors=size(links, 1), process_grid=grid, halo,
        start=:cold, verbose=0)
    for μ in 1:4
        field = U[μ].U
        global_ranges = ntuple(d ->
            (field.coords[d] * field.PN[d] + 1):((field.coords[d] + 1) * field.PN[d]), 4)
        local_ranges = ntuple(d -> (halo + 1):(halo + field.PN[d]), 4)
        block = convert(JACC.array_type(), links[:, :, global_ranges..., μ])
        copyto!(view(field.A, :, :, local_ranges...), block)
        LatticeMatrices.mark_halo_dirty!(field)
    end
    Gaugefields.set_wing_U!(U)
    return U
end

@testset "Density values across JACC and MPI backends" begin
    fixture = TOML.parsefile(joinpath(@__DIR__, "references", "gauge_density_v1.toml"))
    L = Tuple(fixture["lattice"])
    nprocs = MPI.Comm_size(MPI.COMM_WORLD)
    grids = unique([(1, nprocs, 1, 1), (1, 1, 1, nprocs)])
    nprocs == 4 && push!(grids, (1, 2, 1, 2))
    grids = filter(grid -> all(L[d] % grid[d] == 0 && L[d] > grid[d] for d in 1:4), grids)
    @test !isempty(grids)
    layouts = nprocs > 1 ? (false, true) : (false,)
    for colors in (2, 3), grid in grids, legacy in layouts
        label = legacy && colors == 3 ? "legacy_su3" : "lm_su$colors"
        reference = fixture["cases"][label]
        links = density_reference_links(L, colors)
        digest = bytes2hex(sha256(reinterpret(UInt8,
            vec(permutedims(links, (7, 1, 2, 3, 4, 5, 6))))))
        @test digest == reference["input_sha256"]
        U = density_backend_field(links, grid; legacy)
        arrays = legacy ? [link.U for link in U] : [link.U.A for link in U]
        @test legacy ? arrays[1] isa Array : arrays[1] isa JACC.array_type()
        original = [Array(copy(data)) for data in arrays]
        workspace = GaugeDensityWorkspace(U)
        @test workspace.result.A isa JACC.array_type()
        local_size = legacy ? U[1].PN : U[1].U.PN
        @test workspace.result.PN == local_size
        @test size(workspace.result.A)[3:end] == local_size
        output = fill(NaN, L)
        scalar = get_value(measure(TopologicalChargeMeasurement(U;
            TC_methods=["plaquette", "clover"], improved_topological_charge_definition="both"), U))
        for (method, key) in ("plaquette" => "plaquette", "clover" => "clover",
                             "alexandrou" => "clover improved",
                             "bilson_thompson" => "clover improved bilson-thompson")
            expected = reshape(reference["topological_fields"][method], L)
            @test topological_charge_density!(output, U, workspace; method) === output
            @test output isa Array{Float64,4}
            @test size(output) == L
            @test all(isapprox.(output, expected; atol=2e-13, rtol=2e-12))
            @test sum(output) ≈ reference["charges"][key] atol=2e-12 rtol=2e-12
            @test scalar[key] ≈ reference["charges"][key] atol=2e-12 rtol=2e-12
            @test sum(output) ≈ scalar[key] atol=2e-12 rtol=2e-12
            allocated = QCDMeasurements.topological_charge_density(U; method=Symbol(method))
            @test all(isapprox.(allocated, expected; atol=2e-13, rtol=2e-12))
        end
        history = get_value(measure(GradientFlowScaleMeasurement(U;
            number_of_flow_steps=0, measure_topological_charge=false), U))
        for method in ("plaquette", "clover")
            expected = reshape(reference["action_fields"][method], L)
            @test action_density!(output, U, workspace; method) === output
            @test all(isapprox.(output, expected; atol=2e-12, rtol=2e-12))
            @test sum(output)/length(output) ≈ reference["actions"][method] atol=2e-12 rtol=2e-12
            @test history.energy_density[method][1] ≈ reference["actions"][method] atol=2e-12 rtol=2e-12
            @test all(isapprox.(action_density(U; method), expected; atol=2e-12, rtol=2e-12))
        end
        @test get_value(measure(EnergyDensityMeasurement(U), U)) ≈ reference["legacy_energy"] atol=2e-12 rtol=2e-12
        @test all(Array(arrays[μ]) == original[μ] for μ in 1:4)
        # Reuse after a new configuration; this also checks collective output
        # assembly and all methods on unit links on every backend and rank.
        Ucold = legacy ? Initialize_Gaugefields(colors, 0, L...;
            mpi=true, mpiinit=false, PEs=grid, condition="cold", verbose_level=0) :
            gauge_configuration(L; colors, process_grid=grid, start=:cold, verbose=0)
        for method in ("plaquette", "clover", "alexandrou", "bilson_thompson")
            @test all(iszero, topological_charge_density!(output, Ucold, workspace; method))
        end
        for method in ("plaquette", "clover")
            @test all(iszero, action_density!(output, Ucold, workspace; method))
        end
    end
end
