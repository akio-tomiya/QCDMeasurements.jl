using JACC
JACC.@init_backend
using MPI
MPI.Initialized() || MPI.Init()
using LatticeMatrices
using QCDMeasurements
using Test

@testset "QCDMeasurements.jl" begin
    include("dependency_compatibility.jl")
    include("public_api.jl")
    include("basic_measurements.jl")
    include("legacy_measurements.jl")
    include("gradient_flow_scale.jl")
    include("simulateqcd_gradient_flow_reference.jl")
    include("meson_correlator_kernels.jl")
    include("simulateqcd_meson_reference.jl")
    include("grid_meson_reference.jl")
    include("pcac_mass.jl")
    include("domainwall_residual_mass.jl")
    include("wilson_clover_reference.jl")
    include("pion_correlator.jl")
    include("pion_solver_diagnostics.jl")
end

# `dicttest.jl` and `gauge.jl` are retained as legacy integration examples.
# They are intentionally excluded here because they run long scans, write into
# the source tree, and previously completed without making any test assertions.
