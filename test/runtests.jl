using QCDMeasurements
using Test

@testset "QCDMeasurements.jl" begin
    include("pion_solver_diagnostics.jl")
    include("dicttest.jl")
    include("pion_correlator.jl")
    include("gauge.jl")

end
