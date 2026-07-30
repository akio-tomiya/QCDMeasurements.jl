using QCDMeasurements
using Test

@testset "QCDMeasurements.jl" begin
    # Write your tests here.

    include("dicttest.jl")
    include("pion_correlator.jl")
    include("gauge.jl")

end
