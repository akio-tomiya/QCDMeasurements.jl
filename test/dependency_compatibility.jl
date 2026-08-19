using Gaugefields
using LatticeDiracOperators
using LatticeMatrices

@testset "Supported dependency versions" begin
    gaugefields_version = pkgversion(Gaugefields)
    dirac_operators_version = pkgversion(LatticeDiracOperators)
    lattice_matrices_version = pkgversion(LatticeMatrices)

    @test v"1.0.3" <= gaugefields_version < v"2.0.0"
    @test v"1.0.0" <= dirac_operators_version < v"2.0.0"
    @test v"1.1.2" <= lattice_matrices_version < v"2.0.0"
end
