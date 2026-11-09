using Gaugefields
using LatticeDiracOperators
using LatticeMatrices

@testset "Supported dependency versions" begin
    gaugefields_version = pkgversion(Gaugefields)
    dirac_operators_version = pkgversion(LatticeDiracOperators)
    lattice_matrices_version = pkgversion(LatticeMatrices)

    @test v"1.1.12" <= gaugefields_version < v"2.0.0"
    @test v"1.1.5" <= dirac_operators_version < v"2.0.0"
    @test v"1.2.11" <= lattice_matrices_version < v"2.0.0"
end
