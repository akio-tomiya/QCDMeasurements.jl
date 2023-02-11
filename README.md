# QCDMeasurements.jl
Measurements for lattice QCD. 

Lattice QCD is a well-established non-perturbative approach to solving the quantum chromodynamics (QCD) theory of quarks and gluons. 
Gauge field for gluons is treated by [Gaugefields.jl](https://github.com/akio-tomiya/Gaugefields.jl).
Pseudo fermion field for quarks is treated by [LatticeDiracOperators.jl](https://github.com/akio-tomiya/LatticeDiracOperators.jl).
It is important to measure physical observables from gauge fields. 
QCDMeasurements.jl is now the external package for measurements in Lattice QCD. 

<img src="LQCDjl_block.png" width=300> 


This is intended to use in [LatticeQCD.jl](https://github.com/akio-tomiya/LatticeQCD.jl).

# What this package can do:
This package has following functionarities

- Plaquette measurement.
- Poylakov loop measurement.
- Pion correlator measurement.
- Chiral condensate measurement.
- Topological charge measurement.
- Energy density measurement.


# Sample

```julia
using Gaugefields
function test()
    println("SU3test")
    NX = 4
    NY = 4
    NZ = 4
    NT = 4
    Nwing = 0
    Dim = 4
    NC = 3

    U = Initialize_4DGaugefields(NC,Nwing,NX,NY,NZ,NT,condition = "cold")
    #U = Initialize_Gaugefields(NC,Nwing,NX,NY,NZ,NT,condition = "hot",randomnumber="Reproducible")
    filename = "testconf.txt"
    L = [NX,NY,NZ,NT]
    load_BridgeText!(filename,U,L,NC)
    #=
    filename = "./conf_00000008.ildg"
    ildg = ILDG(filename)
    i = 1
    L = [NX,NY,NZ,NT]
    load_gaugefield!(U,i,ildg,L,NC)
    =#

    m_plaq = Plaquette_measurement(U)
    m_poly = Polyakov_measurement(U)

    plaq = get_value(measure(m_plaq,U))
    poly = get_value(measure(m_poly,U))
    println("plaq: $plaq")
    println("poly: $poly")

    m_energy = Energy_density_measurement(U)
    m_topo = Topological_charge_measurement(U)
    energy = get_value(measure(m_energy,U))
    topo = get_value(measure(m_topo,U))
    println("energy: $energy")
    println("topo: $topo")

    m_pion = Pion_correlator_measurement(U)
    m_pion_Staggered = Pion_correlator_measurement(U,fermiontype = "Staggered")
    m_pion_Wilson = Pion_correlator_measurement(U,fermiontype = "Wilson")
    pion = get_value(measure(m_pion,U))
    pion_s = get_value(measure(m_pion_Staggered,U))
    pion_w = get_value(measure(m_pion_Wilson,U))

    println("pion: $pion")
    println("pion correlator with Staggered fermion: $pion_s")
    println("pion correlator with  Wilson fermion: $pion_w")

    m_chiral_Staggered = Chiral_condensate_measurement(U,fermiontype = "Staggered")
    m_chiral_Wilson = Chiral_condensate_measurement(U,fermiontype = "Wilson")
    chiral_s = get_value(measure(m_chiral_Staggered,U))
    chiral_w = get_value(measure(m_chiral_Wilson,U))

    println("Chiral condensate with Staggered fermion: $chiral_s")
    println("Chiral condensatewith  Wilson fermion: $chiral_w")
end
test()
```
