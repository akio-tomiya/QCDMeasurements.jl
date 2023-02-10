# QCDMeasurements.jl
Measurements for lattice QCD. 

This is intended to use in [LatticeQCD.jl](https://github.com/akio-tomiya/LatticeQCD.jl).


# Sample

```julia
    NX = 4
    NY = 4
    NZ = 4
    NT = 4
    Nwing = 0
    Dim = 4
    NC = 3

    U = Initialize_4DGaugefields(NC,Nwing,NX,NY,NZ,NT,condition = "cold")
    
    filename = "./conf_00000008.ildg"
    ildg = ILDG(filename)
    i = 1
    L = [NX,NY,NZ,NT]
    load_gaugefield!(U,i,ildg,L,NC)

    m_plaq = Plaquette_measurement(U)
    m_poly = Polyakov_measurement(U)

    plaq = get_value(measure(m_plaq,U))
    poly = get_value(measure(m_poly,U))
    println("plaq: $plaq")
    println("poly: $poly")
```