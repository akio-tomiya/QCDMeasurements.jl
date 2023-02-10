using Gaugefields
function SU3test()
    println("SU3test")
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
    m_pion = Pion_correlator_measurement(U)
    m_pion_Staggered = Pion_correlator_measurement(U,fermiontype = "Staggered")
    m_pion_Wilson = Pion_correlator_measurement(U,fermiontype = "Wilson")

    plaq = get_value(measure(m_plaq,U))
    poly = get_value(measure(m_poly,U))
    println("plaq: $plaq")
    println("poly: $poly")


    pion = get_value(measure(m_pion,U))
    pion_s = get_value(measure(m_pion_Staggered,U))
    pion_w = get_value(measure(m_pion_Wilson,U))

    println("pion: $pion")
    println("pion correlator with Staggered fermion: $pion_s")
    println("pion correlator with  Wilson fermion: $pion_w")

end
SU3test()