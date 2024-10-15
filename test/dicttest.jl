using Gaugefields
using Wilsonloop
function SU3test()
    println("SU3test")
    NX = 4
    NY = 4
    NZ = 4
    NT = 4
    Nwing = 0
    Dim = 4
    NC = 3

    U = Initialize_4DGaugefields(NC, Nwing, NX, NY, NZ, NT, condition="cold")
    #U = Initialize_Gaugefields(NC,Nwing,NX,NY,NZ,NT,condition = "hot",randomnumber="Reproducible")
    filename = "testconf.txt"
    L = [NX, NY, NZ, NT]
    load_BridgeText!(filename, U, L, NC)

    #=
    filename = "./conf_00000008.ildg"
    ildg = ILDG(filename)
    i = 1
    L = [NX,NY,NZ,NT]
    load_gaugefield!(U,i,ildg,L,NC)
    =#

    method = Dict()
    methodname = "Correlation"
    method["methodname"] = methodname
    loop1 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    method["loop1"] = loop1
    loop2 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    method["loop2"] = loop2
    method["position"] = [0, 0, 0, 2]

    m = prepare_measurement_from_dict(U, method)
    value = get_value(measure(m, U))
    println("$methodname $value")


    method = Dict()
    methodname = "Eigenvalue"
    method["methodname"] = methodname
    method["fermiontype"] = "Wilson"
    κ = 0.141139
    method["hop"] = κ
    method["nev"] = 1 #number of eigenvalues
    m = prepare_measurement_from_dict(U, method)
    value, vectors = get_value(measure(m, U)) #eigenvalues and eigenvectors
    println("$methodname $value")


    method = Dict()
    methodname = "Pion_correlator"
    method["methodname"] = methodname
    method["fermiontype"] = "Staggered"
    method["mass"] = 1
    method["Nf"] = 4
    m = prepare_measurement_from_dict(U, method)
    value = get_value(measure(m, U))
    println("$methodname $value")

    method = Dict()
    methodname = "Pion_correlator"
    method["methodname"] = methodname
    method["fermiontype"] = "Wilson"
    method["hop"] = 1
    m = prepare_measurement_from_dict(U, method)
    value = get_value(measure(m, U))
    println("$methodname $value")


    methodsname = ["Plaquette", "Polyakov_loop", "Topological_charge", "Chiral_condensate",
        "Pion_correlator", "Energy_density", "Wilson_loop", "Eigenvalue"]
    method = Dict()
    for methodname in methodsname
        method["methodname"] = methodname
        m = prepare_measurement_from_dict(U, method)
        value = get_value(measure(m, U))
        if methodname == "Eigenvalue"
            println("$methodname $(value[1])")
        else
            println("$methodname $(value)")
        end
    end

end
SU3test()