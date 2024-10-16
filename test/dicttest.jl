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
    methodname = "Topological_charge_density_correlation"
    method["methodname"] = methodname
    method["kinds_of_topological_charge"] = ["plaquette", "clover"]

    m = prepare_measurement_from_dict(U, method)
    println("$methodname")
    for ix = 0:NX
        value = get_value(measure(m, U, [1, 1, 1, 1], [ix, 0, 0, 0]))
        println("$ix $(value["plaquette"]) $(value["clover"]) $(value["clover improved"])")
    end



    method = Dict()
    methodname = "Guluonic_correlators"
    method["methodname"] = methodname
    loops1 = []
    loop1 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops1, loop1)
    loop1 = [(1, +1), (3, +1), (1, -1), (3, -1)]
    push!(loops1, loop1)
    method["loop1"] = loops1
    loops2 = []
    loop2 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops2, loop2)
    method["loop2"] = loops2

    m = prepare_measurement_from_dict(U, method)
    println("$methodname")
    for ix = 0:NX-1
        value = get_value(measure(m, U, [1, 1, 1, 1], [ix, 0, 0, 0]))
        println("$ix $value")
    end

    method = Dict()
    methodname = "Correlation"
    method["methodname"] = methodname
    loops1 = []
    loop1 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops1, loop1)
    loop1 = [(1, +1), (3, +1), (1, -1), (3, -1)]
    push!(loops1, loop1)
    method["loop1"] = loops1
    loops2 = []
    loop2 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops2, loop2)
    method["loop2"] = loops2
    method["relativeposition"] = [0, 0, 0, 2]

    m = prepare_measurement_from_dict(U, method)
    value = get_value(measure(m, U))
    println("$methodname $value")

    method = Dict()
    methodname = "Correlation"
    method["methodname"] = methodname
    loops1 = []
    loop1 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops1, loop1)
    method["loop1"] = loops1
    loops2 = []
    loop2 = [(1, +1), (2, +1), (1, -1), (2, -1)]
    push!(loops2, loop2)
    method["loop2"] = loops2
    method["relativeposition"] = [0, 0, 0, 2]
    method["loop1position"] = [1, 1, 2, 2]

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