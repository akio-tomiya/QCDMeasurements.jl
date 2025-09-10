

mutable struct Pion_correlator_measurement{Dim,TG,TD,TF,TF_vec,Dim_2,TCov} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    Dim::Int8
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    fermi_action::TF
    _temporary_fermionfields::Vector{TF_vec}
    #Nr::Int64
    Nspinor::Int64
    S::Array{ComplexF64,Dim_2}
    cov_neural_net::TCov#Union{Nothing,CovNeuralnet}

    function Pion_correlator_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        fermiontype="Staggered",
        mass=0.1,
        Nf=2,
        κ=1,
        r=1,
        L5=2,
        M=-1,
        eps_CG=1e-14,
        MaxCGstep=5000,
        BoundaryCondition=nothing,
        cov_neural_net=nothing
    ) where {T}
        NC = U[1].NC

        Dim = length(U)
        if BoundaryCondition == nothing
            if Dim == 4
                boundarycondition = BoundaryCondition_4D_default
            elseif Dim == 2
                boundarycondition = BoundaryCondition_2D_default
            end
        else
            boundarycondition = BoundaryCondition
        end
        #println(boundarycondition)
        params, parameters_action, x, factor = make_fermionparameter_dict(U,
            fermiontype, mass,
            Nf,
            κ,
            r,
            L5,
            M,
        )        #=
        params = Dict()
        parameters_action = Dict()
        if fermiontype == "Staggered"
            x = Initialize_pseudofermion_fields(U[1], "staggered")
            params["Dirac_operator"] = "staggered"
            params["mass"] = mass
            parameters_action["Nf"] = Nf
            Nfbase = 4
            #Nfbase = ifelse( m.fparam.Dirac_operator == "Staggered",4,1)
            factor = Nf / Nfbase

        elseif fermiontype == "Wilson"
            x = Initialize_pseudofermion_fields(U[1], "Wilson", nowing = true)
            params["Dirac_operator"] = "Wilson"
            params["κ"] = κ
            params["r"] = r
            params["faster version"] = true
        elseif fermiontype == "Domainwall"
            params["Dirac_operator"] = "Domainwall"
            params["mass"] = mass
            params["L5"] = L5
            params["M"] = M
            x = Initialize_pseudofermion_fields(U[1], "Domainwall", L5 = L5)
        else
            error(
                "fermion type $fermiontype is not supported in chiral condensate measurement",
            )
        end
        =#

        Nspinor = ifelse(fermiontype == "Staggered", 1, 4)

        _, _, NN... = size(U[1])
        S = zeros(ComplexF64, NN..., Nspinor * NC, Nspinor * NC)


        params["eps_CG"] = eps_CG
        params["verbose_level"] = verbose_level
        params["MaxCGstep"] = MaxCGstep
        params["boundarycondition"] = boundarycondition


        
        D = Dirac_operator(U, x, params)
        fermi_action = FermiAction(D, parameters_action)
        
        
        TD = typeof(D)
        TF = typeof(fermi_action)

        TCov = typeof(cov_neural_net)

        myrank = get_myrank(U)

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end



        numg = 1
        _temporary_gaugefields = Temporalfields(U[1], num=numg)
        #_temporary_gaugefields = Vector{T}(undef, numg)
        #_temporary_gaugefields[1] = similar(U[1])
        #for i = 2:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end

        numf = 2
        TF_vec = typeof(x)
        _temporary_fermionfields = Vector{TF_vec}(undef, numf)
        for i = 1:numf
            _temporary_fermionfields[i] = similar(x)
        end
        Dim_2 = Dim + 2

        return new{Dim,T,TD,TF,TF_vec,Dim_2,TCov}(
            filename, #::String
            _temporary_gaugefields,#::Vector{TG}
            Dim,#::Int8
            verbose_print,#::Union{Verbose_print,Nothing}
            printvalues,#::Bool
            D,#::TD
            fermi_action,#::TF
            _temporary_fermionfields,#::Vector{TF_vec}
            Nspinor,#::Int64
            S,#::Array{ComplexF64,3}
            cov_neural_net,
        )
    end

end

function Pion_correlator_measurement(
    U::Vector{T},
    params::Pion_parameters,
    filename="Pion_correlator.txt",
) where {T}

    #if params.smearing_for_fermion != "nothing"
    #    error("smearing is not implemented in Pion correlator")
    #end

    if params.smearing_for_fermion == "nothing"
        cov_neural_net = nothing
    elseif params.smearing_for_fermion == "stout"
        cov_neural_net = CovNeuralnet(U)
        if params.stout_numlayers == 1
            #st = STOUT_Layer(p.stout_loops, p.stout_ρ, L)
            st = STOUT_Layer(params.stout_loops, params.stout_ρ, U)
            push!(cov_neural_net, st)
        else
            if length(params.stout_ρ) == 1
                @warn "num. of stout layer is $(params.stout_numlayers) but there is only one rho. rho values are all same."
                for ilayer = 1:length(params.stout_ρ)
                    st = STOUT_Layer(params.stout_loops, params.stout_ρ, U)
                    push!(cov_neural_net, st)
                end
            else
                for ilayer = 1:length(params.stout_ρ)
                    st = STOUT_Layer(params.stout_loops, params.stout_ρ[ilayer], U)
                    push!(cov_neural_net, st)
                end
            end
        end
    else
        error("params.smearing_for_fermion = $(params.smearing_for_fermion) is not supported")
    end

    #println(params)

    params_tuple = fermionparameter_params(params)

    fermionparameters = params.fermion_parameters
    if params.fermiontype == "Staggered"
        method = Pion_correlator_measurement(
            U;
            filename=filename, 
            cov_neural_net=cov_neural_net,
            params_tuple...
            #=
            verbose_level = params.verbose_level,
            printvalues = params.printvalues,
            fermiontype = params.fermiontype,
            mass = fermionparameters.mass,
            Nf = fermionparameters.Nf,
            eps_CG = params.eps,
            MaxCGstep = params.MaxCGstep,
            =#
        )
    elseif params.fermiontype == "Wilson" || params.fermiontype == "WilsonClover"
        if fermionparameters.hasclover
            error("WilsonClover is not implemented in Pion measurement")
        end
        method = Pion_correlator_measurement(
            U;
            filename=filename,
            cov_neural_net=cov_neural_net,
            params_tuple...
            #=
            verbose_level = params.verbose_level,
            printvalues = params.printvalues,
            fermiontype = params.fermiontype,
            κ = fermionparameters.hop,
            r = fermionparameters.r,
            eps_CG = params.eps,
            MaxCGstep = params.MaxCGstep,
            =#
        )
    elseif params.fermiontype == "Domainwall"
        # error("Domainwall fermion is not implemented in Pion measurement!")
        @warn "Pion correlator for Domainwall fermion is under construction!!"
        method = Pion_correlator_measurement(
            U;
            filename=filename,
            cov_neural_net=cov_neural_net,
            params_tuple...
            #=
            verbose_level = params.verbose_level,
            printvalues = params.printvalues,
            fermiontype = params.fermiontype,
            mass=params.mass,
            L5 = fermionparameters.N5,
            M = fermionparameters.M,
            eps_CG = params.eps,
            MaxCGstep = params.MaxCGstep
            =#
        )
    else
        error("fermiontype = $(params.fermiontype) is not supported")
    end

    return method
end

@inline function spincolor(ic, is, NC)
    return ic - 1 + (is - 1) * NC + 1
end


function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Pion_correlator_measurement,NC,Dim}
    S = m.S
    S .= 0
    measurestring = ""
    st = "Hadron spectrum started"
    measurestring *= st * "\n"
    println_verbose_level2(U[1], st)
    Nspinor = m.Nspinor
    #D = m.D(U)
    # calculate quark propagators from a point source at he origin
    if m.cov_neural_net === nothing
        propagators, st = calc_quark_propagators_point_source(m, U)
    else
        Uout, Uout_multi, _ = calc_smearedU(U, m.cov_neural_net)
        println("smeared U is used in Pion measurement")
        propagators, st = calc_quark_propagators_point_source(m, Uout)
    end
    measurestring *= st * "\n"
    #=
    #println(propagators)
    for ic=1:NC
        for is=1:Nspinor
            icum = (ic-1)*Nspinor+is
            println("$icum ", dot(propagators[icum],propagators[icum]))
        end
    end
    #error("prop")
    =#

    _, _, NN... = size(U[1])
    #println("NN = $NN")




    #ctr = 0 # a counter
    for ic = 1:NC
        for is = 1:Nspinor
            icum = (ic - 1) * Nspinor + is

            propagator = propagators[icum]
            α0 = spincolor(ic, is, NC) # source(color-spinor) index
            # reconstruction
            if Dim == 4
                @inbounds for t = 1:NN[4]
                    for z = 1:NN[3]
                        for y = 1:NN[2]
                            for x = 1:NN[1]
                                for ic2 = 1:NC
                                    @inbounds @simd for is2 = 1:Nspinor # Nspinor is the number of spinor index in 4d.
                                        β = spincolor(ic2, is2, NC)
                                        S[x, y, z, t, α0, β] +=
                                            propagator[ic, x, y, z, t, is]
                                        #println( propagator[ic,x,y,z,t,is])
                                    end
                                end
                            end
                        end
                    end
                end
            else
                error("Dim = $Dim is not supported")
            end
            # end for the substitution


            #ctr+=1
        end
    end
    #println(sum(S))
    #error("prop")
    # contruction end.
    st = "Hadron spectrum: Reconstruction"
    measurestring *= st * "\n"
    println_verbose_level2(U[1], st)
    #println("Hadron spectrum: Reconstruction")
    Cpi = zeros(NN[end])
    #Cpi = zeros( univ.NT )
    # Construct Pion propagator 
    if Dim == 4
        @inbounds for t = 1:NN[4]
            tmp = 0.0 + 0.0im
            for z = 1:NN[3]
                for y = 1:NN[2]
                    for x = 1:NN[1]
                        for ic = 1:NC
                            for is = 1:Nspinor # Nspinor is the number of spinor index in 4d.
                                α = spincolor(ic, is, NC)
                                for ic2 = 1:NC
                                    for is2 = 1:Nspinor # Nspinor is the number of spinor index in 4d.
                                        β = spincolor(ic2, is2, NC)
                                        tmp += S[x, y, z, t, α, β] * S[x, y, z, t, α, β]'#inner product.
                                        # complex conjugate = g5 S g5.
                                    end
                                end
                                # complex conjugate = g5 S g5.
                            end
                        end
                    end
                end
            end
            # staggered Pion correlator relies on https://itp.uni-frankfurt.de/~philipsen/theses/breitenfelder_ba.pdf (3.33)
            # we adopt ignoreing the staggering factor. See detail above reference.
            ksfact = 1.0 #ifelse( meas.fparam.Dirac_operator == "Staggered" , (-1)^(t-1) * 64, 1)
            Cpi[t] = real(tmp) * ksfact
        end
    end

    #println(typeof(verbose),"\t",verbose)
    st = "Hadron spectrum end"
    measurestring *= st * "\n"
    println_verbose_level2(U[1], st)
    #println("Hadron spectrum end")


    if m.printvalues
        stringcc = " "
        #println_verbose_level1(U[1],"$itrj ")
        #println_verbose_level1(m.verbose_print,"$itrj ")

        for it = 1:length(Cpi)
            cc = Cpi[it]
            #println_verbose_level1(U[1],"$cc ")
            stringcc *= "$cc "
        end
        #println_verbose_level1(U[1],stringcc)
        #measurestring *= stringcc * "\n"
        measurestring *= stringcc * " #pioncorrelator"
        #println_verbose_level1(m.verbose_print, stringcc)
        #println_verbose_level1(U[1],"#pioncorrelator")
        #st = "#pioncorrelator"
        #measurestring *= st * "\n"
        println_verbose_level1(m.verbose_print, stringcc * " #pioncorrelator")
        #println_verbose_level1(m.verbose_print, st)
    end

    output = Measurement_output(Cpi, measurestring)


    return output


end


function calc_quark_propagators_point_source(
    m,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
) where {NC,Dim}
    # D^{-1} for each spin x color element
    D = m.D(U)
    stvec = String[]
    propagators = map(
        i -> calc_quark_propagators_point_source_each(m, U, D, i, stvec),
        1:NC*m.Nspinor,
    )
    st = ""
    for i = 1:NC*m.Nspinor
        st *= stvec[i] * "\n"
    end
    return propagators, st
end

function calc_quark_propagators_point_source_each(m, U, D, i, stvec)
    # calculate D^{-1} for a given source at the origin.
    # Nc*Ns (Ns: dim of spinor, Wilson=4, ks=1) elements has to be gathered.
    # staggered Pion correlator relies on https://itp.uni-frankfurt.de/~philipsen/theses/breitenfelder_ba.pdf (3.33)
    temps_fermi = get_temporary_fermionfields(m)
    measurestring = ""
    b = temps_fermi[1]
    p = temps_fermi[2]
    #b = similar(meas._temporal_fermions[1]) # source is allocated
    #p = similar(b) # sink is allocated (propagator to sink position)
    #k = meas._temporal_fermi2[2]
    #clear_fermion!(b)
    Nspinor = m.Nspinor#ifelse( meas.fparam.Dirac_operator == "Staggered" ,1,4)
    is = ((i - 1) % Nspinor) + 1 # spin index   
    ic = ((i - is) ÷ Nspinor) + 1 # color index
    st = "$ic $is"
    measurestring *= st * "\n"
    println_verbose_level1(U[1], st)
    v = 1
    clear_fermion!(b)
    clear_fermion!(p)
    #b[ic,1,1,1,1,is] = v
    #println(dot(b,b))
    #println("ic = $ic is = $is")
    iorigin = (1, 1, 1, 1)
    setindex_global!(b, v, ic, iorigin..., is)  # source at the origin

    #=
    mul!(p,D,b)
    for it=1:U[1].NT
        for iz=1:U[1].NZ
            for iy=1:U[1].NY
                for ix=1:U[1].NX
                    val = p[ic,ix,iy,iz,it,is]
                    if abs(val) > 1e-16
                        println("$ix $iy $iz $it $val")
                    end
                end
            end

        end
    end
    =#

    #println(p[ic,1,1,1,1,is])
    #println(p[ic,2,1,1,1,is])
    #Z4_distribution_fermi!(b)
    #error("dd")

    @time solve_DinvX!(p, D, b)

    #error("dd")
    #println("norm p ",dot(p,p))
    st = "Hadron spectrum: Inversion $(i)/$(U[1].NC*m.Nspinor) is done"
    measurestring *= st * "\n"
    println_verbose_level1(U[1], st)

    flush(stdout)
    push!(stvec, measurestring)

    return deepcopy(p)
end

# for Domain wall, under construction
function calc_quark_propagators_point_source(
    m::Pion_correlator_measurement{Dim,TG,TD,TF,TF_vec,Dim_2,TCov},
    U::Array{<:AbstractGaugefields{NC,Dim},1},
) where {NC,Dim,TG,TD,TF,TF_vec<:LatticeDiracOperators.Dirac_operators.Abstract_DomainwallFermion_5D,Dim_2,TCov}
    # D^{-1} for each spin x color element
    D = m.D.D5DW(U)
    stvec = String[]
    propagators = map(
        i -> calc_quark_propagators_point_source_each(m, U, D, i, stvec),
        1:NC*m.Nspinor,
    )
    st = ""
    for i = 1:NC*m.Nspinor
        st *= stvec[i] * "\n"
    end
    return propagators, st
end

function calc_quark_propagators_point_source_each(
    m::Pion_correlator_measurement{Dim,TG,TD,TF,TF_vec,Dim_2,TCov},
    U,D,i, stvec
) where {Dim,TG,TD,TF,TF_vec<:LatticeDiracOperators.Dirac_operators.Abstract_DomainwallFermion_5D,Dim_2,TCov}
    temps_fermi = get_temporary_fermionfields(m)
    measurestring = ""
    b = temps_fermi[1]
    p = temps_fermi[2]

    Nspinor = m.Nspinor#ifelse( meas.fparam.Dirac_operator == "Staggered" ,1,4)
    is = ((i - 1) % Nspinor) + 1 # spin index   
    ic = ((i - is) ÷ Nspinor) + 1 # color index
    st = "$ic $is"
    measurestring *= st * "\n"
    println_verbose_level1(U[1], st)
    v = 1
    clear_fermion!(b)
    clear_fermion!(p)

    iorigin = (1, 1, 1, 1)
    setindex!(b, v, ic, iorigin..., is, 1)  # for domain wall fermion

    c = similar(b)
    apply_P!(c,b)
    d = similar(b)
    apply_R!(d,c)
    o = similar(p)
    @time solve_DinvX!(o, D, d)
    apply_P_edge!(p,o)

    st = "Hadron spectrum: Inversion $(i)/$(U[1].NC*m.Nspinor) is done"
    measurestring *= st * "\n"
    println_verbose_level1(U[1], st)

    flush(stdout)
    push!(stvec, measurestring)

    return deepcopy(p.w[1])
    
end
