using LinearAlgebra
mutable struct Chiral_condensate_measurement{Dim,TG,TD,TF,TF_vec,TCov} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    Dim::Int8
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    fermi_action::TF
    _temporary_fermionfields::Vector{TF_vec}
    Nr::Int64
    factor::Float64
    order::Int64
    cov_neural_net::TCov

    function Chiral_condensate_measurement(
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
        b=1,
        c=1,
        bs=[1,1],
        cs=[1,1],
        eps_CG=1e-14,
        MaxCGstep=5000,
        BoundaryCondition=nothing,
        Nr=10,
        order=1,
        cov_neural_net=nothing
    ) where {T}

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

        params, parameters_action, x, factor = make_fermionparameter_dict(
            U,
            fermiontype,
            mass,
            Nf,
            κ,
            r,
            L5,
            M,
            b,
            c,
            bs,
            cs
        )
        #=
        Nfbase = 1
        factor = 1
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
            x = Initialize_pseudofermion_fields(U[1], "Domainwall", L5 = L5)
            params["Dirac_operator"] = "Domainwall"
            params["mass"] = mass
            params["L5"] = L5
            params["M"] = M
        else
            error(
                "fermion type $fermiontype is not supported in chiral condensate measurement",
            )
        end
        =#

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
        #=
        if U[1].mpi == false
            myrank = 0
        else
            myrank = U[1].myrank
        end
        =#
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
        if order > 1
            numf += 1
        end

        TF_vec = typeof(x)
        _temporary_fermionfields = Vector{TF_vec}(undef, numf)
        for i = 1:numf
            _temporary_fermionfields[i] = similar(x)
        end

        return new{Dim,T,TD,TF,TF_vec,TCov}(
            filename,
            _temporary_gaugefields,
            Dim,
            verbose_print,
            printvalues,
            D,#::TD
            fermi_action,#::TF,
            _temporary_fermionfields,
            Nr,
            factor,
            order,
            cov_neural_net
        )

    end



end

function Chiral_condensate_measurement(
    U::Vector{T},
    params::ChiralCondensate_parameters,
    filename="Chiral_condensate.txt",
) where {T}

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

    if params.fermiontype == "Staggered"
        method = Chiral_condensate_measurement(
            U;
            filename=filename,
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            mass=params.mass,
            Nf=params.Nf,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
            Nr=params.Nr,
            cov_neural_net=cov_neural_net,
        )
    elseif params.fermiontype == "Domainwall"
        method = Chiral_condensate_measurement(
            U;
            filename=filename,
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            mass=params.mass,
            L5=params.L5,
            M=params.M,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
            Nr=params.Nr,
        )
    elseif params.fermiontype == "MobiusDomainwall"
        method = Chiral_condensate_measurement(
            U;
            filename=filename,
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            mass=params.mass,
            L5=params.L5,
            M=params.M,
            b=params.b,
            c=params.c,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
            Nr=params.Nr,
        )
    elseif params.fermiontype == "GeneralizedDomainwall"
        method = Chiral_condensate_measurement(
            U;
            filename=filename,
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            mass=params.mass,
            L5=params.L5,
            M=params.M,
            bs=params.bs,
            cs=params.cs,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
            Nr=params.Nr,
        )
    else
        error("$(params.fermiontype) is not supported in Chiral_condensate_measurement")
    end


    return method
end

function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Chiral_condensate_measurement,NC,Dim}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    r = temps_fermi[2]
    if m.cov_neural_net === nothing
        D = m.D(U)
    else
        Uout, Uout_multi, _ = calc_smearedU(U, m.cov_neural_net)
        println("smeared U is used in chiral measurement")
        D = m.D(Uout)
    end
    pbp = 0.0
    #Nr = 100
    Nr = m.Nr
    measurestring = ""
    if m.order != 1
        tmps = zeros(ComplexF64, m.order)
        p2 = temps_fermi[3]
        pbps = zeros(ComplexF64, m.order)
    end

    for ir = 1:Nr
        clear_fermion!(p)
        Z4_distribution_fermi!(r)
        solve_DinvX!(p, D, r)
        tmp = dot(r, p) # hermitian inner product

        if m.order != 1
            tmps[1] = tmp
            for i = 2:m.order
                solve_DinvX!(p2, D, p)
                p, p2 = p2, p
                tmps[i] = dot(r, p)
            end
            pbps .+= tmps
        end

        if m.printvalues
            # println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            measurestring_ir = "# $ir $additional_string $(real(tmp)/U[1].NV) # itrj irand chiralcond"
            if m.order != 1
                measurestring_ir = "# $ir $additional_string"
                for i = 1:m.order
                    measurestring_ir *= " $(real(tmps[i])/U[1].NV) "
                end
                measurestring_ir *= " # itrj irand chiralcond: $(m.order)-th orders"
            end
            println_verbose_level2(m.verbose_print, measurestring_ir)
            measurestring *= measurestring_ir * "\n"
        end
        pbp += tmp
    end

    pbp_value = real(pbp / Nr) / U[1].NV * m.factor
    
    if m.order != 1
        pbp_values = real.(pbps / Nr) / U[1].NV * m.factor
    end

    if m.printvalues
        measurestring_ir = "$pbp_value # pbp Nr=$Nr"
        if m.order != 1
            measurestring_ir = " "
            for i = 1:m.order
                measurestring_ir *= " $(pbp_values[i]) "
            end
            measurestring_ir *= "# pbp Nr=$Nr"
        end
        println_verbose_level1(m.verbose_print, measurestring_ir)
        measurestring *= measurestring_ir * "\n"
        flush(stdout)
    end

    if m.order != 1
        output = Measurement_output(pbp_values, measurestring)
    else
        output = Measurement_output(pbp_value, measurestring)
    end

    return output
end

#=
"""
c-------------------------------------------------c
c     Random number function Z4  Noise
c     https://arxiv.org/pdf/1611.01193.pdf
c-------------------------------------------------c
    """
    function Z4_distribution_fermion!(x::AbstractFermionfields_4D{NC})  where NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f)[6]
        θ = 0.0
        N::Int32 = 4
        Ninv = Float64(1/N)
        for ialpha = 1:n6
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            @inbounds @simd for ic=1:NC
                                θ = Float64(rand(0:N-1))*π*Ninv # r \in [0,π/4,2π/4,3π/4]
                                x[ic,ix,iy,iz,it,ialpha] = cos(θ)+im*sin(θ) 
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermion!(x)

        return
    end

    =#

#for Domainwall
function measure(
    m::Chiral_condensate_measurement{Dim,TG,TD,TF,TF_vec},
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {NC,Dim,TG,TD,TF,TF_vec<:LatticeDiracOperators.Dirac_operators.Abstract_DomainwallFermion_5D}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    r = temps_fermi[2]
    D = m.D.D5DW(U)
    pbp = 0.0
    Nr = m.Nr
    measurestring = ""
    if m.order != 1
        tmps = zeros(ComplexF64, m.order)
        p2 = temps_fermi[3]
        pbps = zeros(ComplexF64, m.order)
    end

    for ir = 1:Nr
        clear_fermion!(p)
        Z4_distribution_fermi!(r)
        r2 = similar(r)
        apply_P!(r2, r)
        apply_R!(p, r2)
        solve_DinvX!(r, D, p)
        tmp = dot(r2, r) # hermitian inner product

        if m.order != 1
            tmps[1] = tmp
            for i = 2:m.order
                solve_DinvX!(p2, D, p)
                p, p2 = p2, p
                tmps[i] = dot(r, p)
            end
            pbps .+= tmps
        end

        if m.printvalues
            # println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            measurestring_ir = "# $ir $additional_string $(real(tmp)/U[1].NV) # itrj irand chiralcond"
            if m.order != 1
                measurestring_ir = "# $ir $additional_string"
                for i = 1:m.order
                    measurestring_ir *= " $(real(tmps[i])/U[1].NV) "
                end
                measurestring_ir *= " # itrj irand chiralcond: $(m.order)-th orders"
            end
            println_verbose_level2(m.verbose_print, measurestring_ir)
            measurestring *= measurestring_ir * "\n"
        end
        pbp += tmp
    end

    pbp_value = real(pbp / Nr) / U[1].NV * m.factor
    
    if m.order != 1
        pbp_values = real.(pbps / Nr) / U[1].NV * m.factor
    end

    if m.printvalues
        measurestring_ir = "$pbp_value # pbp Nr=$Nr"
        if m.order != 1
            measurestring_ir = " "
            for i = 1:m.order
                measurestring_ir *= " $(pbp_values[i]) "
            end
            measurestring_ir *= "# pbp Nr=$Nr"
        end
        println_verbose_level1(m.verbose_print, measurestring_ir)
        measurestring *= measurestring_ir * "\n"
        flush(stdout)
    end

    if m.order != 1
        output = Measurement_output(pbp_values, measurestring)
    else
        output = Measurement_output(pbp_value, measurestring)
    end

    return output
end

#for MobiusDomainwall
function measure(
    m::Chiral_condensate_measurement{Dim,TG,TD,TF,TF_vec},
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {NC,Dim,TG,TD,TF,TF_vec<:LatticeDiracOperators.Dirac_operators.Abstract_MobiusDomainwallFermion_5D}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    r = temps_fermi[2]
    D = m.D.D5DW(U)
    pbp = 0.0
    Nr = m.Nr
    measurestring = ""
    if m.order != 1
        tmps = zeros(ComplexF64, m.order)
        p2 = temps_fermi[3]
        pbps = zeros(ComplexF64, m.order)
    end

    for ir = 1:Nr
        clear_fermion!(p)
        Z4_distribution_fermi!(r)
        r2 = similar(r)
        apply_P!(r2, r)
        apply_R!(p, r2)
        solve_DinvX!(r, D, p)
        tmp = dot(r2, r) # hermitian inner product

        if m.order != 1
            tmps[1] = tmp
            for i = 2:m.order
                solve_DinvX!(p2, D, p)
                p, p2 = p2, p
                tmps[i] = dot(r, p)
            end
            pbps .+= tmps
        end

        if m.printvalues
            # println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            measurestring_ir = "# $ir $additional_string $(real(tmp)/U[1].NV) # itrj irand chiralcond"
            if m.order != 1
                measurestring_ir = "# $ir $additional_string"
                for i = 1:m.order
                    measurestring_ir *= " $(real(tmps[i])/U[1].NV) "
                end
                measurestring_ir *= " # itrj irand chiralcond: $(m.order)-th orders"
            end
            println_verbose_level2(m.verbose_print, measurestring_ir)
            measurestring *= measurestring_ir * "\n"
        end
        pbp += tmp
    end

    pbp_value = real(pbp / Nr) / U[1].NV * m.factor
    
    if m.order != 1
        pbp_values = real.(pbps / Nr) / U[1].NV * m.factor
    end

    if m.printvalues
        measurestring_ir = "$pbp_value # pbp Nr=$Nr"
        if m.order != 1
            measurestring_ir = " "
            for i = 1:m.order
                measurestring_ir *= " $(pbp_values[i]) "
            end
            measurestring_ir *= "# pbp Nr=$Nr"
        end
        println_verbose_level1(m.verbose_print, measurestring_ir)
        measurestring *= measurestring_ir * "\n"
        flush(stdout)
    end

    if m.order != 1
        output = Measurement_output(pbp_values, measurestring)
    else
        output = Measurement_output(pbp_value, measurestring)
    end

    return output
end

#for GeneralizedDomainwall
function measure(
    m::Chiral_condensate_measurement{Dim,TG,TD,TF,TF_vec},
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {NC,Dim,TG,TD,TF,TF_vec<:LatticeDiracOperators.Dirac_operators.Abstract_GeneralizedDomainwallFermion_5D}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    r = temps_fermi[2]
    D = m.D.D5DW(U)
    pbp = 0.0
    Nr = m.Nr
    measurestring = ""
    if m.order != 1
        tmps = zeros(ComplexF64, m.order)
        p2 = temps_fermi[3]
        pbps = zeros(ComplexF64, m.order)
    end

    for ir = 1:Nr
        clear_fermion!(p)
        Z4_distribution_fermi!(r)
        r2 = similar(r)
        apply_P!(r2, r)
        apply_R!(p, r2)
        solve_DinvX!(r, D, p)
        tmp = dot(r2, r) # hermitian inner product

        if m.order != 1
            tmps[1] = tmp
            for i = 2:m.order
                solve_DinvX!(p2, D, p)
                p, p2 = p2, p
                tmps[i] = dot(r, p)
            end
            pbps .+= tmps
        end

        if m.printvalues
            # println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            measurestring_ir = "# $ir $additional_string $(real(tmp)/U[1].NV) # itrj irand chiralcond"
            if m.order != 1
                measurestring_ir = "# $ir $additional_string"
                for i = 1:m.order
                    measurestring_ir *= " $(real(tmps[i])/U[1].NV) "
                end
                measurestring_ir *= " # itrj irand chiralcond: $(m.order)-th orders"
            end
            println_verbose_level2(m.verbose_print, measurestring_ir)
            measurestring *= measurestring_ir * "\n"
        end
        pbp += tmp
    end

    pbp_value = real(pbp / Nr) / U[1].NV * m.factor
    
    if m.order != 1
        pbp_values = real.(pbps / Nr) / U[1].NV * m.factor
    end

    if m.printvalues
        measurestring_ir = "$pbp_value # pbp Nr=$Nr"
        if m.order != 1
            measurestring_ir = " "
            for i = 1:m.order
                measurestring_ir *= " $(pbp_values[i]) "
            end
            measurestring_ir *= "# pbp Nr=$Nr"
        end
        println_verbose_level1(m.verbose_print, measurestring_ir)
        measurestring *= measurestring_ir * "\n"
        flush(stdout)
    end

    if m.order != 1
        output = Measurement_output(pbp_values, measurestring)
    else
        output = Measurement_output(pbp_value, measurestring)
    end

    return output
end