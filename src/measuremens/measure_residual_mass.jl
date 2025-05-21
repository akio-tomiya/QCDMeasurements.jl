using LinearAlgebra
mutable struct Residual_mass_measurement{Dim,TG,TD,TF,TF_vec} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    Dim::Int8
    mass::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    fermi_action::TF
    _temporary_fermionfields::Vector{TF_vec}
    Nr::Int64
    factor::Float64
    order::Int64

    function Residual_mass_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        fermiontype="Domainwall",
        mass=0.1,
        Nf=2,
        κ=1,
        r=1,
        L5=2,
        M=-1,
        b=1,
        c=1,
        eps_CG=1e-14,
        MaxCGstep=5000,
        BoundaryCondition=nothing,
        Nr=10,
        order=1,
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
            c
        )

        params["eps_CG"] = eps_CG
        params["verbose_level"] = verbose_level
        params["MaxCGstep"] = MaxCGstep
        params["boundarycondition"] = boundarycondition

        
        D = Dirac_operator(U, x, params)
        fermi_action = FermiAction(D, parameters_action)


        TD = typeof(D)
        TF = typeof(fermi_action)


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

        return new{Dim,T,TD,TF,TF_vec}(
            filename,
            _temporary_gaugefields,
            Dim,
            mass,
            verbose_print,
            printvalues,
            D,#::TD
            fermi_action,#::TF,
            _temporary_fermionfields,
            Nr,
            factor,
            order
        )

    end



end

function Residual_mass_measurement(
    U::Vector{T},
    params::ResidualMass_parameters,
    filename="Residual_mass.txt",
) where {T}
    if params.fermiontype == "Domainwall"
        method = Residual_mass_measurement(
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
        method = Residual_mass_measurement(
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
    else
        error("$(params.fermiontype) is not supported in Residual_mass_measurement")
    end

    return method
end

function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Residual_mass_measurement,NC,Dim}
    temps_fermi = get_temporary_fermionfields(m)
    p = temps_fermi[1]
    q = temps_fermi[2]

    D = m.D.D5DW(U)
    mass = m.mass

    pbp = 0.0
    Nr = m.Nr
    measurestring = ""

    for ir = 1:Nr
        Z4_distribution_fermi!(p)
        clear_fermion!(q)
        apply_P!(q, p)
        r = similar(p)
        apply_R!(r, q)
        s = similar(p)
        solve_DinvX!(s, D, r)
        t = similar(p)
        apply_P_edge!(t,s)
        den = dot(t,t)

        num = real(dot(q, s))

        #tmp = effective mass = quark mass + residual mass
        tmp = num / den

        if m.printvalues
            # println_verbose_level2(U[1],"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand mres")
            measurestring_ir = "# $ir $additional_string $(real(tmp)) # itrj irand mres"
            if m.order != 1
                measurestring_ir = "# $ir $additional_string"
                for i = 1:m.order
                    measurestring_ir *= " $(real(tmps[i])) "
                end
                measurestring_ir *= " # itrj irand mres: $(m.order)-th orders"
            end
            println_verbose_level2(m.verbose_print, measurestring_ir)
            measurestring *= measurestring_ir * "\n"
        end
        pbp += tmp
    end

    pbp -= mass

    pbp_value = pbp / Nr

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

    output = Measurement_output(pbp_value, measurestring)

    return output
end
