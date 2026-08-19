
struct PionSolverDiagnostic
    source_number::Int
    source_color::Int
    source_spin::Int
    method::Symbol
    iterations::Int
    restart_count::Int
    convergence_branch::Symbol
    recursive_residual_squared::Float64
    target_residual_squared::Float64
    maximum_iterations::Int
    true_relative_residual::Float64
end


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
    S::Union{Nothing,Array{ComplexF64,Dim_2}}
    cov_neural_net::TCov#Union{Nothing,CovNeuralnet}
    solver_diagnostics::Vector{PionSolverDiagnostic}

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
        cSW=1.5612,
        L5=2,
        M=-1,
        eps_CG=1e-14,
        MaxCGstep=5000,
        BoundaryCondition=nothing,
        cov_neural_net=nothing,
        method_CG=nothing,
    ) where {T}
        if fermiontype == "Domainwall"
            throw(ArgumentError(
                "Domain-wall physical boundary-field projection is not yet " *
                "implemented in QCDMeasurements' pion measurement layer",
            ))
        end
        fermiontype in ("Wilson", "WilsonClover", "Staggered") || throw(ArgumentError(
            "PionCorrelatorMeasurement supports Wilson, WilsonClover, and " *
            "Staggered fermions; got $fermiontype",
        ))
        if fermiontype == "WilsonClover" &&
           !(hasproperty(U[1], :mpi) && getproperty(U[1], :mpi))
            throw(ArgumentError(
                "WilsonClover measurements require Gaugefields' MPILattice " *
                "backend so that LatticeMatrices.WilsonDiracCloverOperator4D " *
                "is used",
            ))
        end
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
            cSW=cSW,
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

        S = nothing


        params["eps_CG"] = eps_CG
        params["verbose_level"] = verbose_level
        params["MaxCGstep"] = MaxCGstep
        params["boundarycondition"] = boundarycondition
        if method_CG !== nothing
            method = String(method_CG)
            method in ("bicg", "bicgstab", "preconditiond_bicgstab") ||
                throw(
                    ArgumentError(
                        "method_CG must be bicg, bicgstab, or " *
                        "preconditiond_bicgstab",
                    ),
                )
            params["method_CG"] = method
            if fermiontype in ("Wilson", "WilsonClover") &&
               method == "preconditiond_bicgstab"
                fermiontype == "WilsonClover" && throw(ArgumentError(
                    "preconditiond_bicgstab is not supported for WilsonClover; " *
                    "use bicg or bicgstab",
                ))
                # LatticeDiracOperators currently defines its even-odd
                # wrapper only for the standard Wilson operator.
                params["faster version"] = false
            end
        end

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
            PionSolverDiagnostic[],
        )
    end

end

get_solver_diagnostics(m::Pion_correlator_measurement) =
    copy(m.solver_diagnostics)

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
            method_CG=params.method_CG,
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
        method = Pion_correlator_measurement(
            U;
            filename=filename,
            cov_neural_net=cov_neural_net,
            method_CG=params.method_CG,
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
        throw(ArgumentError(
            "Domain-wall physical boundary-field projection is not yet " *
            "implemented in QCDMeasurements' pion measurement layer",
        ))
    else
        error("fermiontype = $(params.fermiontype) is not supported")
    end

    return method
end

@inline function spincolor(ic, is, NC)
    return ic - 1 + (is - 1) * NC + 1
end

function _measurement_kernel_storage(field)
    if hasproperty(field, :field)
        return getproperty(field, :field)
    elseif hasproperty(field, :f)
        storage = getproperty(field, :f)
        return storage isa AbstractArray ? nothing : storage
    end
    return nothing
end

function _set_fermion_point_source!(
    field,
    value,
    color::Integer,
    spin::Integer,
    source_position::NTuple{D,<:Integer},
) where D
    storage = _measurement_kernel_storage(field)
    if storage !== nothing && applicable(
        set_global_component!, storage, value, color, spin, source_position)
        set_global_component!(storage, value, color, spin, source_position)
    else
        setindex_global!(field, value, color, source_position..., spin)
    end
    return field
end

function _accumulate_pion_correlator!(
    Cpi,
    propagator,
    NC,
    Nspinor,
    NN,
)
    length(NN) == 4 || throw(ArgumentError("only four dimensions are supported"))
    length(Cpi) == NN[4] ||
        throw(DimensionMismatch("expected $(NN[4]) correlator times, got $(length(Cpi))"))
    @inbounds for t = 1:NN[4]
        contribution = 0.0
        for z = 1:NN[3]
            for y = 1:NN[2]
                for x = 1:NN[1]
                    for sink_color = 1:NC
                        @simd for sink_spin = 1:Nspinor
                            contribution += abs2(
                                propagator[
                                    sink_color,
                                    x,
                                    y,
                                    z,
                                    t,
                                    sink_spin,
                                ],
                            )
                        end
                    end
                end
            end
        end
        Cpi[t] += contribution
    end
    return Cpi
end

function _accumulate_pion_propagator_block!(
    Cpi,
    propagators,
    NC,
    Nspinor,
    NN,
)
    length(propagators) == Nspinor || throw(DimensionMismatch(
        "expected $Nspinor source-spin propagators, got $(length(propagators))"))
    storages = ntuple(
        source_spin -> _measurement_kernel_storage(propagators[source_spin]),
        Nspinor,
    )
    spin_identity = Matrix{ComplexF64}(I, Nspinor, Nspinor)
    if all(storage -> storage !== nothing, storages) && applicable(
        projected_bilinear_slices,
        storages,
        storages,
        spin_identity,
        spin_identity,
    )
        Cpi .+= real.(projected_bilinear_slices(
            storages,
            storages,
            spin_identity,
            spin_identity;
            axis=4,
            origin=(1, 1, 1, 1),
            momentum=(0, 0, 0, 0),
            parity_mask=(0, 0, 0, 0),
            coefficient=1,
        ))
    else
        for propagator in propagators
            _accumulate_pion_correlator!(Cpi, propagator, NC, Nspinor, NN)
        end
    end
    return Cpi
end

function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Pion_correlator_measurement,NC,Dim}
    measurestring = ""
    st = "Hadron spectrum started"
    measurestring *= st * "\n"
    println_verbose_level2(U[1], st)
    Nspinor = m.Nspinor
    #D = m.D(U)
    # calculate quark propagators from a point source at he origin
    if m.cov_neural_net === nothing
        Cpi, st = calc_pion_correlator_point_source(m, U)
    else
        Uout, Uout_multi, _ = calc_smearedU(U, m.cov_neural_net)
        println("smeared U is used in Pion measurement")
        Cpi, st = calc_pion_correlator_point_source(m, Uout)
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




    Dim == 4 || error("Dim = $Dim is not supported")
    st = "Hadron spectrum: Contraction"
    measurestring *= st * "\n"
    println_verbose_level2(U[1], st)

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

function calc_pion_correlator_point_source(
    m,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
) where {NC,Dim}
    D = m.D(U)
    empty!(m.solver_diagnostics)
    _, _, NN... = size(U[1])
    Cpi = zeros(NN[end])
    diagnostics = PionSolverDiagnostic[]
    measurestrings = String[]
    try
        for source_color in 1:NC
            results = ntuple(source_spin ->
                calc_quark_propagators_point_source_each(
                    m,
                    U,
                    D,
                    (source_color - 1) * m.Nspinor + source_spin;
                    copy_propagator=true,
                ),
                m.Nspinor,
            )
            propagators = ntuple(
                source_spin -> results[source_spin].propagator,
                m.Nspinor,
            )
            _accumulate_pion_propagator_block!(
                Cpi, propagators, NC, m.Nspinor, NN)
            append!(diagnostics, (result.diagnostic for result in results))
            append!(measurestrings, (result.measurestring for result in results))
        end
    catch
        empty!(m.solver_diagnostics)
        rethrow()
    end
    m.solver_diagnostics = diagnostics
    st = join(measurestrings, "\n") * "\n"
    return Cpi, st
end


function calc_quark_propagators_point_source(
    m,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
) where {NC,Dim}
    # D^{-1} for each spin x color element
    D = m.D(U)
    empty!(m.solver_diagnostics)
    results = map(
        i -> calc_quark_propagators_point_source_each(m, U, D, i),
        1:NC*m.Nspinor,
    )
    propagators = [result.propagator for result in results]
    m.solver_diagnostics = [result.diagnostic for result in results]
    st = join((result.measurestring for result in results), "\n") * "\n"
    return propagators, st
end

function calc_quark_propagators_point_source_each(
    m,
    U,
    D,
    i;
    copy_propagator=true,
    source_position=(1, 1, 1, 1),
)
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
    p#rintln("ic = $ic is = $is")
    length(source_position) == 4 ||
        throw(ArgumentError("source_position must have four entries"))
    source_position_tuple = ntuple(d -> Int(source_position[d]), 4)
    _set_fermion_point_source!(b, v, ic, is, source_position_tuple)

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
    solver_result = solve_DinvX!(p, D, b)
    residual = similar(p)
    clear_fermion!(residual)
    mul!(residual, D, p)
    add_fermion!(residual, -1, b)
    source_norm_squared = real(b ⋅ b)
    source_norm_squared > 0 ||
        error("point source has zero norm for source $i")
    true_relative_residual =
        sqrt(real(residual ⋅ residual) / source_norm_squared)

    required_diagnostic_properties = (
        :method,
        :iterations,
        :restart_count,
        :convergence_branch,
        :recursive_residual_squared,
        :target_residual_squared,
        :maximum_iterations,
    )
    for property in required_diagnostic_properties
        hasproperty(solver_result, property) || error(
            "solver $(D.method_CG) did not report $property for source $i",
        )
    end

    method = Symbol(getproperty(solver_result, :method))
    iterations = Int(getproperty(solver_result, :iterations))
    restart_count = Int(getproperty(solver_result, :restart_count))
    convergence_branch =
        Symbol(getproperty(solver_result, :convergence_branch))
    recursive_residual_squared =
        Float64(getproperty(solver_result, :recursive_residual_squared))
    target_residual_squared =
        Float64(getproperty(solver_result, :target_residual_squared))
    maximum_iterations =
        Int(getproperty(solver_result, :maximum_iterations))

    iterations >= 0 || error(
        "solver $(D.method_CG) reported invalid iteration count $iterations for source $i",
    )
    restart_count >= 0 || error(
        "solver $(D.method_CG) reported invalid restart count $restart_count for source $i",
    )
    convergence_branch !== :unknown || error(
        "solver $(D.method_CG) did not report its convergence branch for source $i",
    )
    isfinite(recursive_residual_squared) || error(
        "solver $(D.method_CG) reported a non-finite recursive residual for source $i",
    )
    diagnostic = PionSolverDiagnostic(
        i,
        ic,
        is,
        method,
        iterations,
        restart_count,
        convergence_branch,
        recursive_residual_squared,
        target_residual_squared,
        maximum_iterations,
        true_relative_residual,
    )
    #error("dd")
    #println("norm p ",dot(p,p))
    st = "Hadron spectrum: Inversion $(i)/$(U[1].NC*m.Nspinor) is done"
    measurestring *= st * "\n"
    println_verbose_level1(U[1], st)

    flush(stdout)
    return (
        propagator=copy_propagator ? deepcopy(p) : p,
        diagnostic,
        measurestring,
    )
end
