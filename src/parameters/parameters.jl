struct2dict(x::T) where {T} =
    Dict{String,Any}(string(fn) => getfield(x, fn) for fn ∈ fieldnames(T))

abstract type Fermion_parameters end

Base.@kwdef mutable struct Quench_parameters <: Fermion_parameters
    Dirac_operator::String = "nothing"
end

Base.@kwdef mutable struct Wilson_parameters <: Fermion_parameters
    Dirac_operator::String = "Wilson"
    hop::Float64 = 0.141139
    r::Float64 = 1
    hasclover::Bool = false
    Clover_coefficient::Float64 = 1.5612
end

Base.@kwdef mutable struct Staggered_parameters <: Fermion_parameters
    Dirac_operator::String = "Staggered"
    mass::Float64 = 0.5 #mass
    Nf::Int64 = 2 #flavor 
end

Base.@kwdef mutable struct Domainwall_parameters <: Fermion_parameters
    Dirac_operator::String = "Domainwall"
    N5::Int64 = 4
    M::Float64 = -1 #mass for Wilson operator which should be negative
    m::Float64 = 0.1 #physical mass
end

function initialize_fermion_parameters(fermion_type)
    if fermion_type == "nothing"
        fermion_parameter = Quench_parameters()
    elseif fermion_type == "Wilson" || fermion_type == "WilsonClover"
        fermion_parameter = Wilson_parameters()
    elseif fermion_type == "Staggered"
        fermion_parameter = Staggered_parameters()
    elseif fermion_type == "Domainwall"
        fermion_parameter = Domainwall_parameters()
    else
        @error "$fermion_type is not implemented in parameters.jl"
    end
    return fermion_parameter
end



abstract type Measurement_parameters end

Base.@kwdef mutable struct Plaq_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Plaquette"
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    measure_every::Int64 = 1
end

Base.@kwdef mutable struct Poly_parameters <: Measurement_parameters
    methodname::String = "Polyakov_loop"
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    measure_every::Int64 = 1
    #common::Measurement_common_parameters = Measurement_common_parameters()
end


Base.@kwdef mutable struct Wilson_loop_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Wilson_loop"
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    Tmax::Int64 = 4
    Rmax::Int64 = 4
    measure_every::Int64 = 10
end

Base.@kwdef mutable struct Pion_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Pion_correlator"
    measure_every::Int64 = 10
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Int64 = 0#Union{Nothing,Int64} = nothing
    stout_ρ::Array{Float64,1} = zeros(1)#Vector{Float64}(undef, 1)#Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Array{String,1} = [""] #Vector{String}(undef, 1)#Union{Nothing,Array{String,1}} = nothing
    #stout_numlayers::Union{Nothing,Int64} = nothing
    #stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    #stout_loops::Union{Nothing,Array{String,1}} = nothing

    #smearing::Smearing_parameters = NoSmearing_parameters()
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct ChiralCondensate_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Chiral_condensate"
    measure_every::Int64 = 10
    fermiontype::String = "Staggered"
    Nf::Int64 = 4
    eps::Float64 = 1e-19
    mass::Float64 = 0.5
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"

    stout_numlayers::Int64 = 0#Union{Nothing,Int64} = nothing
    stout_ρ::Vector{Float64} = zeros(1)#Vector{Float64}(undef, 1)#Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Vector{String} = [""] #Vector{String}(undef, 1)#Union{Nothing,Array{String,1}} = nothing
    #stout_numlayers::Union{Nothing,Int64} = nothing
    #stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    #stout_loops::Union{Nothing,Array{String,1}} = nothing
    verbose_level::Int64 = 2
    printvalues::Bool = true
    Nr = 10
    #smearing::Smearing_parameters = Stout_parameters()
end

Base.@kwdef mutable struct Energy_density_parameters <: Measurement_parameters
    methodname::String = "Energy_density"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    #common::Measurement_common_parameters = Measurement_common_parameters()
end


Base.@kwdef mutable struct Correlation_parameters <: Measurement_parameters
    methodname::String = "Correlation"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    loop1::Vector{Vector{Tuple{Int64,Int64}}} = [Tuple{Int64,Int64}[]]
    loop2::Vector{Vector{Tuple{Int64,Int64}}} = [Tuple{Int64,Int64}[]]
    relativeposition::Vector{Int64} = [0, 0, 0, 0]
    originonly::Bool = true
    loop1position::Vector{Int64} = [1, 1, 1, 1]
    #common::Measurement_common_parameters = Measurement_common_parameters()
end

Base.@kwdef mutable struct Guluonic_correlators_parameters <: Measurement_parameters
    methodname::String = "Correlation"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    loop1::Vector{Vector{Tuple{Int64,Int64}}} = [Tuple{Int64,Int64}[]]
    loop2::Vector{Vector{Tuple{Int64,Int64}}} = [Tuple{Int64,Int64}[]]
end


Base.@kwdef mutable struct TopologicalCharge_parameters <: Measurement_parameters
    methodname::String = "Topological_charge"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    #common::Measurement_common_parameters = Measurement_common_parameters()
    #numflow::Int64 = 1 #number of flows
    #Nflowsteps::Int64 = 1
    #eps_flow::Float64 = 0.01
    verbose_level::Int64 = 2
    printvalues::Bool = true
    kinds_of_topological_charge::Vector{String} = ["plaquette", "clover"]
end

Base.@kwdef mutable struct TopologicalChargeDensityCorrelation_parameters <: Measurement_parameters
    methodname::String = "Topological_charge_density_correlation"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    kinds_of_topological_charge::Vector{String} = ["plaquette", "clover"]
end


Base.@kwdef mutable struct Eigenvalue_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Eigenvalue"
    measure_every::Int64 = 10
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Int64 = 0#Union{Nothing,Int64} = nothing
    stout_ρ::Array{Float64,1} = Vector{Float64}(undef, 1)#Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Array{String,1} = Vector{String}(undef, 1)#Union{Nothing,Array{String,1}} = nothing
    #smearing::Smearing_parameters = NoSmearing_parameters()
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
    nev::Int64 = 10 #num. of eigenvalues
    which::Symbol = :SM # :SM smallest magnitude
end

function initialize_measurement_parameters(methodname)
    if methodname == "Plaquette"
        method = Plaq_parameters()
    elseif methodname == "Polyakov_loop"
        method = Poly_parameters()
    elseif methodname == "Topological_charge"
        method = TopologicalCharge_parameters()
    elseif methodname == "Topological_charge_density_correlation"
        method = TopologicalChargeDensityCorrelation_parameters()
    elseif methodname == "Chiral_condensate"
        method = ChiralCondensate_parameters()
    elseif methodname == "Pion_correlator"
        method = Pion_parameters()
    elseif methodname == "Energy_density"
        method = Energy_density_parameters()
    elseif methodname == "Correlation"
        method = Correlation_parameters()
    elseif methodname == "Guluonic_correlators"
        method = Guluonic_correlators_parameters()
    elseif methodname == "Wilson_loop"
        method = Wilson_loop_parameters()
    elseif methodname == "Eigenvalue"
        method = Eigenvalue_parameters()
    else
        @error "$methodname is not implemented in parameter_structs.jl"
    end
    return method
end

function prepare_measurement_from_dict(U, value_i::Dict, filename="")
    parameters = construct_Measurement_parameters_from_dict(value_i)
    return prepare_measurement(U, parameters, filename)
end


function construct_Measurement_parameters_from_dict(value_i::Dict)
    #println(value)
    @assert haskey(value_i, "methodname") "methodname should be set in measurement."
    methodname = value_i["methodname"]
    method = initialize_measurement_parameters(methodname)
    method_dict = struct2dict(method)
    #println("value_i ", value_i, haskey(value_i, "Dirac_operator"), value_i["Dirac_operator"])
    if haskey(value_i, "Dirac_operator")
        fermiontype = value_i["Dirac_operator"]
    else
        if haskey(value_i, "fermiontype")
            if value_i["fermiontype"] == nothing
                fermiontype = method.fermiontype
                #fermiontype = "nothing"
            else
                fermiontype = value_i["fermiontype"]
            end
        else
            fermiontype = method.fermiontype
            #fermiontype = "nothing"
        end
    end
    method.fermiontype = fermiontype
    #println("fermiontype $fermiontype")
    fermion_parameters = initialize_fermion_parameters(fermiontype)
    #println(fermion_parameters)
    fermion_parameters_dict = struct2dict(fermion_parameters)
    #println("femriontype ",fermiontype)

    for (key_ii, value_ii) in value_i
        #println("$key_ii $value_ii")
        if haskey(method_dict, key_ii)
            if typeof(value_ii) != Nothing
                #println(getfield(method, Symbol(key_ii)), "\t", Symbol(key_ii), "\t", value_ii)
                keytype = typeof(getfield(method, Symbol(key_ii)))
                setfield!(method, Symbol(key_ii), keytype(value_ii))
            end
        else
            if haskey(fermion_parameters_dict, key_ii)
                #println("fermion $key_ii $value_ii")
                keytype = typeof(getfield(fermion_parameters, Symbol(key_ii)))
                setfield!(fermion_parameters, Symbol(key_ii), keytype(value_ii))
            else
                @warn "$key_ii is not found! in $(typeof(method))"
            end
        end
    end

    if haskey(method_dict, "fermion_parameters")
        setfield!(method, Symbol("fermion_parameters"), fermion_parameters)
    end
    value_out = deepcopy(method)
    #println(value_out)

    return value_out
end

function prepare_measurement(U, measurement_parameters::T, filename="") where {T}
    if T == Plaq_parameters
        filename_input = ifelse(filename == "", "Plaquette.txt", filename)
        measurement = Plaquette_measurement(U, measurement_parameters, filename_input)
    elseif T == Poly_parameters
        filename_input = ifelse(filename == "", "Polyakov_loop.txt", filename)
        measurement = Polyakov_measurement(U, measurement_parameters, filename_input)
    elseif T == TopologicalCharge_parameters
        filename_input = ifelse(filename == "", "Topological_charge.txt", filename)
        measurement =
            Topological_charge_measurement(U, measurement_parameters, filename_input)
    elseif T == TopologicalChargeDensityCorrelation_parameters
        filename_input = ifelse(filename == "", "Topological_charge_density_correlation.txt", filename)
        measurement =
            Topological_charge_density_correlation_measurement(U, measurement_parameters, filename_input)
    elseif T == ChiralCondensate_parameters
        filename_input = ifelse(filename == "", "Chiral_condensate.txt", filename)
        measurement =
            Chiral_condensate_measurement(U, measurement_parameters, filename_input)
    elseif T == Pion_parameters
        filename_input = ifelse(filename == "", "Pion_correlator.txt", filename)
        #println(measurement_parameters)
        measurement = Pion_correlator_measurement(U, measurement_parameters, filename_input)
    elseif T == Energy_density_parameters
        filename_input = ifelse(filename == "", "Energy_density.txt", filename)
        measurement = Energy_density_measurement(U, measurement_parameters, filename_input)
    elseif T == Correlation_parameters
        filename_input = ifelse(filename == "", "Correlation.txt", filename)
        measurement = Correlation_measurement(U, measurement_parameters, filename_input)
    elseif T == Guluonic_correlators_parameters
        filename_input = ifelse(filename == "", "Guluonic_correlators.txt", filename)
        measurement = Guluonic_correlators_measurement(U, measurement_parameters, filename_input)
    elseif T == Wilson_loop_parameters
        filename_input = ifelse(filename == "", "Wilson_loop.txt", filename)
        measurement = Wilson_loop_measurement(U, measurement_parameters, filename_input)
    elseif T == Eigenvalue_parameters
        filename_input = ifelse(filename == "", "Eigenvalues.txt", filename)
        measurement = Eigenvalue_measurement(U, measurement_parameters, filename_input)
    else
        error(T, " is not supported in measurements")
    end
    return measurement
end


function make_fermionparameter_dict(U, fermiontype,
    mass,
    Nf,
    κ,
    r,
    L5,
    M,
)
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
        x = Initialize_pseudofermion_fields(U[1], "Wilson", nowing=true)
        params["Dirac_operator"] = "Wilson"
        params["κ"] = κ
        params["r"] = r
        params["faster version"] = true
    elseif fermiontype == "Domainwall"
        params["Dirac_operator"] = "Domainwall"
        params["mass"] = mass
        params["L5"] = L5
        params["M"] = M
        x = Initialize_pseudofermion_fields(U[1], "Domainwall", L5=L5)
    else
        error(
            "fermion type $fermiontype is not supported in chiral condensate measurement",
        )
    end
    return params, parameters_action, x, factor
end

function fermionparameter_params(params)
    fermionparameters = params.fermion_parameters
    #println(fermionparameters)
    #println(params)
    if params.fermiontype == "Staggered"
        params_tuple = (
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            mass=fermionparameters.mass,
            Nf=fermionparameters.Nf,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
        )
    elseif params.fermiontype == "Wilson" || params.fermiontype == "WilsonClover"
        if fermionparameters.hasclover
            error("WilsonClover is not implemented in Pion measurement")
        end
        params_tuple = (
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            κ=fermionparameters.hop,
            r=fermionparameters.r,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
        )
    elseif params.fermiontype == "Domainwall"
        #error("Domainwall fermion is not implemented in Pion measurement!")
        params_tuple = (
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            L5=fermionparameters.N5,
            M=fermionparameters.M,
            eps_CG=params.eps,
            MaxCGstep=params.MaxCGstep,
        )
    else
        error("fermiontype = $(params.fermiontype) is not supported")
    end
    return params_tuple

end