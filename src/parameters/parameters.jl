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

const FERMION_PARAMETER_TYPES = Dict{String,DataType}(
    "nothing" => Quench_parameters,
    "Wilson" => Wilson_parameters,
    "WilsonClover" => Wilson_parameters,
    "Staggered" => Staggered_parameters,
    "Domainwall" => Domainwall_parameters,
)

function initialize_fermion_parameters(fermion_type::Union{AbstractString,Symbol})
    name = String(fermion_type)
    parameter_type = get(FERMION_PARAMETER_TYPES, name, nothing)
    parameter_type === nothing && throw(ArgumentError(
        "fermion type $name is not supported; choose one of " *
        join(sort!(collect(keys(FERMION_PARAMETER_TYPES))), ", "),
    ))
    parameters = parameter_type()
    if name == "WilsonClover"
        parameters.Dirac_operator = "WilsonClover"
        parameters.hasclover = true
    end
    return parameters
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
    method_CG::String = "bicg"
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

Base.@kwdef mutable struct MesonCorrelator_parameters <: Measurement_parameters
    methodname::String = "Meson_correlator"
    measure_every::Int64 = 10
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-14
    MaxCGstep::Int64 = 5000
    method_CG::String = "bicg"
    channels::Vector{String} = ["pseudoscalar"]
    momenta::Vector{Vector{Int64}} = [[0, 0, 0]]
    correlation_axis::Int64 = 4
    source_position::Vector{Int64} = [1, 1, 1, 1]
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct PCACMass_parameters <: Measurement_parameters
    methodname::String = "PCAC_mass"
    measure_every::Int64 = 10
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-14
    MaxCGstep::Int64 = 5000
    method_CG::String = "bicg"
    correlation_axis::Int64 = 4
    source_position::Vector{Int64} = [1, 1, 1, 1]
    improvement_coefficient::Float64 = 0.0
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct DomainWallResidualMass_parameters <: Measurement_parameters
    methodname::String = "Domainwall_residual_mass"
    measure_every::Int64 = 10
    fermiontype::String = "Domainwall"
    eps::Float64 = 1e-14
    MaxCGstep::Int64 = 5000
    method_CG::String = "bicg"
    correlation_axis::Int64 = 4
    source_position::Vector{Int64} = [1, 1, 1, 1]
    momentum::Vector{Int64} = [0, 0, 0, 0]
    BoundaryCondition::Vector{Int64} = [1, 1, 1, -1]
    fermion_parameters::Fermion_parameters = Domainwall_parameters()
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
    hop::Float64 = 0.141139
    r::Float64 = 1.0
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

Base.@kwdef mutable struct GradientFlowScale_parameters <: Measurement_parameters
    methodname::String = "Gradient_flow_scale"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    flow_step_size::Float64 = 0.01
    number_of_flow_steps::Int64 = 100
    flow_measure_every::Int64 = 1
    energy_methods::Vector{String} = ["plaquette", "clover"]
    measure_topological_charge::Bool = true
    kinds_of_topological_charge::Vector{String} = ["clover"]
    improved_topological_charge_definition::String = "alexandrou"
    verbose_level::Int64 = 2
    printvalues::Bool = true
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
    methodname::String = "Guluonic_correlators"
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
    improved_topological_charge_definition::String = "alexandrou"
end

Base.@kwdef mutable struct TopologicalChargeDensityCorrelation_parameters <: Measurement_parameters
    methodname::String = "Topological_charge_density_correlation"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    kinds_of_topological_charge::Vector{String} = ["plaquette", "clover"]
    improved_topological_charge_definition::String = "alexandrou"
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
    isDdagD::Bool = false #measure DdagD
    BoundaryCondition = [1, 1, 1, -1]
    solver = "Arpack"
end

Base.@kwdef mutable struct MdagMspectrum_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "MdagMspectrum"
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
    emin::Float64 = -2
    emax::Float64 = 2
    eta::Float64 = 0.01
    numpoints::Int64 = 3000
    position::Tuple{Int64,Int64} = (1, 1)

end


const MEASUREMENT_PARAMETER_TYPES = Dict{String,DataType}(
    "Plaquette" => Plaq_parameters,
    "Polyakov_loop" => Poly_parameters,
    "Topological_charge" => TopologicalCharge_parameters,
    "Topological_charge_density_correlation" =>
        TopologicalChargeDensityCorrelation_parameters,
    "Chiral_condensate" => ChiralCondensate_parameters,
    "Pion_correlator" => Pion_parameters,
    "Meson_correlator" => MesonCorrelator_parameters,
    "PCAC_mass" => PCACMass_parameters,
    "Domainwall_residual_mass" => DomainWallResidualMass_parameters,
    "Energy_density" => Energy_density_parameters,
    "Gradient_flow_scale" => GradientFlowScale_parameters,
    "Correlation" => Correlation_parameters,
    "Guluonic_correlators" => Guluonic_correlators_parameters,
    "Gluonic_correlators" => Guluonic_correlators_parameters,
    "Wilson_loop" => Wilson_loop_parameters,
    "Eigenvalue" => Eigenvalue_parameters,
    "MdagMspectrum" => MdagMspectrum_parameters,
)

function initialize_measurement_parameters(methodname::Union{AbstractString,Symbol})
    name = String(methodname)
    parameter_type = get(MEASUREMENT_PARAMETER_TYPES, name, nothing)
    parameter_type === nothing && throw(ArgumentError(
        "measurement $name is not supported; choose one of " *
        join(sort!(collect(keys(MEASUREMENT_PARAMETER_TYPES))), ", "),
    ))
    return parameter_type()
end

function _string_keyed_dictionary(values::AbstractDict)
    normalized = Dict{String,Any}()
    for (key, value) in pairs(values)
        normalized[String(key)] = value
    end
    return normalized
end

function _convert_parameter_value(current, value, key)
    value === nothing && return current
    target_type = typeof(current)
    try
        return convert(target_type, value)
    catch convert_exception
        convert_exception isa MethodError ||
            convert_exception isa InexactError ||
            convert_exception isa ArgumentError || rethrow()
    end
    try
        return target_type(value)
    catch constructor_exception
        constructor_exception isa MethodError ||
            constructor_exception isa InexactError ||
            constructor_exception isa ArgumentError || rethrow()
        throw(ArgumentError(
            "parameter $key expects $target_type, but received $(typeof(value))",
        ))
    end
end

function prepare_measurement_from_dict(U, values::AbstractDict, filename="")
    return prepare_measurement(U, values, filename)
end

function prepare_measurement(U, values::AbstractDict, filename="")
    parameters = construct_Measurement_parameters_from_dict(values)
    return prepare_measurement(U, parameters, filename)
end

function construct_Measurement_parameters_from_dict(values::AbstractDict)
    value_i = _string_keyed_dictionary(values)
    haskey(value_i, "methodname") || throw(ArgumentError(
        "methodname must be set in a measurement configuration",
    ))
    methodname = String(value_i["methodname"])
    method = initialize_measurement_parameters(methodname)
    method_dict = struct2dict(method)
    if haskey(value_i, "Dirac_operator")
        fermiontype = String(value_i["Dirac_operator"])
    else
        if haskey(value_i, "fermiontype")
            if value_i["fermiontype"] === nothing
                fermiontype = method.fermiontype
            else
                fermiontype = String(value_i["fermiontype"])
            end
        else
            fermiontype = method.fermiontype
        end
    end
    method.fermiontype = fermiontype
    fermion_parameters = initialize_fermion_parameters(fermiontype)
    fermion_parameters_dict = struct2dict(fermion_parameters)

    for (key_ii, value_ii) in value_i
        if method isa DomainWallResidualMass_parameters && key_ii == "L5"
            key_ii = "N5"
        elseif method isa DomainWallResidualMass_parameters && key_ii == "mass"
            key_ii = "m"
        elseif key_ii == "cSW" && fermiontype == "WilsonClover"
            key_ii = "Clover_coefficient"
        end
        if haskey(method_dict, key_ii)
            field = Symbol(key_ii)
            converted = _convert_parameter_value(
                getfield(method, field), value_ii, key_ii,
            )
            setfield!(method, field, converted)
        elseif haskey(fermion_parameters_dict, key_ii)
            field = Symbol(key_ii)
            converted = _convert_parameter_value(
                getfield(fermion_parameters, field), value_ii, key_ii,
            )
            setfield!(fermion_parameters, field, converted)
        else
            throw(ArgumentError(
                "unknown parameter $key_ii for $(typeof(method))",
            ))
        end
    end

    if haskey(method_dict, "fermion_parameters")
        setfield!(method, :fermion_parameters, fermion_parameters)
    end
    return method
end

_measurement_filename(filename, default) =
    isempty(filename) ? default : String(filename)

prepare_measurement(U, parameters::Plaq_parameters, filename="") =
    Plaquette_measurement(U, parameters, _measurement_filename(filename, "Plaquette.txt"))
prepare_measurement(U, parameters::Poly_parameters, filename="") =
    Polyakov_measurement(U, parameters, _measurement_filename(filename, "Polyakov_loop.txt"))
prepare_measurement(U, parameters::TopologicalCharge_parameters, filename="") =
    Topological_charge_measurement(
        U, parameters, _measurement_filename(filename, "Topological_charge.txt"))
prepare_measurement(
    U, parameters::TopologicalChargeDensityCorrelation_parameters, filename="",
) = Topological_charge_density_correlation_measurement(
    U,
    parameters,
    _measurement_filename(filename, "Topological_charge_density_correlation.txt"),
)
prepare_measurement(U, parameters::ChiralCondensate_parameters, filename="") =
    Chiral_condensate_measurement(
        U, parameters, _measurement_filename(filename, "Chiral_condensate.txt"))
prepare_measurement(U, parameters::Pion_parameters, filename="") =
    Pion_correlator_measurement(
        U, parameters, _measurement_filename(filename, "Pion_correlator.txt"))
prepare_measurement(U, parameters::MesonCorrelator_parameters, filename="") =
    Meson_correlator_measurement(
        U, parameters, _measurement_filename(filename, "Meson_correlator.txt"))
prepare_measurement(U, parameters::PCACMass_parameters, filename="") =
    PCAC_mass_measurement(U, parameters, _measurement_filename(filename, "PCAC_mass.txt"))
prepare_measurement(U, parameters::DomainWallResidualMass_parameters, filename="") =
    Domainwall_residual_mass_measurement(
        U,
        parameters,
        _measurement_filename(filename, "Domainwall_residual_mass.txt"),
    )
prepare_measurement(U, parameters::Energy_density_parameters, filename="") =
    Energy_density_measurement(
        U, parameters, _measurement_filename(filename, "Energy_density.txt"))
prepare_measurement(U, parameters::GradientFlowScale_parameters, filename="") =
    GradientFlowScale_measurement(
        U, parameters, _measurement_filename(filename, "Gradient_flow_scale.txt"))
prepare_measurement(U, parameters::Correlation_parameters, filename="") =
    Correlation_measurement(U, parameters, _measurement_filename(filename, "Correlation.txt"))
prepare_measurement(U, parameters::Guluonic_correlators_parameters, filename="") =
    Guluonic_correlators_measurement(
        U, parameters, _measurement_filename(filename, "Gluonic_correlators.txt"))
prepare_measurement(U, parameters::Wilson_loop_parameters, filename="") =
    Wilson_loop_measurement(U, parameters, _measurement_filename(filename, "Wilson_loop.txt"))
prepare_measurement(U, parameters::Eigenvalue_parameters, filename="") =
    Eigenvalue_measurement(U, parameters, _measurement_filename(filename, "Eigenvalues.txt"))
prepare_measurement(U, parameters::MdagMspectrum_parameters, filename="") =
    MdagMspectrum_measurement(
        U, parameters, _measurement_filename(filename, "MdagMspectrum.txt"))

function prepare_measurement(U, parameters::Measurement_parameters, filename="")
    throw(ArgumentError(
        "measurement parameters $(typeof(parameters)) are not supported",
    ))
end


function make_fermionparameter_dict(U, fermiontype,
    mass,
    Nf,
    κ,
    r,
    L5,
    M;
    cSW=1.5612,
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
    elseif fermiontype == "WilsonClover"
        x = Initialize_pseudofermion_fields(U[1], "Wilson", nowing=true)
        params["Dirac_operator"] = "WilsonClover"
        params["κ"] = κ
        params["r"] = r
        params["cSW"] = cSW
        params["faster version"] = false
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
        params_tuple = (
            verbose_level=params.verbose_level,
            printvalues=params.printvalues,
            fermiontype=params.fermiontype,
            κ=fermionparameters.hop,
            r=fermionparameters.r,
            cSW=fermionparameters.Clover_coefficient,
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
