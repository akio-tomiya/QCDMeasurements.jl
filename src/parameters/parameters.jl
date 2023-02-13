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


abstract type Measurement_parameters end

Base.@kwdef mutable struct Plaq_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Plaquette"
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct Poly_parameters <: Measurement_parameters
    methodname::String = "Polyakov_loop"
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    #common::Measurement_common_parameters = Measurement_common_parameters()
end

Base.@kwdef mutable struct Pion_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Pion_correlator"
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    #smearing::Smearing_parameters = NoSmearing_parameters()
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct ChiralCondensate_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Chiral_condensate"
    #measure_every::Int64 = 10
    fermiontype::String = "Staggered"
    Nf::Int64 = 4
    eps::Float64 = 1e-19
    mass::Float64 = 0.5
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    verbose_level::Int64 = 2
    printvalues::Bool = true
    Nr = 10
    #smearing::Smearing_parameters = Stout_parameters()
end

Base.@kwdef mutable struct Energy_density_parameters <: Measurement_parameters
    methodname::String = "Energy_density"
    #measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    #common::Measurement_common_parameters = Measurement_common_parameters()
end

Base.@kwdef mutable struct TopologicalCharge_parameters <: Measurement_parameters
    methodname::String = "Topological_charge"
    #measure_every::Int64 = 10
    fermiontype::String = "nothing"
    #common::Measurement_common_parameters = Measurement_common_parameters()
    #numflow::Int64 = 1 #number of flows
    #Nflowsteps::Int64 = 1
    #eps_flow::Float64 = 0.01
    verbose_level::Int64 = 2
    printvalues::Bool = true
    kinds_of_topological_charge::Vector{String} = ["plaquette","clover"]
end