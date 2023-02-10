abstract type Measurement_parameters end

Base.@kwdef mutable struct Plaq_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Plaquette"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
end
