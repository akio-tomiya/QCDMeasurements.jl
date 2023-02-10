abstract type AbstractMeasurement end

function get_temporary_gaugefields(m::AbstractMeasurement)
    return m._temporary_gaugefields
end

function get_temporary_fermionfields(m::AbstractMeasurement)
    return m._temporary_fermionfields
end

struct Measurement_output{T}
    value::T
    outputstring::String
    Measurement_output(value,st) = new{typeof(value)}(value,st)
end

function get_value(m::Measurement_output)
    return m.value 
end

function get_string(m::Measurement_output)
    return m.outputstring
end

function measure(measurement::M, itrj, U) where {M<:AbstractMeasurement}
    error("measure with a type $M is not supported")
end



include("measure_plaquette.jl")
include("measure_polyakov.jl")