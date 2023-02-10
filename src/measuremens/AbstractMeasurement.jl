abstract type AbstractMeasurement end

function get_temporary_gaugefields(m::AbstractMeasurement)
    return m._temporary_gaugefields
end

function get_temporary_fermionfields(m::AbstractMeasurement)
    return m._temporary_fermionfields
end

function measure(measurement::M, itrj, U) where {M<:AbstractMeasurement}
    error("measure with a type $M is not supported")
end