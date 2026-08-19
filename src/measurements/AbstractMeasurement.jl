"""Common supertype of all QCDMeasurements observables."""
abstract type AbstractMeasurement end

function get_temporary_fermionfields(m::AbstractMeasurement)
    return m._temporary_fermionfields
end

"""
    MeasurementOutput(value, outputstring)

Result returned by `measure`. `value` is the numerical or typed physics result,
while `outputstring` is the optional legacy text representation.
"""
struct MeasurementOutput{T}
    value::T
    outputstring::String
    MeasurementOutput(value, outputstring::AbstractString) =
        new{typeof(value)}(value, String(outputstring))
end

const Measurement_output = MeasurementOutput

function get_value(m::MeasurementOutput)
    return m.value
end

function get_string(m::MeasurementOutput)
    return m.outputstring
end

function measure(measurement::M, itrj, U) where {M<:AbstractMeasurement}
    throw(MethodError(measure, (measurement, itrj, U)))
end

const BoundaryCondition_4D_default = [1, 1, 1, -1]
const BoundaryCondition_2D_default = [1, -1]


include("measure_plaquette.jl")
include("measure_polyakov.jl")
include("measure_Pion_correlator.jl")
include("measure_Meson_correlator.jl")
include("measure_PCAC_mass.jl")
include("measure_Domainwall_residual_mass.jl")
include("measure_chiral_condensate.jl")
include("measure_energy_density.jl")
include("measure_Guluonic_correlators.jl")
include("measure_correlation.jl")
include("measure_topological_charge.jl")
include("measure_gradient_flow_scale.jl")
include("measure_topological_charge_density_correlation.jl")
include("measure_Wilon_loop.jl")
include("measure_eigenvalues.jl")
include("measure_MdagMspectrum.jl")
