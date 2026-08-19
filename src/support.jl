"""
    supported_fermions(MeasurementType)
    supported_fermions(measurement)

Return the fermion formulations covered by the QCDMeasurements v1 validation
contract for a measurement type. An empty tuple denotes a gauge-only
observable. This reports validated public support, not every experimental code
path that may exist in a dependency.
"""
supported_fermions(::Type{<:AbstractMeasurement}) = ()
supported_fermions(measurement::AbstractMeasurement) =
    supported_fermions(typeof(measurement))

supported_fermions(::Type{<:Pion_correlator_measurement}) =
    (:Wilson, :WilsonClover, :Staggered)
supported_fermions(::Type{<:Meson_correlator_measurement}) =
    (:Wilson, :WilsonClover, :Staggered)
supported_fermions(::Type{<:PCAC_mass_measurement}) =
    (:Wilson, :WilsonClover)
supported_fermions(::Type{<:Domainwall_residual_mass_measurement}) =
    (:Domainwall,)
supported_fermions(::Type{<:Chiral_condensate_measurement}) =
    (:Wilson, :Staggered)
supported_fermions(::Type{<:Eigenvalue_measurement}) = (:Wilson,)
supported_fermions(::Type{<:MdagMspectrum_measurement}) = (:Wilson,)
