module QCDMeasurements
    using Wilsonloop
    using Gaugefields
    using LatticeDiracOperators
    import LatticeDiracOperators.Dirac_operators:clear_fermion!
    

    include("parameters/parameters.jl")
    include("measuremens/AbstractMeasurement.jl")

    export Plaquette_measurement,measure,get_value
    export Polyakov_measurement
    export Pion_correlator_measurement

# Write your package code here.

end
