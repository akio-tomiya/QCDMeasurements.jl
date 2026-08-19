@testset "QCDMeasurements v1 public API" begin
    @test PlaquetteMeasurement === Plaquette_measurement
    @test PolyakovMeasurement === Polyakov_measurement
    @test PionCorrelatorMeasurement === Pion_correlator_measurement
    @test MesonCorrelatorMeasurement === Meson_correlator_measurement
    @test PCACMassMeasurement === PCAC_mass_measurement
    @test DomainWallResidualMassMeasurement ===
          Domainwall_residual_mass_measurement
    @test ChiralCondensateMeasurement === Chiral_condensate_measurement
    @test EnergyDensityMeasurement === Energy_density_measurement
    @test CorrelationMeasurement === Correlation_measurement
    @test GluonicCorrelatorMeasurement === Guluonic_correlators_measurement
    @test TopologicalChargeMeasurement === Topological_charge_measurement
    @test GradientFlowScaleMeasurement === GradientFlowScale_measurement
    @test TopologicalChargeDensityCorrelationMeasurement ===
          Topological_charge_density_correlation_measurement
    @test WilsonLoopMeasurement === Wilson_loop_measurement
    @test EigenvalueMeasurement === Eigenvalue_measurement
    @test MdagMSpectrumMeasurement === MdagMspectrum_measurement
    @test supported_fermions(PlaquetteMeasurement) == ()
    @test supported_fermions(PionCorrelatorMeasurement) ==
          (:Wilson, :WilsonClover, :Staggered)
    @test supported_fermions(MesonCorrelatorMeasurement) ==
          (:Wilson, :WilsonClover, :Staggered)
    @test supported_fermions(PCACMassMeasurement) ==
          (:Wilson, :WilsonClover)
    @test supported_fermions(DomainWallResidualMassMeasurement) ==
          (:Domainwall,)
    @test supported_fermions(ChiralCondensateMeasurement) ==
          (:Wilson, :Staggered)

    @test PlaquetteParameters === QCDMeasurements.Plaq_parameters
    @test MesonCorrelatorParameters ===
          QCDMeasurements.MesonCorrelator_parameters
    @test DomainWallResidualMassParameters ===
          QCDMeasurements.DomainWallResidualMass_parameters
    @test TopologicalChargeParameters ===
          QCDMeasurements.TopologicalCharge_parameters

    result = MeasurementOutput(1.25, "1.25")
    @test !(result isa AbstractMeasurement)
    @test get_value(result) == 1.25
    @test get_string(result) == "1.25"
    @test QCDMeasurements.Measurement_output === MeasurementOutput

    @test QCDMeasurements.initialize_measurement_parameters(:Plaquette) isa
          PlaquetteParameters
    @test QCDMeasurements.initialize_fermion_parameters(:Wilson) isa
          QCDMeasurements.Wilson_parameters
    clover_parameters =
        QCDMeasurements.initialize_fermion_parameters(:WilsonClover)
    @test clover_parameters.Dirac_operator == "WilsonClover"
    @test clover_parameters.hasclover
    @test QCDMeasurements.construct_Measurement_parameters_from_dict(Dict(
        :methodname => :Plaquette,
        :printvalues => false,
    )) isa PlaquetteParameters
    @test_throws ArgumentError QCDMeasurements.construct_Measurement_parameters_from_dict(Dict(
        :methodname => :Plaquette,
        :unknown_option => 1,
    ))
    @test_throws ArgumentError QCDMeasurements.construct_Measurement_parameters_from_dict(Dict(
        :printvalues => false,
    ))
end
