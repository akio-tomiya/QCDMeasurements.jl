module QCDMeasurements
using Wilsonloop
using Gaugefields
using LatticeDiracOperators
using Arpack
using LinearAlgebra
import LatticeDiracOperators: clear_fermion!, Z4_distribution_fermi!
import Gaugefields:
    Temporalfields,
    evaluate_gaugelinks_eachsite!,
    get_temp,
    shift_U,
    unused!
import LatticeMatrices:
    projected_bilinear_slices,
    set_global_component!


include("parameters/parameters.jl")
include("measurements/AbstractMeasurement.jl")
include("support.jl")

# Canonical v1 public names.  The original spellings remain aliases so that
# existing input files and user code can migrate without changing physics.
const PlaquetteMeasurement = Plaquette_measurement
const PolyakovMeasurement = Polyakov_measurement
const PionCorrelatorMeasurement = Pion_correlator_measurement
const MesonCorrelatorMeasurement = Meson_correlator_measurement
const PCACMassMeasurement = PCAC_mass_measurement
const DomainWallResidualMassMeasurement = Domainwall_residual_mass_measurement
const ChiralCondensateMeasurement = Chiral_condensate_measurement
const EnergyDensityMeasurement = Energy_density_measurement
const CorrelationMeasurement = Correlation_measurement
const GluonicCorrelatorMeasurement = Guluonic_correlators_measurement
const TopologicalChargeMeasurement = Topological_charge_measurement
const GradientFlowScaleMeasurement = GradientFlowScale_measurement
const TopologicalChargeDensityCorrelationMeasurement =
    Topological_charge_density_correlation_measurement
const WilsonLoopMeasurement = Wilson_loop_measurement
const EigenvalueMeasurement = Eigenvalue_measurement
const MdagMSpectrumMeasurement = MdagMspectrum_measurement

const PlaquetteParameters = Plaq_parameters
const PolyakovParameters = Poly_parameters
const PionCorrelatorParameters = Pion_parameters
const MesonCorrelatorParameters = MesonCorrelator_parameters
const PCACMassParameters = PCACMass_parameters
const DomainWallResidualMassParameters = DomainWallResidualMass_parameters
const ChiralCondensateParameters = ChiralCondensate_parameters
const EnergyDensityParameters = Energy_density_parameters
const CorrelationParameters = Correlation_parameters
const GluonicCorrelatorParameters = Guluonic_correlators_parameters
const TopologicalChargeParameters = TopologicalCharge_parameters
const GradientFlowScaleParameters = GradientFlowScale_parameters
const TopologicalChargeDensityCorrelationParameters =
    TopologicalChargeDensityCorrelation_parameters
const WilsonLoopParameters = Wilson_loop_parameters
const EigenvalueParameters = Eigenvalue_parameters
const MdagMSpectrumParameters = MdagMspectrum_parameters

export AbstractMeasurement,
    MeasurementOutput,
    measure,
    get_value,
    get_string,
    supported_fermions
export PlaquetteMeasurement,
    PolyakovMeasurement,
    PionCorrelatorMeasurement,
    MesonCorrelatorMeasurement,
    PCACMassMeasurement,
    DomainWallResidualMassMeasurement,
    ChiralCondensateMeasurement,
    EnergyDensityMeasurement,
    CorrelationMeasurement,
    GluonicCorrelatorMeasurement,
    TopologicalChargeMeasurement,
    GradientFlowScaleMeasurement,
    TopologicalChargeDensityCorrelationMeasurement,
    WilsonLoopMeasurement,
    EigenvalueMeasurement,
    MdagMSpectrumMeasurement
export PlaquetteParameters,
    PolyakovParameters,
    PionCorrelatorParameters,
    MesonCorrelatorParameters,
    PCACMassParameters,
    DomainWallResidualMassParameters,
    ChiralCondensateParameters,
    EnergyDensityParameters,
    CorrelationParameters,
    GluonicCorrelatorParameters,
    TopologicalChargeParameters,
    GradientFlowScaleParameters,
    TopologicalChargeDensityCorrelationParameters,
    WilsonLoopParameters,
    EigenvalueParameters,
    MdagMSpectrumParameters

# Compatibility exports retained for the pre-v1 API.
export Measurement_output, Plaquette_measurement
export Polyakov_measurement
export Pion_correlator_measurement,
    PionSolverDiagnostic,
    get_solver_diagnostics
export Meson_correlator_measurement,
    MesonCorrelatorResult,
    LocalMesonChannel,
    StaggeredMesonChannel,
    local_meson_channel,
    standard_local_meson_channels,
    simulateqcd_staggered_channels
export PCAC_mass_measurement, PCACMassResult
export Domainwall_residual_mass_measurement, DomainWallResidualMassResult
export Chiral_condensate_measurement
export Energy_density_measurement
export Correlation_measurement
export Guluonic_correlators_measurement
export Topological_charge_measurement
export GradientFlowScale_measurement,
    GradientFlowHistory,
    GradientFlowScaleEstimate,
    estimate_flow_scales
export Topological_charge_density_correlation_measurement
export Wilson_loop_measurement
export Eigenvalue_measurement, MdagMspectrum_measurement
export prepare_measurement, prepare_measurement_from_dict
#export initialize_fermion_parameters

# Write your package code here.

end
