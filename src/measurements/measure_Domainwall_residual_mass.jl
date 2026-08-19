raw"""
    DomainWallResidualMassResult

Configuration-by-configuration Shamir domain-wall correlators.  The residual
mass estimator is

```math
m_{\rm res}(t)=\frac{C_{J_{5q}P}(t)}{C_{PP}(t)}.
```

The returned vector is the timeslice ratio; choosing a plateau and performing
an ensemble fit are deliberately left to analysis code.
"""
struct DomainWallResidualMassResult
    axis::Int
    source_position::NTuple{4,Int}
    momentum::NTuple{4,Int}
    pseudoscalar_correlator::Vector{ComplexF64}
    midpoint_pseudoscalar_correlator::Vector{ComplexF64}
    residual_mass::Vector{ComplexF64}
end

function Base.getindex(result::DomainWallResidualMassResult, quantity::Symbol)
    quantity in (:pseudoscalar, :PP) &&
        return result.pseudoscalar_correlator
    quantity in (:midpoint_pseudoscalar, :midpoint, :J5qP) &&
        return result.midpoint_pseudoscalar_correlator
    quantity in (:residual_mass, :mres, :ratio) &&
        return result.residual_mass
    throw(KeyError(quantity))
end

mutable struct Domainwall_residual_mass_measurement{TD,TF} <: AbstractMeasurement
    filename::Union{Nothing,String}
    D::TD
    template5::TF
    correlation_axis::Int
    source_position::NTuple{4,Int}
    momentum::NTuple{4,Int}
    printvalues::Bool
    solver_diagnostics::Vector{SolverDiagnostics}
end

raw"""
    Domainwall_residual_mass_measurement(U; kwargs...)

Measure the physical-boundary pseudoscalar correlator ``C_{PP}``, the
mid-plane correlator ``C_{J_{5q}P}``, and their timeslice ratio for Shamir
domain-wall fermions.  This follows the Furman--Shamir residual-mass
definition and Grid's `ContractJ5q` surface convention.

The current LatticeDiracOperators implementation requires a four-dimensional
Gaugefields MPILattice backend (also for one MPI rank), even `L5`, and the
Shamir coefficients ``b=1,c=1``.
"""
function Domainwall_residual_mass_measurement(
    U::Vector;
    fermiontype="Domainwall",
    mass=0.1,
    L5=4,
    M=-1.0,
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
    momentum=(0, 0, 0, 0),
    eps_CG=1e-14,
    MaxCGstep=5000,
    method_CG="bicg",
    BoundaryCondition=nothing,
    verbose_level=2,
    printvalues=false,
    filename=nothing,
)
    length(U) == 4 || throw(ArgumentError(
        "DomainWallResidualMassMeasurement requires a four-dimensional gauge field"))
    fermiontype == "Domainwall" || throw(ArgumentError(
        "DomainWallResidualMassMeasurement supports only fermiontype=\"Domainwall\""))
    hasproperty(U[1], :mpi) && getproperty(U[1], :mpi) || throw(ArgumentError(
        "DomainWallResidualMassMeasurement requires Gaugefields' MPILattice backend"))
    U[1].NC == 3 || throw(ArgumentError(
        "DomainWallResidualMassMeasurement currently requires three colors"))
    1 <= correlation_axis <= 4 || throw(ArgumentError(
        "correlation_axis must be between 1 and 4"))
    L5 >= 2 && iseven(L5) || throw(ArgumentError(
        "L5 must be an even integer greater than or equal to two"))
    length(source_position) == 4 || throw(DimensionMismatch(
        "source_position must have four entries"))
    source = ntuple(d -> Int(source_position[d]), 4)
    lattice_size = ntuple(d -> Int(U[1].U.gsize[d]), 4)
    for d in 1:4
        1 <= source[d] <= lattice_size[d] ||
            throw(BoundsError(1:lattice_size[d], source[d]))
    end
    normalized_momentum = _normalize_meson_momentum(momentum, correlation_axis)
    boundary_condition = BoundaryCondition === nothing ?
        copy(BoundaryCondition_4D_default) : Int.(BoundaryCondition)
    length(boundary_condition) == 4 || throw(DimensionMismatch(
        "BoundaryCondition must have four entries"))
    all(value -> value in (-1, 1), boundary_condition) || throw(ArgumentError(
        "BoundaryCondition entries must be +1 or -1"))
    method = String(method_CG)
    method in ("bicg", "bicgstab") || throw(ArgumentError(
        "method_CG must be bicg or bicgstab for domain-wall propagators"))

    template5 = Initialize_pseudofermion_fields(U[1], "Domainwall"; L5=Int(L5))
    parameters = Dict{String,Any}(
        "Dirac_operator" => "Domainwall",
        "mass" => Float64(mass),
        "L5" => Int(L5),
        "M" => Float64(M),
        "eps_CG" => Float64(eps_CG),
        "MaxCGstep" => Int(MaxCGstep),
        "method_CG" => method,
        "verbose_level" => Int(verbose_level),
        "boundarycondition" => boundary_condition,
    )
    D = Dirac_operator(U, template5, parameters)
    return Domainwall_residual_mass_measurement(
        filename === nothing ? nothing : String(filename),
        D,
        template5,
        Int(correlation_axis),
        source,
        normalized_momentum,
        Bool(printvalues),
        SolverDiagnostics[],
    )
end

function Domainwall_residual_mass_measurement(
    U::Vector,
    params::DomainWallResidualMass_parameters,
    filename="Domainwall_residual_mass.txt",
)
    fermion = params.fermion_parameters
    fermion isa Domainwall_parameters || throw(ArgumentError(
        "DomainWallResidualMassParameters requires Domainwall_parameters"))
    return Domainwall_residual_mass_measurement(
        U;
        filename,
        fermiontype=params.fermiontype,
        mass=fermion.m,
        L5=fermion.N5,
        M=fermion.M,
        correlation_axis=params.correlation_axis,
        source_position=Tuple(params.source_position),
        momentum=Tuple(params.momentum),
        eps_CG=params.eps,
        MaxCGstep=params.MaxCGstep,
        method_CG=params.method_CG,
        BoundaryCondition=params.BoundaryCondition,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
    )
end

get_solver_diagnostics(measurement::Domainwall_residual_mass_measurement) =
    copy(measurement.solver_diagnostics)

function measure(
    measurement::Domainwall_residual_mass_measurement,
    U::Array{<:AbstractGaugefields{NC,4},1};
    additional_string="",
) where NC
    current_D = measurement.D(U)
    propagators = domainwall_physical_point_propagators(
        current_D,
        measurement.template5;
        source_position=measurement.source_position,
    )
    correlators = domainwall_residual_mass_correlator(
        propagators.five_dimensional;
        axis=measurement.correlation_axis,
        origin=measurement.source_position,
        momentum=measurement.momentum,
    )
    measurement.solver_diagnostics = collect(propagators.diagnostics)
    result = DomainWallResidualMassResult(
        measurement.correlation_axis,
        measurement.source_position,
        measurement.momentum,
        ComplexF64.(collect(correlators.PP)),
        ComplexF64.(collect(correlators.J5qP)),
        ComplexF64.(collect(correlators.ratio)),
    )

    lines = String[
        "$additional_string #domainwall_residual_mass " *
        "axis=$(result.axis) source=$(result.source_position) " *
        "momentum=$(result.momentum)",
    ]
    if measurement.printvalues
        push!(lines, "PP " * join(result.pseudoscalar_correlator, " "))
        push!(lines, "J5qP " * join(result.midpoint_pseudoscalar_correlator, " "))
        push!(lines, "mres " * join(result.residual_mass, " "))
    end
    return MeasurementOutput(result, join(lines, "\n") * "\n")
end
