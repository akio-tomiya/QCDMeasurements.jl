raw"""
    PCACMassResult

Configuration-by-configuration Wilson PCAC data in lattice units.  With
periodic indexing along the correlation axis,

```math
am_{\rm PCAC}(t) =
\frac{\widetilde\partial C_{A_aP}(t)
      + c_A\,\partial^*\partial C_{PP}(t)}{2 C_{PP}(t)},
```

where `symmetric_axial_derivative` stores
`(C_AP(t+1) - C_AP(t-1))/2` and `pseudoscalar_laplacian` stores
`C_PP(t+1) - 2C_PP(t) + C_PP(t-1)`.  The result is bare: current
renormalization and ensemble plateau fits are analysis-stage operations.
"""
struct PCACMassResult
    axis::Int
    source_position::NTuple{4,Int}
    improvement_coefficient::Float64
    pseudoscalar_correlator::Vector{ComplexF64}
    axial_pseudoscalar_correlator::Vector{ComplexF64}
    symmetric_axial_derivative::Vector{ComplexF64}
    pseudoscalar_laplacian::Vector{ComplexF64}
    mass::Vector{ComplexF64}
end

function Base.getindex(result::PCACMassResult, quantity::Symbol)
    quantity in (:pseudoscalar, :PP) && return result.pseudoscalar_correlator
    quantity in (:axial_pseudoscalar, :AP) &&
        return result.axial_pseudoscalar_correlator
    quantity in (:derivative, :symmetric_axial_derivative) &&
        return result.symmetric_axial_derivative
    quantity in (:laplacian, :pseudoscalar_laplacian) &&
        return result.pseudoscalar_laplacian
    quantity in (:mass, :pcac_mass) && return result.mass
    throw(KeyError(quantity))
end

mutable struct PCAC_mass_measurement{M} <: AbstractMeasurement
    meson_measurement::M
    improvement_coefficient::Float64
    printvalues::Bool
end

function _periodic_symmetric_derivative(values::AbstractVector)
    length(values) >= 3 || throw(ArgumentError(
        "the correlation axis must contain at least three sites"))
    return ComplexF64[
        (values[mod1(t + 1, length(values))] -
         values[mod1(t - 1, length(values))]) / 2
        for t in eachindex(values)
    ]
end

function _periodic_laplacian(values::AbstractVector)
    length(values) >= 3 || throw(ArgumentError(
        "the correlation axis must contain at least three sites"))
    return ComplexF64[
        values[mod1(t + 1, length(values))] - 2values[t] +
        values[mod1(t - 1, length(values))]
        for t in eachindex(values)
    ]
end

raw"""
    PCAC_mass_measurement(U; kwargs...)

Measure the zero-momentum local Wilson `PP` and `A_a P` correlators required
for the bare PCAC mass.  `correlation_axis=a` also selects the axial current
`A_a = bar(psi) gamma_a gamma_5 psi`.  `improvement_coefficient` is the
axial-current coefficient `c_A`; its default of zero gives the unimproved
nearest-neighbour definition.

This measurement only composes QCDMeasurements' generalized meson
measurement and does not require a new Dirac-operator implementation.
"""
function PCAC_mass_measurement(
    U::Vector;
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
    improvement_coefficient=0.0,
    fermiontype="Wilson",
    κ=0.141139,
    κ2=nothing,
    kappa=nothing,
    kappa2=nothing,
    r=1,
    cSW=1.5612,
    eps_CG=1e-14,
    MaxCGstep=5000,
    BoundaryCondition=nothing,
    method_CG=nothing,
    verbose_level=2,
    printvalues=false,
    filename=nothing,
)
    fermiontype in ("Wilson", "WilsonClover") || throw(ArgumentError(
        "PCAC_mass_measurement currently supports Wilson and WilsonClover fermions"))
    1 <= correlation_axis <= 4 || throw(ArgumentError(
        "correlation_axis must be between 1 and 4"))
    isfinite(improvement_coefficient) || throw(ArgumentError(
        "improvement_coefficient must be finite"))
    first_kappa = kappa === nothing ? κ : kappa
    second_kappa = kappa2 === nothing ? κ2 : kappa2
    isfinite(first_kappa) || throw(ArgumentError("kappa must be finite"))
    (second_kappa === nothing || isfinite(second_kappa)) ||
        throw(ArgumentError("kappa2 must be finite"))

    gamma5 = _euclidean_gamma5()
    gamma_axis = _euclidean_gamma_matrices()[correlation_axis]
    channels = LocalMesonChannel[
        LocalMesonChannel(
            :pseudoscalar_pseudoscalar,
            gamma5;
            source_matrix=gamma5,
        ),
        LocalMesonChannel(
            :axial_pseudoscalar,
            gamma_axis * gamma5;
            source_matrix=gamma5,
        ),
    ]
    meson = Meson_correlator_measurement(
        U;
        channels,
        momenta=[(0, 0, 0)],
        correlation_axis,
        source_position,
        fermiontype,
        κ=first_kappa,
        κ2=second_kappa,
        r,
        cSW,
        eps_CG,
        MaxCGstep,
        BoundaryCondition,
        method_CG,
        verbose_level,
        printvalues=false,
        filename,
    )
    axis_length = size(U[1])[correlation_axis + 2]
    axis_length >= 3 || throw(ArgumentError(
        "the correlation axis must contain at least three sites"))
    return PCAC_mass_measurement(
        meson,
        Float64(improvement_coefficient),
        Bool(printvalues),
    )
end

function PCAC_mass_measurement(
    U::Vector,
    params::PCACMass_parameters,
    filename="PCAC_mass.txt",
)
    params_tuple = fermionparameter_params(params)
    return PCAC_mass_measurement(
        U;
        filename,
        correlation_axis=params.correlation_axis,
        source_position=Tuple(params.source_position),
        improvement_coefficient=params.improvement_coefficient,
        method_CG=params.method_CG,
        params_tuple...,
    )
end

get_solver_diagnostics(m::PCAC_mass_measurement) =
    get_solver_diagnostics(m.meson_measurement)

function measure(
    measurement::PCAC_mass_measurement,
    U::Array{<:AbstractGaugefields{NC,4},1};
    additional_string="",
) where NC
    meson_output = measure(
        measurement.meson_measurement,
        U;
        additional_string,
    )
    meson_result = get_value(meson_output)
    pp = vec(copy(meson_result[:pseudoscalar_pseudoscalar][:, 1]))
    ap = vec(copy(meson_result[:axial_pseudoscalar][:, 1]))
    derivative = _periodic_symmetric_derivative(ap)
    laplacian = _periodic_laplacian(pp)
    numerator = derivative .+
        measurement.improvement_coefficient .* laplacian
    mass = numerator ./ (2 .* pp)
    result = PCACMassResult(
        meson_result.axis,
        meson_result.source_position,
        measurement.improvement_coefficient,
        pp,
        ap,
        derivative,
        laplacian,
        mass,
    )

    output_lines = String[
        chomp(get_string(meson_output)),
        "$additional_string #pcac_mass axis=$(result.axis) " *
        "source=$(result.source_position) c_A=$(result.improvement_coefficient)",
    ]
    if measurement.printvalues
        push!(output_lines, "PP " * join(pp, " "))
        push!(output_lines, "AP " * join(ap, " "))
        push!(output_lines, "am_PCAC " * join(mass, " "))
    end
    return Measurement_output(result, join(output_lines, "\n") * "\n")
end
