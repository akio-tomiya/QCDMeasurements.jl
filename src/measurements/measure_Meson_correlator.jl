abstract type AbstractMesonChannel end

raw"""
    LocalMesonChannel(name, sink_matrix; source_matrix=sink_matrix, coefficient=-1)

A connected local Wilson-like meson channel. The coefficient defaults to the
fermionic Wick-contraction sign in

```math
-\operatorname{Tr}[\Gamma_{\rm sink} S_2(x,0)
\bar\Gamma_{\rm source} S_1(0,x)].
```
"""
struct LocalMesonChannel <: AbstractMesonChannel
    name::Symbol
    sink_matrix::Matrix{ComplexF64}
    source_matrix::Matrix{ComplexF64}
    coefficient::ComplexF64

    function LocalMesonChannel(
        name,
        sink_matrix;
        source_matrix=sink_matrix,
        coefficient=-1,
    )
        size(sink_matrix) == (4, 4) || throw(DimensionMismatch(
            "sink_matrix must be 4 by 4, got $(size(sink_matrix))"))
        size(source_matrix) == (4, 4) || throw(DimensionMismatch(
            "source_matrix must be 4 by 4, got $(size(source_matrix))"))
        return new(
            Symbol(name),
            ComplexF64.(sink_matrix),
            ComplexF64.(source_matrix),
            ComplexF64(coefficient),
        )
    end
end

"""A local staggered meson channel represented by a coordinate-parity mask."""
struct StaggeredMesonChannel <: AbstractMesonChannel
    name::Symbol
    parity_mask::NTuple{4,Int}
    coefficient::ComplexF64

    function StaggeredMesonChannel(name, parity_mask; coefficient=1)
        length(parity_mask) == 4 ||
            throw(DimensionMismatch("parity_mask must have four entries"))
        mask = ntuple(d -> begin
            value = Int(parity_mask[d])
            value in (0, 1) || throw(ArgumentError(
                "parity-mask entries must be zero or one, got $value"))
            value
        end, 4)
        return new(Symbol(name), mask, ComplexF64(coefficient))
    end
end

function _euclidean_gamma_matrices()
    return ntuple(mu -> ComplexF64.((γ1, γ2, γ3, γ4)[mu]), 4)
end

function _euclidean_gamma5()
    gammas = _euclidean_gamma_matrices()
    return gammas[1] * gammas[2] * gammas[3] * gammas[4]
end

"""
    local_meson_channel(kind; direction=nothing, coefficient=-1)

Construct a standard local scalar, pseudoscalar, vector, axial-vector, or
tensor channel in the gamma-matrix convention used by LatticeDiracOperators.
Directions are one-based Euclidean indices.
"""
function local_meson_channel(kind::Symbol; direction=nothing, coefficient=-1)
    gammas = _euclidean_gamma_matrices()
    gamma5 = gammas[1] * gammas[2] * gammas[3] * gammas[4]
    matrix, name = if kind in (:scalar, :S)
        Matrix{ComplexF64}(I, 4, 4), :scalar
    elseif kind in (:pseudoscalar, :P, :pion)
        gamma5, :pseudoscalar
    elseif kind in (:vector, :V)
        direction isa Integer || throw(ArgumentError(
            "a vector channel requires direction=1, 2, 3, or 4"))
        1 <= direction <= 4 || throw(ArgumentError("invalid vector direction $direction"))
        gammas[direction], Symbol("vector_", direction)
    elseif kind in (:axial, :A, :axialvector)
        direction isa Integer || throw(ArgumentError(
            "an axial channel requires direction=1, 2, 3, or 4"))
        1 <= direction <= 4 || throw(ArgumentError("invalid axial direction $direction"))
        gammas[direction] * gamma5, Symbol("axial_", direction)
    elseif kind in (:tensor, :T)
        direction isa Tuple && length(direction) == 2 || throw(ArgumentError(
            "a tensor channel requires direction=(mu, nu)"))
        mu, nu = direction
        1 <= mu < nu <= 4 || throw(ArgumentError(
            "tensor directions must satisfy 1 <= mu < nu <= 4"))
        (gammas[mu] * gammas[nu] - gammas[nu] * gammas[mu]) / 2,
            Symbol("tensor_", mu, nu)
    else
        throw(ArgumentError("unsupported local meson channel $kind"))
    end
    return LocalMesonChannel(name, matrix; coefficient)
end

function standard_local_meson_channels()
    channels = LocalMesonChannel[
        local_meson_channel(:scalar),
        local_meson_channel(:pseudoscalar),
    ]
    append!(channels, [local_meson_channel(:vector; direction=mu) for mu in 1:4])
    append!(channels, [local_meson_channel(:axial; direction=mu) for mu in 1:4])
    append!(channels, [
        local_meson_channel(:tensor; direction=(mu, nu))
        for mu in 1:3 for nu in (mu + 1):4
    ])
    return channels
end

"""
    simulateqcd_staggered_channels(axis=4)

Return the eight local staggered M1--M8 parity phases implemented by
SIMULATeQCD's `measureHadrons` module. The ordering of the three transverse
directions follows increasing lattice-axis number.
"""
function simulateqcd_staggered_channels(axis::Integer=4)
    1 <= axis <= 4 || throw(ArgumentError("axis must be between 1 and 4"))
    transverse = Tuple(d for d in 1:4 if d != axis)
    transverse_masks = (
        (1, 1, 1),
        (0, 0, 0),
        (0, 1, 1),
        (1, 0, 1),
        (1, 1, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    )
    names = (:M1_scalar, :M2_pseudoscalar, :M3_axial_1, :M4_axial_2,
        :M5_axial_3, :M6_vector_1, :M7_vector_2, :M8_vector_3)
    return [StaggeredMesonChannel(names[i], ntuple(d -> begin
        location = findfirst(==(d), transverse)
        location === nothing ? 0 : transverse_masks[i][location]
    end, 4)) for i in 1:8]
end

struct MesonCorrelatorResult
    axis::Int
    source_position::NTuple{4,Int}
    momenta::Vector{NTuple{4,Int}}
    correlators::Dict{Symbol,Matrix{ComplexF64}}
end

Base.getindex(result::MesonCorrelatorResult, channel::Symbol) =
    result.correlators[channel]

mutable struct Meson_correlator_measurement{B1,B2,C} <: AbstractMeasurement
    first_quark::B1
    second_quark::B2
    same_quarks::Bool
    channels::Vector{C}
    momenta::Vector{NTuple{4,Int}}
    correlation_axis::Int
    source_position::NTuple{4,Int}
    printvalues::Bool
    solver_diagnostics::Vector{PionSolverDiagnostic}
end

function _normalize_meson_momentum(momentum, axis)
    if length(momentum) == 4
        result = ntuple(d -> Int(momentum[d]), 4)
    elseif length(momentum) == 3
        result = ntuple(d -> begin
            if d == axis
                0
            elseif d < axis
                Int(momentum[d])
            else
                Int(momentum[d - 1])
            end
        end, 4)
    else
        throw(DimensionMismatch("each momentum must have three or four entries"))
    end
    iszero(result[axis]) || throw(ArgumentError(
        "momentum along correlation_axis must be zero"))
    return result
end

function _normalize_local_channels(channels)
    return LocalMesonChannel[
        channel isa LocalMesonChannel ? channel : _parse_local_meson_channel(channel)
        for channel in channels
    ]
end

function _parse_local_meson_channel(channel)
    name = lowercase(String(channel))
    name in ("scalar", "s") && return local_meson_channel(:scalar)
    name in ("pseudoscalar", "p", "pion") &&
        return local_meson_channel(:pseudoscalar)
    match_result = match(r"^(vector|axial)_(\d)$", name)
    if match_result !== nothing
        kind = Symbol(match_result.captures[1])
        return local_meson_channel(kind; direction=parse(Int, match_result.captures[2]))
    end
    match_result = match(r"^tensor_([1-4])_?([1-4])$", name)
    if match_result !== nothing
        directions = parse.(Int, match_result.captures)
        return local_meson_channel(:tensor; direction=Tuple(directions))
    end
    throw(ArgumentError(
        "unknown local meson channel $channel; use scalar, pseudoscalar, " *
        "vector_1--vector_4, axial_1--axial_4, or tensor_12--tensor_34"))
end

function _normalize_staggered_channels(channels, axis)
    defaults = simulateqcd_staggered_channels(axis)
    by_name = Dict(channel.name => channel for channel in defaults)
    for (index, channel) in pairs(defaults)
        by_name[Symbol("M", index)] = channel
        by_name[Symbol("m", index)] = channel
    end
    return StaggeredMesonChannel[
        if channel isa StaggeredMesonChannel
            channel
        elseif channel isa Integer
            1 <= channel <= 8 || throw(ArgumentError("staggered channel must be M1--M8"))
            defaults[channel]
        else
            name = Symbol(channel)
            haskey(by_name, name) || throw(ArgumentError(
                "unknown staggered channel $channel; pass a StaggeredMesonChannel for a custom phase"))
            by_name[name]
        end
        for channel in channels
    ]
end

function Meson_correlator_measurement(
    U::Vector;
    channels=nothing,
    momenta=[(0, 0, 0)],
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
    fermiontype="Wilson",
    mass=0.1,
    mass2=nothing,
    Nf=2,
    κ=0.141139,
    κ2=nothing,
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
    length(U) == 4 || throw(ArgumentError("only four-dimensional meson correlators are supported"))
    1 <= correlation_axis <= 4 ||
        throw(ArgumentError("correlation_axis must be between 1 and 4"))
    _, _, lattice_size... = size(U[1])
    length(source_position) == 4 || throw(DimensionMismatch(
        "source_position must have four entries"))
    source = ntuple(d -> Int(source_position[d]), 4)
    for d in 1:4
        1 <= source[d] <= lattice_size[d] ||
            throw(BoundsError(1:lattice_size[d], source[d]))
    end
    normalized_momenta = [
        _normalize_meson_momentum(momentum, correlation_axis) for momentum in momenta
    ]
    isempty(normalized_momenta) && throw(ArgumentError("at least one momentum is required"))

    base_keywords = (
        fermiontype=fermiontype,
        Nf=Nf,
        r=r,
        cSW=cSW,
        eps_CG=eps_CG,
        MaxCGstep=MaxCGstep,
        BoundaryCondition=BoundaryCondition,
        method_CG=method_CG,
        verbose_level=verbose_level,
        printvalues=false,
        filename=filename,
    )
    first = Pion_correlator_measurement(U; mass, κ, base_keywords...)
    second_mass = mass2 === nothing ? mass : mass2
    second_kappa = κ2 === nothing ? κ : κ2
    same_quarks = fermiontype in ("Wilson", "WilsonClover") ? second_kappa == κ :
        second_mass == mass
    second = same_quarks ? first : Pion_correlator_measurement(
        U; mass=second_mass, κ=second_kappa, base_keywords...)

    normalized_channels = if fermiontype in ("Wilson", "WilsonClover")
        _normalize_local_channels(channels === nothing ? (:pseudoscalar,) : channels)
    elseif fermiontype == "Staggered"
        selected_channels = channels === nothing ? (1:8) : channels
        _normalize_staggered_channels(selected_channels, correlation_axis)
    else
        throw(ArgumentError(
            "Meson_correlator_measurement currently supports Wilson, " *
            "WilsonClover, and Staggered fermions"))
    end
    isempty(normalized_channels) && throw(ArgumentError("at least one channel is required"))
    names = getfield.(normalized_channels, :name)
    length(unique(names)) == length(names) ||
        throw(ArgumentError("meson channel names must be unique"))

    return Meson_correlator_measurement(
        first,
        second,
        same_quarks,
        normalized_channels,
        normalized_momenta,
        Int(correlation_axis),
        source,
        printvalues,
        PionSolverDiagnostic[],
    )
end

function Meson_correlator_measurement(
    U::Vector,
    params::MesonCorrelator_parameters,
    filename="Meson_correlator.txt",
)
    params_tuple = fermionparameter_params(params)
    return Meson_correlator_measurement(
        U;
        filename,
        channels=params.channels,
        momenta=params.momenta,
        correlation_axis=params.correlation_axis,
        source_position=params.source_position,
        method_CG=params.method_CG,
        params_tuple...,
    )
end

get_solver_diagnostics(m::Meson_correlator_measurement) =
    copy(m.solver_diagnostics)

function _wilson_contraction_matrices(channel::LocalMesonChannel)
    gamma4 = ComplexF64.(γ4)
    gamma5 = _euclidean_gamma5()
    source_adjoint = gamma4 * channel.source_matrix' * gamma4
    left = gamma5 * channel.sink_matrix
    right = source_adjoint * gamma5
    return left, right
end

function _reference_projected_bilinear_slices(
    propagators1,
    propagators2,
    left,
    right;
    axis,
    origin,
    momentum,
    parity_mask,
    coefficient,
)
    NS = length(propagators1)
    length(propagators2) == NS || throw(DimensionMismatch(
        "propagator blocks must have the same source-spin size"))
    field_size = size(propagators1[1])
    length(field_size) == 6 || throw(DimensionMismatch(
        "four-dimensional fermion fields must have six logical dimensions"))
    NC = field_size[1]
    lattice_size = field_size[2:5]
    sink_spins = field_size[6]
    sink_spins == NS || throw(DimensionMismatch(
        "propagator sink-spin size $sink_spins does not match source-spin size $NS"))
    result = zeros(ComplexF64, lattice_size[axis])
    for site in CartesianIndices(Tuple(lattice_size))
        position = Tuple(site)
        relative = ntuple(d -> position[d] - origin[d], 4)
        separation = mod(relative[axis], lattice_size[axis]) + 1
        parity = sum(parity_mask[d] * relative[d] for d in 1:4)
        angle = -2pi * sum(momentum[d] * relative[d] / lattice_size[d] for d in 1:4)
        weight = (iseven(parity) ? 1 : -1) * cis(angle)
        contraction = 0.0 + 0.0im
        for color in 1:NC
            P = Matrix{ComplexF64}(undef, NS, NS)
            Q = Matrix{ComplexF64}(undef, NS, NS)
            for source_spin in 1:NS, sink_spin in 1:NS
                P[sink_spin, source_spin] =
                    propagators1[source_spin][color, position..., sink_spin]
                Q[sink_spin, source_spin] =
                    propagators2[source_spin][color, position..., sink_spin]
            end
            contraction += sum((left * P * right) .* conj.(Q))
        end
        result[separation] += coefficient * weight * contraction
    end
    return result
end

function _contract_meson_block(
    propagators1,
    propagators2,
    left,
    right;
    axis,
    origin,
    momentum,
    parity_mask,
    coefficient,
)
    storages1 = ntuple(i -> _measurement_kernel_storage(propagators1[i]), length(propagators1))
    storages2 = ntuple(i -> _measurement_kernel_storage(propagators2[i]), length(propagators2))
    use_extension = all(x -> x !== nothing, storages1) &&
        all(x -> x !== nothing, storages2) &&
        applicable(projected_bilinear_slices, storages1, storages2, left, right)
    if use_extension
        return projected_bilinear_slices(
            storages1,
            storages2,
            left,
            right;
            axis,
            origin,
            momentum,
            parity_mask,
            coefficient,
        )
    end
    return _reference_projected_bilinear_slices(
        propagators1,
        propagators2,
        left,
        right;
        axis,
        origin,
        momentum,
        parity_mask,
        coefficient,
    )
end

function _solve_meson_source_color(base, U, D, source_color, source_position)
    NS = base.Nspinor
    results = ntuple(source_spin -> calc_quark_propagators_point_source_each(
        base,
        U,
        D,
        (source_color - 1) * NS + source_spin;
        copy_propagator=true,
        source_position,
    ), NS)
    propagators = ntuple(source_spin -> results[source_spin].propagator, NS)
    diagnostics = [result.diagnostic for result in results]
    messages = [result.measurestring for result in results]
    return propagators, diagnostics, messages
end

function measure(
    measurement::Meson_correlator_measurement,
    U::Array{<:AbstractGaugefields{NC,4},1};
    additional_string="",
) where NC
    first_D = measurement.first_quark.D(U)
    second_D = measurement.same_quarks ? first_D : measurement.second_quark.D(U)
    axis_length = size(U[1])[measurement.correlation_axis + 2]
    correlators = Dict(
        channel.name => zeros(ComplexF64, axis_length, length(measurement.momenta))
        for channel in measurement.channels
    )
    diagnostics = PionSolverDiagnostic[]
    messages = String[]

    for source_color in 1:NC
        propagators1, diagnostics1, messages1 = _solve_meson_source_color(
            measurement.first_quark, U, first_D, source_color, measurement.source_position)
        append!(diagnostics, diagnostics1)
        append!(messages, messages1)
        if measurement.same_quarks
            propagators2 = propagators1
        else
            propagators2, diagnostics2, messages2 = _solve_meson_source_color(
                measurement.second_quark,
                U,
                second_D,
                source_color,
                measurement.source_position,
            )
            append!(diagnostics, diagnostics2)
            append!(messages, messages2)
        end

        for channel in measurement.channels
            if channel isa LocalMesonChannel
                left, right = _wilson_contraction_matrices(channel)
                parity_mask = (0, 0, 0, 0)
            else
                left = ones(ComplexF64, 1, 1)
                right = ones(ComplexF64, 1, 1)
                parity_mask = channel.parity_mask
            end
            for (momentum_index, momentum) in pairs(measurement.momenta)
                correlators[channel.name][:, momentum_index] .+= _contract_meson_block(
                    propagators1,
                    propagators2,
                    left,
                    right;
                    axis=measurement.correlation_axis,
                    origin=measurement.source_position,
                    momentum,
                    parity_mask,
                    coefficient=channel.coefficient,
                )
            end
        end
    end
    measurement.solver_diagnostics = diagnostics
    measurement.first_quark.solver_diagnostics = diagnostics

    result = MesonCorrelatorResult(
        measurement.correlation_axis,
        measurement.source_position,
        copy(measurement.momenta),
        correlators,
    )
    output_lines = String[
        "$additional_string #meson_correlator axis=$(measurement.correlation_axis) " *
        "source=$(measurement.source_position)",
    ]
    append!(output_lines, messages)
    if measurement.printvalues
        for channel in measurement.channels, momentum_index in eachindex(measurement.momenta)
            values = correlators[channel.name][:, momentum_index]
            push!(output_lines,
                "$(channel.name) momentum=$(measurement.momenta[momentum_index]) " *
                join(("$(real(value)),$(imag(value))" for value in values), " "))
        end
    end
    return Measurement_output(result, join(output_lines, "\n") * "\n")
end
