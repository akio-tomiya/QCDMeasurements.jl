# Gaugefields owns loop construction and halo exchange; JACC owns the local
# contractions. Only the final scalar field is transferred/gathered to the host.
function _density_lattice_size(U::Vector{<:AbstractGaugefields{NC,Dim}}) where {NC,Dim}
    Dim == 4 && length(U) == 4 || throw(ArgumentError(
        "gauge densities require four-dimensional fields and four link directions",
    ))
    lattice_size = size(U[1])[3:end]
    for link in U
        size(link) == (NC, NC, lattice_size...) || throw(DimensionMismatch(
            "all link directions must have the same lattice dimensions",
        ))
    end
    return lattice_size
end

"""
    GaugeDensityWorkspace(U)

Reusable loop fields and scratch storage for [`topological_charge_density!`](@ref)
and [`action_density!`](@ref). Local contractions run through JACC on CPU or
GPU. LatticeMatrices handles MPI collection, including GPU-backed distributed
fields. Construction and use are collective on the input communicator. Reuse
with the same field type, dimensions, halo width, communicator, and process
decomposition. Input links are not retained or modified. A workspace must not
be used concurrently.
"""
struct GaugeDensityWorkspace{TG,TT,TL,TR,TS}
    tensor::Matrix{TG}
    temps::TT
    loops::TL
    result::TR
    scratch::TS
    lattice_size::NTuple{4,Int}
    halo_width::Int
end

function GaugeDensityWorkspace(U::Vector{T}) where {T<:AbstractGaugefields}
    lattice_size = _density_lattice_size(U)
    tensor = [similar(U[1]) for μ in 1:4, ν in 1:4]
    loops = (
        plaquette=first(calc_loopset_μν_name("plaq", 4)),
        clover=first(calc_loopset_μν_name("clover", 4)),
        rectangle=first(calc_loopset_μν_name("rect", 4)),
    )
    result = if hasproperty(U[1], :myrank_xyzt) && hasproperty(U[1], :comm)
        _density_legacy_mpi_result(U[1], lattice_size)
    else
        _density_result_field(U[1].U, lattice_size)
    end
    return GaugeDensityWorkspace(tensor, Temporalfields(U[1]; num=6), loops,
        result, similar(_density_array(result)), lattice_size, U[1].NDW)
end

function _check_density_workspace(density, U, workspace::GaugeDensityWorkspace{T}) where {T}
    lattice_size = _density_lattice_size(U)
    size(density) == lattice_size == workspace.lattice_size || throw(DimensionMismatch(
        "density, workspace, and input lattice dimensions must match",
    ))
    eltype(U) === T || throw(ArgumentError("workspace and input field types must match"))
    all(link -> link.NDW == workspace.halo_width, U) || throw(ArgumentError(
        "workspace and input halo widths must match",
    ))
    _check_density_layout(U[1].U, workspace.tensor[1, 2].U)
    if hasproperty(U[1], :myrank_xyzt) && hasproperty(U[1], :comm)
        reference = workspace.tensor[1, 2]
        U[1].comm == reference.comm && U[1].PEs == reference.PEs &&
            U[1].myrank_xyzt == reference.myrank_xyzt || throw(ArgumentError(
                "workspace and input communicator and decomposition must match",
            ))
    end
    return nothing
end

# Store only the 24 nonzero epsilon terms; use the existing sign convention.
const _DENSITY_EPSILON_TERMS = Tuple(
    (μ, ν, ρ, σ, epsilon_tensor(μ, ν, ρ, σ))
    for μ in 1:4 for ν in 1:4 for ρ in 1:4 for σ in 1:4
    if length(Set((μ, ν, ρ, σ))) == 4
)

_density_array(field::LatticeMatrices.LatticeMatrix) = field.A
_density_array(field::AbstractArray) = field

# Implemented by the MPI extension for the legacy rank ordering (x fastest).
function _density_legacy_mpi_result end

function _density_result_field(source::LatticeMatrices.LatticeMatrix, lattice_size)
    return LatticeMatrices.LatticeMatrix(1, 1, 4, lattice_size,
        LatticeMatrices.get_PEs(source); nw=0, elementtype=Float64,
        comm0=source.comm, device_mapping=:current)
end

_density_result_field(source::AbstractArray, lattice_size) =
    JACC.zeros(Float64, 1, 1, lattice_size...)

_check_density_layout(source::AbstractArray, reference::AbstractArray) = nothing
function _check_density_layout(source::LatticeMatrices.LatticeMatrix,
                               reference::LatticeMatrices.LatticeMatrix)
    source.comm == reference.comm && source.dims == reference.dims &&
        source.coords == reference.coords || throw(ArgumentError(
            "workspace and input communicator and decomposition must match",
        ))
    return nothing
end

@inline function _density_trace_product(A, B, indices, ::Val{NC}) where {NC}
    value = 0.0 + 0.0im
    @inbounds for j in 1:NC, i in 1:NC
        value += A[i, j, indices...] * B[j, i, indices...]
    end
    return real(value)
end

@inline _density_plane(μ, ν) = 3 * (μ - 1) + ν - (ν > μ)

function _density_tensor_arrays(workspace)
    return ntuple(12) do plane
        μ = (plane - 1) ÷ 3 + 1
        ν = (plane - 1) % 3 + 1
        ν += ν >= μ
        field = _density_array(workspace.tensor[μ, ν].U)
        # No copy for native arrays; legacy host fields can also use JACC.
        convert(JACC.array_type(), field)
    end
end

@inline function _density_indices(i, local_size, halo)
    x = (i - 1) % local_size[1] + 1
    j = (i - 1) ÷ local_size[1]
    y = j % local_size[2] + 1
    j ÷= local_size[2]
    z = j % local_size[3] + 1
    t = j ÷ local_size[3] + 1
    return (x + halo, y + halo, z + halo, t + halo)
end

@inline function _kernel_topological_density!(i, output, fields, local_size, halo, nc, factor)
    indices = _density_indices(i, local_size, halo)
    value = 0.0
    @inbounds for (μ, ν, ρ, σ, ε) in _DENSITY_EPSILON_TERMS
        value += ε * _density_trace_product(
            fields[_density_plane(μ, ν)], fields[_density_plane(ρ, σ)], indices, nc)
    end
    @inbounds output[i] = factor * value
    return nothing
end

function _contract_topological_density!(output, workspace::GaugeDensityWorkspace{TG}, nloops) where {NC,TG<:AbstractGaugefields{NC,4}}
    fields = _density_tensor_arrays(workspace)
    local_size = size(output)[3:end]
    JACC.parallel_for(prod(local_size), _kernel_topological_density!, output, fields,
        local_size, workspace.halo_width, Val(NC), -1 / (32 * π^2 * nloops^2))
    return output
end

@inline function _kernel_improve_density!(i, density, scratch)
    @inbounds density[i] = (5 / 3) * density[i] - (1 / 6) * scratch[i]
    return nothing
end

_density_to_host(result::LatticeMatrices.LatticeMatrix) =
    LatticeMatrices.gather_and_bcast_matrix(result)
_density_to_host(result::AbstractArray) = Array(result)

function _copy_density_to_host!(density, workspace)
    values = _density_to_host(workspace.result)
    copyto!(vec(density), vec(values))
    return density
end

function _improve_density_tensor!(workspace, U)
    # Existing Bilson--Thompson normalization, in raw loop sums:
    # A_imp = (5/3) A_clover - (1/6) A_rectangle; contraction divides by 4^2.
    # Construct rectangles one plane at a time to avoid a second tensor field.
    temps = workspace.temps
    path_temps, path_indices = get_temp(temps, 4)
    evaluated, evaluated_index = get_temp(temps)
    rectangle, rectangle_index = get_temp(temps)
    try
        for μ in 1:4, ν in 1:4
            μ == ν && continue
            evaluate_gaugelinks!(evaluated, workspace.loops.rectangle[μ, ν], U, path_temps)
            Traceless_antihermitian!(rectangle, evaluated)
            clover = workspace.tensor[μ, ν]
            clear_U!(evaluated)
            add_U!(evaluated, 5 / 3, clover)
            add_U!(evaluated, -1 / 6, rectangle)
            substitute_U!(clover, evaluated)
        end
    finally
        unused!(temps, path_indices)
        unused!(temps, evaluated_index)
        unused!(temps, rectangle_index)
    end
    return nothing
end

function _topological_density_method(method)
    name = string(method)
    name in ("plaquette", "clover", "alexandrou", "bilson_thompson") ||
        throw(ArgumentError(
            "topological density method must be plaquette, clover, alexandrou, or bilson_thompson; got $name",
        ))
    return name
end

"""
    QCDMeasurements.topological_charge_density(U; method="clover")

Return `Array{Float64,4}` in `(X,Y,Z,T)` order, without changing `U`. Each
element is the dimensionless charge in one lattice cell (`a^4 q_phys`);
`sum(q)` agrees, up to summation roundoff, with `TopologicalChargeMeasurement`.
No volume average or additional powers of the lattice spacing are applied.

Methods (strings or symbols) and matching measurement dictionary keys:
- `"plaquette"`: `"plaquette"`.
- `"clover"`: `"clover"` (only the unimproved density is calculated).
- `"alexandrou"`: `"clover improved"` (direct density improvement).
- `"bilson_thompson"`: `"clover improved bilson-thompson"` (field-strength
  improvement, including the clover/rectangle cross terms).

Local contractions use JACC on CPU or GPU. For MPI fields this is a collective
operation: all ranks must call it, and every rank receives the complete global
array. Only the final scalar density is copied to host memory and gathered.
Qualify the function with `QCDMeasurements` when
also using Gaugefields, which exports its own function with the same name.
For repeated measurements use [`GaugeDensityWorkspace`](@ref) and
[`topological_charge_density!`](@ref).
"""
function topological_charge_density(U::Vector{<:AbstractGaugefields}; method="clover")
    name = _topological_density_method(method)
    workspace = GaugeDensityWorkspace(U)
    density = Array{Float64,4}(undef, workspace.lattice_size)
    return topological_charge_density!(density, U, workspace; method=name)
end

"""
    topological_charge_density!(q::Array{Float64,4}, U, workspace; method="clover")

Overwrite and return `q`, reusing a [`GaugeDensityWorkspace`](@ref).
Definitions and supported inputs are those of [`topological_charge_density`](@ref).
"""
function topological_charge_density!(
    density::Array{Float64,4}, U::Vector{<:AbstractGaugefields},
    workspace::GaugeDensityWorkspace; method="clover",
)
    name = _topological_density_method(method)
    _check_density_workspace(density, U, workspace)
    output = _density_array(workspace.result)
    if name == "plaquette"
        calc_UμνTA!(workspace.tensor, workspace.loops.plaquette, U, workspace.temps)
        _contract_topological_density!(output, workspace, 1)
        return _copy_density_to_host!(density, workspace)
    end

    calc_UμνTA!(workspace.tensor, workspace.loops.clover, U, workspace.temps)
    if name == "bilson_thompson"
        _improve_density_tensor!(workspace, U)
    end
    _contract_topological_density!(output, workspace, 4)
    if name == "alexandrou"
        calc_UμνTA!(workspace.tensor, workspace.loops.rectangle, U, workspace.temps)
        _contract_topological_density!(workspace.scratch, workspace, 8)
        # Q_rectangle in the existing Alexandrou definition is twice calc_Q.
        JACC.parallel_for(length(output), _kernel_improve_density!, output, workspace.scratch)
    end
    return _copy_density_to_host!(density, workspace)
end

function _action_density_method(method)
    name = string(method)
    name in ("plaquette", "clover") || throw(ArgumentError(
        "action density method must be plaquette or clover; got $name",
    ))
    return name
end

"""
    action_density(U; method="clover")

Return `Array{Float64,4}` in `(X,Y,Z,T)` order, using the dimensionless
Wilson-flow action-density convention of `GradientFlowScaleMeasurement`:

- `"clover"`: `e(x) = -sum(real(tr(F[μ,ν](x)^2)), μ<ν)`, with the
  traceless anti-Hermitian clover tensor `F = TA(sum of 4 loops)/4`.
- `"plaquette"`: `e(x) = 2sum(Nc - real(tr(P[μ,ν](x))), μ<ν)`.

Thus `sum(e)/length(e)` agrees with the corresponding gradient-flow energy
observable. Unit links give zero. No factor of `β`, `1/g^2`, or `a^-4` is
included. This is distinct from the historical `EnergyDensityMeasurement`,
which returns one for unit links. Strings and symbols are accepted.

Uses the same CPU/GPU JACC kernels and MPI collective return convention as
[`topological_charge_density`](@ref) and does not modify `U`.
"""
function action_density(U::Vector{<:AbstractGaugefields}; method="clover")
    name = _action_density_method(method)
    workspace = GaugeDensityWorkspace(U)
    return action_density!(Array{Float64,4}(undef, workspace.lattice_size), U, workspace; method=name)
end

"""
    action_density!(e::Array{Float64,4}, U, workspace; method="clover")

Overwrite and return `e`, reusing a [`GaugeDensityWorkspace`](@ref).
Definitions and supported inputs are those of [`action_density`](@ref).
"""
function action_density!(
    density::Array{Float64,4}, U::Vector{<:AbstractGaugefields{NC,Dim}},
    workspace::GaugeDensityWorkspace; method="clover",
) where {NC,Dim}
    name = _action_density_method(method)
    _check_density_workspace(density, U, workspace)
    if name == "clover"
        calc_UμνTA!(workspace.tensor, workspace.loops.clover, U, workspace.temps)
    else
        calc_large_wilson_loop!(workspace.tensor, workspace.loops.plaquette, U, workspace.temps)
    end
    output = _density_array(workspace.result)
    fields = _density_tensor_arrays(workspace)
    local_size = size(output)[3:end]
    JACC.parallel_for(prod(local_size), _kernel_action_density!, output, fields,
        local_size, workspace.halo_width, Val(NC), Val(name == "clover"))
    return _copy_density_to_host!(density, workspace)
end

@inline function _kernel_action_density!(i, output, fields, local_size, halo,
                                        nc::Val{NC}, ::Val{Clover}) where {NC,Clover}
    indices = _density_indices(i, local_size, halo)
    value = 0.0
    @inbounds for μ in 1:4, ν in (μ + 1):4
        field = fields[_density_plane(μ, ν)]
        if Clover
            value -= _density_trace_product(field, field, indices, nc) / 16
        else
            trace = 0.0
            for c in 1:NC
                trace += real(field[c, c, indices...])
            end
            value += 2 * (NC - trace)
        end
    end
    @inbounds output[i] = value
    return nothing
end
