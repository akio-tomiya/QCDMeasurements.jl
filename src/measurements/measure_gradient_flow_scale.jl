"""
    GradientFlowHistory

Per-configuration history produced by [`GradientFlowScale_measurement`](@ref).
`energy_density` contains the plaquette and/or clover discretizations of
the Wilson-flow action density.  Scale setting must be performed on an
ensemble with [`estimate_flow_scales`](@ref), not configuration by
configuration.
"""
struct GradientFlowHistory
    flow_time::Vector{Float64}
    plaquette::Vector{Float64}
    energy_density::Dict{String,Vector{Float64}}
    t2_energy::Dict{String,Vector{Float64}}
    topological_charge::Dict{String,Vector{Float64}}
    volume::Int
end

function GradientFlowHistory(
    flow_time::AbstractVector,
    energy_density::AbstractDict;
    plaquette=fill(NaN, length(flow_time)),
    topological_charge=Dict{String,Vector{Float64}}(),
    volume=1,
)
    times = Float64.(flow_time)
    energies = Dict{String,Vector{Float64}}(
        string(key) => Float64.(value) for (key, value) in energy_density
    )
    t2_energy = Dict(
        key => times .^ 2 .* value for (key, value) in energies
    )
    topology = Dict{String,Vector{Float64}}(
        string(key) => Float64.(value) for (key, value) in topological_charge
    )
    return GradientFlowHistory(
        times,
        Float64.(plaquette),
        energies,
        t2_energy,
        topology,
        Int(volume),
    )
end

"""
    GradientFlowScaleEstimate

Ensemble estimates of the standard gradient-flow scales.  For each action
density discretization, `t0` solves `t^2 * <E(t)> = c`, while `w0` solves
`t * d/dt(t^2 * <E(t)>) = c` at `t = w0^2`.  Errors are leave-one-out
jackknife errors and are `missing` for a one-configuration ensemble or when
the requested crossing is outside the measured flow-time range.
"""
struct GradientFlowScaleEstimate
    c::Float64
    flow_time::Vector{Float64}
    mean_energy_density::Dict{String,Vector{Float64}}
    t2_energy::Dict{String,Vector{Float64}}
    w_observable::Dict{String,Vector{Float64}}
    t0::Dict{String,Union{Missing,Float64}}
    t0_error::Dict{String,Union{Missing,Float64}}
    w0::Dict{String,Union{Missing,Float64}}
    w0_error::Dict{String,Union{Missing,Float64}}
    number_of_configurations::Int
end

"""
    GradientFlowScale_measurement(U; kwargs...)

Apply fixed-step Wilson flow to a private copy of `U` and record the
plaquette, plaquette/clover action density, and optional topological charge.
The input gauge field is not modified.

Important keywords are `flow_step_size`, `number_of_flow_steps`, and
`flow_measure_every`.  Supported `energy_methods` are `"plaquette"` and
`"clover"`; supported `topological_charge_methods` are those of
[`Topological_charge_measurement`](@ref).  The optional topology calculation
uses its three-way `improved_topological_charge_definition` selection.
"""
mutable struct GradientFlowScale_measurement{TG,TF,TP,TT} <: AbstractMeasurement
    filename::Union{Nothing,String}
    Uflow::Vector{TG}
    gradient_flow::TF
    plaquette_measurement::TP
    topological_measurement::TT
    flow_step_size::Float64
    number_of_flow_steps::Int
    flow_measure_every::Int
    energy_methods::Vector{String}
    measure_topological_charge::Bool
    topological_charge_methods::Vector{String}
    improved_topological_charge_definition::String
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool

    function GradientFlowScale_measurement(
        U::Vector{T};
        filename=nothing,
        flow_step_size=0.01,
        number_of_flow_steps=100,
        flow_measure_every=1,
        energy_methods=["plaquette", "clover"],
        measure_topological_charge=true,
        topological_charge_methods=["clover"],
        improved_topological_charge_definition="alexandrou",
        verbose_level=2,
        printvalues=false,
    ) where {T}
        length(U) == 4 || throw(ArgumentError(
            "GradientFlowScale_measurement supports four dimensions only",
        ))
        flow_step_size > 0 || throw(ArgumentError("flow_step_size must be positive"))
        number_of_flow_steps >= 0 || throw(ArgumentError(
            "number_of_flow_steps must be nonnegative",
        ))
        flow_measure_every >= 1 || throw(ArgumentError(
            "flow_measure_every must be at least one",
        ))

        methods = String.(energy_methods)
        isempty(methods) && throw(ArgumentError("energy_methods must not be empty"))
        length(unique(methods)) == length(methods) || throw(ArgumentError(
            "energy_methods must not contain duplicates",
        ))
        for method in methods
            method in ("plaquette", "clover") || throw(ArgumentError(
                "energy method $method is not supported",
            ))
        end

        topology_methods = String.(topological_charge_methods)
        improved_definition = _validate_improved_topological_charge_definition(
            improved_topological_charge_definition,
        )
        if measure_topological_charge
            isempty(topology_methods) && throw(ArgumentError(
                "topological_charge_methods must not be empty when topology is enabled",
            ))
            for method in topology_methods
                method in ("plaquette", "clover") || throw(ArgumentError(
                    "topological-charge method $method is not supported",
                ))
            end
        end

        Uflow = similar(U)
        for μ in eachindex(U)
            substitute_U!(Uflow[μ], U[μ])
        end
        gradient_flow = Gradientflow(Uflow; Nflow=1, eps=flow_step_size)
        plaquette_measurement = Plaquette_measurement(Uflow)
        topological_measurement = Topological_charge_measurement(
            Uflow;
            TC_methods=topology_methods,
            improved_topological_charge_definition=improved_definition,
        )

        verbose_print = if printvalues
            Verbose_print(verbose_level; myid=get_myrank(U), filename)
        else
            nothing
        end

        return new{T,typeof(gradient_flow),typeof(plaquette_measurement),
            typeof(topological_measurement)}(
            filename,
            Uflow,
            gradient_flow,
            plaquette_measurement,
            topological_measurement,
            Float64(flow_step_size),
            Int(number_of_flow_steps),
            Int(flow_measure_every),
            methods,
            Bool(measure_topological_charge),
            topology_methods,
            improved_definition,
            verbose_print,
            Bool(printvalues),
        )
    end
end

function GradientFlowScale_measurement(
    U::Vector{T},
    params::GradientFlowScale_parameters,
    filename="Gradient_flow_scale.txt",
) where {T}
    return GradientFlowScale_measurement(
        U;
        filename,
        flow_step_size=params.flow_step_size,
        number_of_flow_steps=params.number_of_flow_steps,
        flow_measure_every=params.flow_measure_every,
        energy_methods=params.energy_methods,
        measure_topological_charge=params.measure_topological_charge,
        topological_charge_methods=params.kinds_of_topological_charge,
        improved_topological_charge_definition=
            params.improved_topological_charge_definition,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
    )
end

function calculate_clover_action_density(
    U::Array{<:AbstractGaugefields{NC,4},1},
    temp_UμνTA,
    temps,
) where {NC}
    numofloops = calc_UμνTA!(temp_UμνTA, "clover", U, temps)
    value = 0.0
    for μ in 1:4
        for ν in (μ + 1):4
            value -= real(tr(temp_UμνTA[μ, ν], temp_UμνTA[μ, ν]))
        end
    end
    return value / (numofloops^2 * U[1].NV)
end

function _gradient_flow_observables(m::GradientFlowScale_measurement, U)
    plaquette = get_value(measure(m.plaquette_measurement, U))
    energies = Dict{String,Float64}()
    if "plaquette" in m.energy_methods
        number_of_planes = length(U) * (length(U) - 1) ÷ 2
        energies["plaquette"] = 2 * U[1].NC * number_of_planes * (1 - plaquette)
    end
    if "clover" in m.energy_methods
        topological = m.topological_measurement
        energies["clover"] = calculate_clover_action_density(
            U, topological.temp_UμνTA, topological._temporary_gaugefields,
        )
    end

    topology = if m.measure_topological_charge
        get_value(measure(m.topological_measurement, U))
    else
        Dict{String,Float64}()
    end
    return plaquette, energies, topology
end

function measure(
    m::GradientFlowScale_measurement,
    U::Array{<:AbstractGaugefields{NC,4},1};
    additional_string="",
) where {NC}
    for μ in eachindex(U)
        substitute_U!(m.Uflow[μ], U[μ])
    end

    steps = collect(0:m.flow_measure_every:m.number_of_flow_steps)
    if isempty(steps) || last(steps) != m.number_of_flow_steps
        push!(steps, m.number_of_flow_steps)
    end
    flow_time = Float64[]
    plaquette = Float64[]
    energy_density = Dict(method => Float64[] for method in m.energy_methods)
    topological_charge = Dict{String,Vector{Float64}}()
    measurestring = ""

    next_measurement = 1
    for step in 0:m.number_of_flow_steps
        if step == steps[next_measurement]
            time = step * m.flow_step_size
            plaq, energies, topology = _gradient_flow_observables(m, m.Uflow)
            push!(flow_time, time)
            push!(plaquette, plaq)
            for method in m.energy_methods
                push!(energy_density[method], energies[method])
            end
            for (method, value) in topology
                values = get!(topological_charge, method, Float64[])
                push!(values, value)
            end

            if m.printvalues
                fields = String[additional_string, string(time), string(plaq)]
                append!(fields, string(energies[method]) for method in m.energy_methods)
                append!(fields, string(topology[key]) for key in sort!(collect(keys(topology))))
                line = join(filter(field -> !isempty(field), fields), " ")
                measurestring *= line * "\n"
                println_verbose_level2(m.verbose_print, line)
            end

            next_measurement += 1
            next_measurement > length(steps) && break
        end
        flow!(m.Uflow, m.gradient_flow)
    end

    t2_energy = Dict(
        method => flow_time .^ 2 .* values
        for (method, values) in energy_density
    )
    history = GradientFlowHistory(
        flow_time,
        plaquette,
        energy_density,
        t2_energy,
        topological_charge,
        U[1].NV,
    )
    return Measurement_output(history, measurestring)
end

function _validate_flow_histories(histories::AbstractVector{<:GradientFlowHistory})
    isempty(histories) && throw(ArgumentError("at least one flow history is required"))
    reference = first(histories)
    length(reference.flow_time) >= 2 || throw(ArgumentError(
        "at least two flow-time samples are required",
    ))
    all(
        reference.flow_time[i] > reference.flow_time[i - 1]
        for i in 2:length(reference.flow_time)
    ) || throw(ArgumentError(
        "flow times must be strictly increasing",
    ))
    methods = sort!(collect(keys(reference.energy_density)))
    isempty(methods) && throw(ArgumentError("flow histories contain no energy density"))
    for history in histories
        history.flow_time == reference.flow_time || throw(ArgumentError(
            "all flow histories must use the same flow-time grid",
        ))
        sort!(collect(keys(history.energy_density))) == methods || throw(ArgumentError(
            "all flow histories must contain the same energy methods",
        ))
        for method in methods
            length(history.energy_density[method]) == length(reference.flow_time) ||
                throw(ArgumentError("energy history length does not match flow-time grid"))
            all(isfinite, history.energy_density[method]) || throw(ArgumentError(
                "energy histories must contain finite values",
            ))
        end
    end
    return reference.flow_time, methods
end

function _flow_derivative(times, values)
    n = length(times)
    derivative = similar(values, Float64)
    derivative[1] = (values[2] - values[1]) / (times[2] - times[1])
    for i in 2:(n - 1)
        derivative[i] = (values[i + 1] - values[i - 1]) /
                        (times[i + 1] - times[i - 1])
    end
    derivative[n] = (values[n] - values[n - 1]) / (times[n] - times[n - 1])
    return derivative
end

function _first_upward_crossing(times, values, target)
    for i in eachindex(times)
        values[i] == target && return Float64(times[i])
        if i > firstindex(times) && values[i - 1] < target < values[i]
            fraction = (target - values[i - 1]) / (values[i] - values[i - 1])
            return Float64(times[i - 1] + fraction * (times[i] - times[i - 1]))
        elseif i > firstindex(times) && values[i - 1] < target == values[i]
            return Float64(times[i])
        end
    end
    return missing
end

function _flow_scale_from_mean(times, mean_energy, c)
    t2_energy = times .^ 2 .* mean_energy
    w_observable = times .* _flow_derivative(times, t2_energy)
    t0 = _first_upward_crossing(times, t2_energy, c)
    w0_squared = _first_upward_crossing(times, w_observable, c)
    w0 = ismissing(w0_squared) ? missing : sqrt(w0_squared)
    return t2_energy, w_observable, t0, w0
end

function _jackknife_error(values)
    any(ismissing, values) && return missing
    n = length(values)
    n >= 2 || return missing
    numeric_values = Float64.(values)
    center = sum(numeric_values) / n
    return sqrt((n - 1) / n * sum((value - center)^2 for value in numeric_values))
end

"""
    estimate_flow_scales(histories; c=0.3)

Average `E(t)` over configurations first, then determine `t0` and `w0` by
linear interpolation.  This ordering implements the standard ensemble
definitions `t0^2 * <E(t0)> = c` and
`t*d/dt(t^2*<E(t)>)|t=w0^2 = c`.

The definitions follow M. Lüscher, JHEP 08 (2010) 071
(`arXiv:1006.4518`) and S. Borsanyi et al., JHEP 09 (2012) 010
(`arXiv:1203.4469`).
"""
function estimate_flow_scales(
    histories::AbstractVector{<:GradientFlowHistory};
    c=0.3,
)
    c > 0 || throw(ArgumentError("c must be positive"))
    times, methods = _validate_flow_histories(histories)
    nconfigurations = length(histories)
    mean_energy_density = Dict{String,Vector{Float64}}()
    t2_energy = Dict{String,Vector{Float64}}()
    w_observable = Dict{String,Vector{Float64}}()
    t0 = Dict{String,Union{Missing,Float64}}()
    t0_error = Dict{String,Union{Missing,Float64}}()
    w0 = Dict{String,Union{Missing,Float64}}()
    w0_error = Dict{String,Union{Missing,Float64}}()

    for method in methods
        total = zeros(Float64, length(times))
        for history in histories
            total .+= history.energy_density[method]
        end
        mean_energy = total ./ nconfigurations
        mean_energy_density[method] = mean_energy
        t2, w, t0_value, w0_value = _flow_scale_from_mean(times, mean_energy, c)
        t2_energy[method] = t2
        w_observable[method] = w
        t0[method] = t0_value
        w0[method] = w0_value

        if nconfigurations >= 2
            jackknife_t0 = Union{Missing,Float64}[]
            jackknife_w0 = Union{Missing,Float64}[]
            for omitted in eachindex(histories)
                jackknife_mean = (
                    total .- histories[omitted].energy_density[method]
                ) ./ (nconfigurations - 1)
                _, _, t0_jackknife, w0_jackknife = _flow_scale_from_mean(
                    times, jackknife_mean, c,
                )
                push!(jackknife_t0, t0_jackknife)
                push!(jackknife_w0, w0_jackknife)
            end
            t0_error[method] = _jackknife_error(jackknife_t0)
            w0_error[method] = _jackknife_error(jackknife_w0)
        else
            t0_error[method] = missing
            w0_error[method] = missing
        end
    end

    return GradientFlowScaleEstimate(
        Float64(c),
        copy(times),
        mean_energy_density,
        t2_energy,
        w_observable,
        t0,
        t0_error,
        w0,
        w0_error,
        nconfigurations,
    )
end
