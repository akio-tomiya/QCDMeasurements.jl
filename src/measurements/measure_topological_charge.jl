const _TOPOLOGICAL_CHARGE_IMPROVED_DEFINITIONS = (
    "alexandrou",
    "bilson_thompson",
    "both",
)
const _ALEXANDROU_IMPROVED_KEY = "clover improved"
const _BILSON_THOMPSON_IMPROVED_KEY = "clover improved bilson-thompson"

function _validate_improved_topological_charge_definition(definition)
    value = string(definition)
    value in _TOPOLOGICAL_CHARGE_IMPROVED_DEFINITIONS || throw(ArgumentError(
        "improved_topological_charge_definition must be one of " *
        join(_TOPOLOGICAL_CHARGE_IMPROVED_DEFINITIONS, ", ") *
        "; got $value",
    ))
    return value
end

"""
    Topological_charge_measurement(U; TC_methods=["plaquette"],
        improved_topological_charge_definition="alexandrou", kwargs...)

Measure gluonic topological charge on one gauge configuration.  When
`TC_methods` contains `"clover"`, `improved_topological_charge_definition`
selects the improved operator:

- `"alexandrou"` (default): directly improved density of Alexandrou,
  Athenodorou, and Jansen, Phys. Rev. D 92, 125014 (2015), Eqs. (15)-(17).
- `"bilson_thompson"`: improve the field-strength tensor first, following the
  construction of Bilson-Thompson, Leinweber, and Williams, Ann. Phys. 304,
  1 (2003), using the 1x1+1x2 coefficients also used by SIMULATeQCD.
- `"both"`: compute both improved operators.

The Alexandrou result retains the output key `"clover improved"`; the
field-strength-improved result uses `"clover improved bilson-thompson"`.
"""
mutable struct Topological_charge_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    temp_UμνTA::Matrix{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    TC_methods::Vector{String}
    improved_topological_charge_definition::String

    function Topological_charge_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        TC_methods=["plaquette"],
        improved_topological_charge_definition="alexandrou",
    ) where {T}
        myrank = get_myrank(U)
        improved_definition = _validate_improved_topological_charge_definition(
            improved_topological_charge_definition,
        )

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        temp_UμνTA = Array{T,2}(undef, Dim, Dim)

        for μ = 1:Dim
            for ν = 1:Dim
                temp_UμνTA[ν, μ] = similar(U[1])
            end
        end


        # Four path-evaluation temporaries plus fields for the evaluated loop
        # and its traceless anti-Hermitian part.  Six are needed by the
        # rectangle/clover cross term in the improved definition.
        numg = 6
        _temporary_gaugefields = Temporalfields(U[1], num=numg)
        #_temporary_gaugefields = Vector{T}(undef, numg)
        #_temporary_gaugefields[1] = similar(U[1])
        #for i = 2:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end



        return new{Dim,T}(
            filename,
            _temporary_gaugefields,
            temp_UμνTA,
            Dim,
            verbose_print,
            printvalues,
            TC_methods,
            improved_definition,
        )

    end



end

function Topological_charge_measurement(
    U::Vector{T},
    params::TopologicalCharge_parameters,
    filename="Topological_charge.txt",
) where {T}

    return Topological_charge_measurement(
        U,
        filename=filename,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
        TC_methods=params.kinds_of_topological_charge, #["plaquette"]
        improved_topological_charge_definition=
            params.improved_topological_charge_definition,
    )

end

function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Topological_charge_measurement,NC,Dim}
    temps = m._temporary_gaugefields# get_temporary_gaugefields(m)
    #temp1 = temps[1]
    #temp2 = temps[2]
    measurestring = ""

    nummethod = length(m.TC_methods)
    values = Float64[]
    valuedic = Dict{String,Float64}()
    output_labels = String[]
    printstring = " " * additional_string
    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            Qplaq = calculate_topological_charge_plaq(U, m.temp_UμνTA, temps)
            push!(values, real(Qplaq))
            valuedic["plaquette"] = real(Qplaq)
            push!(output_labels, "Qplaq")
        elseif methodname == "clover"
            Qclover = calculate_topological_charge_clover(U, m.temp_UμνTA, temps)
            push!(values, real(Qclover))
            valuedic["clover"] = real(Qclover)
            push!(output_labels, "Qclover")
            definition = m.improved_topological_charge_definition
            if definition in ("alexandrou", "both")
                Qimproved = calculate_topological_charge_improved(
                    U, m.temp_UμνTA, Qclover, temps,
                )
                push!(values, real(Qimproved))
                valuedic[_ALEXANDROU_IMPROVED_KEY] = real(Qimproved)
                push!(output_labels, "Qimproved_alexandrou")
            end
            if definition in ("bilson_thompson", "both")
                if definition == "both"
                    # The Alexandrou calculation overwrites temp_UμνTA
                    # with rectangles, so reconstruct the clover tensor.
                    Qclover = calculate_topological_charge_clover(
                        U, m.temp_UμνTA, temps,
                    )
                end
                Qimproved = calculate_topological_charge_bilson_thompson(
                    U, m.temp_UμνTA, Qclover, temps,
                )
                push!(values, real(Qimproved))
                valuedic[_BILSON_THOMPSON_IMPROVED_KEY] = real(Qimproved)
                push!(output_labels, "Qimproved_bilson_thompson")
            end
        else
            error("method $methodname is not supported in topological charge measurement")
        end
        #printstring *= "$(values[i]) "
    end
    for value in values
        printstring *= "$(value) "
    end
    printstring *= "#  " * join(output_labels, " ") * " "

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        measurestring = printstring
        println_verbose_level2(m.verbose_print, printstring)
        #println_verbose_level2(U[1],"-----------------")
    end

    output = Measurement_output(valuedic, measurestring)

    return output
end

function calculate_topological_charge_plaq(U::Array{T,1}, temp_UμνTA, temps) where {T}
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA, "plaq", U, temps)
    Q = calc_Q(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_clover(U::Array{T,1}, temp_UμνTA, temps) where {T}
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA, "clover", U, temps)
    Q = calc_Q(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_improved(
    U::Array{T,1},
    temp_UμνTA,
    Qclover,
    temps,
) where {T}
    # Alexandrou, Athenodorou, and Jansen, Phys. Rev. D 92, 125014
    # (2015), Eqs. (15)-(17): improve the composite density directly.
    rectangle_loops = calc_UμνTA!(temp_UμνTA, "rect", U, temps)
    Qrectangle = 2 * calc_Q(temp_UμνTA, rectangle_loops, U)
    return (5 / 3) * Qclover - (1 / 12) * Qrectangle
end

function calculate_topological_charge_bilson_thompson(
    U::Array{T,1},
    temp_UμνTA,
    Qclover,
    temps,
) where {T}
    # Bilson-Thompson-style construction: improve the field-strength tensor
    # first and then insert it into the quadratic topological density.  The
    # 1x1+1x2 tree-level coefficients here are also those used by
    # SIMULATeQCD's FieldStrengthTensor_imp.
    #
    #   F_imp = (5/3) F_1x1 - (1/3) F_1x2,
    #
    # where the rectangular clover carries twice the area normalization.
    # In terms of the raw traceless anti-Hermitian loop sums used here this is
    # A_imp = (5/3) A_plaq - (1/6) A_rect.
    clover_rectangle_cross = calculate_topological_charge_cross_rect_clover(
        U, temp_UμνTA, temps,
    )
    rectangle_loops = calc_UμνTA!(temp_UμνTA, "rect", U, temps)
    Qrectangle = calc_Q(temp_UμνTA, rectangle_loops, U)

    clover_coefficient = 5 / 3
    rectangle_coefficient = -1 / 6
    return clover_coefficient^2 * Qclover +
           2 * clover_coefficient * rectangle_coefficient *
           clover_rectangle_cross +
           4 * rectangle_coefficient^2 * Qrectangle
end

function calculate_topological_charge_cross_rect_clover(
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    clover_UμνTA,
    temps_g,
) where {NC,Dim}
    Dim == 4 || error("Dimension $Dim is not supported")
    rectangle_loops, _ = calc_loopset_μν_name("rect", Dim)
    path_temps, path_temp_indices = get_temp(temps_g, 4)
    evaluated_loop, evaluated_loop_index = get_temp(temps_g)
    rectangle_UμνTA, rectangle_index = get_temp(temps_g)
    cross_term = 0.0

    for μ in 1:Dim
        for ν in 1:Dim
            μ == ν && continue
            evaluate_gaugelinks!(
                evaluated_loop, rectangle_loops[μ, ν], U, path_temps,
            )
            Traceless_antihermitian!(rectangle_UμνTA, evaluated_loop)
            for ρ in 1:Dim
                for σ in 1:Dim
                    ρ == σ && continue
                    cross_term += epsilon_tensor(μ, ν, ρ, σ) *
                                  tr(rectangle_UμνTA, clover_UμνTA[ρ, σ])
                end
            end
        end
    end

    unused!(temps_g, path_temp_indices)
    unused!(temps_g, evaluated_loop_index)
    unused!(temps_g, rectangle_index)
    return -cross_term / (32 * π^2 * 4^2)
end

function calc_UμνTA!(
    temp_UμνTA,
    name::String,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    loops_μν, numofloops = calc_loopset_μν_name(name, Dim)
    calc_UμνTA!(temp_UμνTA, loops_μν, U, temps)
    return numofloops
end


function calc_UμνTA!(
    temp_UμνTA,
    loops_μν,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps_g,
) where {NC,Dim}
    UμνTA = temp_UμνTA
    temps, its_temps = get_temp(temps_g, 4)
    temp1, it_temp1 = get_temp(temps_g)
    for μ = 1:Dim
        for ν = 1:Dim
            if ν == μ
                continue
            end

            evaluate_gaugelinks!(temp1, loops_μν[μ, ν], U, temps)
            Traceless_antihermitian!(UμνTA[μ, ν], temp1)
            #loopset = Loops(U,loops_μν[μ,ν])
            #UμνTA[μ,ν] = evaluate_loops(loopset,U)

            #UμνTA[μ,ν] = Traceless_antihermitian(UμνTA[μ,ν])
        end
    end
    unused!(temps_g, its_temps)
    unused!(temps_g, it_temp1)
    return
end


#=
implementation of topological charge is based on
https://arxiv.org/abs/1509.04259
=#
function calc_Q(UμνTA, numofloops, U::Array{<:AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
    Q = 0.0
    if Dim == 4
        ε(μ, ν, ρ, σ) = epsilon_tensor(μ, ν, ρ, σ)
    else
        error("Dimension $Dim is not supported")
    end
    for μ = 1:Dim
        for ν = 1:Dim
            if ν == μ
                continue
            end
            Uμν = UμνTA[μ, ν]
            for ρ = 1:Dim
                for σ = 1:Dim
                    if ρ == σ
                        continue
                    end
                    Uρσ = UμνTA[ρ, σ]
                    s = tr(Uμν, Uρσ)
                    Q += ε(μ, ν, ρ, σ) * s / numofloops^2
                end
            end
        end
    end

    return -Q / (32 * (π^2))
end




#topological charge
function epsilon_tensor(mu::Int, nu::Int, rho::Int, sigma::Int)
    sign = 1 # (3) 1710.09474 extended epsilon tensor
    if mu < 0
        sign *= -1
        mu = -mu
    end
    if nu < 0
        sign *= -1
        nu = -nu
    end
    if rho < 0
        sign *= -1
        rho = -rho
    end
    if sigma < 0
        sign *= -1
        sigma = -sigma
    end
    epsilon = zeros(Int, 4, 4, 4, 4)
    epsilon[1, 2, 3, 4] = 1
    epsilon[1, 2, 4, 3] = -1
    epsilon[1, 3, 2, 4] = -1
    epsilon[1, 3, 4, 2] = 1
    epsilon[1, 4, 2, 3] = 1
    epsilon[1, 4, 3, 2] = -1
    epsilon[2, 1, 3, 4] = -1
    epsilon[2, 1, 4, 3] = 1
    epsilon[2, 3, 1, 4] = 1
    epsilon[2, 3, 4, 1] = -1
    epsilon[2, 4, 1, 3] = -1
    epsilon[2, 4, 3, 1] = 1
    epsilon[3, 1, 2, 4] = 1
    epsilon[3, 1, 4, 2] = -1
    epsilon[3, 2, 1, 4] = -1
    epsilon[3, 2, 4, 1] = 1
    epsilon[3, 4, 1, 2] = 1
    epsilon[3, 4, 2, 1] = -1
    epsilon[4, 1, 2, 3] = -1
    epsilon[4, 1, 3, 2] = 1
    epsilon[4, 2, 1, 3] = 1
    epsilon[4, 2, 3, 1] = -1
    epsilon[4, 3, 1, 2] = -1
    epsilon[4, 3, 2, 1] = 1
    return epsilon[mu, nu, rho, sigma] * sign
end



function calc_loopset_μν_name(name, Dim)
    loops_μν = Array{Vector{Wilsonline{Dim}},2}(undef, Dim, Dim)
    if name == "plaq"
        numofloops = 1
        for μ = 1:Dim
            for ν = 1:Dim
                loops_μν[μ, ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                plaq = make_plaq(μ, ν, Dim=Dim)
                push!(loops_μν[μ, ν], plaq)
            end
        end
    elseif name == "clover"
        numofloops = 4
        for μ = 1:Dim
            for ν = 1:Dim
                loops_μν[μ, ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                loops_μν[μ, ν] = make_cloverloops_topo(μ, ν, Dim=Dim)
            end
        end
    elseif name == "rect"
        numofloops = 8
        for μ = 1:4
            for ν = 1:4
                if ν == μ
                    continue
                end
                loops = Wilsonline{Dim}[]
                loop_righttop = Wilsonline([(μ, 2), (ν, 1), (μ, -2), (ν, -1)])
                loop_lefttop = Wilsonline([(ν, 1), (μ, -2), (ν, -1), (μ, 2)])
                loop_rightbottom = Wilsonline([(ν, -1), (μ, 2), (ν, 1), (μ, -2)])
                loop_leftbottom = Wilsonline([(μ, -2), (ν, -1), (μ, 2), (ν, 1)])
                push!(loops, loop_righttop)
                push!(loops, loop_lefttop)
                push!(loops, loop_rightbottom)
                push!(loops, loop_leftbottom)

                loop_righttop = Wilsonline([(μ, 1), (ν, 2), (μ, -1), (ν, -2)])
                loop_lefttop = Wilsonline([(ν, 2), (μ, -1), (ν, -2), (μ, 1)])
                loop_rightbottom = Wilsonline([(ν, -2), (μ, 1), (ν, 2), (μ, -1)])
                loop_leftbottom = Wilsonline([(μ, -1), (ν, -2), (μ, 1), (ν, 2)])
                push!(loops, loop_righttop)
                push!(loops, loop_lefttop)
                push!(loops, loop_rightbottom)
                push!(loops, loop_leftbottom)

                loops_μν[μ, ν] = loops
            end
        end
    else
        error("$name is not supported")
    end
    return loops_μν, numofloops
end


function make_cloverloops_topo(μ, ν; Dim=4)
    loops = Wilsonline{Dim}[]
    loop_righttop = Wilsonline([(μ, 1), (ν, 1), (μ, -1), (ν, -1)])
    loop_lefttop = Wilsonline([(ν, 1), (μ, -1), (ν, -1), (μ, 1)])
    loop_rightbottom = Wilsonline([(ν, -1), (μ, 1), (ν, 1), (μ, -1)])
    loop_leftbottom = Wilsonline([(μ, -1), (ν, -1), (μ, 1), (ν, 1)])
    push!(loops, loop_righttop)
    push!(loops, loop_lefttop)
    push!(loops, loop_rightbottom)
    push!(loops, loop_leftbottom)
    return loops
end
