mutable struct Topological_charge_density_correlation_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    #_temporary_gaugefields::Vector{TG}
    _temporary_matrices::Vector{Matrix{ComplexF64}}
    temp_UμνTA::Matrix{Matrix{ComplexF64}}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    TC_methods::Vector{String}
    improved_topological_charge_definition::String

    function Topological_charge_density_correlation_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        TC_methods=["plaquette"],
        improved_topological_charge_definition="alexandrou",
    ) where {T}
        myrank = get_myrank(U)
        NC = U[1].NC
        improved_definition = _validate_improved_topological_charge_definition(
            improved_topological_charge_definition,
        )

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        temp_UμνTA = Array{Matrix{ComplexF64},2}(undef, Dim, Dim)


        for μ = 1:Dim
            for ν = 1:Dim
                temp_UμνTA[ν, μ] = zeros(ComplexF64, NC, NC)
            end
        end

        numg = 5
        _temporary_matrices = Vector{Matrix{ComplexF64}}(undef, numg)
        for i = 1:numg
            _temporary_matrices[i] = zeros(ComplexF64, NC, NC)
        end




        return new{Dim,T}(
            filename,
            _temporary_matrices,
            temp_UμνTA,
            Dim,
            verbose_print,
            printvalues,
            TC_methods,
            improved_definition,
        )

    end



end

function Topological_charge_density_correlation_measurement(
    U::Vector{T},
    params::TopologicalChargeDensityCorrelation_parameters,
    filename="Topological_charge_density_correlation.txt",
) where {T}

    return Topological_charge_density_correlation_measurement(
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
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    loop1position,
    relativeposition;
    additional_string="",
) where {M<:Topological_charge_density_correlation_measurement,NC,Dim}
    temps = m._temporary_matrices
    measurestring = ""

    nummethod = length(m.TC_methods)
    values = Float64[]
    valuedic = Dict{String,Float64}()
    output_labels = String[]
    printstring = " " * additional_string
    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            Qplaq1 = calculate_topological_charge_plaq(U, loop1position,
                m.temp_UμνTA, temps)
            Qplaq2 = calculate_topological_charge_plaq(U, loop1position .+ relativeposition,
                m.temp_UμνTA, temps)
            QQ = Qplaq1 * Qplaq2
            push!(values, real(QQ))
            valuedic["plaquette"] = real(QQ)
            push!(output_labels, "Qplaq")
        elseif methodname == "clover"
            Qclover1 = calculate_topological_charge_clover(U, loop1position,
                m.temp_UμνTA, temps)
            Qclover2 = calculate_topological_charge_clover(U, loop1position .+ relativeposition,
                m.temp_UμνTA, temps)
            QQ = Qclover1 * Qclover2
            push!(values, real(QQ))
            valuedic["clover"] = real(QQ)
            push!(output_labels, "Qclover")
            definition = m.improved_topological_charge_definition
            if definition in ("alexandrou", "both")
                Qimproved1 = calculate_topological_charge_improved(
                    U, loop1position, m.temp_UμνTA, Qclover1, temps,
                )
                Qimproved2 = calculate_topological_charge_improved(
                    U, loop1position .+ relativeposition,
                    m.temp_UμνTA, Qclover2, temps,
                )
                QQ = Qimproved1 * Qimproved2
                push!(values, real(QQ))
                valuedic[_ALEXANDROU_IMPROVED_KEY] = real(QQ)
                push!(output_labels, "Qimproved_alexandrou")
            end
            if definition in ("bilson_thompson", "both")
                # The field-strength construction needs the clover tensor at
                # the same site as each rectangle tensor, so evaluate each
                # clover/improved pair consecutively.
                Qclover1 = calculate_topological_charge_clover(
                    U, loop1position, m.temp_UμνTA, temps,
                )
                Qimproved1 = calculate_topological_charge_bilson_thompson(
                    U, loop1position, m.temp_UμνTA, Qclover1, temps,
                )
                second_position = loop1position .+ relativeposition
                Qclover2 = calculate_topological_charge_clover(
                    U, second_position, m.temp_UμνTA, temps,
                )
                Qimproved2 = calculate_topological_charge_bilson_thompson(
                    U, second_position, m.temp_UμνTA, Qclover2, temps,
                )
                QQ = Qimproved1 * Qimproved2
                push!(values, real(QQ))
                valuedic[_BILSON_THOMPSON_IMPROVED_KEY] = real(QQ)
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

function calculate_topological_charge_plaq(U::Array{T,1}, position,
    temp_UμνTA,
    temps) where {T}

    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA,
        position,
        "plaq", U, temps)
    Q = calc_Q_each(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_clover(U::Array{T,1}, position,
    temp_UμνTA, temps) where {T}
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA,
        position,
        "clover", U, temps)
    Q = calc_Q_each(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_improved(
    U::Array{T,1},
    position,
    temp_UμνTA,
    Qclover,
    temps,
) where {T}
    rectangle_loops = calc_UμνTA!(
        temp_UμνTA, position, "rect", U, temps,
    )
    Qrectangle = 2 * calc_Q_each(temp_UμνTA, rectangle_loops, U)
    return (5 / 3) * Qclover - (1 / 12) * Qrectangle
end

function calculate_topological_charge_bilson_thompson(
    U::Array{T,1},
    position,
    temp_UμνTA,
    Qclover,
    temps,
) where {T}
    clover_rectangle_cross =
        calculate_topological_charge_cross_rect_clover_eachsite(
            U, position, temp_UμνTA, temps,
        )
    rectangle_loops = calc_UμνTA!(
        temp_UμνTA, position, "rect", U, temps,
    )
    Qrectangle = calc_Q_each(temp_UμνTA, rectangle_loops, U)

    clover_coefficient = 5 / 3
    rectangle_coefficient = -1 / 6
    return clover_coefficient^2 * Qclover +
           2 * clover_coefficient * rectangle_coefficient *
           clover_rectangle_cross +
           4 * rectangle_coefficient^2 * Qrectangle
end

function calculate_topological_charge_cross_rect_clover_eachsite(
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    position,
    clover_UμνTA,
    temps,
) where {NC,Dim}
    Dim == 4 || error("Dimension $Dim is not supported")
    rectangle_loops, _ = calc_loopset_μν_name("rect", Dim)
    indices = Tuple(position)
    evaluated_loop = temps[1]
    rectangle_UμνTA = temps[2]
    cross_term = 0.0

    for μ in 1:Dim
        for ν in 1:Dim
            μ == ν && continue
            evaluate_gaugelinks_eachsite!(
                evaluated_loop,
                rectangle_loops[μ, ν],
                U,
                temps[2:end],
                indices...,
            )
            _traceless_antihermitian_matrix!(rectangle_UμνTA, evaluated_loop)
            for ρ in 1:Dim
                for σ in 1:Dim
                    ρ == σ && continue
                    cross_term += epsilon_tensor(μ, ν, ρ, σ) *
                                  tr(rectangle_UμνTA * clover_UμνTA[ρ, σ])
                end
            end
        end
    end

    return -cross_term / (32 * π^2 * 4^2)
end

function _traceless_antihermitian_matrix!(output, input)
    NC = size(input, 1)
    diagonal_imaginary_part = sum(imag(input[k, k]) for k in 1:NC) / NC
    for k in 1:NC
        output[k, k] = (imag(input[k, k]) - diagonal_imaginary_part) * im
    end
    for column in 1:NC
        for row in (column + 1):NC
            value = 0.5 * (input[column, row] - conj(input[row, column]))
            output[column, row] = value
            output[row, column] = -conj(value)
        end
    end
    return output
end

function calc_UμνTA!(
    temp_UμνTA,
    position,
    name::String,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    loops_μν, numofloops = calc_loopset_μν_name(name, Dim)
    calc_UμνTA!(temp_UμνTA, position,
        loops_μν, U, temps)
    return numofloops
end


function calc_UμνTA!(
    temp_UμνTA,
    position,
    loops_μν,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    indices = Tuple(position)#(1, 1, 1, 1)
    fac1N = 1 / NC


    UμνTA = temp_UμνTA
    for μ = 1:Dim
        for ν = 1:Dim
            if ν == μ
                continue
            end
            V1 = temps[1]

            evaluate_gaugelinks_eachsite!(
                V1,
                loops_μν[μ, ν],
                U,
                temps[2:end],
                indices...,
            )
            tri = 0.0
            @simd for k = 1:NC
                tri += imag(V1[k, k])
            end
            tri *= fac1N
            @simd for k = 1:NC
                UμνTA[μ, ν][k, k] =
                    (imag(V1[k, k]) - tri) * im
            end

            for k1 = 1:NC
                @simd for k2 = k1+1:NC
                    vv =
                        0.5 * (
                            V1[k1, k2] -
                            conj(V1[k2, k1])
                        )
                    UμνTA[μ, ν][k1, k2] = vv
                    UμνTA[μ, ν][k2, k1] = -conj(vv)
                end
            end
        end
    end
    return
end

function calc_Q_each(UμνTA, numofloops, U::Array{<:AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
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
                    s = tr(Uμν * Uρσ)
                    Q += ε(μ, ν, ρ, σ) * s / numofloops^2
                end
            end
        end
    end

    return -Q / (32 * (π^2))
end
