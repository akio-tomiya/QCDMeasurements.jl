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

    function Topological_charge_density_correlation_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        TC_methods=["plaquette"],
    ) where {T}
        myrank = get_myrank(U)
        NC = U[1].NC

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
    temp1 = temps[1]
    temp2 = temps[2]
    measurestring = ""

    nummethod = length(m.TC_methods)
    values = Float64[]
    valuedic = Dict{String,Float64}()
    printstring = " " * additional_string
    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            Qplaq1 = calculate_topological_charge_plaq(U, loop1position,
                m.temp_UμνTA, temps)
            Qplaq2 = calculate_topological_charge_plaq(U, loop1position .+ relativeposition,
                m.temp_UμνTA, temps)
            QQ = Qplaq1 * Qplaq2
            push!(values, QQ)
            valuedic["plaquette"] = QQ
        elseif methodname == "clover"
            Qclover1 = calculate_topological_charge_clover(U, loop1position,
                m.temp_UμνTA, temps)
            Qclover2 = calculate_topological_charge_clover(U, loop1position .+ relativeposition,
                m.temp_UμνTA, temps)
            QQ = Qclover1 * Qclover2
            push!(values, QQ)
            valuedic["clover"] = QQ
            Qimproved1 =
                calculate_topological_charge_improved(U, loop1position,
                    m.temp_UμνTA, Qclover1, temps)
            Qimproved2 =
                calculate_topological_charge_improved(U, loop1position .+ relativeposition,
                    m.temp_UμνTA, Qclover2, temps)
            QQ = Qimproved1 * Qimproved2
            push!(values, QQ)
            valuedic["clover improved"] = QQ
        else
            error("method $methodname is not supported in topological charge measurement")
        end
        #printstring *= "$(values[i]) "
    end
    for value in values
        printstring *= "$(value) "
    end
    printstring *= "#  "

    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            printstring *= "Qplaq "
        elseif methodname == "clover"
            printstring *= "Qclover Qimproved "
        else
            error("method $methodname is not supported in topological charge measurement")
        end
    end

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
    UμνTA = temp_UμνTA
    #numofloops = calc_UμνTA!(UμνTA,"clover",U)
    #Qclover = calc_Q(UμνTA,numofloops,U)

    numofloops = calc_UμνTA!(UμνTA,
        position,
        "rect", U, temps)

    #numofloops = calc_UμνTA!(UμνTA, "rect", U, temps)

    Qrect = 2 * calc_Q_each(UμνTA, numofloops, U)
    c1 = -1 / 12
    c0 = 5 / 3
    Q = c0 * Qclover + c1 * Qrect
    return Q
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

            Gaugefields.AbstractGaugefields_module.evaluate_gaugelinks_eachsite!(
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


