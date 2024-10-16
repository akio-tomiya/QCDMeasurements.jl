mutable struct Guluonic_correlators_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    Dim::Int8
    loop1::Vector{Wilsonline{Dim}}
    loop2::Vector{Wilsonline{Dim}}
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    temporal_matrix::Vector{Matrix{ComplexF64}}

    function Guluonic_correlators_measurement(
        U::Vector{T},
        loop1,
        loop2;
        filename=nothing,
        verbose_level=2,
        printvalues=false
    ) where {T}
        myrank = get_myrank(U)
        NC = U[1].NC

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)



        wloop1 = []
        for loop in loop1
            push!(wloop1, Wilsonline(loop))
        end
        wloop2 = []
        for loop in loop2
            push!(wloop2, Wilsonline(loop))
        end

        temporal_matrix = Matrix{ComplexF64}[]
        for i = 1:5
            push!(temporal_matrix, zeros(ComplexF64, NC, NC))
        end

        return new{Dim,T}(
            filename,
            Dim,
            wloop1,#Wilsonline(loop1),
            wloop2,#Wilsonline(loop2),
            verbose_print,
            printvalues,
            temporal_matrix
        )

    end

end

function Guluonic_correlators_measurement(U::Vector{TG}, params::T, filename::String="Guluonic_correlators.txt") where {TG,T<:Guluonic_correlators_parameters}

    return Guluonic_correlators_measurement(
        U,
        params.loop1,
        params.loop2,
        filename=filename,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues
    )
end


function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    loop1position,
    relativeposition;
    additional_string="",
) where {M<:Guluonic_correlators_measurement,NC,Dim}


    indices = Tuple(loop1position)#(1, 1, 1, 1)
    indices2 = Tuple(loop1position .+ relativeposition)

    V1 = zeros(ComplexF64, NC, NC)
    V2 = zeros(ComplexF64, NC, NC)

    Gaugefields.AbstractGaugefields_module.evaluate_gaugelinks_eachsite!(
        V1,
        m.loop1,
        U,
        m.temporal_matrix,
        indices...,
    )
    Gaugefields.AbstractGaugefields_module.evaluate_gaugelinks_eachsite!(
        V2,
        m.loop2,
        U,
        m.temporal_matrix,
        indices2...,
    )
    value = tr(V1) * tr(V2)

    measurestring = ""

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        measurestring = "$additional_string $(real(value)) $(imag(value)) # Guluonic_correlators"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end


    output = Measurement_output(value, measurestring)
    return output

end

