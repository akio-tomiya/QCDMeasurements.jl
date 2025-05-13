

mutable struct Correlation_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    temp_g1g2::Vector{TG}
    Dim::Int8
    loop1::Vector{Wilsonline{Dim}}
    loop2::Vector{Wilsonline{Dim}}
    relativeposition::Vector{Int64}
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    originonly::Bool
    loop1position::Vector{Int64}

    function Correlation_measurement(
        U::Vector{T},
        loop1,
        loop2,
        relativeposition, ;
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        originonly=true,
        loop1position=Int64[1, 1, 1, 1]
    ) where {T}
        myrank = get_myrank(U)

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        temp_g1g2 = Array{T,1}(undef, 2)
        for i = 1:2
            temp_g1g2[i] = similar(U[1])
        end

        numg = 4
        _temporary_gaugefields = Temporalfields(U[1], num=numg)

        #_temporary_gaugefields = Vector{T}(undef, numg)
        #_temporary_gaugefields[1] = similar(U[1])
        #for i = 2:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end

        wloop1 = []
        for loop in loop1
            push!(wloop1, Wilsonline(loop))
        end
        wloop2 = []
        for loop in loop2
            push!(wloop2, Wilsonline(loop))
        end
        #println("originonly = $originonly")

        return new{Dim,T}(
            filename,
            _temporary_gaugefields,
            temp_g1g2,
            Dim,
            wloop1,#Wilsonline(loop1),
            wloop2,#Wilsonline(loop2),
            relativeposition,
            verbose_print,
            printvalues,
            originonly,
            loop1position
        )

    end

end

function Correlation_measurement(U, params, filename="Correlation.txt")

    return Correlation_measurement(
        U,
        params.loop1,
        params.loop2,
        params.relativeposition;
        filename=filename,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
        originonly=params.originonly,
        loop1position=params.loop1position
    )
end


function measure(
    m::M,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string="",
) where {M<:Correlation_measurement,NC,Dim}


    function site_trace_product!(B, A)
        c = tr(B) * tr(A)
        B .= c / NC
    end

    if m.originonly
        #if Dim == 4
        indices = Tuple(m.loop1position)#(1, 1, 1, 1)
        indices2 = Tuple(m.loop1position .+ m.relativeposition)
        #value = NC * g1[1, 1, 1, 1, 1, 1]
        #elseif Dim == 2
        #    indices = (1, 1)
        #    indices2 = (1 + m.position[1], 1 + m.relativeposition[2])
        #    #value = NC * g1[1, 1, 1, 1]
        #end
        mat_temps = Matrix{ComplexF64}[]
        for i = 1:5
            push!(mat_temps, zeros(ComplexF64, NC, NC))
        end

        V1 = zeros(ComplexF64, NC, NC)
        V2 = zeros(ComplexF64, NC, NC)
        Gaugefields.AbstractGaugefields_module.evaluate_gaugelinks_eachsite!(
            V1,
            m.loop1,
            U,
            mat_temps,
            indices...,
        )
        Gaugefields.AbstractGaugefields_module.evaluate_gaugelinks_eachsite!(
            V2,
            m.loop2,
            U,
            mat_temps,
            indices2...,
        )
        value = tr(V1) * tr(V2)
        #println(value)


        #if Dim == 4
        #    value = NC * g1[1, 1, 1, 1, 1, 1]
        #elseif Dim == 2
        #    value = NC * g1[1, 1, 1, 1]
        #end
        #println(value)
    else
        #temps = get_temporary_gaugefields(m)
        temps, its_temps = get_temp(m._temporary_gaugefields, 4)
        g1 = m.temp_g1g2[1]
        g2 = m.temp_g1g2[2]

        evaluate_gaugelinks!(g1, m.loop1, U, temps)
        evaluate_gaugelinks!(g2, m.loop2, U, temps)

        g2shifted = Gaugefields.AbstractGaugefields_module.shift_U(g2, Tuple(m.relativeposition))
        substitute_U!(temps[1], g2shifted)

        map_U!(
            g1,
            site_trace_product!,
            temps[1])

        value = tr(g1)
        unused!(m._temporary_gaugefields, its_temps)
    end

    measurestring = ""

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        measurestring = "$additional_string $(real(value)) $(imag(value)) # two-point correlation"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end


    output = Measurement_output(value, measurestring)
    return output

end

#=
# = = = calc energy density = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
function calculate_energy_density(U::Array{T,1}, Wmat, temps) where {T<:AbstractGaugefields}
    # Making a ( Ls × Lt) Wilson loop operator for potential calculations
    WL = 0.0 + 0.0im
    NV = U[1].NV
    NC = U[1].NC
    #Wmat = Array{T,2}(undef,4,4)
    #
    make_energy_density!(Wmat, U, temps) # make wilon loop operator and evaluate as a field, not traced.
    WL = make_energy_density_core(Wmat, U, NV) # tracing over color and average over spacetime and x,y,z.
    NDir = 4.0 * 3.0 / 2 # choice of 2 axis from 4.
    return real(WL) / NV / NDir / NC / 8
end

function make_energy_density!(
    Wmat,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    W_operator, numofloops = calc_loopset_μν_name("clover", Dim)#make_Wilson_loop(Lt,Ls)
    calc_large_wilson_loop!(Wmat, W_operator, U, temps)
    return
end

function make_energy_density_core(
    Wmat::Array{<:AbstractGaugefields{NC,Dim},2},
    U::Array{T,1},
    NV,
) where {T<:AbstractGaugefields,NC,Dim}
    @assert Dim == 4

    W = 0.0 + 0.0im
    for μ = 1:Dim # all directions
        for ν = 1:Dim
            if μ == ν
                continue
            end
            W += tr(Wmat[μ, ν], Wmat[μ, ν]) / 4
        end
    end
    return W
end



function calc_large_wilson_loop!(
    temp_Wmat::Array{<:AbstractGaugefields{NC,Dim},2},
    W_operator,
    U::Array{T,1},
    temps,
) where {T<:AbstractGaugefields,NC,Dim}
    W = temp_Wmat
    for μ = 1:Dim
        for ν = 1:Dim
            if μ == ν
                continue
            end
            #println(typeof(μ)," ",typeof(ν))
            #exit()
            #loopset = Loops(U,W_operator[μ,ν])
            evaluate_gaugelinks!(W[μ, ν], W_operator[μ, ν], U, temps)
            #W[μ,ν] = evaluate_loops(loopset,U)
        end
    end
    return
end
=#