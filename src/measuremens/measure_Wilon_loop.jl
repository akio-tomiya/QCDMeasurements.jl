mutable struct Wilson_loop_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    _temporary_gaugefields_mat::Matrix{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    Tmax::Int64
    Rmax::Int64
    outputvalues::Matrix{Float64}
    wilsonloops::Matrix{Matrix{Vector{Wilsonline{Dim}}}}


    function Wilson_loop_measurement(
        U::Vector{TG};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        Tmax=4,
        Rmax=4,
    ) where {TG}
        myrank = get_myrank(U)

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)


        numg = 4
        _temporary_gaugefields = Temporalfields(U[1], num=numg)
        #_temporary_gaugefields = Vector{TG}(undef, numg)
        #for i = 1:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end

        _temporary_gaugefields_mat = Matrix{TG}(undef, Dim, Dim)
        for μ = 1:Dim
            for ν = 1:Dim
                _temporary_gaugefields_mat[ν, μ] = similar(U[1])
            end
        end
        wilsonloops = Matrix{Matrix{Vector{Wilsonline{Dim}}}}(undef, Tmax, Rmax)

        for T = 1:Tmax
            for R = 1:Rmax
                Wloops = make_Wilson_loop(T, R, Dim)
                wilsonloops[T, R] = Wloops
            end
        end

        outputvalues = zeros(Float64, Tmax, Rmax)

        return new{Dim,TG}(
            filename,
            _temporary_gaugefields,
            _temporary_gaugefields_mat,
            Dim,
            verbose_print,
            printvalues,
            Tmax,
            Rmax,
            outputvalues,
            wilsonloops,
        )

    end
end

function Wilson_loop_measurement(
    U::Vector{T},
    params::Wilson_loop_parameters,
    filename="Wilson_loop.txt",
) where {T}
    return Wilson_loop_measurement(
        U,
        filename=filename,
        Tmax=params.Tmax,
        Rmax=params.Rmax,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
    )
end




function measure(
    m::Wilson_loop_measurement{Dim,TG},
    U;
    additional_string="",
) where {Dim,TG}
    temps, its_temps = get_temp(m._temporary_gaugefields, 4)#get_temporary_gaugefields(m)
    NC, _, NN... = size(U[1])
    NV = prod(NN)
    measurestring = ""


    for T = 1:m.Tmax
        for R = 1:m.Rmax
            Wmat = m._temporary_gaugefields_mat
            Wloops = m.wilsonloops[T, R]
            calc_large_wiloson_loop!(Wmat, Wloops, U, temps, Dim)
            W = 0.0im
            for μ = 1:Dim-1 # spatial directions
                ν = Dim  # T-direction is not summed over
                W += tr(Wmat[μ, ν])
            end
            NDir = 3.0 # in 4 diemension, 3 associated staples. t-x plane, t-y plane, t-z plane
            WL = real(W) / NV / NDir / NC
            m.outputvalues[T, R] = WL

            if m.printvalues
                measurestring = " $additional_string $T $R $WL # Wilson_loops # T R W(T,R)"
                println_verbose_level2(m.verbose_print, measurestring)
            end
        end
    end
    unused!(temps, its_temps)


    output = Measurement_output(m.outputvalues, measurestring)

    return output
end

#=
function calc_Wilson_loop(U::Array{T,1},Lt,Ls) where T <: GaugeFields
    # Making a ( Ls × Lt) Wilson loop operator for potential calculations
    WL = 0.0+0.0im
    NV = U[1].NV
    NC = U[1].NC
    Wmat = Array{GaugeFields_1d,2}(undef,4,4)
    #
    calc_large_wiloson_loop!(Wmat,Lt,Ls,U) # make wilon loop operator and evaluate as a field, not traced.
    WL = calc_Wilson_loop_core(Wmat,U,NV) # tracing over color and average over spacetime and x,y,z.
    NDir = 3.0 # in 4 diemension, 3 associated staples. t-x plane, t-y plane, t-z plane
    return real(WL)/NV/NDir/NC
end
=#

#=
function calc_Wilson_loop_core(Wmat, U::Array{GaugeFields{S},1} ,NV) where S <: SUn
    if S == SU3
        NC = 3
    elseif S == SU2
        NC = 2
    else
        NC = U[1].NC
        #error("NC != 2,3 is not supported")
    end
    W = 0.0 + 0.0im
    for n=1:NV
        for μ=1:3 # spatial directions
            ν=4  # T-direction is not summed over
            W += tr(Wmat[μ,ν][:,:,n])
        end
    end
    return W
end
function calc_large_wiloson_loop!(Wmat,Lt,Ls,U)
    W_operator = make_Wilson_loop(Lt,Ls)
    calc_large_wiloson_loop!(Wmat,W_operator,U)
    return 
end
=#
function make_Wilson_loop(Lt, Ls, Dim)
    #= Making a Wilson loop operator for potential calculations
        Ls × Lt
        ν=4
        ↑
        +--+ 
        |  |
        |  |
        |  |
        +--+ → μ=1,2,3 (averaged)
    =#
    Wmatset = Array{Vector{Wilsonline{Dim}},2}(undef, 4, 4)
    for μ = 1:Dim-1 # spatial directions
        ν = Dim # T-direction is not summed over
        loops = Wilsonline{Dim}[]
        #loops = Wilson_loop_set()
        loop = Wilsonline([(μ, Ls), (ν, Lt), (μ, -Ls), (ν, -Lt)])
        push!(loops, loop)
        Wmatset[μ, ν] = loops
    end
    return Wmatset
end
function calc_large_wiloson_loop!(temp_Wmat, loops_μν, U, temps, Dim)
    W = temp_Wmat
    for μ = 1:Dim-1 # spatial directions
        ν = Dim # T-direction is not summed over
        evaluate_gaugelinks!(W[μ, ν], loops_μν[μ, ν], U, temps)
    end
    return
end
