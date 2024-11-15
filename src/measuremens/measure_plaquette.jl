"""
Struct for measuring the plaquette values. 

# Constructor

    Plaquette_measurement(
        U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = false,
    )

## Arguments
- `U<:AbstractGaugefields`: The gauge fields 
- `filename::Union{String,Nothing}=nothing`: The filename of the output data if ```printvalues = true```.
- `verbose_level<:Int = 2` :Verbose level for printing
- `printvalues::Bool = false` : The values are printed during the measurement if ```printvalue = true```.


## Examples
```julia
plaq_m = Plaquette_measurement(U,filename="test.txt",printvalue=true)
```

# Functions
    measure(m::Plaquette_measurement, U; additional_string = "")

```additional_string``` is a header for the printed values if ```m.printvalue = true```. 

"""
mutable struct Plaquette_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    Dim::Int8
    factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool

    function Plaquette_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
    ) where {T}
        myrank = get_myrank(U)
        #=
        if U[1].mpi == false
            myrank = 0
        else
            myrank = U[1].myrank
        end
        =#
        if printvalues
            verbose_print = Verbose_print(verbose_level, myid=myrank, filename=filename)
            #println(verbose_print)
            #error("v")
        else
            verbose_print = nothing
        end
        Dim = length(U)

        if Dim == 4
            comb = 6
        elseif Dim == 2
            comb = 1
        else
            error("Dim = $Dim is not supported in Plaquette_measurement")
        end
        factor = 1 / (comb * U[1].NV * U[1].NC)

        numg = 2
        _temporary_gaugefields = Temporalfields(U[1], num=numg)
        #_temporary_gaugefields = Vector{T}(undef, numg)
        #_temporary_gaugefields[1] = similar(U[1])
        #for i = 2:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end

        return new{Dim,T}(
            filename,
            _temporary_gaugefields,
            Dim,
            factor,
            verbose_print,
            printvalues,
        )

    end
end


function Plaquette_measurement(
    U::Vector{T},
    params::Plaq_parameters,
    filename="Plaquette.txt",
) where {T}
    return Plaquette_measurement(
        U,
        filename=filename,
        verbose_level=params.verbose_level,
        printvalues=params.printvalues,
    )
end






function measure(m::M, U; additional_string="") where {M<:Plaquette_measurement}
    temps = m._temporary_gaugefields#get_temporary_gaugefields(m)
    temp1, it_temp1 = get_temp(temps)
    temp2, it_temp2 = get_temp(temp2)
    plaq = real(calculate_Plaquette(U, temp1, temp2) * m.factor)
    unused!(temps, it_temp1)
    unused!(temps, it_temp2)

    measurestring = ""

    if m.printvalues
        #println("m.verbose_print ",m.verbose_print)
        #println_verbose_level2(U[1],"-----------------")
        measurestring = "$additional_string $plaq # plaq"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end
    output = Measurement_output(plaq, measurestring)

    return output
end
