using RSCG
mutable struct MdagMspectrum_measurement{Dim,TG,TD} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Temporalfields{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    σ::Vector{ComplexF64}
    position::Tuple{Int64,Int64}



    function MdagMspectrum_measurement(
        U::Vector{T};
        filename=nothing,
        verbose_level=2,
        printvalues=false,
        fermiontype="Wilson",
        mass=0.1,
        Nf=2,
        κ=1,
        r=1,
        L5=2,
        M=-1,
        eps_CG=1e-14,
        MaxCGstep=3000,
        BoundaryCondition=nothing,
        emin=-3,
        emax=3,
        eta=0.01,
        numdatapoints=3000,
        position=(1, 1)
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
        else
            verbose_print = nothing
        end
        Dim = length(U)
        if BoundaryCondition == nothing
            if Dim == 4
                boundarycondition = BoundaryCondition_4D_default
            elseif Dim == 2
                boundarycondition = BoundaryCondition_2D_default
            end
        else
            boundarycondition = BoundaryCondition
        end


        params, parameters_action, x, factor = make_fermionparameter_dict(U,
            fermiontype, mass,
            Nf,
            κ,
            r,
            L5,
            M,
        )

        params["eps_CG"] = eps_CG
        params["verbose_level"] = verbose_level
        params["MaxCGstep"] = MaxCGstep
        params["boundarycondition"] = boundarycondition

        D = Dirac_operator(U, x, params)
        TD = typeof(D)



        numg = 2
        _temporary_gaugefields = Temporalfields(U[1], num=numg)
        #_temporary_gaugefields = Vector{T}(undef, numg)
        #_temporary_gaugefields[1] = similar(U[1])
        #for i = 2:numg
        #    _temporary_gaugefields[i] = similar(U[1])
        #end

        σ = zeros(ComplexF64, numdatapoints)
        σmin = emin + im * eta
        σmax = emax + im * eta
        for i = 1:numdatapoints
            σ[i] = (i - 1) * (σmax - σmin) / (numdatapoints - 1) + σmin
        end

        return new{Dim,T,TD}(filename, _temporary_gaugefields,
            Dim, verbose_print, printvalues,
            D,
            σ,
            position)

    end
end

#=
const  MdagMspectrum_keys = ["verbose_level","printvalue","fermiontype","mass",
        "Nf",
        "κ",
        "r",
        "L5",
        "M",
]
=#

function MdagMspectrum_measurement(
    U::Vector{T},
    params::MdagMspectrum_parameters,
    filename="MdagMspectrums.txt",
) where {T}

    params_tuple = fermionparameter_params(params)

    method = MdagMspectrum_measurement(
        U;
        filename=filename,
        emin=params.emin,
        emax=params.emax,
        eta=params.eta,
        numdatapoints=params.numpoints,
        position=params.position,
        params_tuple...
    )
    return method

end




function measure(m::M, U; additional_string="", maxiter=3000) where {M<:MdagMspectrum_measurement}
    #temps = get_temporary_gaugefields(m)
    #poly = calculate_Polyakov_loop(U, temps[1], temps[2])
    #println_verbose_level2(m.verbose_print,"constructing sparse matrix...")
    Du = m.D(U)
    Ds = construct_sparsematrix(Du)
    DdagD = Ds' * Ds
    n, _ = size(DdagD)
    i, j = m.position
    Gij1 = greensfunctions(i, j, m.σ, DdagD)
    ρ = imag.(Gij1) / (-1 / (π))
    vals = ρ

    #println_verbose_level2(m.verbose_print,"done...")
    measurestring = ""

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        for i = 1:length(vals)
            measurestring *= "$additional_string $(real.(m.σ)[i]) $(vals[i]) \n "
        end

        #measurestring = "$additional_string $(real(poly)) $(imag(poly)) # poly"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end
    output = Measurement_output((real.(m.σ), vals), measurestring)

    return output
end
