mutable struct Eigenvalue_measurement{Dim,TG,TD} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    D::TD
    nev::Int64
    which::Symbol


    function Eigenvalue_measurement(
        U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = false,
        fermiontype = "Wilson",
        nev = 10,
        which = :SM,
        mass = 0.1,
        Nf = 2,
        κ = 1,
        r = 1,
        L5 = 2,
        M = -1,
        eps_CG = 1e-14,
        MaxCGstep = 3000,
        BoundaryCondition = nothing,
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
            verbose_print = Verbose_print(verbose_level, myid = myrank, filename = filename)
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


        params,parameters_action,x,factor = make_fermionparameter_dict(U,
            fermiontype,mass,
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
        _temporary_gaugefields = Vector{T}(undef, numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i = 2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        return new{Dim,T,TD}(filename, _temporary_gaugefields, 
            Dim, verbose_print, printvalues,
            D,
            nev,which)

    end
end

#=
const eigenvalue_keys = ["verbose_level","printvalue","fermiontype","mass",
        "Nf",
        "κ",
        "r",
        "L5",
        "M",
]
=#

function Eigenvalue_measurement(
    U::Vector{T},
    params::Eigenvalue_parameters,
    filename = "Eigenvalues.txt",
) where {T}

    params_tuple = fermionparameter_params(params)
    method = Eigenvalue_measurement(
            U;
            filename = filename,
            nev = params.nev,
            which = params.which,
            params_tuple...
    )
    return method

end




function measure(m::M, U; additional_string = "",maxiter=3000) where {M<:Eigenvalue_measurement}
    #temps = get_temporary_gaugefields(m)
    #poly = calculate_Polyakov_loop(U, temps[1], temps[2])
    #println_verbose_level2(m.verbose_print,"constructing sparse matrix...")
    Du = m.D(U)
    Ds = construct_sparsematrix(Du)
    n,_ = size(Ds)
    #=
    e,v = eigen(Matrix(Ds))
    fp = open("testwilson.txt","w")
    for ei in e
        println(fp,"$(real(ei)) $(imag(ei))")
    end
    close(fp)
    display(Ds)
    =#
    vals, vectors = eigs(Ds, nev = m.nev, which=m.which, maxiter=maxiter);
    #println_verbose_level2(m.verbose_print,"done...")
    measurestring = ""

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        for i=1:length(vals)
            measurestring *= "$additional_string $i $(real(vals[i])) \t $(imag(vals[i])) \n"
        end
        measurestring *= "#eigenvectors\n"
        n1,n2 = size(vectors)
        for i=1:n2
            for j=1:n1
                v = vectors[j,i]
                measurestring *= "$additional_string $j $i $(real(v)) \t $(imag(v)) \n"
            end
        end
        #measurestring = "$additional_string $(real(poly)) $(imag(poly)) # poly"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end
    output = Measurement_output((vals, vectors), measurestring)

    return output
end
