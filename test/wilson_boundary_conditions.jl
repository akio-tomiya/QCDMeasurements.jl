module MeasurementBoundaryTests
using Test, LinearAlgebra, Gaugefields, LatticeDiracOperators, LatticeMatrices, QCDMeasurements, JACC
JACC.@init_backend
const snapshots = Dict{String,Vector{Float64}}()
# The package runner initializes MPI before including this file; standalone
# backend checks opt in through LDO_TEST_MPI.
const distributed = get(ENV,"LDO_TEST_MPI","false") == "true" ||
    (isdefined(Main, :MPI) && Main.MPI.Initialized())
if distributed
    using MPI
    MPI.Initialized() || MPI.Init()
end

# Free-field Fourier propagator: independent of the iterative solver and halos.
function free_pion(dims, bc; kappa=0.1, colors=3)
    any(iszero,bc) && return wall_pion(dims,bc;kappa,colors)
    gamma = LatticeDiracOperators.Dirac_operators.mk_gamma(1.0)[1]
    propagator = zeros(ComplexF64,4,4,dims...)
    for mode in CartesianIndices(dims)
        p = [(2pi*(mode[d]-1)+(bc[d]==1 ? 0 : pi))/dims[d] for d in 1:4]
        Dp = (1-2kappa*sum(cos.(p)))*Matrix{ComplexF64}(I,4,4)
        for d in 1:4
            Dp .+= 2im*kappa*sin(p[d]).*gamma[:,:,d]
        end
        inverse = inv(Dp)/prod(dims)
        for site in CartesianIndices(dims)
            propagator[:,:,Tuple(site)...] .+= exp(im*sum(p[d]*(site[d]-1) for d in 1:4)).*inverse
        end
    end
    [colors*sum(abs2,@view propagator[:,:,:,:,:,t]) for t in 1:dims[4]]
end

# Independent dense Dirac matrix for a free field with an absorbing wall.
function wall_pion(dims,bc;kappa=0.1,colors=3)
    gamma=LatticeDiracOperators.Dirac_operators.mk_gamma(1.0)[1]
    D=Matrix{ComplexF64}(I,4prod(dims),4prod(dims))
    indices=LinearIndices(dims)
    block(pos)=4*(indices[pos...]-1).+(1:4)
    for site in CartesianIndices(dims),mu in 1:4,sign in (-1,1)
        pos=Tuple(site)
        raw=pos[mu]+sign
        factor=1<=raw<=dims[mu] ? 1 : bc[mu]
        src=ntuple(d->d==mu ? mod1(raw,dims[d]) : pos[d],4)
        D[block(pos),block(src)].-=kappa*factor.*(I-sign*gamma[:,:,mu])
    end
    source=zeros(ComplexF64,4prod(dims),4);source[1:4,:].=Matrix{ComplexF64}(I,4,4)
    propagator=D\source
    [colors*sum(abs2,view(propagator,4prod(dims[1:3])*(t-1)+1:4prod(dims[1:3])*t,:)) for t in 1:dims[4]]
end

@testset "Wilson measurement boundary conditions" begin
    nprocs = distributed ? MPI.Comm_size(MPI.COMM_WORLD) : 1
    dims = (4,2,2,4)
    U = Initialize_Gaugefields(3,1,dims...;condition="cold",isMPILattice=true,
        PEs=(1,1,1,nprocs),verbose_level=0)
    for bc in ([1,1,1,-1],[1,1,1,1],[-1,1,1,-1],[1,1,1,0],[0,1,1,0])
        expected = free_pion(dims,bc)
        for family in ("Wilson","WilsonClover")
            kwargs = (;fermiontype=family,κ=0.1,cSW=1.17,
                BoundaryCondition=bc,eps_CG=1e-24,verbose_level=0)
            pion = PionCorrelatorMeasurement(U;kwargs...)
            @test pion.D.boundarycondition == bc
            @test pion._temporary_fermionfields[1].f.phases == bc
            actual = get_value(measure(pion,U))
            @test actual ≈ expected rtol=1e-10 atol=1e-12
            meson = MesonCorrelatorMeasurement(U;channels=["pseudoscalar"],kwargs...)
            observed = get_value(measure(meson,U))[:pseudoscalar][:,1]
            @test observed ≈ expected rtol=1e-10 atol=1e-12
            @test observed ≈ actual rtol=1e-11 atol=1e-12
            snapshots[family*string(bc)] = collect(real.(observed))
        end
    end
    @test snapshots["Wilson[1, 1, 1, 1]"] != snapshots["Wilson[1, 1, 1, -1]"]
end
end
