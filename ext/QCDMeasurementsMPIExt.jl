module QCDMeasurementsMPIExt

using MPI, LatticeMatrices
import QCDMeasurements: _density_legacy_mpi_result

function _density_legacy_mpi_result(link, lattice_size)
    # Legacy gauge fields number ranks with x fastest; MPI Cartesian grids
    # number ranks with the last coordinate fastest. Reorder the communicator
    # so the scalar lattice's physical block matches the input's local block.
    key = 0
    for d in 1:4
        key = key * link.PEs[d] + link.myrank_xyzt[d]
    end
    ordered = MPI.Comm_split(link.comm, 0, key)
    return LatticeMatrix(1, 1, 4, lattice_size, link.PEs;
        nw=0, elementtype=Float64, comm0=ordered, device_mapping=:current)
end

end
