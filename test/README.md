# Test layout

`runtests.jl` contains deterministic, assertion-based tests that are safe to
run in CI. The tests use small lattices and do not write generated results into
the package source tree.

- `basic_measurements.jl` checks gluonic observables on a cold field.
- `legacy_measurements.jl` gives the previously example-only Wilson-loop,
  loop-correlation, gluonic-correlation, topological-density-correlation,
  chiral-condensate, eigenvalue, and local spectral-density measurements
  deterministic regression coverage.  It uses exact cold-field limits, an
  analytic nontrivial Wilson loop, and a dense-resolvent reference on a small
  lattice.
- `grid_meson_reference.jl` compares the generalized Wilson pseudoscalar
  correlator with independent [Grid](https://github.com/paboyle/Grid)
  calculations on both cold and non-cold `4^4` lattices. The compact
  non-cold field changes one link to `diag(exp(1.1im), exp(-1.1im), 1)`, so it
  can be reconstructed exactly in both codes without a binary fixture. The
  pinned Grid commit, raw values, plaquette, boundary condition, solver
  tolerance conversion, and Wilson-operator normalization are recorded in
  the test.
- `pcac_mass.jl` compares both Wilson `PP` and `A_4P` correlators with Grid
  on the same cold and localized non-cold fields, then checks the symmetric
  PCAC derivative and Grid's final `c_A=0` ratio, optional `c_A` term,
  dictionary interface, and validation.
- `wilson_clover_reference.jl` runs the LatticeMatrices Wilson-clover operator
  through QCDMeasurements' MPILattice API and compares `PP`, `A_4P`, and the
  final `c_A=0` PCAC ratio with Grid's `WilsonCloverFermionD` on the same
  reconstructible non-cold field at `cSW=1.2`. It also checks dictionary
  construction, solver validation, and rejection of the legacy serial backend.
  The same frozen global-field values are checked on one and two MPI ranks.
- `domainwall_residual_mass.jl` runs the public
  `DomainWallResidualMassMeasurement` on the reconstructible non-cold field
  and compares `PP`, midpoint `J5qP`, and `m_res(t)=J5qP(t)/PP(t)` with Grid.
  It also checks all twelve solver outcomes, the dictionary interface, and
  MPILattice/`L5` validation.
- `references/grid/qcdm_domainwall_mres_reference.cc` is the executable Grid
  reference used by that test. It uses
  `DomainWallFermionD` physical source/solution projection and `ContractJ5q`
  to record `PP`, midpoint `J5q`, and their timeslice ratio on the same
  reconstructible one-link field. The frozen values are in the adjacent
  `.txt` file.
- `gradient_flow_scale.jl` checks Wilson-flow histories, `t0`/`w0` ensemble
  reduction, input immutability, and dictionary construction.
- `simulateqcd_gradient_flow_reference.jl` compares plaquette, clover energy,
  clover charge, and the Bilson--Thompson-style improved-field-strength charge
  against frozen SIMULATeQCD v1.2.0 values for the same deterministic
  non-Abelian field and fixed-step Wilson flow. It also checks the three-way
  improved-topology selection and the Alexandrou formula directly.
- `dependency_compatibility.jl` checks the supported Gaugefields,
  LatticeDiracOperators, and LatticeMatrices version ranges.
- `pion_correlator.jl` compares the Wilson pion measurement with a direct
  color-spin contraction and guards the sink-index regression.
- `pion_solver_diagnostics.jl` checks solver outcomes, failure handling,
  streaming storage, and agreement of normal and even-odd solves.

`dicttest.jl` and `gauge.jl` are retained as legacy, long-running integration
examples. They write output files and are not included by `runtests.jl`.

## H100/CUDA validation

`gpu/runtests.jl` is an opt-in accelerator test. It uses the sibling local v1
work trees through `gpu/Project.toml` and selects JACC's CUDA backend through
`gpu/LocalPreferences.toml`. It checks that the gauge storage is a `CuArray`,
then compares cold-field gauge observables and Wilson pion/generalized-meson
correlators with the same frozen CPU/Grid-normalized values used by the normal
suite. It also constructs the non-cold Grid field with the LatticeMatrices
JACC setter and checks a Wilson-flow trajectory against the CPU baseline and
the Wilson and Wilson-clover final PCAC results against Grid's frozen values.
Finally it runs the Shamir domain-wall physical/midpoint measurement and checks
`PP`, `J5qP`, and `m_res(t)` against the independent Grid output.

Run it on one GPU with:

```sh
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=0 \
    julia --project=test/gpu test/gpu/runtests.jl
```

## Local two-rank validation

`mpi/Project.toml` is the CPU/MPI counterpart of the GPU environment. It also
resolves the four sibling v1 work trees and runs the normal suite on every
rank, including the distributed LatticeMatrices pion/domain-wall contractions
and the Wilson-clover/domain-wall Grid references. From the repository root:

```sh
julia --project=test/mpi -e '
    using MPI
    run(`$(MPI.mpiexec()) -n 2 $(Base.julia_cmd()) \
        --startup-file=no --project=test/mpi test/mpi/runtests.jl`)
'
```

## SIMULATeQCD reference provenance

The frozen gradient-flow values were produced with SIMULATeQCD v1.2.0 at git
commit `767a1b1` on an NVIDIA H100.  The deterministic field construction is in
`simulateqcd_gradient_flow_reference.jl`; it was exported as an ILDG field and
run with the Wilson force, fixed-step RK3, `start_step_size = 0.01`, and required
flow times `0 0.01 0.02 0.03 0.04 0.05`.  SIMULATeQCD and QCDMeasurements use
opposite orientations for the sign of the topological charge, so the reference
test reverses that sign explicitly.  Plaquette, clover action density, and the
magnitudes of clover and improved topological charge are otherwise compared
directly.
