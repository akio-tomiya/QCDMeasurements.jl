# Test layout

`runtests.jl` contains deterministic, assertion-based tests that are safe to
run in CI. The tests use small lattices and do not write generated results into
the package source tree.

- `basic_measurements.jl` checks gluonic observables on a cold field.
- `gauge_density.jl` checks the four public topological density definitions on
  unit links and deterministic non-Abelian, noncubic SU(2)/SU(3) fields. It
  compares sums with the scalar measurements and selected sites with the
  existing pointwise loop implementation. Action-density averages are checked
  against the gradient-flow observables. It also tests workspace reuse, input
  immutability, CPU backend parity, and argument validation.
- `gauge_density_backends.jl` checks CPU, GPU, MPI, and GPU+MPI calculations
  against `references/gauge_density_v1.toml`, generated using the unmodified
  pre-change commit `5576729`. It compares all 360 sites of every density on
  noncubic SU(2)/SU(3) fields, scalar charges and action averages, both allocating
  and reusable APIs, input immutability, and reuse on unit links. Multiple MPI
  decompositions and legacy MPI rank ordering are covered. Input links are
  checked against a SHA-256 digest of the baseline configuration.
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

## Density API validation (2026-10-21)

The density change was validated on Linux with Julia 1.11.8,
Gaugefields 1.1.11, LatticeMatrices 1.2.8, LatticeDiracOperators 1.1.3,
JACC 1.4.0, MPI 0.20.27, and CUDA.jl 6.4.0.
CPU and GPU use the same JACC kernels; MPI collects the
local scalar fields into the complete host array on every rank.

| Run | Result |
|---|---|
| Complete CPU suite, 4 Julia threads, one MPI rank | 575/575 passed |
| Complete CPU suite, 2 MPI ranks, 2 Julia threads per rank | 641/641 passed on each rank |
| Density regression, CPU, 4 MPI ranks | 221/221 passed on each rank |
| Complete GPU suite, NVIDIA H100 NVL | 136/136 passed |
| Density regression, H100, 2 MPI ranks sharing one GPU | 441/441 passed on each rank |
| Density regression, RTX PRO 6000 Blackwell, 4 MPI ranks sharing one GPU | 221/221 passed on each rank |
| Density regression, 2 MPI ranks on separate H100 and Blackwell GPUs | 441/441 passed on each rank |
| Independent pre-change/post-change comparison without MPI loaded | 48/48 passed |

The baseline is the unmodified commit
`5576729f81c07955df01bee1f13d1ae12017f2e1`, run in a separate environment with
the same dependencies and deterministic input links. The three cases are
noncubic `3×4×5×6` SU(3) LatticeMatrices, SU(3) legacy, and SU(2)
LatticeMatrices fields. An independent before/after run passed **48/48**
checks: the existing four scalar charges, both action averages, and historical
energy observable were bitwise unchanged at four Julia threads. Against the
old pointwise implementation, the largest absolute differences were
`2.78e-17` for topological density and `8.89e-16` for action density. Differences
in the density sums/averages arise from floating-point reduction order.

`references/gauge_density_v1.toml` freezes every lattice site, scalar
observables, and input SHA-256 digests. To regenerate it, activate an environment
pointing to the old commit (not this working tree), with the versions above:

```sh
JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 julia --project=/path/to/baseline-env \
    test/references/generate_gauge_density_v1.jl /tmp/gauge_density_v1.toml
```

The generator checks the old source's SHA-256 and the thread count. Its output
was independently regenerated and compared exactly with the frozen fixture.
Backend regressions use tolerances because device arithmetic and MPI/thread
reduction orders can differ. They compare every site, both allocating and
reusable calls, all scalar observables, input immutability, and reuse after
switching to unit links. The two-rank tests split Y and T separately; four-rank
tests split both Y and T (`1×2×1×2`), including the legacy rank permutation.

Run the standalone density regression with CPU or GPU preferences and any
supported rank count, for example:

```sh
julia --project=test/gpu -e '
    using MPI
    run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd())
        --startup-file=no --project=test/gpu test/run_gauge_density_backends.jl`)
'
```

For JACC 1.4, select `default_backend = "cuda"` and `backends = ["cuda"]`
in the environment's JACC preferences. The legacy `backend` preference alone
does not select CUDA in that version. A GPU UUID in `CUDA_VISIBLE_DEVICES`
unambiguously selects the device. GPU+MPI was tested both with ranks sharing
one visible GPU and with two ranks assigned to separate H100 and Blackwell
GPUs. The latter run confirmed rank 0 on H100 and rank 1 on Blackwell. These
runs use host-staged MPI communication.

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
