# QCDMeasurements v1 audit and implementation status

Date: 2026-08-19

This is the current release-oriented inventory for the independent
`/home/nagai/JuliaQCD/QCDMeasurements-v1` work tree.  It records what is now
implemented and tested, as well as the remaining release gates.

## Current result

- QCDMeasurements is version `1.0.0` on branch `qcdmeasurements-v1`.
- It resolves the requested local v1 stack and loads without dependency-
  internal imports.
- The normal CPU suite passes **200/200** tests.
- The H100 CUDA parity suite passes **25/25** tests on an NVIDIA H100 NVL.
- The two-rank MPI suite passes **199/199 on each rank**.
- Gaugefields 1.0.3's complete package suite passes with LatticeMatrices 1.1.3,
  including High-level API, Enzyme MD, flow, HMC, and heatbath coverage.
- Shamir domain-wall `PP`, `J5qP`, and their pointwise residual-mass ratio now
  agree with the independent Grid reference on CPU, two MPI ranks, and H100.
- An isolated minimum-version environment using registered LatticeMatrices
  1.1.2 and registered LatticeDiracOperators 1.0.0 passes **200/200** tests.

These results establish a local compatibility baseline.  A remotely
reproducible CI result still depends on publishing or otherwise making the
exact Gaugefields and LatticeDiracOperators v1 revisions available to the CI
runner.

LatticeMatrices 1.1.3 and LatticeDiracOperators 1.0.0 are registered. The only
remaining dependency publication-order constraint introduced here is that
Gaugefields 1.0.3 must be published before QCDMeasurements 1.0 remote CI can
resolve its corrected lower bound.

## Local dependency snapshot

| Package | Version | Local path | Branch | State relevant to this audit |
|---|---:|---|---|---|
| QCDMeasurements | 1.0.0 | `/home/nagai/JuliaQCD/QCDMeasurements-v1` | `qcdmeasurements-v1` | Active v1 work tree |
| LatticeMatrices | 1.1.3 | `/home/nagai/JuliaQCD/LatticeMatrices-v1.1.0` | `release/v1.1.3` | Domain-wall and generalized contraction kernels |
| Gaugefields | 1.0.3 | `/home/nagai/JuliaQCD/Gaugefields-v1` | `gaugefields-v1` | Pkg/TOML empty-vector dispatch fix |
| LatticeDiracOperators | 1.0.0 | `/home/nagai/JuliaQCD/LatticeDiracOperators-v1` | `ldo-v1` | Pre-existing uncommitted v1 migration work |

Gaugefields 1.0.2 promoted the gauge-link evaluation and temporary-field
operations used by QCDMeasurements to top-level exports. Gaugefields 1.0.3
adds an exact `Vector{Union{}}` `similar` method so Pkg/TOML empty arrays do not
ambiguously match all gauge-field-vector methods. QCDMeasurements therefore
has the lower bound `Gaugefields = "1.0.3"`.

LatticeDiracOperators 1.0.0 requires LatticeMatrices 1.1.2 for the public
domain-wall physical-source API, so QCDMeasurements uses the matching lower
bound `LatticeMatrices = "1.1.2"`. The compatibility CI checks both 1.1.2 and
1.1.3.

The package `Manifest.toml` files are intentionally ignored.  The main,
`test/gpu`, and `test/mpi` local environments resolve the sibling work trees
through path sources.

## V1 public contract now in place

### Names and results

- Canonical exported measurement names use Julia type casing, for example
  `PlaquetteMeasurement`, `MesonCorrelatorMeasurement`, and
  `TopologicalChargeMeasurement`.
- Canonical parameter aliases are exported, such as `PlaquetteParameters` and
  `MesonCorrelatorParameters`.
- The original underscore-style names remain aliases for source compatibility.
- `AbstractMeasurement` and the canonical `MeasurementOutput{T}` are public.
  `get_value` returns the physics value and `get_string` returns the legacy
  printable representation.
- Measurements which existed but were not exported (gluonic correlation,
  topological-density correlation, eigenvalue, and `MdagM` spectrum) now have
  explicit public canonical names.

### Construction

- `prepare_measurement(U, configuration)` is the canonical configuration API;
  `prepare_measurement_from_dict` remains a compatibility alias.
- Configurations accept `AbstractDict` with string or symbol keys and values.
- Measurement and fermion selection use registries plus multiple dispatch
  instead of a single exact-type `if/elseif` factory.
- Missing method names, unknown options, and conversion failures raise clear
  `ArgumentError`s rather than being silently ignored.
- Wilson chiral-condensate configuration now exposes the required `hop` and
  `r` parameters and is covered by an exact test.

### Dependency boundary

QCDMeasurements imports the required functionality only from the public
top-level modules:

- `clear_fermion!` and `Z4_distribution_fermi!` from
  LatticeDiracOperators;
- `Temporalfields`, `get_temp`, `unused!`, `shift_U`, and
  `evaluate_gaugelinks_eachsite!` from Gaugefields; and
- `projected_bilinear_slices` and `set_global_component!` from
  LatticeMatrices 1.1.

The source tree contains no remaining imports or qualified calls through
`Gaugefields.AbstractGaugefields_module`, `Gaugefields.Temporalfields_module`,
or `LatticeDiracOperators.Dirac_operators`.

## Supported and validated measurements

| Area | Measurement and fermion support | Validation |
|---|---|---|
| Basic gauge | plaquette, Polyakov loop, Wilson loops, energy density | Exact cold limits and analytic nontrivial Wilson loop |
| Gauge correlations | loop, gluonic, topological-density correlations | Deterministic regression and exact cold checks |
| Topology | plaquette/clover charge; Alexandrou and Bilson--Thompson improved choices | Formula checks and frozen SIMULATeQCD v1.2.0 values |
| Wilson flow | per-configuration history; ensemble `t0` and `w0` estimator | Input/ensemble tests and frozen SIMULATeQCD flow values |
| Legacy pion | Wilson, Wilson-clover, and staggered | Direct contractions, LM kernel path, hot normal/even-odd agreement; Grid clover reference |
| Generalized meson | Wilson/Wilson-clover local S/P/V/A/T; staggered M1--M8 | Direct kernels, SIMULATeQCD, and Grid cold/non-cold references |
| PCAC mass | Wilson and Wilson-clover | Grid cold/non-cold correlators and final PCAC values |
| Domain-wall residual mass | Shamir `PP`, midpoint `J5qP`, and pointwise ratio | Grid non-cold reference on CPU, two MPI ranks, and H100 |
| Chiral condensate | Wilson and staggered | Exact small-field tests for both public construction paths |
| Spectral | Wilson eigenvalue and local `MdagM` spectrum | Dense exact small-lattice references |

The machine-readable support contract is `supported_fermions`:

| Measurement | Result |
|---|---|
| pion, generalized meson | `(:Wilson, :WilsonClover, :Staggered)` |
| PCAC | `(:Wilson, :WilsonClover)` |
| domain-wall residual mass | `(:Domainwall,)` |
| chiral condensate | `(:Wilson, :Staggered)` |
| eigenvalue, `MdagM` spectrum | `(:Wilson,)` |
| gauge-only measurements | `()` |

Wilson-clover pion, generalized meson, and PCAC are validated on the
LatticeMatrices MPILattice backend; legacy serial clover is not in the v1
contract.

`DomainWallResidualMassMeasurement` composes LDO's twelve physical point-source
solves with LM's physical-surface and midpoint projections. It returns a typed
result containing `PP`, `J5qP`, and `m_res(t)=J5qP(t)/PP(t)`. A plateau fit is
an ensemble analysis and is intentionally not performed per configuration.
The compiled Grid driver, frozen output, and public-QCDMeasurements regression
are stored under `test/references/grid` and `test/domainwall_residual_mass.jl`.

## Test matrix and provenance

### CPU

```text
QCDMeasurements.jl | 200 passed / 200 total | 1m28.9s
```

This uses Julia 1.11, the Threads JACC backend, and one MPI rank.

### H100 CUDA

```text
QCDMeasurements H100 CUDA parity | 25 passed / 25 total | 1m33.7s
```

The run used an NVIDIA H100 NVL (compute capability 9.0).  It verifies CUDA
storage and cold plaquette, Polyakov loop, energy density, topology, Wilson
pion, and generalized Wilson pseudoscalar results.  It additionally constructs
the reconstructible non-cold one-link field with the LatticeMatrices JACC
setter, checks a Wilson-flow trajectory against the CPU result, and compares
Wilson and Wilson-clover `PP`, `AP`, and final PCAC mass directly with frozen
Grid values. It also compares Shamir domain-wall `PP`, `J5qP`, and the
residual-mass ratio with Grid.

### MPI

```text
rank 0: QCDMeasurements.jl | 199 passed / 199 total | 1m27.4s
rank 1: QCDMeasurements.jl | 199 passed / 199 total | 1m27.4s
```

This runs the normal deterministic suite with two MPI ranks, including the
distributed LatticeMatrices contraction path and the Wilson-clover `PP`,
`A_4P`, and final PCAC comparison with the same global Grid reference field,
plus the distributed Shamir domain-wall Grid comparison.
The per-rank count differs from the single-rank count because one assertion is
conditional on the lattice decomposition.

The repository records independent Grid values for Wilson and Wilson-clover
generalized-meson and PCAC observables on reconstructible non-cold fields (and
cold Wilson fields). It records
SIMULATeQCD v1.2.0 values for staggered M1--M8 channels, gradient flow, energy
density, and topology.  Commits, conventions, normalizations, raw values, and
sign choices are stored beside the relevant tests.

## Remaining work in priority order

### P0: release blockers

1. Commit/publish Gaugefields 1.0.3, then make the remote Julia 1.11/1.12, OS,
   LM 1.1.2/1.1.3, and two-rank MPI CI matrix authoritative. Local path
   dependencies cannot prove remote CI reproducibility.
2. Decide whether every public measurement needs a dedicated typed value
   struct.  `MeasurementOutput{T}` is stable, but some `value`s are still
   dictionaries, tuples, or vectors.
3. Separate remaining verbose/file formatting side effects from numerical
   computation where legacy implementations still interleave them.

### P1: broaden validated backends

1. Extend H100 parity to the full non-Abelian SIMULATeQCD flow fixture and add
   an accelerator allocation/performance guard; the current non-cold flow and
   PCAC numerical paths are covered.
2. Retain frozen Grid/SIMULATeQCD inputs and values as regression fixtures.

### P2: maintenance

1. Rename legacy misspelled source filenames and implementation types only
   after compatibility/deprecation policy is settled.
2. Relax unnecessary concrete container annotations to suitable abstract
   interfaces.
3. Split the large legacy parameter and pion implementations into smaller
   setup, solver, and measurement units.

## Proposed v1 release gate

Tag QCDMeasurements 1.0.0 only when:

- dependency lower bounds resolve to published reproducible revisions;
- remote Julia 1.11/1.12 and two-rank MPI CI is green;
- README examples use the canonical v1 API;
- every advertised measurement/fermion/backend combination has a recorded
  exact or independent reference; and
- the chosen public value/result contracts are documented and frozen.
