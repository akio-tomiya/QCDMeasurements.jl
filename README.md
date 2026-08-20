# QCDMeasurements.jl

QCDMeasurements.jl calculates lattice-QCD observables from gauge
configurations provided by
[Gaugefields.jl](https://github.com/akio-tomiya/Gaugefields.jl). Fermionic
measurements use
[LatticeDiracOperators.jl](https://github.com/akio-tomiya/LatticeDiracOperators.jl)
and [LatticeMatrices.jl](https://github.com/JuliaQCD/LatticeMatrices.jl).

<img src="LQCDjl_block.png" width=300> 

## Breaking changes in v1.0

Compared with QCDMeasurements v0.3, v1.0 has the following compatibility
changes:

- The supported dependency stack is now Gaugefields 1.0.3 or later within
  major version 1, LatticeDiracOperators 1.x, and LatticeMatrices 1.1.2 or
  later within major version 1. Environments pinned to Gaugefields 0.7,
  LatticeDiracOperators 0.6, or an earlier LatticeMatrices release must update
  those packages together.
- Measurement configuration dictionaries are now validated strictly. Missing
  `methodname` entries, unknown option names, invalid value conversions, and
  unsupported fermion/measurement combinations raise `ArgumentError` instead
  of allowing an incomplete or silently ignored configuration.

The underscore-style measurement names and `prepare_measurement_from_dict`
remain compatibility aliases, so existing callers can migrate gradually.
Gaugefields' `LegacyBackend()` has not been removed; the examples below use
the Gaugefields v1 default LatticeMatrices/JACC backend.

## Contents

- [Breaking changes in v1.0](#breaking-changes-in-v10)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Supported measurements](#supported-measurements)
- [Gauge observables](#gauge-observables)
- [Generalized local meson correlators](#generalized-local-meson-correlators)
- [Independent validation](#independent-validation)
- [Shamir domain-wall residual mass](#shamir-domain-wall-residual-mass)
- [Wilson PCAC mass](#wilson-pcac-mass)
- [Gradient-flow scales](#gradient-flow-scales-t0-and-w0)
- [Improved topological-charge definitions](#improved-topological-charge-definitions)
- [Citation](#citation)

## Requirements

QCDMeasurements 1.0 requires Gaugefields 1.0.3, LatticeDiracOperators 1.x, and
LatticeMatrices 1.1.2 or later. Gaugefields 1.0.2 introduced the public
top-level gauge-link evaluation and temporary-field APIs used here, while
1.0.3 fixes an ambiguity between gauge-field-vector `similar` methods and
Pkg/TOML's empty `Vector{Union{}}`. LatticeMatrices 1.1.2 introduced the
physical-surface and midpoint kernels required by the domain-wall measurement.
QCDMeasurements no longer imports dependency-internal modules.
Gaugefields v1's default `gauge_configuration` backend uses JACC, so add JACC
to the application environment and initialize the selected backend once at
startup.

## Quick start

New code should use the canonical v1 names.  The original underscore-style
names remain compatibility aliases during migration.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

U = gauge_configuration(
    (4, 4, 4, 4);
    colors=3,
    start=:cold,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)

measurement = PlaquetteMeasurement(U; printvalues=false)
output = measure(measurement, U)       # MeasurementOutput
plaquette = get_value(output)          # numerical or typed physics result
text = get_string(output)              # empty here because printvalues=false

# String or Symbol keys and values are accepted by the configuration API.
configured = prepare_measurement(U, Dict(
    :methodname => :Plaquette,
    :printvalues => false,
))
```

Unknown configuration keys and unsupported fermion combinations raise an
`ArgumentError`.  `supported_fermions(MeasurementType)` reports the fermion
formulations covered by the v1 validation contract.

## Supported measurements

| Measurement | Fermion formulation |
|---|---|
| `PionCorrelatorMeasurement` | Wilson, Wilson-clover, staggered |
| `MesonCorrelatorMeasurement` | Wilson, Wilson-clover, staggered |
| `PCACMassMeasurement` | Wilson, Wilson-clover |
| `DomainWallResidualMassMeasurement` | Shamir domain-wall |
| `ChiralCondensateMeasurement` | Wilson, staggered |
| `EigenvalueMeasurement`, `MdagMSpectrumMeasurement` | Wilson |
| `PlaquetteMeasurement` | — (no fermion operator is used) |
| `PolyakovMeasurement` | — (no fermion operator is used) |
| `EnergyDensityMeasurement` | — (no fermion operator is used) |
| `TopologicalChargeMeasurement` | — (no fermion operator is used) |
| `TopologicalChargeDensityCorrelationMeasurement` | — (no fermion operator is used) |
| `GradientFlowScaleMeasurement` | — (no fermion operator is used) |
| `WilsonLoopMeasurement` | — (no fermion operator is used) |
| `CorrelationMeasurement`, `GluonicCorrelatorMeasurement` | — (no fermion operator is used) |

Gauge-only measurements include the plaquette, Polyakov loop, Wilson loops,
energy density, topological charge and its density correlation, and
gradient-flow histories and scales. Fermionic measurements additionally
include pion and generalized local-meson correlators, the Wilson PCAC mass,
the Shamir domain-wall residual mass, the chiral condensate, and spectra.

Wilson-clover pion, generalized local-meson, and PCAC measurements are part of
the v1 contract with the default LatticeMatrices/JACC backend returned by
`gauge_configuration`. This selects LatticeMatrices'
`WilsonDiracCloverOperator4D`; pass the clover coefficient as `cSW`.

`DomainWallResidualMassMeasurement` uses LatticeDiracOperators' Shamir
physical-source solve and LatticeMatrices' surface/midpoint contractions. It
returns `PP`, `J5qP`, and the timeslice estimator
`m_res(t) = J5qP(t)/PP(t)`. The legacy `PionCorrelatorMeasurement` remains a
Wilson/Wilson-clover/staggered API; domain-wall physics is exposed through the
dedicated typed measurement instead.

## Gauge observables

The following measurements use only the input gauge configuration. This block
continues from the `U` constructed in [Quick start](#quick-start).

```julia
plaquette = get_value(measure(
    PlaquetteMeasurement(U; printvalues=false), U))

polyakov_loop = get_value(measure(
    PolyakovMeasurement(U; printvalues=false), U))

energy_density = get_value(measure(
    EnergyDensityMeasurement(U; printvalues=false), U))

topological_charge = get_value(measure(
    TopologicalChargeMeasurement(
        U; TC_methods=["plaquette", "clover"], printvalues=false),
    U,
))

wilson_loops = get_value(measure(
    WilsonLoopMeasurement(
        U; Tmax=2, Rmax=2, printvalues=false),
    U,
))
```

On a cold SU(3) field, `plaquette == 1`, `polyakov_loop == 3 + 0im`, and all
reported topological charges vanish. The later
[topological-charge section](#improved-topological-charge-definitions) shows
how to select either improved definition or calculate both.

## Generalized local meson correlators

`MesonCorrelatorMeasurement` computes connected point-to-all correlators
from one input gauge configuration.  For Wilson fermions it implements the
local scalar, pseudoscalar, vector, axial-vector, and tensor bilinears.  The
connected contraction convention is

```math
C_{\Gamma_2\Gamma_1}(\Delta,p) =
-\sum_{x:\,x_a-x_{0,a}=\Delta} e^{-ip\cdot(x-x_0)}
\mathrm{tr}_{c,s}
[\Gamma_2 S_2(x,x_0)\bar\Gamma_1 S_1(x_0,x)],
```

where `a` is `correlation_axis` and
`bar(Γ) = gamma[4] * adjoint(Γ) * gamma[4]`.  Gamma5 hermiticity is used so
that both propagators can be solved forward from the same point source.  Only
the connected Wick contraction is included; flavor-singlet disconnected
diagrams are not.

More explicitly, this new measurement provides:

- connected local meson two-point functions from a point source;
- momentum projection in the three directions transverse to
  `correlation_axis`;
- arbitrary source position and any of the four correlation axes;
- for Wilson fermions, all 16 local spin bilinears: scalar, pseudoscalar,
  four vectors, four axial vectors, and six antisymmetric tensors;
- for staggered fermions, the SIMULATeQCD-compatible M1--M8 local phase
  channels; and
- equal-mass or non-degenerate valence quarks, selected with `κ2` for Wilson
  fermions or `mass2` for staggered fermions.

For each requested channel, the returned value is a matrix indexed by
separation along `correlation_axis` and by the requested momentum. It is a
configuration-by-configuration observable: ensemble averages and subsequent
mass/effective-mass fits are intentionally left to the analysis stage.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

U = gauge_configuration(
    (4, 4, 4, 4);               # small runnable example
    colors=3,
    start=:cold,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)

measurement = MesonCorrelatorMeasurement(
    U;
    fermiontype="Wilson",
    κ=0.141139,
    channels=["pseudoscalar", "vector_1", "axial_1", "tensor_12"],
    momenta=[(0, 0, 0), (1, 0, 0)],
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
)
result = get_value(measure(measurement, U))

# Each channel is an axis_length × number_of_momenta matrix.
pion_at_rest = result[:pseudoscalar][:, 1]
vector_px = result[:vector_1][:, 2]
diagnostics = get_solver_diagnostics(measurement)
```

Use `mass2` (staggered) or `κ2` (Wilson) for a non-degenerate second quark.
Passing a `LocalMesonChannel` permits a custom pair of 4 by 4 sink/source spin
matrices.  `standard_local_meson_channels()` returns the complete set of 16
local S/P/V/A/T matrices.

For staggered fermions, the default channels are the M1--M8 phase projections
used by SIMULATeQCD's `measureHadrons` module.  The phase coordinates are
ordered transverse to the chosen correlation axis.

```julia
staggered_measurement = MesonCorrelatorMeasurement(
    U;
    fermiontype="Staggered",
    mass=0.05,
    channels=["M1", "M2", "M6"],
    correlation_axis=4,
)
staggered_result = get_value(measure(staggered_measurement, U))
scalar = staggered_result[:M1_scalar][:, 1]
pseudoscalar = staggered_result[:M2_pseudoscalar][:, 1]
```

The dictionary interface is also available:

```julia
measurement = prepare_measurement(
    U,
    Dict(
        "methodname" => "Meson_correlator",
        "fermiontype" => "Wilson",
        "hop" => 0.141139,
        "channels" => ["pseudoscalar", "vector_1"],
        "momenta" => [[0, 0, 0]],
        "correlation_axis" => 4,
        "source_position" => [1, 1, 1, 1],
        "printvalues" => false,
    ),
)
```

The accelerator/distributed contraction methods live in LatticeMatrices 1.1.2,
which is a direct QCDMeasurements dependency. QCDMeasurements imports its
JACC/MPI implementation directly; it performs one `MPI.Allreduce` per projected
correlator. QCDMeasurements supports LatticeMatrices 1.1.2 and later.

`gauge_configuration` selects the LatticeMatrices backend by default. The
same constructor is used for CPU, GPU, and MPI execution; select the JACC
backend before constructing `U` and adjust `process_grid` for MPI jobs. The
measurement API itself is unchanged:

```julia
U = gauge_configuration(
    (4, 4, 4, 4);               # small runnable example
    colors=3,
    halo=1,
    start=:cold,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)
measurement = MesonCorrelatorMeasurement(
    U; fermiontype="Wilson", κ=0.141139,
    channels=["pseudoscalar", "vector_1"])
result = get_value(measure(measurement, U))
```

Wilson-clover uses the same Gaugefields v1 interface and an explicit `cSW`:

```julia
U = gauge_configuration(
    (4, 4, 4, 4);               # small runnable example
    colors=3,
    halo=1,
    start=:cold,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)
clover_measurement = MesonCorrelatorMeasurement(
    U;
    fermiontype="WilsonClover",
    κ=0.141139,
    cSW=1.2,
    channels=["pseudoscalar", "axial_4"],
    BoundaryCondition=[1, 1, 1, -1],
)
clover_result = get_value(measure(clover_measurement, U))
```

The Gaugefields v1 default keeps CPU, MPI, and GPU execution on the same
LatticeMatrices/JACC clover operator.

## Independent validation

The v1 numerical contract uses frozen outputs from independent codes in
addition to internal unit tests. The executable drivers, raw values, and
normalization notes are kept under [`test/references`](test/references), while
[`test/README.md`](test/README.md) describes how each comparison is exercised.

| Observable | Independent reference | Fields and execution paths |
|---|---|---|
| Wilson `PP` and `A_4P` | Grid | cold, deterministic hot, reconstructible non-cold; CPU/MPI/H100 |
| Wilson-clover `PP`, `A_4P`, PCAC | Grid; Lattice-Tool-Kit source cross-check | non-cold field; CPU/MPI/H100 |
| Shamir domain-wall `PP`, `J5qP`, `m_res` | Grid `ContractJ5q` | non-cold field; CPU/MPI/H100 |
| Staggered M1--M8 phases | SIMULATeQCD and published spin-taste assignments | source-level convention tests |

### Wilson meson and PCAC

The quantity measured in the Grid comparison was the connected Wilson
zero-momentum pseudoscalar two-point function (the pion channel),

```math
C_{PP}(t,\mathbf{0}) =
\sum_{\mathbf{x}} \mathrm{tr}_{c,s}
[S(\mathbf{x},t;0) S(\mathbf{x},t;0)^\dagger],
```

from a point source at the origin. This checks the new measurement's Wilson
pseudoscalar channel, its point-source solve, spin-color contraction,
zero-momentum projection, time-slice reduction, and GPU/MPI reduction path.
The other Wilson S/V/A/T channels and nonzero momenta are covered by internal
contraction tests. The PCAC measurement described below additionally has an
independent Grid comparison for the local `A_4P` cross channel.

The result was compared with an independent calculation using
[Grid](https://github.com/paboyle/Grid) at commit
[`0ac72cb6a30ccdc41d664e7e0759f0c8833078f1`](https://github.com/paboyle/Grid/commit/0ac72cb6a30ccdc41d664e7e0759f0c8833078f1).
The Grid contractions follow channels 0 (`PP`) and 3 (`A_4P`) of its
[`Example_wall_wall_spectrum.cc`](https://github.com/paboyle/Grid/blob/0ac72cb6a30ccdc41d664e7e0759f0c8833078f1/examples/Example_wall_wall_spectrum.cc).
The comparison used a `4^4` lattice, a point source at the origin,
`κ = 0.1`, zero momentum, and periodic spatial/antiperiodic temporal fermion
boundary conditions. Grid's mass was set to
`m = 1/(2κ) - 4 = 1`.

Grid uses `D_Grid = (4 + m) - H/2`, while QCDMeasurements/LDO uses
`D_QCDM = 1 - κH`. Therefore the correlators are related by
`C_QCDM = C_Grid/(4κ^2)`; this normalization was applied before comparison.

Both cold and non-cold fields were checked. For the full hot-field check, Grid
generated a deterministic SU(3) hot configuration with fixed seed
`{45, 12, 81, 9}`. Its complete double-precision link field was loaded into
the QCDMeasurements LatticeMatrices/JACC path, so the two programs did not use
separate random configurations. The plaquette agreed to `8.7e-19` absolute
precision.
On an NVIDIA H100 NVL, the legacy `PionCorrelatorMeasurement` and the new
`MesonCorrelatorMeasurement` produced identical results, and both agreed
with Grid with maximum absolute and relative differences of `1.1e-14` and
`7.9e-15`, respectively.

To keep the automated test compact, the non-cold CI reference starts from a
cold field and changes one link to
`diag(exp(1.1im), exp(-1.1im), 1)`. This exactly reconstructible field checks
both its Grid plaquette and meson correlator without storing a large binary
gauge fixture. The pinned raw values and complete provenance are recorded in
[`test/grid_meson_reference.jl`](test/grid_meson_reference.jl).
The corresponding frozen `PP`/`A_4P` values used by the PCAC measurement are
in [`test/pcac_mass.jl`](test/pcac_mass.jl).

### Wilson-clover

Wilson-clover was checked separately on the same reconstructible non-cold
field, where the clover field strength is nonzero. The comparison used
`cSW=1.2`; QCDMeasurements' LatticeMatrices/JACC results for `PP`, `A_4P`, and
the final unimproved PCAC ratio agree with Grid's `WilsonCloverFermionD` results to the
recorded `2e-11` relative tolerance on one CPU rank, the same field split over
two MPI ranks, and an NVIDIA H100 NVL. The raw Grid values and normalization
are frozen in
[`test/wilson_clover_reference.jl`](test/wilson_clover_reference.jl), and the
independent Grid driver is
[`test/references/grid/qcdm_wilson_clover_reference.cc`](test/references/grid/qcdm_wilson_clover_reference.cc).

The adjacent
[Lattice-Tool-Kit](https://github.com/cometscome/Lattice-Tool-Kit/tree/1bfddc93da3f29e2e8811b65abf459715f6f98dc/FERMIONwithClover)
was also checked as an independent Wilson-clover implementation candidate:
its `FERMIONwithClover/propa.f` performs the same 12 point-source solves and
the connected pion spin-color norm (with an additional spatial-volume
average). The automated v1 regression uses Grid because its deterministic
link field can be reconstructed exactly in QCDMeasurements without a file
format conversion.

### Shamir domain-wall residual mass

The domain-wall comparison uses Grid's
`DomainWallFermionD::ImportPhysicalFermionSource`,
`ExportPhysicalFermionSolution`, and `ContractJ5q` as the independent
implementation of the physical surface and mid-plane quantities in
`m_res(t)=C_{J_{5q}P}(t)/C_{PP}(t)`.

A reproducible Grid driver and its frozen output are already recorded in
[`test/references/grid/qcdm_domainwall_mres_reference.cc`](test/references/grid/qcdm_domainwall_mres_reference.cc)
and
[`test/references/grid/qcdm_domainwall_mres_reference.txt`](test/references/grid/qcdm_domainwall_mres_reference.txt).
They use the same reconstructible non-cold one-link `4^4` field, Shamir mass
`0.1`, Grid `M5=1` (the LM/LDO convention `M=-1`), and `Ls=4`.
QCDMeasurements agrees with the frozen Grid `PP`, `J5qP`, and pointwise ratio
at the recorded relative tolerances (`2e-11`, `2e-11`, and `3e-11`). The
comparison is exercised on one CPU rank, the two-rank MPI decomposition, and
the H100 CUDA path.

The midpoint term and axial Ward identity follow V. Furman and Y. Shamir,
*Nucl. Phys. B* **439**, 54--78 (1995),
[DOI: 10.1016/0550-3213(95)00031-M](https://doi.org/10.1016/0550-3213(95)00031-M),
[arXiv:hep-lat/9405004](https://arxiv.org/abs/hep-lat/9405004). Grid's pinned
[`ContractJ5q` implementation](https://github.com/paboyle/Grid/blob/0ac72cb6a30ccdc41d664e7e0759f0c8833078f1/Grid/qcd/action/fermion/implementation/CayleyFermion5DImplementation.h#L587-L625)
uses the two central fifth-direction slices in this construction.

Grid's solver tolerance is a residual norm, whereas LDO's `eps_CG` is compared
with the squared residual. Thus the Grid setting `1e-14` corresponds to
`eps_CG=1e-28` for this precision-level comparison.

### Channel conventions and literature

The spin-color contraction follows the standard meson two-point construction
discussed by T. DeGrand and S. Schaefer, *Comput. Phys. Commun.* **159**,
185--191 (2004),
[DOI: 10.1016/j.cpc.2004.02.006](https://doi.org/10.1016/j.cpc.2004.02.006),
[arXiv:hep-lat/0401011](https://arxiv.org/abs/hep-lat/0401011).
The clover action is the on-shell improved Wilson action of B. Sheikholeslami
and R. Wohlert, *Nucl. Phys. B* **259**, 572--596 (1985),
[DOI: 10.1016/0550-3213(85)90002-1](https://doi.org/10.1016/0550-3213(85)90002-1).
The staggered channel set is checked against the M1--M8 phases in the
[SIMULATeQCD `measureHadrons` source](https://github.com/LatticeQCD/SIMULATeQCD/tree/main/src/modules/measureHadrons),
whose code and HISQ measurement capabilities are described by L. Mazur et al.,
*Comput. Phys. Commun.* **300**, 109164 (2024),
[DOI: 10.1016/j.cpc.2024.109164](https://doi.org/10.1016/j.cpc.2024.109164),
[arXiv:2306.01098](https://arxiv.org/abs/2306.01098).
For the local staggered screening channels and their spin-taste assignments,
see A. Bazavov et al., *Phys. Rev. D* **100**, 094510 (2019),
[DOI: 10.1103/PhysRevD.100.094510](https://doi.org/10.1103/PhysRevD.100.094510),
[arXiv:1908.09552](https://arxiv.org/abs/1908.09552).

## Shamir domain-wall residual mass

`DomainWallResidualMassMeasurement` consumes one Gaugefields v1 configuration
using the default LatticeMatrices/JACC backend and solves the twelve physical
point-source columns. Its typed
result exposes `result[:PP]`, `result[:J5qP]`, and `result[:mres]`. The last is
a timeslice estimator; selecting a plateau and fitting an ensemble are not
performed configuration by configuration.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

U = gauge_configuration(
    (4, 4, 4, 4);               # small runnable example
    colors=3,
    halo=1,
    start=:cold,                 # replace with the configuration to measure
    process_grid=(1, 1, 1, 1),
    verbose=0,
)

measurement = DomainWallResidualMassMeasurement(
    U;
    mass=0.1,
    L5=4,
    M=-1,
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
    momentum=(0, 0, 0, 0),
    BoundaryCondition=[1, 1, 1, -1],
    eps_CG=1e-14,
    MaxCGstep=5000,
    method_CG="bicg",
)
result = get_value(measure(measurement, U))

PP = result[:PP]
J5qP = result[:J5qP]
mres_t = result[:mres]
diagnostics = get_solver_diagnostics(measurement)
```

The same measurement can be constructed from a configuration dictionary;
`mass`/`L5` are accepted as aliases of the historical fermion-parameter fields
`m`/`N5`:

```julia
measurement = prepare_measurement(U, Dict(
    "methodname" => "Domainwall_residual_mass",
    "mass" => 0.1,
    "L5" => 4,
    "M" => -1.0,
    "printvalues" => false,
))
```

The definition follows the axial Ward identity of V. Furman and Y. Shamir,
*Nucl. Phys. B* **439**, 54--78 (1995),
[DOI: 10.1016/0550-3213(95)00031-M](https://doi.org/10.1016/0550-3213(95)00031-M),
and the implementation is checked against Grid's pinned
[`ContractJ5q`](https://github.com/paboyle/Grid/blob/0ac72cb6a30ccdc41d664e7e0759f0c8833078f1/Grid/qcd/action/fermion/implementation/CayleyFermion5DImplementation.h#L587-L625)
calculation described above.

## Wilson PCAC mass

`PCACMassMeasurement` is a configuration-by-configuration observable built
from the generalized Wilson/Wilson-clover meson measurement. It computes the zero-momentum
local correlators

```math
C_{PP}(t)=\sum_{\mathbf{x}}\langle P(\mathbf{x},t)P(0)\rangle,
\qquad
C_{A_4P}(t)=\sum_{\mathbf{x}}\langle A_4(\mathbf{x},t)P(0)\rangle,
```

and returns the bare PCAC mass in lattice units,

```math
am_{\mathrm{PCAC}}(t)=
\frac{\widetilde\partial_4 C_{A_4P}(t)
      +c_A\,\partial_4^*\partial_4 C_{PP}(t)}{2C_{PP}(t)}.
```

Here $\widetilde\partial$ is the nearest-neighbour symmetric derivative. The
default `improvement_coefficient=0` is the unimproved definition; setting it
to a known axial-current coefficient supplies the optional `c_A` term. The
returned quantity is not renormalized: multiply an ensemble result by the
appropriate `Z_A/Z_P` where required. Contact points and points without a
PCAC plateau should be excluded during the ensemble analysis.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

U = gauge_configuration(
    (4, 4, 4, 4);               # small runnable example
    colors=3,
    start=:cold,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)

measurement = PCACMassMeasurement(
    U;
    kappa=0.141139,              # `κ=0.141139` is equivalent
    correlation_axis=4,
    source_position=(1, 1, 1, 1),
    improvement_coefficient=0.0,
    BoundaryCondition=[1, 1, 1, -1],
)
result = get_value(measure(measurement, U))

cpp = result[:PP]
cap = result[:AP]
am_pcac = result[:mass]
```

Use `kappa2`/`κ2` in the direct constructor for a non-degenerate second
valence quark. As with `t0`, production analysis should first collect the raw
correlators over an ensemble, then form a correlated ratio and fit a plateau;
an average of per-configuration ratios is generally not the preferred final
estimator.

The dictionary interface uses `hop` for the Wilson hopping parameter:

```julia
measurement = prepare_measurement(
    U,
    Dict(
        "methodname" => "PCAC_mass",
        "fermiontype" => "Wilson",
        "hop" => 0.141139,
        "correlation_axis" => 4,
        "source_position" => [1, 1, 1, 1],
        "improvement_coefficient" => 0.0,
        "printvalues" => false,
    ),
)
```

The definition follows the axial Ward-identity construction of M. Lüscher
et al., *Nucl. Phys. B* **491**, 323--343 (1997),
[DOI: 10.1016/S0550-3213(97)00080-1](https://doi.org/10.1016/S0550-3213(97)00080-1),
[arXiv:hep-lat/9609035](https://arxiv.org/abs/hep-lat/9609035).
The symmetric-derivative formula and its `c_A`-improved form are given
explicitly in J. Rolf and S. Sint (ALPHA Collaboration), *JHEP* **12** (2002)
007, Eqs. (3.24) and (3.37),
[arXiv:hep-ph/0209255](https://arxiv.org/abs/hep-ph/0209255).

For independent numerical validation, QCDMeasurements' `PP` and `A_4P`
correlators and final `c_A=0` ratio were compared with Grid's `Gamma5` and
`GammaTGamma5` channels. Both cold and reconstructible non-cold fields agree;
the raw values, normalization, and covered execution paths are documented in
[Independent validation](#independent-validation).

## Gradient-flow scales (`t0` and `w0`)

The standard scale `t0` is an ensemble quantity defined by
`t0^2 * <E(t0)> = c`, usually with `c = 0.3`.  QCDMeasurements therefore uses
two stages: `measure` maps each input configuration to a `GradientFlowHistory`,
and `estimate_flow_scales` averages the histories before finding the crossing.
The input gauge configuration is not modified.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

# Replace these small hot-start examples with an ensemble of stored
# configurations having the same lattice geometry.
configurations = [
    gauge_configuration(
        (4, 4, 4, 4);
        colors=3,
        start=:hot,
        seed=seed,
        process_grid=(1, 1, 1, 1),
        verbose=0,
    ) for seed in (1234, 5678)
]

flow_measurement = GradientFlowScaleMeasurement(
    first(configurations);
    flow_step_size=0.01,
    number_of_flow_steps=200, # choose a range that contains the c=0.3 crossing
    flow_measure_every=2,
    energy_methods=["clover"],
    measure_topological_charge=false,
)

histories = GradientFlowHistory[]
for U in configurations
    push!(histories, get_value(measure(flow_measurement, U)))
end

scales = estimate_flow_scales(histories; c=0.3)
t0_over_a2 = scales.t0["clover"]
t0_error = scales.t0_error["clover"]
w0_over_a = scales.w0["clover"]
w0_error = scales.w0_error["clover"]
```

Do not determine a crossing separately on every configuration and average the
result: that is not the standard definition because the ensemble average
`<E(t)>` must be formed first.  The implementation follows
[Lüscher, JHEP 08 (2010) 071](https://arxiv.org/abs/1006.4518) and
[Borsanyi et al., JHEP 09 (2012) 010](https://arxiv.org/abs/1203.4469).

## Improved topological-charge definitions

With `TC_methods=["clover"]`, the option
`improved_topological_charge_definition` has exactly three choices:

- `"alexandrou"` (default) computes the directly improved density
  `q_imp = (5/3)q_clover - (1/12)q_rectangle`.  Its output key is
  `"clover improved"`.
- `"bilson_thompson"` first constructs the improved field strength
  `F_imp = (5/3)F_clover - (1/3)F_rectangle` and then computes `q[F_imp]`.
  Its output key is `"clover improved bilson-thompson"`.
- `"both"` computes both results.

The following is a complete example.  Replace the hot-start field by a field
loaded from a configuration file for a production measurement.

```julia
import JACC
using Gaugefields
using QCDMeasurements

JACC.@init_backend

U = gauge_configuration(
    (4, 4, 4, 4);
    colors=3,
    start=:hot,
    seed=1234,
    process_grid=(1, 1, 1, 1),
    verbose=0,
)

# Default: the directly improved Alexandrou definition.
m_default = TopologicalChargeMeasurement(
    U;
    TC_methods=["clover"],
)
default_result = get_value(measure(m_default, U))
q_alexandrou = default_result["clover improved"]

# Alternative: improve F_{mu nu} first, then construct q[F].
m_field_strength = TopologicalChargeMeasurement(
    U;
    TC_methods=["clover"],
    improved_topological_charge_definition="bilson_thompson",
)
field_strength_result = get_value(measure(m_field_strength, U))
q_bilson_thompson =
    field_strength_result["clover improved bilson-thompson"]

# Calculate both definitions in one measurement.
m_both = TopologicalChargeMeasurement(
    U;
    TC_methods=["clover"],
    improved_topological_charge_definition="both",
)
both_results = get_value(measure(m_both, U))
q_alexandrou_both = both_results["clover improved"]
q_bilson_thompson_both =
    both_results["clover improved bilson-thompson"]
```

The same selection is available through the dictionary API:

```julia
measurement = prepare_measurement(
    U,
    Dict(
        "methodname" => "Topological_charge",
        "kinds_of_topological_charge" => ["clover"],
        "improved_topological_charge_definition" => "both",
        "printvalues" => false,
    ),
)
values = get_value(measure(measurement, U))
```

It is also accepted by `GradientFlowScaleMeasurement` and the
topological-charge-density correlation measurement.  For continuum physics,
measure these observables after a specified smoothing or gradient-flow time
and report both the flow prescription and the chosen charge definition.

### References for the topological-charge definitions

The default `"alexandrou"` implementation uses the directly improved density
of Eqs. (15)--(17) in C. Alexandrou, A. Athenodorou, and K. Jansen,
"Topological charge using cooling and the gradient flow," *Phys. Rev. D* **92**,
125014 (2015),
[DOI: 10.1103/PhysRevD.92.125014](https://doi.org/10.1103/PhysRevD.92.125014),
[arXiv:1509.04259](https://arxiv.org/abs/1509.04259).

The `"bilson_thompson"` option follows the construction principle of first
improving the lattice field-strength tensor and then forming the topological
charge, from S. O. Bilson-Thompson, D. B. Leinweber, and A. G. Williams,
"Highly-improved lattice field-strength tensor," *Annals of Physics* **304**,
1--21 (2003),
[DOI: 10.1016/S0003-4916(03)00009-5](https://doi.org/10.1016/S0003-4916(03)00009-5),
[arXiv:hep-lat/0203008](https://arxiv.org/abs/hep-lat/0203008).
The present implementation specializes that construction to the tree-level
`1x1 + 1x2` improved tensor.  This two-loop tensor and its use in the
topological charge are also described by P. T. Jahn, G. D. Moore, and
D. Robaina, "Estimating $\chi_{\mathrm{top}}$ lattice artifacts from flowed SU(2)
calorons," *Eur. Phys. J. C* **79**, 510 (2019),
[DOI: 10.1140/epjc/s10052-019-7008-9](https://doi.org/10.1140/epjc/s10052-019-7008-9),
[arXiv:1805.11511](https://arxiv.org/abs/1805.11511).

## Citation

If you use this package in a paper, please cite:

```
@article{Nagai:2024yaf,
    author = "Nagai, Yuki and Tomiya, Akio",
    title = "{JuliaQCD: Portable lattice QCD package in Julia language}",
    eprint = "2409.03030",
    archivePrefix = "arXiv",
    primaryClass = "hep-lat",
    month = "9",
    year = "2024"
}
```
and the paper is [arXiv:2409.03030](https://arxiv.org/abs/2409.03030).
