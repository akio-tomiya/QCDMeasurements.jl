# Changes

## 1.0.1 — 2026-10-21

### Density APIs

- Add `QCDMeasurements.topological_charge_density` and
  `topological_charge_density!` for plaquette, clover, Alexandrou, and
  Bilson–Thompson definitions, including the cross terms in the improved
  field-strength definition.
- Add `action_density` and `action_density!` for plaquette and clover
  definitions. Their volume averages match the gradient-flow energy
  observables; unit links give zero. The historical `EnergyDensityMeasurement`
  convention remains unchanged.
- Return `Array{Float64,4}` in `(X,Y,Z,T)` order. Topological-density sums
  match the existing scalar charges up to floating-point reduction roundoff.
- Add `GaugeDensityWorkspace` to reuse loop paths, temporary fields, and
  scalar storage without retaining or modifying the input configuration.

### CPU, GPU, and MPI execution

- Use JACC for local contractions on CPU and GPU, and Gaugefields and
  LatticeMatrices for loop construction, halo exchange, and MPI collection.
  Distributed GPU fields are supported. MPI calls are collective and return
  the complete global host array on every rank.
- Support legacy CPU and MPI fields, including conversion between legacy
  rank ordering and the Cartesian ordering used for scalar-field collection.
- Preserve the selected GPU when allocating density result storage. Validate
  workspace field type, dimensions, halo width, communicator, and process
  decomposition.
- Make JACC a direct dependency and add an optional MPI extension for legacy
  distributed fields. Serial use does not require loading MPI.
- Skip zero Levi-Civita terms before global traces in scalar topological
  measurements, reducing the ordered plane-pair contractions from 144 to 24.

### Numerical regression and validation

- Freeze scalar observables and all 360 sites of each density from unmodified
  commit `5576729f81c07955df01bee1f13d1ae12017f2e1` on deterministic noncubic
  SU(2)/SU(3) configurations. Include input SHA-256 digests and a generator
  that checks the baseline source and thread count.
- Verify both allocating and reusable APIs, input immutability, workspace
  reuse on unit links, multiple MPI decompositions, and legacy rank ordering.
- The independent before/after comparison passed 48/48 checks. Existing
  scalar observables were bitwise unchanged under the same four-thread CPU
  setup. Maximum pointwise absolute differences from the old implementation
  were `2.78e-17` for topological density and `8.89e-16` for action density.
- Full suites passed: CPU 575/575, GPU 136/136, and two-rank MPI 641/641 per
  rank. Density regressions passed on four CPU ranks (221/221 per rank),
  two GPU/MPI ranks (441/441 per rank), and four GPU/MPI ranks (221/221 per
  rank). Two ranks on separate H100 and Blackwell GPUs also passed 441/441
  per rank using host-staged MPI communication.
- Record dependency versions and reproduction instructions in
  [the validation record](test/README.md#density-api-validation-2026-10-21).
