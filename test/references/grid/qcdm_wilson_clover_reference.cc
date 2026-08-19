// Independent Wilson--clover PP/A4P reference for QCDMeasurements.jl.
// Build this against Grid commit 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1.

#include <Grid/Grid.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace Grid;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);
  const Coordinate lattice({4, 4, 4, 4});
  const Coordinate mpi({1, 1, 1, 1});
  const Coordinate simd = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto *grid = SpaceTimeGrid::makeFourDimGrid(lattice, simd, mpi);
  auto *rb_grid = SpaceTimeGrid::makeFourDimRedBlackGrid(grid);

  LatticeGaugeField gauge(grid);
  SU<Nc>::ColdConfiguration(gauge);
  LorentzColourMatrixD origin_links;
  const Coordinate origin({0, 0, 0, 0});
  peekSite(origin_links, gauge, origin);
  constexpr RealD theta = 1.1;
  origin_links(0)()(0, 0) = ComplexD(std::cos(theta), std::sin(theta));
  origin_links(0)()(1, 1) = ComplexD(std::cos(theta), -std::sin(theta));
  origin_links(0)()(2, 2) = ComplexD(1.0, 0.0);
  pokeSite(origin_links, gauge, origin);

  constexpr RealD kappa = 0.1;
  constexpr RealD mass = 1.0 / (2.0 * kappa) - 4.0;
  constexpr RealD csw = 1.2;
  WilsonCloverFermionD::ImplParams params;
  params.boundary_phases[Nd - 1] = -1.0;
  WilsonAnisotropyCoefficients anisotropy;
  WilsonCloverFermionD dirac(
      gauge, *grid, *rb_grid, mass, csw, csw, anisotropy, params);

  LatticePropagator source(grid);
  LatticePropagator propagator(grid);
  source = Zero();
  propagator = Zero();
  SpinColourMatrix unit;
  unit = 1.0;
  pokeSite(unit, source, origin);

  ConjugateGradient<LatticeFermion> cg(1.0e-14, 100000);
  SchurRedBlackDiagTwoSolve<LatticeFermion> solve(cg);
  ZeroGuesser<LatticeFermion> guess;
  for (int spin = 0; spin < Ns; ++spin) {
    for (int colour = 0; colour < Nc; ++colour) {
      LatticeFermion rhs(grid);
      LatticeFermion solution(grid);
      PropToFerm<WilsonCloverFermionD>(rhs, source, spin, colour);
      solution = Zero();
      solve(dirac, rhs, solution, guess);
      FermToProp<WilsonCloverFermionD>(propagator, solution, spin, colour);
    }
  }

  Gamma gamma5(Gamma::Algebra::Gamma5);
  Gamma gamma_t_gamma5(Gamma::Algebra::GammaTGamma5);
  LatticeComplex pp =
      trace(gamma5 * adj(propagator) * gamma5 * gamma5 * propagator *
            adj(gamma5));
  LatticeComplex ap =
      trace(gamma5 * adj(propagator) * gamma5 * gamma_t_gamma5 * propagator *
            adj(gamma5));
  std::vector<TComplex> pp_timeslices;
  std::vector<TComplex> ap_timeslices;
  sliceSum(pp, pp_timeslices, Nd - 1);
  sliceSum(ap, ap_timeslices, Nd - 1);

  std::cout << std::setprecision(17);
  std::cout << "GRID_COMMIT 0ac72cb6a30ccdc41d664e7e0759f0c8833078f1\n";
  std::cout << "GAUGE localized theta=" << theta << "\n";
  std::cout << "PLAQUETTE "
            << WilsonLoops<PeriodicGimplR>::avgPlaquette(gauge) << "\n";
  std::cout << "KAPPA " << kappa << " MASS " << mass
            << " CSW " << csw << "\n";
  std::cout << "BOUNDARY 1 1 1 -1\n";
  for (int t = 0; t < static_cast<int>(pp_timeslices.size()); ++t) {
    const Complex pp_raw = TensorRemove(pp_timeslices[t]);
    const Complex ap_raw = TensorRemove(ap_timeslices[t]);
    std::cout << "MESON t=" << t
              << " pp_grid_raw=" << pp_raw
              << " pp_qcdm_normalized=" << pp_raw / (4.0 * kappa * kappa)
              << " ap_grid_raw=" << ap_raw
              << " ap_qcdm_normalized=" << ap_raw / (4.0 * kappa * kappa)
              << "\n";
  }
  for (int t = 0; t < static_cast<int>(pp_timeslices.size()); ++t) {
    const int nt = static_cast<int>(pp_timeslices.size());
    const Complex pp_raw = TensorRemove(pp_timeslices[t]);
    const Complex ap_forward = TensorRemove(ap_timeslices[(t + 1) % nt]);
    const Complex ap_backward = TensorRemove(ap_timeslices[(t + nt - 1) % nt]);
    const Complex derivative = (ap_forward - ap_backward) / 2.0;
    std::cout << "PCAC t=" << t
              << " c_A=0 am_pcac=" << derivative / (2.0 * pp_raw)
              << "\n";
  }

  delete rb_grid;
  delete grid;
  Grid_finalize();
  return 0;
}
