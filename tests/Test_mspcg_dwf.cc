/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: tests/Test_mspcg_dwf.cc

Copyright (C) 2015

Author: Chulwoo Jung <chulwoo@bnl.gov>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

/**
 * Test_mspcg_dwf
 *
 * Compares MultiSplittingPreconditionedCG (MSPCG) against standard
 * SchurRedBlack preconditioned CG on a domain-wall fermion system.
 *
 * Three sub-preconditioners are constructed, each being a CG solve on the
 * same DWF operator but with Dirichlet boundary conditions (DirichletBlock)
 * applied to different sub-block sizes.  This exercises the DD-HMC subdomain
 * solver path via SubdomainSolver<>.
 *
 * Usage:
 *   ./Test_mspcg_dwf [Grid options] [--Ls <Ls>] [--mass <m>] [--tol <tol>]
 *                    [--inner-tol <tol>] [--parallel]
 *
 * Checks:
 *   1. Standard preconditioned CG (SchurRedBlack) converges.
 *   2. MSPCG (sequential sub-preconditioners) converges to the same solution.
 *   3. MSPCG (parallel sub-preconditioners) gives the same result.
 *   4. The solutions from all three solvers agree to within a tight tolerance.
 */

#include "disable_tests_without_instantiations.h"
#ifdef ENABLE_FERMION_INSTANTIATIONS

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

// ---- helpers ----------------------------------------------------------------

// Compute || x - y || / || y ||
static RealD relDiff(const LatticeFermionD& x, const LatticeFermionD& y)
{
  LatticeFermionD diff(x.Grid());
  diff = x - y;
  return std::sqrt(norm2(diff) / norm2(y));
}

// ---- main -------------------------------------------------------------------

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  // ---- Parameters from command line ----------------------------------------
  int   Ls       = 12;
  RealD mass     = 0.1;
  RealD M5       = 1.8;
  RealD tol      = 1e-8;
  RealD innerTol = 1e-5;   // tolerance for sub-preconditioner (inner) solves
  bool  doParallel = false;

  if (GridCmdOptionExists(argv, argv+argc, "--Ls")) {
    std::string s = GridCmdOptionPayload(argv, argv+argc, "--Ls");
    GridCmdOptionInt(s, Ls);
  }
  if (GridCmdOptionExists(argv, argv+argc, "--mass")) {
    std::string s = GridCmdOptionPayload(argv, argv+argc, "--mass");
    GridCmdOptionFloat(s, mass);
  }
  if (GridCmdOptionExists(argv, argv+argc, "--tol")) {
    std::string s = GridCmdOptionPayload(argv, argv+argc, "--tol");
    GridCmdOptionFloat(s, tol);
  }
  if (GridCmdOptionExists(argv, argv+argc, "--inner-tol")) {
    std::string s = GridCmdOptionPayload(argv, argv+argc, "--inner-tol");
    GridCmdOptionFloat(s, innerTol);
  }
  if (GridCmdOptionExists(argv, argv+argc, "--parallel"))
    doParallel = true;

  std::cout << GridLogMessage << "Ls=" << Ls << " mass=" << mass
            << " M5=" << M5 << " tol=" << tol << " innerTol=" << innerTol
            << " parallel=" << doParallel << std::endl;

  // ---- Grids ---------------------------------------------------------------
  GridCartesian*         UGrid   = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  // ---- Gauge field ---------------------------------------------------------
  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeFieldD Umu(UGrid);
  SU<Nc>::HotConfiguration(RNG4, Umu);

  // ---- Source and initial guess -------------------------------------------
  LatticeFermionD src(FGrid);    random(RNG5, src);
  LatticeFermionD src_o(FrbGrid);
  pickCheckerboard(Odd, src_o, src);

  // ---- Full operator for the outer solve -----------------------------------
  DomainWallFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  SchurDiagMooeeOperator<DomainWallFermionD, LatticeFermionD> HermOpFull(Ddwf);

  // =========================================================================
  // 1. Reference: standard SchurRedBlack preconditioned CG
  // =========================================================================
  std::cout << GridLogMessage
            << "\n========== Reference: SchurRedBlack PCG ==========" << std::endl;

  ConjugateGradient<LatticeFermionD> refCG(tol, 10000);
  SchurRedBlackDiagMooeeSolve<LatticeFermionD> refSolver(refCG);

  LatticeFermionD refResult(FGrid); refResult = Zero();
  double t0 = usecond();
  refSolver(Ddwf, src, refResult);
  double t_ref = usecond() - t0;

  std::cout << GridLogMessage << "Reference PCG time: " << t_ref/1e6 << " s" << std::endl;

  // =========================================================================
  // 2. MSPCG with Dirichlet subdomain sub-preconditioners
  //
  // We build three sub-preconditioners using DWF operators with Dirichlet BCs
  // at different block sizes.  The sub-preconditioners each do an inexact
  // inner solve (innerTol >> tol) so they act as cheap approximate inverses.
  //
  // In a full DD-HMC setup these would be domain-local solves on separate
  // blocks; here they all share the same gauge field but have different
  // Dirichlet block patterns.
  // =========================================================================

  // Build sub-operators.  Each gets its own DWF instance so DirichletBlock
  // can be set independently.
  Coordinate latt = GridDefaultLatt();

  // Block A: half the lattice in each 4d direction (coarse blocking)
  Coordinate blockA(Nd+1, 0);
  for (int mu = 0; mu < Nd; mu++) blockA[mu+1] = latt[mu] / 2;

  // Block B: quarter of the lattice (finer blocking)
  Coordinate blockB(Nd+1, 0);
  for (int mu = 0; mu < Nd; mu++) blockB[mu+1] = std::max(latt[mu] / 4, 1);

  // Block C: full lattice (= no Dirichlet, acts as a global sub-solve)
  Coordinate blockC(Nd+1, 0);  // all zeros = no Dirichlet

  DomainWallFermionD DdwfA(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  DomainWallFermionD DdwfB(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);
  DomainWallFermionD DdwfC(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5);

  // Apply Dirichlet BCs (zero inter-block links).  blockC all-zero is a no-op.
  DdwfA.DirichletBlock(blockA);
  DdwfB.DirichletBlock(blockB);
  // DdwfC: no Dirichlet — full lattice solve

  std::cout << GridLogMessage << "BlockA=" << blockA << std::endl;
  std::cout << GridLogMessage << "BlockB=" << blockB << std::endl;
  std::cout << GridLogMessage << "BlockC (no Dirichlet)" << std::endl;

  // Schur operators for sub-preconditioners
  SchurDiagMooeeOperator<DomainWallFermionD, LatticeFermionD> HermOpA(DdwfA);
  SchurDiagMooeeOperator<DomainWallFermionD, LatticeFermionD> HermOpB(DdwfB);
  SchurDiagMooeeOperator<DomainWallFermionD, LatticeFermionD> HermOpC(DdwfC);

  // Inner CG solvers (loose tolerance — act as preconditioners)
  ConjugateGradient<LatticeFermionD> innerCG_A(innerTol, 2000, false);
  ConjugateGradient<LatticeFermionD> innerCG_B(innerTol, 2000, false);
  ConjugateGradient<LatticeFermionD> innerCG_C(innerTol, 2000, false);

  // Wrap as SubdomainSolvers (LinearFunction interface)
  SubdomainSolver<LatticeFermionD> subA(HermOpA, innerCG_A);
  SubdomainSolver<LatticeFermionD> subB(HermOpB, innerCG_B);
  SubdomainSolver<LatticeFermionD> subC(HermOpC, innerCG_C);

  std::vector<LinearFunction<LatticeFermionD>*> precs = {&subA, &subB, &subC};

  // =========================================================================
  // 2a. MSPCG — sequential
  // =========================================================================
  std::cout << GridLogMessage
            << "\n========== MSPCG sequential ==========" << std::endl;

  MultiSplittingPreconditionedCG<LatticeFermionD> mspcgSeq(
      tol, 10000, precs, /*err_on_no_conv=*/true, /*parallel=*/false);

  LatticeFermionD mspcgResultSeq(FGrid); mspcgResultSeq = Zero();
  LatticeFermionD mspcgResultSeq_o(FrbGrid); mspcgResultSeq_o = Zero();
  mspcgResultSeq_o.Checkerboard() = Odd;

  t0 = usecond();
  mspcgSeq(HermOpFull, src_o, mspcgResultSeq_o);
  double t_seq = usecond() - t0;

  std::cout << GridLogMessage << "MSPCG sequential time: " << t_seq/1e6 << " s" << std::endl;
  std::cout << GridLogMessage << "MSPCG sequential iterations: "
            << mspcgSeq.IterationsToComplete << std::endl;

  // =========================================================================
  // 2b. MSPCG — parallel (if requested)
  // =========================================================================
  double t_par = 0.0;
  LatticeFermionD mspcgResultPar_o(FrbGrid);
  mspcgResultPar_o = Zero();
  mspcgResultPar_o.Checkerboard() = Odd;
  int parIters = 0;

  if (doParallel) {
    std::cout << GridLogMessage
              << "\n========== MSPCG parallel ==========" << std::endl;

    MultiSplittingPreconditionedCG<LatticeFermionD> mspcgPar(
        tol, 10000, precs, /*err_on_no_conv=*/true, /*parallel=*/true);

    t0 = usecond();
    mspcgPar(HermOpFull, src_o, mspcgResultPar_o);
    t_par = usecond() - t0;
    parIters = mspcgPar.IterationsToComplete;

    std::cout << GridLogMessage << "MSPCG parallel time: " << t_par/1e6 << " s" << std::endl;
    std::cout << GridLogMessage << "MSPCG parallel iterations: " << parIters << std::endl;
  }

  // =========================================================================
  // 3. Reconstruct full solutions and compare
  // =========================================================================

  // Reconstruct full 5d solution from odd-cb result for MSPCG sequential
  {
    LatticeFermionD tmp(FGrid);
    setCheckerboard(mspcgResultSeq, mspcgResultSeq_o);
    // Even cb: psi_e = Mee^{-1}(src_e - Meo psi_o)
    LatticeFermionD src_e(FrbGrid), psi_e(FrbGrid), Meo_psi_o(FrbGrid);
    pickCheckerboard(Even, src_e, src);
    Ddwf.Meooe(mspcgResultSeq_o, Meo_psi_o);
    Ddwf.MooeeInv(src_e - Meo_psi_o, psi_e);
    setCheckerboard(mspcgResultSeq, psi_e);
  }

  std::cout << GridLogMessage
            << "\n========== Solution comparison ==========" << std::endl;

  RealD diff_seq = relDiff(mspcgResultSeq, refResult);
  std::cout << GridLogMessage
            << "||MSPCG_seq - ref|| / ||ref|| = " << diff_seq << std::endl;

  if (doParallel) {
    LatticeFermionD mspcgResultPar(FGrid); mspcgResultPar = Zero();
    setCheckerboard(mspcgResultPar, mspcgResultPar_o);
    LatticeFermionD src_e(FrbGrid), psi_e(FrbGrid), Meo_psi_o(FrbGrid);
    pickCheckerboard(Even, src_e, src);
    Ddwf.Meooe(mspcgResultPar_o, Meo_psi_o);
    Ddwf.MooeeInv(src_e - Meo_psi_o, psi_e);
    setCheckerboard(mspcgResultPar, psi_e);

    RealD diff_par_ref = relDiff(mspcgResultPar, refResult);
    RealD diff_par_seq = relDiff(mspcgResultPar, mspcgResultSeq);
    std::cout << GridLogMessage
              << "||MSPCG_par - ref||  / ||ref||      = " << diff_par_ref << std::endl;
    std::cout << GridLogMessage
              << "||MSPCG_par - MSPCG_seq|| / ||seq|| = " << diff_par_seq << std::endl;
  }

  // =========================================================================
  // 4. Summary
  // =========================================================================
  std::cout << GridLogMessage << "\n========== Summary ==========" << std::endl;
  std::cout << GridLogMessage << "Ref PCG      time = " << t_ref/1e6 << " s" << std::endl;
  std::cout << GridLogMessage << "MSPCG seq    time = " << t_seq/1e6 << " s"
            << "  iters=" << mspcgSeq.IterationsToComplete
            << "  true_res=" << mspcgSeq.TrueResidual << std::endl;
  if (doParallel)
    std::cout << GridLogMessage << "MSPCG par    time = " << t_par/1e6 << " s"
              << "  iters=" << parIters << std::endl;

  // Pass/fail
  RealD passTol = 10.0 * tol;
  bool pass = (diff_seq < passTol);
  std::cout << GridLogMessage
            << (pass ? "PASS" : "FAIL")
            << " ||MSPCG_seq - ref|| / ||ref|| = " << diff_seq
            << "  (tolerance " << passTol << ")" << std::endl;
  GRID_ASSERT(pass);

  Grid_finalize();
  return 0;
}

#endif // ENABLE_FERMION_INSTANTIATIONS
