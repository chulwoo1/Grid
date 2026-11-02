/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./Grid/algorithms/iterative/MultiSplittingPreconditionedCG.h

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
#ifndef GRID_MULTI_SPLITTING_PRECONDITIONED_CG_H
#define GRID_MULTI_SPLITTING_PRECONDITIONED_CG_H

#include <thread>

NAMESPACE_BEGIN(Grid);

//--------------------------------------------------------------------
// SubdomainSolver
//
// Wraps a (LinearOperatorBase, OperatorFunction) pair as a
// LinearFunction so it can be used as a sub-preconditioner in
// MultiSplittingPreconditionedCG.
//
// For DD-HMC usage: configure the fermion operator with
//   fermOp.DirichletBlock(block)
// before constructing the SubdomainSolver, or call
//   solver.setBlock(block)
// to (re-)apply Dirichlet BCs on demand.
//
// Example:
//   MobiusFermionD subOp(Umu, *FGrid, *FrbGrid, mass, M5, b, c);
//   subOp.DirichletBlock(block);
//   SchurRedBlackDiagMooeeSolve<FermionField> schur(cg);
//   SubdomainSolver<FermionField> prec(subOp, schur);
//--------------------------------------------------------------------
template<class Field>
class SubdomainSolver : public LinearFunction<Field> {
public:
  LinearOperatorBase<Field>& Op;
  OperatorFunction<Field>&   Solver;

  SubdomainSolver(LinearOperatorBase<Field>& op, OperatorFunction<Field>& solver)
    : Op(op), Solver(solver) {}

  void operator()(const Field& src, Field& psi) override {
    Solver(Op, src, psi);
  }
};

/**
 * Multi-splitting Preconditioned Conjugate Gradient
 *
 * Solves  A x = b  where A is Hermitian positive definite.
 *
 * Preconditioner:
 *   Given k sub-preconditioners M_1^{-1}, ..., M_k^{-1} corresponding to
 *   splittings A = M_i - N_i, the composite preconditioner is
 *
 *       M^{-1} v = sum_{i=1}^{k} w_i * M_i^{-1} v
 *
 *   where w_i >= 0 and sum w_i = 1 (default: uniform w_i = 1/k).
 *   If every M_i^{-1} is SPD the composite is SPD, so standard PCG applies.
 *
 * The sub-preconditioners are provided as LinearFunction<Field> objects.
 * Typical choices:
 *   - Domain-restricted solves (additive Schwarz / DD-HMC subdomain CG)
 *   - Block-Jacobi applied to disjoint subsets of unknowns
 *   - Red-black / multi-colour sub-solves
 *   - Inexact sub-solves (e.g., a few steps of CG on a sub-domain)
 *
 * Parallel execution (parallelApply = true):
 *   All sub-preconditioner applies are launched as std::threads and run
 *   concurrently.  Each thread operates on its own pre-allocated scratch
 *   field so there is no shared mutable state.  This is most beneficial
 *   when each sub-preconditioner is a Dirichlet (DD-HMC) subdomain solve
 *   whose inter-block coupling is zeroed — the solves are then truly
 *   independent and can overlap in time.
 *
 *   CAUTION: Grid's internal OpenMP parallelism is active inside each
 *   sub-solve, so enabling parallelApply launches k*nThreads OS threads
 *   simultaneously.  On GPU builds this is less relevant; on CPU builds
 *   consider setting OMP_NUM_THREADS accordingly.
 *
 * PCG iteration:
 *   r    = b - A x
 *   z    = M^{-1} r          (multi-splitting apply)
 *   p    = z
 *   rz   = (r, z)
 *   loop:
 *     q     = A p
 *     alpha = rz / (p, q)
 *     x    += alpha * p
 *     r    -= alpha * q
 *     z     = M^{-1} r
 *     rz'   = (r, z)
 *     beta  = rz' / rz
 *     p     = z + beta * p
 *     rz    = rz'
 *
 * References:
 *   O'Leary & White, SIAM J. Algebraic Discrete Methods 6(4):630-640, 1985.
 *   White, Linear Algebra Appl. 131:55-65, 1990.
 */
template <class Field>
class MultiSplittingPreconditionedCG : public OperatorFunction<Field> {
public:
  using OperatorFunction<Field>::operator();

  RealD   Tolerance;
  int     MaxIterations;
  bool    ErrorOnNoConverge;
  bool    parallelApply;      // run sub-preconditioners concurrently via std::thread
  int     IterationsToComplete;
  RealD   TrueResidual;

  // Sub-preconditioners M_i^{-1} and their weights w_i.
  // weights are normalised to sum to 1 on construction.
  std::vector<LinearFunction<Field>*> Subpreconditioners;
  std::vector<RealD>                  Weights;

  GridStopWatch MatrixTimer;
  GridStopWatch PrecTimer;
  GridStopWatch LinalgTimer;
  GridStopWatch SolverTimer;

  // Uniform weights
  MultiSplittingPreconditionedCG(RealD tol, int maxit,
                                 std::vector<LinearFunction<Field>*> precs,
                                 bool err_on_no_conv = true,
                                 bool parallel = false)
    : Tolerance(tol)
    , MaxIterations(maxit)
    , ErrorOnNoConverge(err_on_no_conv)
    , parallelApply(parallel)
    , IterationsToComplete(0)
    , TrueResidual(0.0)
    , Subpreconditioners(precs)
  {
    int k = (int)precs.size();
    GRID_ASSERT(k > 0);
    Weights.assign(k, 1.0 / k);
  }

  // User-supplied weights (normalised internally)
  MultiSplittingPreconditionedCG(RealD tol, int maxit,
                                 std::vector<LinearFunction<Field>*> precs,
                                 std::vector<RealD> weights,
                                 bool err_on_no_conv = true,
                                 bool parallel = false)
    : Tolerance(tol)
    , MaxIterations(maxit)
    , ErrorOnNoConverge(err_on_no_conv)
    , parallelApply(parallel)
    , IterationsToComplete(0)
    , TrueResidual(0.0)
    , Subpreconditioners(precs)
    , Weights(weights)
  {
    GRID_ASSERT(precs.size() == weights.size() && precs.size() > 0);
    RealD wsum = 0.0;
    for (auto w : Weights) { GRID_ASSERT(w >= 0.0); wsum += w; }
    GRID_ASSERT(wsum > 0.0);
    for (auto& w : Weights) w /= wsum;
  }

  void operator()(LinearOperatorBase<Field>& LinOp, const Field& src, Field& psi)
  {
    psi.Checkerboard() = src.Checkerboard();
    conformable(psi, src);

    RealD guess = norm2(psi);
    GRID_ASSERT(std::isnan(guess) == 0);

    RealD ssq = norm2(src);
    if (ssq == 0.0) { psi = Zero(); IterationsToComplete = 0; TrueResidual = 0.0; return; }
    RealD rsq = Tolerance * Tolerance * ssq;

    std::cout << GridLogIterative << std::setprecision(8)
              << "MultiSplittingPCG: guess  " << guess << std::endl;
    std::cout << GridLogIterative << std::setprecision(8)
              << "MultiSplittingPCG: src    " << ssq   << std::endl;
    std::cout << GridLogIterative << std::setprecision(8)
              << "MultiSplittingPCG: target " << rsq   << std::endl;
    std::cout << GridLogIterative
              << "MultiSplittingPCG: " << Subpreconditioners.size()
              << " sub-preconditioners"
              << (parallelApply ? " [parallel]" : " [sequential]") << std::endl;

    MatrixTimer.Reset(); PrecTimer.Reset();
    LinalgTimer.Reset(); SolverTimer.Reset();
    SolverTimer.Start();

    // Pre-allocate one scratch field per sub-preconditioner so parallel
    // threads have no shared mutable state.
    int k = (int)Subpreconditioners.size();
    std::vector<Field> scratch(k, Field(src.Grid()));

    Field r(src.Grid());
    Field z(src.Grid());
    Field p(src.Grid());
    Field q(src.Grid());

    // ---- Initial residual ----
    MatrixTimer.Start();
    if (guess == 0.0) {
      r = src;
    } else {
      LinOp.HermOp(psi, q);
      r = src - q;
    }
    MatrixTimer.Stop();

    // ---- Initial preconditioner apply: z = M^{-1} r ----
    applyPreconditioner(r, z, scratch);

    LinalgTimer.Start();
    p = z;
    RealD rz = real(innerProduct(r, z));
    RealD cp = norm2(r);
    LinalgTimer.Stop();

    for (int iter = 1; iter <= MaxIterations; ++iter) {

      // ---- q = A p ----
      MatrixTimer.Start();
      LinOp.HermOp(p, q);
      MatrixTimer.Stop();

      // ---- alpha = (r,z) / (p, Ap) ----
      LinalgTimer.Start();
      RealD pq    = real(innerProduct(p, q));
      RealD alpha = rz / pq;

      psi = psi + alpha * p;
      r   = r   - alpha * q;
      cp  = norm2(r);
      LinalgTimer.Stop();

      std::cout << GridLogIterative
                << "MultiSplittingPCG: iter=" << iter
                << "  res^2=" << cp
                << "  target=" << rsq << std::endl;

      if (cp <= rsq) {
        SolverTimer.Stop();

        MatrixTimer.Start();
        LinOp.HermOp(psi, q);
        MatrixTimer.Stop();
        r = src - q;
        TrueResidual      = std::sqrt(norm2(r) / ssq);
        IterationsToComplete = iter;

        std::cout << GridLogMessage
                  << "MultiSplittingPCG converged in " << iter << " iterations"
                  << "  computed res " << std::sqrt(cp / ssq)
                  << "  true res "     << TrueResidual
                  << "  target "       << Tolerance << std::endl;
        printTimers();
        return;
      }

      // ---- z = M^{-1} r ----
      applyPreconditioner(r, z, scratch);

      // ---- beta = (r',z') / (r,z),  p = z + beta p ----
      LinalgTimer.Start();
      RealD rz_new = real(innerProduct(r, z));
      RealD beta   = rz_new / rz;
      p  = z + beta * p;
      rz = rz_new;
      LinalgTimer.Stop();
    }

    SolverTimer.Stop();

    MatrixTimer.Start();
    LinOp.HermOp(psi, q);
    MatrixTimer.Stop();
    r = src - q;
    TrueResidual      = std::sqrt(norm2(r) / ssq);
    IterationsToComplete = MaxIterations;

    std::cout << GridLogMessage
              << "MultiSplittingPCG did NOT converge after " << MaxIterations << " iterations"
              << "  true res " << TrueResidual << std::endl;
    printTimers();
    if (ErrorOnNoConverge) GRID_ASSERT(0);
  }

private:
  // Apply composite preconditioner: out = sum_i w_i * M_i^{-1}(in).
  // scratch[i] is dedicated to sub-preconditioner i — no sharing.
  void applyPreconditioner(const Field& in, Field& out, std::vector<Field>& scratch)
  {
    PrecTimer.Start();
    int k = (int)Subpreconditioners.size();

    if (parallelApply) {
      // Launch each sub-solve in its own thread.  The threads share only
      // read-only data (in, Subpreconditioners, Weights); each writes
      // exclusively to scratch[i].
      std::vector<std::thread> threads(k);
      for (int i = 0; i < k; i++) {
        threads[i] = std::thread([this, &in, &scratch, i]() {
          (*Subpreconditioners[i])(in, scratch[i]);
        });
      }
      for (auto& t : threads) t.join();
    } else {
      for (int i = 0; i < k; i++)
        (*Subpreconditioners[i])(in, scratch[i]);
    }

    out = Zero();
    for (int i = 0; i < k; i++)
      out = out + Weights[i] * scratch[i];

    PrecTimer.Stop();
  }

  void printTimers()
  {
    std::cout << GridLogMessage << "MultiSplittingPCG Time: Total  " << SolverTimer.Elapsed() << std::endl;
    std::cout << GridLogMessage << "MultiSplittingPCG Time: Matrix " << MatrixTimer.Elapsed() << std::endl;
    std::cout << GridLogMessage << "MultiSplittingPCG Time: Prec   " << PrecTimer.Elapsed()   << std::endl;
    std::cout << GridLogMessage << "MultiSplittingPCG Time: Linalg " << LinalgTimer.Elapsed() << std::endl;
  }
};

NAMESPACE_END(Grid);

#endif // GRID_MULTI_SPLITTING_PRECONDITIONED_CG_H
