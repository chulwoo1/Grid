/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/CountingLinearOperator.h

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

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_COUNTING_LINEAR_OPERATOR_H
#define GRID_COUNTING_LINEAR_OPERATOR_H

NAMESPACE_BEGIN(Grid);

/**
 * Pass-through LinearOperator that counts how many times it is applied.
 *
 * Wrap the operator, hand the wrapper to the solver, and read the counters
 * afterwards.  Counting at the operator rather than inside each solver is the
 * point: it catches every application wherever it happens, including those
 * buried inside an inner solve (a shifted GMRES/CG inside FEAST, a PV solve
 * inside PV^{-1}M), which no solver-level iteration count reports.
 *
 * Why bother when wall time is right there: wall time conflates algorithmic
 * cost with machine conditions.  A GPU shared with another job has produced a
 * 100x swing in this tree for an unchanged calculation, and inner-solver
 * tuning moves wall time without moving the algorithm.  The operator count is
 * reproducible and is the number to compare across methods and machines.
 *
 * costPerCall converts one Op()/AdjOp() into applications of the underlying
 * sparse matrix, so the reported matvecs() is the quantity that costs time:
 * 1 for a plain Dirac operator, 2 for MdagM or PV^dag M (two applications per
 * call).  HermOp is counted as two calls on the assumption HermOp = AdjOp(Op),
 * which holds for the Mdag M style operators; for an operator where it does
 * not, read nOp/nAdjOp/nHermOp directly instead of calls().
 */
template<class Field>
class CountingLinearOperator : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field>& _A;
  int _costPerCall;
public:
  mutable long nOp = 0, nAdjOp = 0, nHermOp = 0;

  CountingLinearOperator(LinearOperatorBase<Field>& A, int costPerCall = 1)
    : _A(A), _costPerCall(costPerCall) {}

  void reset() { nOp = 0; nAdjOp = 0; nHermOp = 0; }

  /// Calls to the wrapped operator; the natural unit for comparing solvers.
  long calls() const { return nOp + nAdjOp + 2*nHermOp; }
  /// Sparse-matrix applications; the unit that costs wall time.
  long matvecs() const { return _costPerCall * calls(); }
  int  costPerCall() const { return _costPerCall; }

  void Op    (const Field& in, Field& out) { nOp++;     _A.Op(in, out); }
  void AdjOp (const Field& in, Field& out) { nAdjOp++;  _A.AdjOp(in, out); }
  void HermOp(const Field& in, Field& out) { nHermOp++; _A.HermOp(in, out); }

  void OpDiag(const Field& in, Field& out)                   { _A.OpDiag(in, out); }
  void OpDir (const Field& in, Field& out, int d, int disp)   { _A.OpDir(in, out, d, disp); }
  void OpDirAll(const Field& in, std::vector<Field>& out)     { _A.OpDirAll(in, out); }
  void HermOpAndNorm(const Field& in, Field& out, RealD& n1, RealD& n2)
  { nHermOp++; _A.HermOpAndNorm(in, out, n1, n2); }
};

NAMESPACE_END(Grid);

#endif
