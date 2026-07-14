/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/algorithms/iterative/HarmonicBlockKrylovSchur.h

Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
#define GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H

#include <iomanip>
#include <numeric>

NAMESPACE_BEGIN(Grid);

/**
 * Block shift-targeted Krylov-Schur eigensolver.
 *
 * This is BlockKrylovSchur with a single algorithmic change: the thick-restart
 * rotation targets eigenvalues near a shift sigma instead of the extremal Ritz
 * values.  Everything else (block Arnoldi, truncation, convergence check,
 * verification, parity/gamma5 seed expansion) is inherited unchanged from
 * BlockKrylovSchur; only restartRotation() is overridden.
 *
 * Algorithm
 * ---------
 * The block Arnoldi factorisation is
 *
 *   A V = V H + F B^dag                                         (1)
 *
 * with V orthonormal (Nm columns), H the Nm x Nm block upper-Hessenberg
 * Rayleigh quotient, F the Nblock residual vectors and B the Nm x Nblock
 * coupling matrix.
 *
 * Shift-targeted thick restart
 * ----------------------------
 * To target eigenvalues near shift sigma, the Schur decomposition is computed
 * for the shifted Rayleigh quotient:
 *
 *   (H - sigma I) = Q^dag S Q                                   (2)
 *
 * Sorting the Schur values of (H - sigma I) by smallest |S(i,i)| = |lambda - sigma|
 * and retaining the leading Nk is equivalent to selecting the Ritz values of H
 * closest to sigma.  Since Q diagonalises (H - sigma I) (and hence H itself),
 * the rotated Rayleigh quotient is exactly upper triangular:
 *
 *   H_new = Q H Q^dag = S + sigma I  (upper triangular)         (3)
 *
 * so truncation to Nk is exact, exactly as in the base class.
 *
 * Parameters
 * ----------
 * shift    : target shift sigma (default 0.0); Schur values sorted by |lambda - sigma|
 *
 * (all other parameters as in BlockKrylovSchur)
 *
 * Usage
 * -----
 *   HarmonicBlockKrylovSchur<Field> hbks(LinOp, Grid, tol, shift, EvalNormSmall);
 *   std::vector<Field> v0(Nblock, Field(Grid));
 *   // fill v0 with random starting vectors
 *   hbks(v0, maxIter, Nm, Nk, Nstop, Nblock);
 *   auto evals = hbks.getEvals();
 *   auto evecs = hbks.getEvecs();
 */
template<class Field>
class HarmonicBlockKrylovSchur : public BlockKrylovSchur<Field> {

protected:
  typedef BlockKrylovSchur<Field> Base;
  typedef typename Base::CMat CMat;
  typedef typename Base::CVec CVec;

  // Dependent-base members referenced in restartRotation()
  using Base::H;
  using Base::Nm;
  using Base::Nk;
  using Base::ritzFilter;
  using Base::className;

  ComplexD shift;       // target shift sigma

public:
  //--------------------------------------------------------------------
  // Constructor
  //--------------------------------------------------------------------
  HarmonicBlockKrylovSchur(LinearOperatorBase<Field>& _Linop, GridBase* _Grid,
                             RealD _Tolerance, ComplexD _shift = 0.0,
                             RitzFilter _rf = EvalNormSmall)
    : Base(_Linop, _Grid, _Tolerance, _rf), shift(_shift)
  {
    className = "HarmonicBlockKrylovSchur";
  }

protected:
  //--------------------------------------------------------------------
  // Restart rotation: shifted Schur decomposition
  //--------------------------------------------------------------------
  /**
   * Schur-decompose (H - sigma I), reorder so the Nk Ritz values nearest sigma
   * lead, and return Q with Hnew = Q H Q^dag = S + sigma I (upper triangular).
   */
  void restartRotation(CMat& Q, CMat& Hnew) override
  {
    int N = Nm;

    CMat Hshift = H - shift * CMat::Identity(N, N);
    ComplexSchurDecomposition schur(Hshift, false, ritzFilter);
    schur.schurReorder(Nk);

    CMat S = schur.getMatrixS();
    std::cout << GridLogMessage << className
              << ": Ritz values nearest shift (first Nk):" << std::endl;
    for (int i = 0; i < Nk; i++)
      std::cout << GridLogMessage << "  [" << i << "] " << S(i, i) + shift << std::endl;

    Q    = schur.getMatrixQ();
    Hnew = S + shift * CMat::Identity(N, N);
  }
};

NAMESPACE_END(Grid);

#endif  // GRID_HARMONIC_BLOCKED_KRYLOV_SCHUR_H
