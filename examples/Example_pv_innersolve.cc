/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pv_innersolve.cc

    Will FEAST's inner solve for the generalised problem  M x = lambda PV x
    actually converge, or stall?

    FEAST's cost is entirely the Ne*M0 shifted solves per iteration.  For the
    generalised form the quadrature needs

        (z PV - M)^{-1}

    at each contour node z.  This probe measures exactly that, node by node,
    with no eigensolver around it -- iterations, wall time, achieved residual.

    The specific worry: on a small contour near the origin |z| is small, so
    z PV - M  ~  -M, i.e. the inner solve degenerates into solving the *light*
    Dirac operator, which is precisely what one would be deflating in the
    first place.  Measure rather than assume.

    For reference the same probe is run for the standard FEAST operator
    (z - PV^dag M) on its own contour.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

namespace Grid {

struct ProbeParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(ProbeParameters,
    RealD,   mass,
    RealD,   M5,
    Integer, Ls,
    Integer, Ne,
    RealD,   EmidRe,
    RealD,   EmidIm,
    RealD,   radius,
    RealD,   innerTol,
    Integer, innerMaxIt,
    Integer, gmresRestart,
    Integer, doPVdagM)

  ProbeParameters()
    : mass(0.01), M5(1.8), Ls(8), Ne(8),
      EmidRe(0.05), EmidIm(0.0), radius(0.15),
      innerTol(1e-5), innerMaxIt(3000), gmresRestart(50), doPVdagM(1)
  {}
};

} // namespace Grid

typedef DomainWallFermionD  DwfOp;
typedef DwfOp::FermionField FermionField;

/// PV^dag M, for the standard-problem reference probe.
template<class Matrix, class Field>
class PVdagMOp : public LinearOperatorBase<Field> {
  Matrix &_M, &_PV;
public:
  PVdagMOp(Matrix& M, Matrix& PV) : _M(M), _PV(PV) {}
  void Op(const Field& in, Field& out) {
    Field t(in.Grid()); _M.M(in, t); _PV.Mdag(t, out);
  }
  void AdjOp(const Field& in, Field& out) {
    Field t(in.Grid()); _PV.M(in, t); _M.Mdag(t, out);
  }
  void OpDiag(const Field&, Field&)                        { assert(0); }
  void OpDir (const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)         { assert(0); }
  void HermOp(const Field&, Field&)                        { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&) { assert(0); }
};


/**
 * Normal-equations wrapper: presents S^dag S as a Hermitian positive
 * operator so ConjugateGradient can be used on a LinearOperatorBase that
 * only offers Op/AdjOp.  Grid's MdagMLinearOperator needs a
 * SparseMatrixBase (M/Mdag/MdagM), which a shifted operator is not.
 */
template<class Field>
class NormalEqOp : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field>& S_;
public:
  NormalEqOp(LinearOperatorBase<Field>& S) : S_(S) {}
  void HermOp(const Field& in, Field& out) {
    Field t(in.Grid());
    S_.Op(in, t);
    S_.AdjOp(t, out);
  }
  void HermOpAndNorm(const Field& in, Field& out, RealD& n1, RealD& n2) {
    HermOp(in, out);
    ComplexD dot = innerProduct(in, out);
    n1 = real(dot);
    n2 = norm2(out);
  }
  void Op   (const Field& in, Field& out) { HermOp(in, out); }
  void AdjOp(const Field& in, Field& out) { HermOp(in, out); }
  void OpDiag(const Field&, Field&)                { assert(0); }
  void OpDir (const Field&, Field&, int, int)      { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&) { assert(0); }
};

/**
 * AA^dag wrapper, for CGNE (Craig's method): solve (S S^dag) y = b, x = S^dag y.
 * Minimises the error ||x - x*|| rather than the residual.
 */
template<class Field>
class NormalEqOpAAdag : public LinearOperatorBase<Field> {
  LinearOperatorBase<Field>& S_;
public:
  NormalEqOpAAdag(LinearOperatorBase<Field>& S) : S_(S) {}
  void HermOp(const Field& in, Field& out) {
    Field t(in.Grid());
    S_.AdjOp(in, t);
    S_.Op(t, out);
  }
  void HermOpAndNorm(const Field& in, Field& out, RealD& n1, RealD& n2) {
    HermOp(in, out);
    ComplexD dot = innerProduct(in, out);
    n1 = real(dot); n2 = norm2(out);
  }
  void Op   (const Field& in, Field& out) { HermOp(in, out); }
  void AdjOp(const Field& in, Field& out) { HermOp(in, out); }
  void OpDiag(const Field&, Field&)                { assert(0); }
  void OpDir (const Field&, Field&, int, int)      { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&) { assert(0); }
};

/**
 * One shifted solve, reporting what actually happened.
 * op must already have had setShift() called.
 */
static void probeSolve(LinearOperatorBase<FermionField>& shifted,
                       const FermionField& rhs, GridBase* grid,
                       RealD tol, int maxit, int restart,
                       const std::string& tag)
{
  FermionField x(grid); x = Zero();
  GeneralisedMinimalResidual<FermionField> gmres(tol, maxit, restart,
                                                 /*err_on_no_conv=*/false);
  double t0 = usecond();
  gmres(shifted, rhs, x);
  double t1 = usecond();

  // true relative residual, independent of what the solver claims
  FermionField r(grid);
  shifted.Op(x, r);
  r = r - rhs;
  RealD relres = std::sqrt(norm2(r) / norm2(rhs));

  std::cout << GridLogMessage << "    " << tag
            << "  wall = " << (t1 - t0)/1e6 << " s"
            << "   true rel.res = " << relres
            << "   ||x|| = " << std::sqrt(norm2(x)) << std::endl;
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  ProbeParameters P;
  { XmlReader rd("ProbeParams.xml"); read(rd, "ProbeParameters", P); }
  std::cout << GridLogMessage << P << std::endl;

  const int Ls = P.Ls;
  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, std::string("config"));

  DwfOp Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, P.mass, P.M5);
  DwfOp Dpv (Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, 1.0,    P.M5);

  NonHermitianLinearOperator<DwfOp, FermionField> Mop(Ddwf);
  NonHermitianLinearOperator<DwfOp, FermionField> PVop(Dpv);

  FermionField rhs(FGrid); random(RNG5, rhs);

  // Contour nodes: build a FeastNonHermitian purely to reuse its quadrature.
  std::complex<double> Emid(P.EmidRe, P.EmidIm);
  FeastNonHermitian<FermionField> F(Mop, FGrid, Emid, P.radius, P.Ne, 1,
                                    1e-7, 1, 1.0, 0.0, FeastGaussLegendre);
  auto nodes = F.getNodes();

  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " Inner-solve probe   Ls=" << Ls
            << "  mass=" << P.mass << "  M5=" << P.M5 << std::endl;
  std::cout << GridLogMessage << " contour Emid=" << Emid << " radius=" << P.radius
            << " Ne=" << P.Ne << std::endl;
  std::cout << GridLogMessage << " GMRES tol=" << P.innerTol
            << " maxit=" << P.innerMaxIt << " restart=" << P.gmresRestart << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;

  // ---- generalised: (z PV - M) ----
  std::cout << GridLogMessage << "--- (z PV - M)  [generalised, gives PV^{-1}M] ---" << std::endl;
  {
    FermionField scratch(FGrid);
    FeastShiftedOperator<FermionField> shifted(Mop, &PVop, &scratch);
    for (int k = 0; k < (int)nodes.size(); k++) {
      ComplexD z(nodes[k].real(), nodes[k].imag());
      std::cout << GridLogMessage << "  node " << k << "  z = " << nodes[k] << std::endl;
      shifted.setShift(z, false);
      probeSolve(shifted, rhs, FGrid, P.innerTol, P.innerMaxIt, P.gmresRestart, "right (zPV-M)   ");
      shifted.setShift(z, true);
      probeSolve(shifted, rhs, FGrid, P.innerTol, P.innerMaxIt, P.gmresRestart, "left  (z*PV-M)^H");
    }
  }

  // ---- reference: plain M, the |z|->0 limit of the above ----
  std::cout << GridLogMessage << "--- reference: plain M solve (the |z|->0 limit) ---" << std::endl;
  {
    FermionField scratch(FGrid);
    FeastShiftedOperator<FermionField> shifted(Mop, nullptr, &scratch);
    shifted.setShift(ComplexD(0.0, 0.0), false);   // represents -M
    probeSolve(shifted, rhs, FGrid, P.innerTol, P.innerMaxIt, P.gmresRestart, "(-M)            ");
  }

  // ---- standard problem reference: (z - PV^dag M) ----
  if (P.doPVdagM) {
    std::cout << GridLogMessage << "--- (z - PV^dag M)  [standard problem] ---" << std::endl;
    PVdagMOp<DwfOp, FermionField> PVdagM(Ddwf, Dpv);
    FermionField scratch(FGrid);
    FeastShiftedOperator<FermionField> shifted(PVdagM, nullptr, &scratch);
    for (int k = 0; k < (int)nodes.size(); k++) {
      ComplexD z(nodes[k].real(), nodes[k].imag());
      std::cout << GridLogMessage << "  node " << k << "  z = " << nodes[k] << std::endl;
      shifted.setShift(z, false);
      probeSolve(shifted, rhs, FGrid, P.innerTol, P.innerMaxIt, P.gmresRestart, "right (z-PVdagM)");
    }
  }


  // ---- structural check:  z PV - M  ==  (z-1) D(m_eff),  m_eff=(z-m)/(z-1) ----
  //
  // DWF is linear in the mass, D(m) = A + mB, so
  //   z D(1) - D(m) = (z-1)A + (z-m)B = (z-1)[A + m_eff B] = (z-1) D(m_eff).
  // If true, the FEAST inner solve is just a DWF solve at an effective
  // (generally complex) mass -- which brings even-odd preconditioning and
  // the rest of the standard DWF toolchain into play.
  // Tested with REAL z so that m_eff is real and a stock DomainWallFermion
  // can represent D(m_eff); linearity in m then gives the complex case.
  {
    std::cout << GridLogMessage << "--- identity check: zPV - M == (z-1) D(m_eff) ---" << std::endl;
    FermionField x(FGrid); random(RNG5, x);
    FermionField lhs(FGrid), rhs2(FGrid), t(FGrid);
    for (RealD zr : {2.0, 5.0, -1.0, 0.5}) {
      RealD meff = (zr - P.mass) / (zr - 1.0);
      DwfOp Deff(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, meff, P.M5);
      // lhs = z PV x - M x
      Dpv.M(x, t); Ddwf.M(x, lhs); lhs = zr * t - lhs;
      // rhs = (z-1) D(meff) x
      Deff.M(x, rhs2); rhs2 = (zr - 1.0) * rhs2;
      FermionField d = lhs - rhs2;
      std::cout << GridLogMessage << "    z = " << zr
                << "   m_eff = " << meff
                << "   ||lhs-rhs||/||lhs|| = " << std::sqrt(norm2(d)/norm2(lhs)) << std::endl;
    }
  }

  // ---- CGNE on the stalled nodes ----
  // GMRES(restart) stagnates flat at ~0.77; CG on the normal equations
  // cannot break down, at the cost of squaring the condition number.
  {
    std::cout << GridLogMessage << "--- CGNR  (S^dag S x = S^dag b, minimises residual) ---" << std::endl;
    FermionField scratch(FGrid);
    FeastShiftedOperator<FermionField> shifted(Mop, &PVop, &scratch);
    for (int k = 0; k < (int)nodes.size(); k++) {
      ComplexD z(nodes[k].real(), nodes[k].imag());
      shifted.setShift(z, false);
      NormalEqOp<FermionField> NE(shifted);
      FermionField b2(FGrid), xx(FGrid);
      shifted.AdjOp(rhs, b2);          // b2 = S^dag rhs
      xx = Zero();
      ConjugateGradient<FermionField> cg(P.innerTol, P.innerMaxIt, false);
      double t0 = usecond();
      cg(NE, b2, xx);
      double t1 = usecond();
      FermionField r(FGrid);
      shifted.Op(xx, r); r = r - rhs;
      std::cout << GridLogMessage << "    node " << k << " z=" << nodes[k]
                << "  CGNR wall = " << (t1-t0)/1e6 << " s"
                << "   true rel.res = " << std::sqrt(norm2(r)/norm2(rhs)) << std::endl;
    }
  }


  // ---- CGNE (Craig): (S S^dag) y = b, x = S^dag y.  Minimises the error. ----
  {
    std::cout << GridLogMessage << "--- CGNE  (S S^dag y = b, x = S^dag y, minimises error) ---" << std::endl;
    FermionField scratch(FGrid);
    FeastShiftedOperator<FermionField> shifted(Mop, &PVop, &scratch);
    for (int k = 0; k < (int)nodes.size(); k++) {
      ComplexD z(nodes[k].real(), nodes[k].imag());
      shifted.setShift(z, false);
      NormalEqOpAAdag<FermionField> NE(shifted);
      FermionField y(FGrid), xx(FGrid);
      y = Zero();
      ConjugateGradient<FermionField> cg(P.innerTol, P.innerMaxIt, false);
      double t0 = usecond();
      cg(NE, rhs, y);
      shifted.AdjOp(y, xx);        // x = S^dag y
      double t1 = usecond();
      FermionField r(FGrid);
      shifted.Op(xx, r); r = r - rhs;
      std::cout << GridLogMessage << "    node " << k << " z=" << nodes[k]
                << "  CGNE wall = " << (t1-t0)/1e6 << " s"
                << "   true rel.res = " << std::sqrt(norm2(r)/norm2(rhs)) << std::endl;
    }
  }

  std::cout << GridLogMessage << std::endl << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
