/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_pv_compare.cc

    Compare eigenvector generation for  PV^dag M  versus  PV^{-1} M
    for domain wall fermions.

    Background
    ----------
    The physically interesting operator is PV^{-1} M (Pauli-Villars
    preconditioned DWF), whose spectrum relates to the 4D effective Dirac
    operator.  Applying it costs a PV solve per matvec.

    PV^dag M is far cheaper -- two sparse applications -- and is sometimes
    used in its place.  The two are similar ONLY if PV^dag PV is proportional
    to the identity, i.e. PV proportional to unitary:

        PV^dag M = (PV^dag PV) PV^{-1} M

    So the spectral spread of PV^dag PV bounds how good the surrogate can
    possibly be.  Stage 1 below measures exactly that, which is far cheaper
    than running two eigensolvers and may settle the question outright.

    Stages (select with <stage> in PVParams.xml):
      1 : PV^dag PV spectral range + operator sanity checks  (cheap)
      2 : stage 1, plus Krylov-Schur on PV^dag M and on PV^{-1} M
      3 : stage 2, plus FEAST on PV^dag M (standard) and on (M, PV)
          as a generalised problem  M x = lambda PV x

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See the full license in the file "LICENSE" in the top level distribution
    directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>
#include <sstream>

using namespace std;
using namespace Grid;

namespace Grid {

struct PVParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(PVParameters,
    RealD,   mass,
    RealD,   M5,
    Integer, Ls,
    Integer, stage,
    RealD,   pvSolverTol,   // tolerance of the inner PV solve for PV^{-1}M
    Integer, pvSolverMaxIt,
    Integer, powerIters,
    // --- Krylov-Schur on PV^dag M ---
    Integer, ksNk,
    Integer, ksNp,
    Integer, ksNstop,
    Integer, ksMaxIter,
    RealD,   ksResid,
    // --- FEAST ---
    Integer, feastNe,
    Integer, feastM0,
    Integer, feastMaxIter,
    Integer, gmresRestart,
    Integer, gmresMaxIt,
    // contour for the generalised run (physical PV^{-1}M low modes)
    RealD,   genEmidRe,
    RealD,   genEmidIm,
    RealD,   genRadius,
    // --- head to head ---
    RealD,   h2hResTarget,  // ||A v - lambda v||/||v|| counted as converged
    Integer, feastInner,    // 0 = GMRES, 1 = BiCGSTAB, 2 = CGNR
    // The inner stopping condition and iteration cap are the dominant FEAST
    // tuning knobs: the contour nodes sit close to the spectrum, so (z-A) is
    // near-singular and nearly all the time is spent in these solves.  Both
    // are swept as a cross product so the trade-off is visible in one run.
    // A cap that bites is not a failure -- FEAST tolerates inexact solves
    // (IFEAST); the question is whether the outer iteration still converges.
    std::vector<RealD>,   feastInnerTols,
    std::vector<Integer>, feastInnerMaxIts,
    // Contour radius = feastRadiusFactor * (max distance of a KS eigenvalue
    // from their centroid).  This trades two things directly against each
    // other: a larger factor pushes the quadrature nodes away from the
    // spectrum, so (z-A) is better conditioned and the inner solves are much
    // cheaper, but it encloses more eigenvalues and so needs a larger M0.
    RealD,   feastRadiusFactor)

  PVParameters()
    : mass(0.01), M5(1.8), Ls(8), stage(1),
      pvSolverTol(1e-10), pvSolverMaxIt(10000), powerIters(200),
      ksNk(32), ksNp(32), ksNstop(16), ksMaxIter(100), ksResid(1e-8),
      feastNe(8), feastM0(12), feastMaxIter(3),
      gmresRestart(50), gmresMaxIt(3000),
      genEmidRe(0.05), genEmidIm(0.0), genRadius(0.15),
      h2hResTarget(1e-6), feastInner(2),
      feastInnerTols({1e-2, 1e-3}), feastInnerMaxIts({200, 1000}),
      feastRadiusFactor(1.5)
  {}
};

} // namespace Grid

typedef DomainWallFermionD           DwfOp;
typedef DwfOp::FermionField          FermionField;

/**
 * PV^dag M : two sparse applications, no solve.
 */
template<class Matrix, class Field>
class PVdagMOp : public LinearOperatorBase<Field> {
  Matrix &_M, &_PV;
public:
  PVdagMOp(Matrix& M, Matrix& PV) : _M(M), _PV(PV) {}
  void Op(const Field& in, Field& out) {
    Field t(in.Grid());
    _M.M(in, t);
    _PV.Mdag(t, out);
  }
  void AdjOp(const Field& in, Field& out) {
    // (PV^dag M)^dag = M^dag PV
    Field t(in.Grid());
    _PV.M(in, t);
    _M.Mdag(t, out);
  }
  void OpDiag(const Field&, Field&)                        { assert(0); }
  void OpDir (const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)         { assert(0); }
  void HermOp(const Field&, Field&)                        { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&) { assert(0); }
};

/**
 * PV^{-1} M : a PV solve per matvec, via CG on the normal equations
 * (PV^dag PV) x = PV^dag b.  This is the expensive reference.
 */
template<class Matrix, class Field>
class PVinvMOp : public LinearOperatorBase<Field> {
  Matrix &_M, &_PV;
  RealD   _tol;
  Integer _maxit;
public:
  mutable long nSolves = 0;
  mutable long nIters  = 0;
  PVinvMOp(Matrix& M, Matrix& PV, RealD tol, Integer maxit)
    : _M(M), _PV(PV), _tol(tol), _maxit(maxit) {}

  void applyPVinv(const Field& in, Field& out) const {
    Field rhs(in.Grid());
    _PV.Mdag(in, rhs);
    MdagMLinearOperator<Matrix, Field> HermOp(_PV);
    ConjugateGradient<Field> CG(_tol, _maxit, /*err_on_no_conv=*/false);
    out = Zero();
    CG(HermOp, rhs, out);
    nSolves++;
  }
  void Op(const Field& in, Field& out) {
    Field t(in.Grid());
    _M.M(in, t);
    applyPVinv(t, out);
  }
  void AdjOp(const Field& in, Field& out) {
    // (PV^{-1} M)^dag = M^dag PV^{-dag};  PV^{-dag} b solves PV^dag x = b
    Field t(in.Grid());
    MdagMLinearOperator<Matrix, Field> HermOp(_PV);
    ConjugateGradient<Field> CG(_tol, _maxit, false);
    Field y(in.Grid()); y = Zero();
    CG(HermOp, in, y);       // y = (PV^dag PV)^{-1} in
    _PV.M(y, t);             // t = PV y = PV^{-dag} in
    _M.Mdag(t, out);
    nSolves++;
  }
  void OpDiag(const Field&, Field&)                        { assert(0); }
  void OpDir (const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)         { assert(0); }
  void HermOp(const Field&, Field&)                        { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&) { assert(0); }
};


// The counting wrapper now lives in the library as CountingLinearOperator
// (Grid/algorithms/CountingLinearOperator.h), shared with Example_reuse.cc.
template<class Field> using CountingOp = CountingLinearOperator<Field>;

/**
 * One row of the head-to-head table.  Residuals are recomputed here with a
 * single definition, ||A v - lambda v|| / ||v||, rather than trusting each
 * solver's internal estimate -- KS reports a Ritz estimate and FEAST its own
 * block residual, and they are not the same quantity.
 */
struct CostRow {
  std::string method;
  double      wall    = 0.0;   // seconds
  long        calls   = 0;     // operator applications
  long        matvecs = 0;     // sparse applications
  int         found   = 0;     // eigenpairs returned
  int         conv    = 0;     // of those, meeting resTarget
  double      bestRes = -1.0;
  double      worstRes = -1.0;
};

static void printCostTable(const std::vector<CostRow>& rows, double resTarget)
{
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "================ head to head ================" << std::endl;
  std::cout << GridLogMessage << "  converged means ||A v - lambda v||/||v|| < "
            << resTarget << ", recomputed identically for both methods" << std::endl;
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage
            << "  method                        wall_s     calls   matvecs  found conv"
            << "   mv/pair   best_res  worst_res" << std::endl;
  for (const auto& r : rows) {
    std::cout << GridLogMessage << "  " << std::left << std::setw(28) << r.method
              << std::right
              << std::setw(9)  << std::fixed << std::setprecision(2) << r.wall
              << std::setw(10) << r.calls
              << std::setw(10) << r.matvecs
              << std::setw(7)  << r.found
              << std::setw(5)  << r.conv
              << std::setw(10) << (r.conv > 0 ? r.matvecs / r.conv : 0)
              << std::setw(11) << std::scientific << std::setprecision(2) << r.bestRes
              << std::setw(11) << r.worstRes
              << std::endl;
  }
  std::cout << GridLogMessage << "==============================================" << std::endl;
}

/**
 * Recompute residuals for a set of eigenpairs under one definition and fill in
 * the convergence columns of a row.
 */
template<class Field>
static void scorePairs(LinearOperatorBase<Field>& A,
                       const Eigen::VectorXcd& evals,
                       const std::vector<Field>& evecs,
                       double resTarget, CostRow& row,
                       const std::string& label)
{
  int n = std::min((int)evals.size(), (int)evecs.size());
  row.found = n;
  row.conv  = 0;
  row.bestRes = 1e30; row.worstRes = 0.0;
  if (n == 0) { row.bestRes = -1.0; row.worstRes = -1.0; return; }

  Field w(evecs[0].Grid());
  for (int i = 0; i < n; i++) {
    RealD nv = std::sqrt(norm2(evecs[i]));
    if (nv < 1e-30) continue;
    A.Op(evecs[i], w);
    ComplexD lam(evals(i).real(), evals(i).imag());
    w = w - lam * evecs[i];
    RealD r = std::sqrt(norm2(w)) / nv;
    row.bestRes  = std::min(row.bestRes,  (double)r);
    row.worstRes = std::max(row.worstRes, (double)r);
    if (r < resTarget) row.conv++;
    std::cout << GridLogMessage << "    " << label << "[" << std::setw(2) << i
              << "]  lambda = " << evals(i) << "   |res| = " << r
              << (r < resTarget ? "" : "   (not converged)") << std::endl;
  }
}

/**
 * Cross-match two eigenvalue lists.  Reports how many of each side's values
 * the other side also found, which is what decides whether the cost figures
 * are comparable at all: a cheap run that found half the spectrum is not
 * cheaper, it answered a smaller question.
 */
static void matchEvals(const Eigen::VectorXcd& a, const std::string& la,
                       const Eigen::VectorXcd& b, const std::string& lb,
                       double rtol)
{
  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "--- eigenvalue cross-match, relative tol " << rtol
            << " ---" << std::endl;
  auto countIn = [&](const Eigen::VectorXcd& x, const Eigen::VectorXcd& y) {
    int hit = 0;
    for (int i = 0; i < x.size(); i++) {
      double scale = std::max(1e-30, std::abs(x(i)));
      for (int j = 0; j < y.size(); j++)
        if (std::abs(x(i) - y(j)) / scale < rtol) { hit++; break; }
    }
    return hit;
  };
  int aInB = countIn(a, b), bInA = countIn(b, a);
  std::cout << GridLogMessage << "  " << la << " values also found by " << lb
            << " : " << aInB << " / " << a.size() << std::endl;
  std::cout << GridLogMessage << "  " << lb << " values also found by " << la
            << " : " << bInA << " / " << b.size() << std::endl;
  if (aInB < a.size() || bInA < b.size())
    std::cout << GridLogMessage
              << "  NOTE: the two did not answer the same question; compare costs with care"
              << std::endl;
}

/// L2 modified Gram-Schmidt in place; returns the rank actually retained.
template<class Field>
int orthonormalise(std::vector<Field>& v)
{
  int keep = 0;
  for (size_t j = 0; j < v.size(); j++) {
    for (int i = 0; i < keep; i++) {
      ComplexD c = innerProduct(v[i], v[j]);
      v[j] = v[j] - c * v[i];
    }
    RealD n = std::sqrt(norm2(v[j]));
    if (n > 1e-10) { v[j] = v[j] * (1.0/n); std::swap(v[keep], v[j]); keep++; }
  }
  return keep;
}

/**
 * Principal angles between span(U) and span(V): singular values of U^dag V.
 * Values near 1 are shared directions.  This, not an eigenvalue comparison,
 * is what decides whether one operator's low modes can deflate the other's.
 */
template<class Field>
void subspaceOverlap(std::vector<Field> U, std::vector<Field> V,
                     const std::string& label)
{
  int nu = orthonormalise(U);
  int nv = orthonormalise(V);
  int k  = std::min(nu, nv);
  if (k == 0) { std::cout << GridLogMessage << label << ": empty" << std::endl; return; }

  Eigen::MatrixXcd M(nu, nv);
  for (int i = 0; i < nu; i++)
    for (int j = 0; j < nv; j++)
      M(i, j) = toStdCmplx(innerProduct(U[i], V[j]));

  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M);
  auto sv = svd.singularValues();

  std::cout << GridLogMessage << "=== " << label << " ===" << std::endl;
  std::cout << GridLogMessage << "  dims: " << nu << " x " << nv << std::endl;
  std::cout << GridLogMessage << "  cos(principal angles), 1 = shared direction:" << std::endl;
  int nshared = 0;
  for (int i = 0; i < sv.size(); i++) {
    if (sv(i) > 0.9) nshared++;
    std::cout << GridLogMessage << "    sigma[" << i << "] = " << sv(i)
              << "   angle = " << std::acos(std::min(1.0, (double)sv(i))) * 180.0/M_PI
              << " deg" << std::endl;
  }
  std::cout << GridLogMessage << "  directions with cos > 0.9 : " << nshared
            << " / " << k << std::endl;
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  PVParameters P;
  { XmlReader rd("PVParams.xml"); read(rd, "PVParameters", P); }
  std::cout << GridLogMessage << P << std::endl;

  const int Ls = P.Ls;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);

  std::vector<int> seeds4({1,2,3,4}), seeds5({5,6,7,8});
  GridParallelRNG RNG4(UGrid); RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  NerscIO::readConfiguration(Umu, header, std::string("config"));
  std::cout << GridLogMessage << "Loaded gauge configuration" << std::endl;

  // M at the physical mass, PV at mass = 1
  DwfOp Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, P.mass, P.M5);
  DwfOp Dpv (Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, 1.0,    P.M5);

  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " PV^dag M  vs  PV^{-1} M      Ls=" << Ls
            << "  mass=" << P.mass << "  M5=" << P.M5 << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;

  // ================= Stage 1 : how far is PV from unitary? =================
  //
  // PV^dag M = (PV^dag PV) PV^{-1} M, so PV^dag M is a faithful stand-in for
  // PV^{-1} M only to the extent that PV^dag PV ~ const * I.  Measure the
  // spectral range of the Hermitian positive-definite operator PV^dag PV.
  {
    MdagMLinearOperator<DwfOp, FermionField> PVdagPV(Dpv);
    FermionField src(FGrid); random(RNG5, src);

    PowerMethod<FermionField> PM;
    RealD lmax = PM(PVdagPV, src);

    // smallest eigenvalue by inverse iteration: one CG solve per step
    FermionField v(FGrid), w(FGrid);
    random(RNG5, v);
    RealD nv = std::sqrt(norm2(v)); v = v * (1.0/nv);
    RealD lmin = 0.0;
    ConjugateGradient<FermionField> CG(1e-10, P.pvSolverMaxIt, false);
    for (int it = 0; it < 20; it++) {
      w = Zero();
      CG(PVdagPV, v, w);                       // w = (PV^dag PV)^{-1} v
      RealD nw = std::sqrt(norm2(w));
      v = w * (1.0/nw);
      lmin = 1.0/nw;                           // Rayleigh-ish estimate
    }

    std::cout << GridLogMessage << "--- Stage 1: PV^dag PV spectral range ---" << std::endl;
    std::cout << GridLogMessage << "  lambda_max(PV^dag PV) = " << lmax << std::endl;
    std::cout << GridLogMessage << "  lambda_min(PV^dag PV) = " << lmin << std::endl;
    std::cout << GridLogMessage << "  spread (max/min)      = "
              << (lmin > 0 ? lmax/lmin : -1.0) << std::endl;
    std::cout << GridLogMessage
              << "  => PV^dag M is a faithful stand-in for PV^{-1}M only if this"
              << " spread is ~1." << std::endl;

    // sanity: PV^{-1}(PV x) == x
    PVinvMOp<DwfOp, FermionField> PVinvM(Ddwf, Dpv, P.pvSolverTol, P.pvSolverMaxIt);
    FermionField x(FGrid), px(FGrid), y(FGrid);
    random(RNG5, x);
    Dpv.M(x, px);
    PVinvM.applyPVinv(px, y);
    y = y - x;
    std::cout << GridLogMessage << "  sanity ||PV^{-1}(PV x) - x||/||x|| = "
              << std::sqrt(norm2(y)/norm2(x)) << std::endl;

    // how different are the two operators on a random vector?
    PVdagMOp<DwfOp, FermionField> PVdagM(Ddwf, Dpv);
    FermionField a(FGrid), b(FGrid);
    PVdagM.Op(x, a);
    PVinvM.Op(x, b);
    RealD na = std::sqrt(norm2(a)), nb = std::sqrt(norm2(b));
    FermionField d = a*(1.0/na) - b*(1.0/nb);
    std::cout << GridLogMessage << "  ||PV^dag M x||/||PV^{-1}M x|| = " << na/nb << std::endl;
    std::cout << GridLogMessage << "  direction difference (normalised) = "
              << std::sqrt(norm2(d)) << std::endl;
  }


  // ================= Stage 2 : Krylov-Schur on PV^dag M =================
  // One counted PV^dag M, shared by stage 2 and stage 3a so that both sides
  // of the head-to-head are measured through the same instrument.  costPerOp
  // = 2: one Op() is M followed by PV^dag.
  PVdagMOp<DwfOp, FermionField> PVdagMbare(Ddwf, Dpv);
  CountingOp<FermionField>      PVdagM(PVdagMbare, 2);

  std::vector<FermionField> ksEvecs;
  Eigen::VectorXcd          ksEvals;
  std::vector<CostRow>      rows;
  const double resTarget = P.h2hResTarget;

  if (P.stage >= 2) {
    std::cout << GridLogMessage << "--- Stage 2: Krylov-Schur on PV^dag M ---" << std::endl;
    std::cout << GridLogMessage << "  Nk=" << P.ksNk << " Np=" << P.ksNp
              << " Nstop=" << P.ksNstop << " maxIter=" << P.ksMaxIter
              << " resid=" << P.ksResid << "  filter=EvalNormSmall" << std::endl;

    FermionField src(FGrid); random(RNG5, src);

    PVdagM.reset();
    double t0 = usecond();
    KrylovSchur<FermionField> KS(PVdagM, FGrid, P.ksResid, EvalNormSmall);
    KS.doEvalCheck = true;
    KS(src, P.ksMaxIter, P.ksNk + P.ksNp, P.ksNk, P.ksNstop);
    double t1 = usecond();

    ksEvecs = KS.evecs;
    ksEvals = KS.getEvals();   // evals itself is private

    CostRow r;
    r.method  = "KrylovSchur";
    r.wall    = (t1-t0)/1e6;
    r.calls   = PVdagM.calls();
    r.matvecs = PVdagM.matvecs();
    // Score only the Nstop that KS was asked to converge; the tail of the
    // basis is not a claim.
    {
      int nk = std::min((int)ksEvals.size(), (int)ksEvecs.size());
      nk = std::min(nk, (int)P.ksNstop);
      Eigen::VectorXcd ev = ksEvals.head(nk);
      std::vector<FermionField> vv(ksEvecs.begin(), ksEvecs.begin() + nk);
      scorePairs(PVdagMbare, ev, vv, resTarget, r, "ks");
      ksEvals = ev;
    }
    rows.push_back(r);

    std::cout << GridLogMessage << "  KS(PV^dag M) wall = " << r.wall
              << " s,  calls = " << r.calls << ",  matvecs = " << r.matvecs
              << ",  converged " << r.conv << "/" << r.found << std::endl;
  }

  // ================= Stage 3 : FEAST, both operators =================
  if (P.stage >= 3) {
    // Inner-tolerance sweep.  FEAST (IFEAST) tolerates inexact solves; the
    // question is how loose we can go before the outer rate degrades.
    // Everything is recorded, loosest-first.
    std::vector<RealD>   tols    = P.feastInnerTols;
    std::vector<Integer> maxIts  = P.feastInnerMaxIts;
    if (tols.empty())   tols   = {1e-3};
    if (maxIts.empty()) maxIts = {(Integer)P.gmresMaxIt};

    // --- contour for PV^dag M, taken from the measured KS spectrum ---
    std::complex<double> cen(0.0, 0.0);
    double rad = 1.0;
    if (ksEvals.size() > 0) {
      int n = std::min((int)ksEvals.size(), (int)P.ksNstop);
      for (int i = 0; i < n; i++) cen += std::complex<double>(ksEvals(i).real(), ksEvals(i).imag());
      cen /= (double)n;
      rad = 0.0;
      for (int i = 0; i < n; i++)
        rad = std::max(rad, std::abs(std::complex<double>(ksEvals(i).real(), ksEvals(i).imag()) - cen));
      rad *= P.feastRadiusFactor;
    }
    std::cout << GridLogMessage << "--- Stage 3a: FEAST on PV^dag M ---" << std::endl;
    std::cout << GridLogMessage << "  contour from KS: Emid=" << cen << " radius=" << rad
              << "  (factor " << P.feastRadiusFactor << " on the KS spread)" << std::endl;

    std::vector<FermionField> feastPVdagM;
    Eigen::VectorXcd          feastBestEvals;
    for (auto tol : tols) {
     for (auto mit : maxIts) {
      std::cout << GridLogMessage << "  [sweep] PV^dag M, innerTol = " << tol
                << "  innerMaxIt = " << mit << std::endl;
      FeastNonHermitian<FermionField> F(PVdagM, FGrid, cen, rad,
                                        P.feastNe, P.feastM0,
                                        1e-7, P.feastMaxIter, 1.0, 0.0,
                                        FeastGaussLegendre);
      std::string innerName;
      if      (P.feastInner == 2) { F.useCGNR(tol, mit);                      innerName = "CGNR"; }
      else if (P.feastInner == 1) { F.useBiCGSTAB(tol, mit);                  innerName = "BiCGSTAB"; }
      else                        { F.useGMRES(tol, mit, P.gmresRestart);     innerName = "GMRES"; }
      std::vector<FermionField> Y0(P.feastM0, FermionField(FGrid));
      for (int m = 0; m < (int)P.feastM0; m++) random(RNG5, Y0[m]);
      PVdagM.reset();
      F.resetCounters();
      double t0 = usecond();
      F(Y0);
      double t1 = usecond();
      auto ev = F.getEvals();

      CostRow r;
      std::ostringstream nm;
      nm << "FEAST " << innerName << " t=" << std::scientific << std::setprecision(0) << tol
         << " n=" << mit;
      r.method  = nm.str();
      r.wall    = (t1-t0)/1e6;
      r.calls   = PVdagM.calls();
      r.matvecs = PVdagM.matvecs();
      scorePairs(PVdagMbare, ev, F.getEvecsR(), resTarget, r, "feast");
      rows.push_back(r);

      long ns = F.nInnerSolves;
      std::cout << GridLogMessage << "  [sweep] PV^dag M tol=" << tol << " maxIt=" << mit
                << "  wall=" << r.wall << " s  outer=" << F.getIterations()
                << "  innerSolves=" << ns
                << "  calls/solve=" << (ns ? r.calls/ns : 0)
                << "  matvecs=" << r.matvecs
                << "  converged " << r.conv << "/" << r.found << std::endl;

      if (!F.getEvecsR().empty()) { feastPVdagM = F.getEvecsR(); feastBestEvals = ev; }
     }
    }

    // --- the head-to-head on PV^dag M: same operator, same residual test ---
    if (!rows.empty()) {
      printCostTable(rows, resTarget);
      if (ksEvals.size() > 0 && feastBestEvals.size() > 0)
        matchEvals(ksEvals, "KS", feastBestEvals, "FEAST", 1e-3);
    }

    // --- generalised:  M x = lambda PV x,  i.e. spectrum of PV^{-1} M ---
    std::cout << GridLogMessage << "--- Stage 3b: FEAST generalised (M, PV) = PV^{-1}M ---" << std::endl;
    std::complex<double> gcen(P.genEmidRe, P.genEmidIm);
    std::cout << GridLogMessage << "  contour: Emid=" << gcen
              << " radius=" << P.genRadius
              << "   (DWF PV^{-1}M is overlap-like: circle through 0 and 1,"
              << " physical low modes near the origin)" << std::endl;

    NonHermitianLinearOperator<DwfOp, FermionField> MopBare(Ddwf);
    NonHermitianLinearOperator<DwfOp, FermionField> PVopBare(Dpv);
    // Count M and PV separately: in the generalised form z PV - M, each
    // shifted application costs one of each, so the two counters should
    // track one another and their sum is the Dirac-application cost.
    CountingOp<FermionField> Mop (MopBare,  1);
    CountingOp<FermionField> PVop(PVopBare, 1);

    std::vector<FermionField> feastPVinvM;
    std::vector<CostRow>      genRows;
    for (auto tol : tols) {
     for (auto mit : maxIts) {
      std::cout << GridLogMessage << "  [sweep] generalised, innerTol = " << tol
                << "  innerMaxIt = " << mit << std::endl;
      FeastNonHermitian<FermionField> G(Mop, FGrid, gcen, P.genRadius,
                                        P.feastNe, P.feastM0,
                                        1e-7, P.feastMaxIter, 1.0, 0.0,
                                        FeastGaussLegendre);
      G.setB(PVop);                      // must precede the solver factory
      std::string innerName;
      if      (P.feastInner == 2) { G.useCGNR(tol, mit);                   innerName = "CGNR"; }
      else if (P.feastInner == 1) { G.useBiCGSTAB(tol, mit);               innerName = "BiCGSTAB"; }
      else                        { G.useGMRES(tol, mit, P.gmresRestart);  innerName = "GMRES"; }
      std::vector<FermionField> Y0(P.feastM0, FermionField(FGrid));
      for (int m = 0; m < (int)P.feastM0; m++) random(RNG5, Y0[m]);
      Mop.reset(); PVop.reset();
      G.resetCounters();
      double t0 = usecond();
      G(Y0);
      double t1 = usecond();
      auto ev = G.getEvals(); auto rs = G.getResiduals();

      CostRow r;
      std::ostringstream nm;
      nm << "FEAST gen " << innerName << " t=" << std::scientific << std::setprecision(0)
         << tol << " n=" << mit;
      r.method  = nm.str();
      r.wall    = (t1-t0)/1e6;
      r.calls   = Mop.calls() + PVop.calls();
      r.matvecs = Mop.matvecs() + PVop.matvecs();
      r.found   = (int)ev.size();
      r.conv    = 0;
      r.bestRes = rs.empty() ? -1.0 : 1e30; r.worstRes = rs.empty() ? -1.0 : 0.0;
      for (auto x : rs) {
        r.bestRes  = std::min(r.bestRes,  (double)x);
        r.worstRes = std::max(r.worstRes, (double)x);
        if (x < resTarget) r.conv++;
      }
      genRows.push_back(r);

      long ns = G.nInnerSolves;
      std::cout << GridLogMessage << "  [sweep] generalised tol=" << tol << " maxIt=" << mit
                << "  wall=" << r.wall << " s  outer=" << G.getIterations()
                << "  innerSolves=" << ns
                << "  calls/solve=" << (ns ? r.calls/ns : 0)
                << "  M=" << Mop.calls() << " PV=" << PVop.calls()
                << "  found=" << r.found << "  conv=" << r.conv << std::endl;
      for (int i = 0; i < (int)ev.size(); i++)
        std::cout << GridLogMessage << "      lambda=" << ev(i) << "  |res|=" << rs[i] << std::endl;
      if (!G.getEvecsR().empty()) feastPVinvM = G.getEvecsR();
     }
    }

    if (!genRows.empty()) {
      std::cout << GridLogMessage << std::endl;
      std::cout << GridLogMessage
                << "  NOTE: the table below is FEAST on PV^{-1}M, a DIFFERENT operator"
                << " from the PV^dag M table above. Costs are not comparable across"
                << " the two tables; residuals are FEAST's own." << std::endl;
      printCostTable(genRows, resTarget);
    }

    // --- the comparison that actually matters ---
    if (!ksEvecs.empty() && !feastPVinvM.empty()) {
      std::vector<FermionField> A(ksEvecs.begin(),
                                  ksEvecs.begin() + std::min((size_t)P.ksNstop, ksEvecs.size()));
      subspaceOverlap(A, feastPVinvM,
                      "subspace overlap: KS low modes of PV^dag M  vs  FEAST low modes of PV^{-1}M");
    }
    if (!feastPVdagM.empty() && !feastPVinvM.empty()) {
      subspaceOverlap(feastPVdagM, feastPVinvM,
                      "subspace overlap: FEAST PV^dag M  vs  FEAST PV^{-1}M (same method both sides)");
    }
  }

  std::cout << GridLogMessage << std::endl << "Done (stage " << P.stage << ")" << std::endl;

  Grid_finalize();
  return 0;
}
