/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_reff.cc

    Measure the effective rank r_eff of the restart residual between two
    gauge configurations.

    Motivation
    ----------
    Reusing a deflation / eigenvector subspace V across HMC trajectories
    requires restarting an eigensolver from an arbitrary orthonormal V.
    That is legitimate provided the exact block relation is carried:

        A V = V H + F B,    H = V^dag A V,   F orthonormal, F perp V,
                            B = F^dag (A V - V H)

    The cost of carrying it is the rank of the residual R = (I - V V^dag) A V,
    so the question that decides whether reuse is worth anything is:

        how fast do the singular values of R decay?

    Three regimes:
      uniformly small  -> accept V as is, no extension needed at all
      uniformly large  -> rebuild from scratch, nothing gained
      decaying         -> truncate at r_eff, block machinery pays off

    What this program does
    ----------------------
      1. read gauge config 1, build A1 = MdagM (Hermitian positive)
      2. compute the Nk lowest eigenvectors V of A1 with IRL
      3. read gauge config 2, build A2 the same way
      4. form W = A2 V, H = V^dag A2 V, R = W - V H
      5. singular values of R via the Nk x Nk Gram matrix R^dag R
         (eigenvalues are sigma^2; far cheaper than an SVD of field vectors)
      6. report the normalised spectrum sigma_i / ||A2 V|| and r_eff(eps)

    Also reported, as the sanity baseline, is the same measurement with
    config 2 = config 1, where R must vanish to eigensolver accuracy.

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

struct ReffParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(ReffParameters,
    std::string, config1,
    std::string, config2,
    RealD,       mass,
    Integer,     Nk,          // wanted eigenvectors = size of V
    Integer,     Nm,          // Krylov basis size, Nm > Nk
    Integer,     Nstop,       // how many must converge; V is truncated to this
    RealD,       lanTol,
    Integer,     lanMaxIter,
    Integer,     Ls,          // NON-HERMITIAN MODE ONLY (Mobius DWF)
    RealD,       M5,          // NON-HERMITIAN MODE ONLY
    RealD,       mobiusB,     // NON-HERMITIAN MODE ONLY
    RealD,       mobiusC,     // NON-HERMITIAN MODE ONLY
    RealD,       chebLo,      // HERMITIAN MODE ONLY: Chebyshev window
    RealD,       chebHi,
    Integer,     chebOrd,
    Integer,     selfTest,    // 1 = also run config2 == config1 baseline
    Integer,     driftScan,   // 1 = ignore config2, scan synthetic HMC-like drift
    RealD,       driftMin,    // smallest eps in U -> exp(i eps T) U
    RealD,       driftMax,
    Integer,     driftPts)

  ReffParameters()
    : config1("config1"), config2("config2"), mass(0.01),
      Nk(16), Nm(48), Nstop(16), lanTol(1e-8), lanMaxIter(10000),
      Ls(12), M5(1.8), mobiusB(1.5), mobiusC(0.5),
      chebLo(0.01), chebHi(64.0), chebOrd(40), selfTest(1),
      driftScan(0), driftMin(1e-3), driftMax(1e-1), driftPts(5)
  {}
};

} // namespace Grid

/// ComplexD is thrust::complex in CUDA builds; Eigen uses std::complex.
static inline ComplexD toGridCmplx(const std::complex<double>& c)
{ return ComplexD(c.real(), c.imag()); }

typedef WilsonFermionD                        WilsonOp;
typedef typename WilsonFermionD::FermionField FermionField;

/**
 * PV^dag M, lifted from Example_spec_kryschur.cc.  That file defines the
 * operator inline; it is not in the library, which is why the references to
 * PVdagMLinearOperator in Example_krylov_schur.cc are commented out.
 */
template<class Matrix, class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix& Mat, Matrix& PV) : _Mat(Mat), _PV(PV) {}
  void Op(const Field& in, Field& out) {
    Field tmp(in.Grid());
    _Mat.M(in, tmp);
    _PV.Mdag(tmp, out);
  }
  void AdjOp(const Field& in, Field& out) {
    Field tmp(in.Grid());
    _PV.M(in, tmp);
    _Mat.Mdag(tmp, out);
  }
  void HermOp(const Field& in, Field& out) {
    Field tmp(in.Grid());
    Op(in, tmp);
    AdjOp(tmp, out);
  }
  void OpDiag(const Field&, Field&)                        { assert(0); }
  void OpDir (const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)         { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&) { assert(0); }
};

// ===== MODE DIFFERENCE 1 and 2: action, operator, and how it is applied ====
//
// NB the modes differ in FIELD TYPE and GRID as well:
//   Hermitian      4D Wilson MdagM, fields live on UGrid
//   non-Hermitian  5D Mobius DWF PV^dag M, fields live on FGrid = Ls x UGrid
// Downstream code is written against FermType and a workGrid pointer so the
// measurement itself is untouched.
#ifdef REFF_NONHERMITIAN
typedef MobiusFermionD                          ActionD;
typedef LatticeFermionD                         FermType;
typedef PVdagMLinearOperator<ActionD, FermType> ReffOp;
static const char* kModeName =
  "NON-HERMITIAN (Mobius DWF, PV^dag M, KrylovSchur) - follows Example_spec_kryschur";
#else
typedef WilsonOp                                ActionD;
typedef FermionField                            FermType;
typedef MdagMLinearOperator<ActionD, FermType>  ReffOp;
static const char* kModeName =
  "HERMITIAN (Wilson MdagM, Implicitly Restarted Lanczos)";
#endif

/// Applying A to a field; mode independent from here on.
typedef std::function<void(const FermType&, FermType&)> ApplyOp;

static inline ApplyOp makeApply(ReffOp& A)
{
#ifdef REFF_NONHERMITIAN
  return [&A](const FermType& in, FermType& out) { A.Op(in, out); };
#else
  return [&A](const FermType& in, FermType& out) { A.HermOp(in, out); };
#endif
}

/**
 * Owns the fermion action(s) and the operator for one gauge field, so the
 * call sites are identical in both modes.  The actions must outlive the
 * operator, hence holding them by value here.
 */
struct OpHolder {
#ifdef REFF_NONHERMITIAN
  ActionD Ddwf, Dpv;
  ReffOp  A;
  OpHolder(LatticeGaugeField& U,
           GridCartesian* UG, GridRedBlackCartesian* Urb,
           GridCartesian* FG, GridRedBlackCartesian* Frb,
           const ReffParameters& P)
    : Ddwf(U, *FG, *Frb, *UG, *Urb, P.mass, P.M5, P.mobiusB, P.mobiusC),
      Dpv (U, *FG, *Frb, *UG, *Urb, 1.0,    P.M5, P.mobiusB, P.mobiusC),
      A(Ddwf, Dpv) {}
#else
  std::vector<Complex>  boundary;
  WilsonOp::ImplParams  ip;
  ActionD               Dw;
  ReffOp                A;
  OpHolder(LatticeGaugeField& U,
           GridCartesian* UG, GridRedBlackCartesian* Urb,
           GridCartesian* FG, GridRedBlackCartesian* Frb,
           const ReffParameters& P)
    : boundary({1,1,1,-1}), ip(boundary),
      Dw(U, *UG, *Urb, P.mass, ip),
      A(Dw) {}
#endif
};

/**
 * Low Nk modes of the chosen operator.  MODE DIFFERENCES 3 and 4 live
 * entirely inside this function.
 */
static void lowModes(LatticeGaugeField& Umu,
                     GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid,
                     GridCartesian* FGrid, GridRedBlackCartesian* FrbGrid,
                     GridParallelRNG& RNG,   // must match workGrid
                     const ReffParameters& P,
                     std::vector<FermType>& evec,
                     std::vector<std::complex<double> >& eval)
{
  OpHolder H1(Umu, UGrid, UrbGrid, FGrid, FrbGrid, P);
  ReffOp& A = H1.A;
#ifdef REFF_NONHERMITIAN
  GridCartesian* workGrid = FGrid;
#else
  GridCartesian* workGrid = UGrid;
#endif

  FermType src(workGrid);
  gaussian(RNG, src);

  // NB: std::vector<Lattice>::resize needs a default ctor, which Lattice
  // does not have.  Use the fill constructor and erase() to shrink.

#ifdef REFF_NONHERMITIAN
  // ---- MODE DIFFERENCE 3: KrylovSchur; RitzFilter, no Chebyshev.
  //      chebLo / chebHi / chebOrd are ignored in this mode.
  KrylovSchur<FermType> KS(A, workGrid, P.lanTol, EvalNormSmall);
  KS(src, P.lanMaxIter, P.Nm, P.Nk, P.Nstop);

  // Keep only the Nstop vectors that were required to converge; including
  // unconverged ones pollutes the baseline, where R must vanish.
  std::vector<FermType> ev = KS.evecs;
  int keep = std::min((int)P.Nstop, (int)ev.size());
  ev.erase(ev.begin() + keep, ev.end());
  evec.swap(ev);
  std::cout << GridLogMessage << "  KrylovSchur returned " << evec.size()
            << " vectors (Nstop=" << P.Nstop << ", Nk=" << P.Nk << ")" << std::endl;

  // ---- MODE DIFFERENCE 4: KrylovSchur::evals is private, and the
  //      eigenvalues are complex; recover them as Rayleigh quotients.
  {
    ApplyOp apply = makeApply(A);
    FermType w(workGrid);
    eval.clear();
    for (size_t i = 0; i < evec.size(); i++) {
      RealD n = std::sqrt(norm2(evec[i]));
      FermType v = evec[i] * (n > 1e-30 ? 1.0/n : 1.0);
      apply(v, w);
      eval.push_back(toStdCmplx(innerProduct(v, w)));
    }
  }
#else
  // ---- MODE DIFFERENCE 3: IRL, needs a Chebyshev accelerator and the
  //      PlainHermOp / FunctionHermOp wrappers.
  PlainHermOp<FermType>    hermop_wrap(A);
  Chebyshev<FermType>      Cheb(P.chebLo, P.chebHi, P.chebOrd + 1);
  FunctionHermOp<FermType> Cheb_wrap(Cheb, A);

  ImplicitlyRestartedLanczos<FermType> IRL(Cheb_wrap, hermop_wrap,
                                               P.Nstop, P.Nk, P.Nm,
                                               P.lanTol, P.lanMaxIter);
  std::vector<RealD> rval(P.Nm);
  std::vector<FermType> ev(P.Nm, FermType(workGrid));

  int Nconv = 0;
  IRL.calc(rval, ev, src, Nconv);
  std::cout << GridLogMessage << "  IRL converged " << Nconv
            << " of " << P.Nstop << " wanted (Nk=" << P.Nk << ")" << std::endl;

  int keep = std::min((int)P.Nstop, std::max(Nconv, 0));
  ev.erase(ev.begin() + keep, ev.end());
  evec.swap(ev);

  // ---- MODE DIFFERENCE 4: IRL returns real eigenvalues directly.
  eval.clear();
  for (int i = 0; i < keep; i++) eval.push_back(std::complex<double>(rval[i], 0.0));
#endif

  for (size_t i = 0; i < eval.size(); i++)
    std::cout << GridLogMessage << "    eval[" << i << "] = " << eval[i] << std::endl;
}

/**
 * HMC-like link update:  U2_mu = exp(i eps T_mu) U1_mu  with random
 * Lie-algebra T.  This is exactly the integrator's update_field, so eps
 * plays the role of an MD step size and lets us scan the drift that the
 * subspace-reuse argument depends on.
 */
static void driftConfig(const LatticeGaugeField& U1, LatticeGaugeField& U2,
                        GridParallelRNG& RNG, RealD eps)
{
  LatticeGaugeField P(U1.Grid());
  LatticeColourMatrix Pmu(U1.Grid());
  for (int mu = 0; mu < Nd; mu++) {
    SU<Nc>::GaussianFundamentalLieAlgebraMatrix(RNG, Pmu);
    PokeIndex<LorentzIndex>(P, Pmu, mu);
  }
  U2 = U1;
  PeriodicGimplD::update_field(P, U2, eps);
}

/**
 * The measurement.  V is assumed L2-orthonormal.
 */
static void measureReff(std::vector<FermType>& V,
                        ApplyOp apply,
                        GridCartesian* workGrid,
                        const std::string& label)
{
  int k = V.size();
  if (k == 0) { std::cout << GridLogMessage << label << ": empty V" << std::endl; return; }

  // Build R column by column, converting A2 V into R in place.  Holding V,
  // W and R at once would be 3k fields; at 32^4 a fermion field is ~201 MB,
  // so this matters.  Peak here is V (k) + R (k) + 1 scratch.
  std::vector<FermType> R(k, FermType(workGrid));
  Eigen::MatrixXcd H(k, k);
  RealD normAV2 = 0.0, normR2 = 0.0;
  for (int j = 0; j < k; j++) {
    apply(V[j], R[j]);                       // R[j] = (A2 V)_j for now
    normAV2 += norm2(R[j]);
    for (int i = 0; i < k; i++)                  // column j of H
      H(i, j) = toStdCmplx(innerProduct(V[i], R[j]));
    for (int i = 0; i < k; i++)                  // subtract in place
      R[j] = R[j] - V[i] * toGridCmplx(H(i, j));
    normR2 += norm2(R[j]);
  }

  // singular values of R from the Gram matrix R^dag R (eigenvalues = sigma^2)
  Eigen::MatrixXcd G(k, k);
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      G(i, j) = toStdCmplx(innerProduct(R[i], R[j]));

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(G);
  Eigen::VectorXd  ev = es.eigenvalues();     // ascending

  RealD scale = std::sqrt(normAV2);
  std::cout << GridLogMessage << "=== " << label << " ===" << std::endl;
  std::cout << GridLogMessage << "  k = " << k
            << "   ||A2 V||_F = " << scale
            << "   ||R||_F / ||A2 V||_F = " << std::sqrt(normR2)/scale << std::endl;
  std::cout << GridLogMessage << "  normalised singular values of R, descending:" << std::endl;
  for (int i = k - 1; i >= 0; i--) {
    RealD s = std::sqrt(std::max(0.0, ev(i))) / scale;
    std::cout << GridLogMessage << "    sigma[" << (k-1-i) << "] = " << s << std::endl;
  }
  for (RealD eps : {1e-1, 1e-2, 1e-3, 1e-4, 1e-6}) {
    int cnt = 0;
    for (int i = 0; i < k; i++)
      if (std::sqrt(std::max(0.0, ev(i)))/scale > eps) cnt++;
    std::cout << GridLogMessage << "  r_eff(eps=" << eps << ") = " << cnt
              << " / " << k << std::endl;
  }
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  ReffParameters P;
  { XmlReader rd("ReffParams.xml"); read(rd, "ReffParameters", P); }
  std::cout << GridLogMessage << P << std::endl;
  std::cout << GridLogMessage << "MODE: " << kModeName << std::endl;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(P.Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(P.Ls, UGrid);

  std::vector<int> seeds4({1,2,3,4}), seeds5({5,6,7,8});
  GridParallelRNG RNG(UGrid);  RNG.SeedFixedIntegers(seeds4);   // gauge drift
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);

  // the RNG whose grid matches the fermion fields in this mode
#ifdef REFF_NONHERMITIAN
  GridParallelRNG& RNGwork = RNG5;
  GridCartesian*   workGrid = FGrid;
#else
  GridParallelRNG& RNGwork = RNG;
  GridCartesian*   workGrid = UGrid;
#endif

  FieldMetaData h1, h2;
  LatticeGaugeField U1(UGrid), U2(UGrid);
  NerscIO::readConfiguration(U1, h1, P.config1);
  NerscIO::readConfiguration(U2, h2, P.config2);
  std::cout << GridLogMessage << "config1 = " << P.config1
            << "  plaq = " << h1.plaquette << std::endl;
  std::cout << GridLogMessage << "config2 = " << P.config2
            << "  plaq = " << h2.plaquette << std::endl;

  // gauge-field distance, the natural x-axis for this measurement
  {
    LatticeGaugeField dU(UGrid);
    dU = U2 - U1;
    std::cout << GridLogMessage << "  ||U2 - U1||_F / ||U1||_F = "
              << std::sqrt(norm2(dU)/norm2(U1)) << std::endl;
  }

  // --- V = low modes of config 1 ---
  std::cout << GridLogMessage << "--- low modes on config1 ---" << std::endl;
  std::vector<FermType> V;
  std::vector<std::complex<double> > eval1;
  lowModes(U1, UGrid, UrbGrid, FGrid, FrbGrid, RNGwork, P, V, eval1);

  // IRL output should already be orthonormal; enforce it so H and R are clean
  for (size_t j = 0; j < V.size(); j++) {
    for (size_t i = 0; i < j; i++)
      V[j] = V[j] - V[i] * toGridCmplx(toStdCmplx(innerProduct(V[i], V[j])));
    RealD n = std::sqrt(norm2(V[j]));
    if (n > 1e-12) V[j] = V[j] * (1.0/n);
  }

  // --- baseline: same config, R must vanish.  Check this first: if it is
  //     not tiny, V is not converged and nothing below means anything.
  if (P.selfTest) {
    OpHolder H1(U1, UGrid, UrbGrid, FGrid, FrbGrid, P);
    measureReff(V, makeApply(H1.A), workGrid,
                "baseline: config2 == config1 (R should vanish)");
  }

  // --- the measurement ---
  if (!P.driftScan) {
    OpHolder H2(U2, UGrid, UrbGrid, FGrid, FrbGrid, P);
    measureReff(V, makeApply(H2.A), workGrid,
                "config1 low modes vs config2 operator");
  } else {
    // scan synthetic drift: this is the regime the HMC reuse argument needs,
    // which two configs ten trajectories apart do not reach.
    for (int i = 0; i < (int)P.driftPts; i++) {
      RealD f = (P.driftPts == 1) ? 0.0 : (double)i/(double)(P.driftPts - 1);
      RealD eps = P.driftMin * std::pow(P.driftMax/P.driftMin, f);
      LatticeGaugeField Ud(UGrid);
      driftConfig(U1, Ud, RNG, eps);
      LatticeGaugeField dU(UGrid); dU = Ud - U1;
      RealD dnorm = std::sqrt(norm2(dU)/norm2(U1));
      std::cout << GridLogMessage << "drift eps = " << eps
                << "   ||dU||/||U|| = " << dnorm << std::endl;
      OpHolder Hd(Ud, UGrid, UrbGrid, FGrid, FrbGrid, P);
      measureReff(V, makeApply(Hd.A), workGrid,
                  "drift eps=" + std::to_string(eps)
                  + "  ||dU||/||U||=" + std::to_string(dnorm));
    }
  }

  std::cout << GridLogMessage << std::endl << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
