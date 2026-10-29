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
    Integer,     Nm,          // Lanczos basis size, Nm > Nk
    RealD,       lanTol,
    Integer,     lanMaxIter,
    RealD,       chebLo,      // Chebyshev window for IRL
    RealD,       chebHi,
    Integer,     chebOrd,
    Integer,     selfTest,    // 1 = also run config2 == config1 baseline
    Integer,     driftScan,   // 1 = ignore config2, scan synthetic HMC-like drift
    RealD,       driftMin,    // smallest eps in U -> exp(i eps T) U
    RealD,       driftMax,
    Integer,     driftPts)

  ReffParameters()
    : config1("config1"), config2("config2"), mass(0.01),
      Nk(16), Nm(48), lanTol(1e-8), lanMaxIter(10000),
      chebLo(0.01), chebHi(64.0), chebOrd(40), selfTest(1),
      driftScan(0), driftMin(1e-3), driftMax(1e-1), driftPts(5)
  {}
};

} // namespace Grid

/// ComplexD is thrust::complex in CUDA builds; Eigen uses std::complex.
static inline ComplexD toGridCmplx(const std::complex<double>& c)
{ return ComplexD(c.real(), c.imag()); }

typedef WilsonFermionD                      WilsonOp;
typedef typename WilsonFermionD::FermionField FermionField;

/// Applying A to a field; mode independent from here on.
typedef std::function<void(const FermionField&, FermionField&)> ApplyOp;

// ===== MODE DIFFERENCE 1 and 2: operator wrapper, and how it is applied =====
#ifdef REFF_NONHERMITIAN
typedef NonHermitianLinearOperator<WilsonOp, FermionField> ReffOp;
static inline ApplyOp makeApply(ReffOp& A)
{ return [&A](const FermionField& in, FermionField& out) { A.Op(in, out); }; }
static const char* kModeName = "NON-HERMITIAN (bare Dirac, KrylovSchur)";
#else
typedef MdagMLinearOperator<WilsonOp, FermionField>        ReffOp;
static inline ApplyOp makeApply(ReffOp& A)
{ return [&A](const FermionField& in, FermionField& out) { A.HermOp(in, out); }; }
static const char* kModeName = "HERMITIAN (MdagM, Implicitly Restarted Lanczos)";
#endif

/**
 * Low Nk modes of the chosen operator.  MODE DIFFERENCES 3 and 4 live
 * entirely inside this function.
 */
static void lowModes(LatticeGaugeField& Umu,
                     GridCartesian* UGrid, GridRedBlackCartesian* UrbGrid,
                     GridParallelRNG& RNG,
                     const ReffParameters& P,
                     std::vector<FermionField>& evec,
                     std::vector<std::complex<double> >& eval)
{
  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonOp::ImplParams ip(boundary);
  WilsonOp Dw(Umu, *UGrid, *UrbGrid, P.mass, ip);
  ReffOp A(Dw);

  FermionField src(UGrid);
  gaussian(RNG, src);

  // NB: std::vector<Lattice>::resize needs a default ctor, which Lattice
  // does not have.  Use the fill constructor and erase() to shrink.

#ifdef REFF_NONHERMITIAN
  // ---- MODE DIFFERENCE 3: KrylovSchur; RitzFilter, no Chebyshev.
  //      chebLo / chebHi / chebOrd are ignored in this mode.
  KrylovSchur<FermionField> KS(A, UGrid, P.lanTol, EvalNormSmall);
  KS(src, P.lanMaxIter, P.Nm, P.Nk, P.Nk);

  std::vector<FermionField> ev = KS.evecs;
  int keep = std::min((int)P.Nk, (int)ev.size());
  ev.erase(ev.begin() + keep, ev.end());
  evec.swap(ev);
  std::cout << GridLogMessage << "  KrylovSchur returned " << evec.size()
            << " vectors (wanted " << P.Nk << ")" << std::endl;

  // ---- MODE DIFFERENCE 4: KrylovSchur::evals is private, and the
  //      eigenvalues are complex; recover them as Rayleigh quotients.
  {
    ApplyOp apply = makeApply(A);
    FermionField w(UGrid);
    eval.clear();
    for (size_t i = 0; i < evec.size(); i++) {
      RealD n = std::sqrt(norm2(evec[i]));
      FermionField v = evec[i] * (n > 1e-30 ? 1.0/n : 1.0);
      apply(v, w);
      eval.push_back(toStdCmplx(innerProduct(v, w)));
    }
  }
#else
  // ---- MODE DIFFERENCE 3: IRL, needs a Chebyshev accelerator and the
  //      PlainHermOp / FunctionHermOp wrappers.
  PlainHermOp<FermionField>    hermop_wrap(A);
  Chebyshev<FermionField>      Cheb(P.chebLo, P.chebHi, P.chebOrd + 1);
  FunctionHermOp<FermionField> Cheb_wrap(Cheb, A);

  ImplicitlyRestartedLanczos<FermionField> IRL(Cheb_wrap, hermop_wrap,
                                               P.Nk, P.Nk, P.Nm,
                                               P.lanTol, P.lanMaxIter);
  std::vector<RealD> rval(P.Nm);
  std::vector<FermionField> ev(P.Nm, FermionField(UGrid));

  int Nconv = 0;
  IRL.calc(rval, ev, src, Nconv);
  std::cout << GridLogMessage << "  IRL converged " << Nconv
            << " of " << P.Nk << " wanted" << std::endl;

  int keep = std::min((int)P.Nk, std::max(Nconv, 0));
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
static void measureReff(std::vector<FermionField>& V,
                        ApplyOp apply,
                        GridCartesian* UGrid,
                        const std::string& label)
{
  int k = V.size();
  if (k == 0) { std::cout << GridLogMessage << label << ": empty V" << std::endl; return; }

  // Build R column by column, converting A2 V into R in place.  Holding V,
  // W and R at once would be 3k fields; at 32^4 a fermion field is ~201 MB,
  // so this matters.  Peak here is V (k) + R (k) + 1 scratch.
  std::vector<FermionField> R(k, FermionField(UGrid));
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

  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG RNG(UGrid); RNG.SeedFixedIntegers(seeds);

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
  std::vector<FermionField> V;
  std::vector<std::complex<double> > eval1;
  lowModes(U1, UGrid, UrbGrid, RNG, P, V, eval1);

  // IRL output should already be orthonormal; enforce it so H and R are clean
  for (size_t j = 0; j < V.size(); j++) {
    for (size_t i = 0; i < j; i++)
      V[j] = V[j] - V[i] * toGridCmplx(toStdCmplx(innerProduct(V[i], V[j])));
    RealD n = std::sqrt(norm2(V[j]));
    if (n > 1e-12) V[j] = V[j] * (1.0/n);
  }

  std::vector<Complex> boundary = {1,1,1,-1};
  WilsonOp::ImplParams ip(boundary);

  // --- baseline: same config, R must vanish ---
  if (P.selfTest) {
    WilsonOp Dw1(U1, *UGrid, *UrbGrid, P.mass, ip);
    ReffOp A1(Dw1);
    measureReff(V, makeApply(A1), UGrid, "baseline: config2 == config1 (R should vanish)");
  }

  // --- the measurement ---
  if (!P.driftScan) {
    WilsonOp Dw2(U2, *UGrid, *UrbGrid, P.mass, ip);
    ReffOp A2(Dw2);
    measureReff(V, makeApply(A2), UGrid, "config1 low modes vs config2 operator");
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
      WilsonOp Dwd(Ud, *UGrid, *UrbGrid, P.mass, ip);
      ReffOp Ad(Dwd);
      measureReff(V, makeApply(Ad), UGrid,
                  "drift eps=" + std::to_string(eps)
                  + "  ||dU||/||U||=" + std::to_string(dnorm));
    }
  }

  std::cout << GridLogMessage << std::endl << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
