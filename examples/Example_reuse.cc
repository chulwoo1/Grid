/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_reuse.cc

    Does reusing config1's eigenvectors actually help on config2?

    Example_reff.cc measures whether the subspace is stale (the singular
    value spectrum of R).  This programme asks the operational question:
    given V from config1, is anything gained on config2?

    Four phases:

      A  V = Nstop low modes of config1                        (timed)
      B  Rayleigh-Ritz of V against config2, no new solve:
             H = V^dag A2 V,  Ritz pairs (lambda_i, V y_i),
             EXACT residual  ||A2 (V y) - lambda (V y)|| = ||R y||
         This is free apart from k applications of A2, and tells you how
         many modes are already good enough to deflate with.
      C  fresh eigensolve on config2 from a random start       (timed)
      D  seeded eigensolve on config2, starting vector built from V (timed)

    B answers "can I skip the solve entirely".  C versus D answers "if I do
    have to solve, does the old subspace make it cheaper".  Both matter:
    they are the two ways subspace reuse can pay off in HMC.

    Note on D: IRL and KrylovSchur both take a single starting vector, so
    the seed is the normalised sum of the columns of V.  That is a weaker
    form of reuse than restarting from the full subspace -- doing the
    latter properly needs the block relation A V = V H + F B and a block
    eigensolver, which this tree does not yet have (see r_eff.tex).

    Mode switch, identical to Example_reff.cc:
      default                 HERMITIAN      Mobius MdagM,   IRL
      -DREFF_NONHERMITIAN     NON-HERMITIAN  Mobius PV^dag M, KrylovSchur

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

struct ReuseParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(ReuseParameters,
    std::string, config1,
    std::string, config2,
    RealD,       mass,
    Integer,     Ls,
    RealD,       M5,
    RealD,       mobiusB,
    RealD,       mobiusC,
    Integer,     Nk,
    Integer,     Nm,
    Integer,     Nstop,
    RealD,       lanTol,
    Integer,     lanMaxIter,
    RealD,       chebLo,      // HERMITIAN MODE ONLY
    RealD,       chebHi,
    Integer,     chebOrd,
    std::string, checkpointPrefix,  // both modes; "" disables checkpointing
    Integer,     checkpointInterval,
    Integer,     resume,
    Integer,     compressMode,// 0 = off, 1 = grouped sum, 2 = SVD best-preserved
    Integer,     Nreduced,    // target block size after compression, in BOTH
                              // modes.  Ignored if 0 or >= k.
    Integer,     doFresh,     // run phase C
    Integer,     doSeeded)    // run phase D

  ReuseParameters()
    : config1("config1"), config2("config2"), mass(0.01),
      Ls(8), M5(1.8), mobiusB(1.5), mobiusC(0.5),
      Nk(6), Nm(24), Nstop(4), lanTol(1e-6), lanMaxIter(20000),
      chebLo(1.0), chebHi(70.0), chebOrd(24),
      checkpointPrefix(""), checkpointInterval(1), resume(0),
      compressMode(0), Nreduced(2), doFresh(1), doSeeded(1)
  {}
};

} // namespace Grid


template <class T> void writeFile(T& in, std::string const fname) {
#ifdef HAVE_LIME
  std::cout << Grid::GridLogMessage << "Writes to: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacWriter WR(in.Grid()->IsBoss());
  WR.open(fname);
  WR.writeScidacFieldRecord(in, record, 0);
  WR.close();
#endif
}

template <class T> void readFile(T& out, std::string const fname) {
  std::cout << Grid::GridLogMessage << "Reads from: " << fname << std::endl;
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out, record);
  RD.close();
}

/**
 * Checkpointing works in both modes, with the same parameters.
 *
 * KrylovSchur, BlockKrylovSchur, TrueHarmonicBlockKrylovSchur and (newly)
 * ImplicitlyRestartedLanczos all expose the same four controls --
 * checkpointPrefix, checkpointInterval, resumeFromCheckpoint and the
 * fieldWrite/fieldRead hooks -- so this helper is one mode-independent
 * template rather than two branches.
 *
 * What is checkpointed differs, because the solvers' resumable states differ:
 * KrylovSchur stores the basis, the residual vector u, the Rayleigh quotient
 * and b; IRL stores evec[0:Nk] and the real tridiagonal alpha/beta.  Both
 * resume at a restart boundary, and both can change Nm and Nstop across a
 * resume but not Nk.
 *
 * Note this only checkpoints the eigensolve itself.  The phases that follow it
 * -- the config2 Rayleigh-Ritz, the compression, the fresh and seeded solves --
 * are not checkpointed, so a resume re-enters the phase A solve and then runs
 * the rest from there.
 */
template<class Solver>
static void setupCheckpoint(Solver& S, const ReuseParameters& P,
                            const std::string& phase)
{
  if (P.checkpointPrefix.empty()) return;
  // Each phase gets its own slot pair.  lowModes() runs up to three separate
  // eigensolves (A on config1, C fresh on config2, D seeded on config2) and
  // they must not share a prefix, or phase C would resume from phase A's
  // state on a different gauge field.
  S.checkpointPrefix     = P.checkpointPrefix + "." + phase;
  S.checkpointInterval   = (int)P.checkpointInterval;
  S.resumeFromCheckpoint = (P.resume != 0);
  // Generic lambdas so this helper does not depend on FermType being
  // declared yet; they convert to the solver's std::function hooks.
  S.fieldWrite = [](auto& f, const std::string& fn){ writeFile(f, fn); };
  S.fieldRead  = [](auto& f, const std::string& fn){ readFile(f, fn); };
  std::cout << GridLogMessage << "checkpointing: prefix=" << S.checkpointPrefix
            << "  interval=" << P.checkpointInterval
            << "  resume=" << P.resume << std::endl;
}

static inline ComplexD toGridCmplx(const std::complex<double>& c)
{ return ComplexD(c.real(), c.imag()); }

/// PV^dag M, as in Example_spec_kryschur.cc (not in the library).
template<class Matrix, class Field>
class PVdagMLinearOperator : public LinearOperatorBase<Field> {
  Matrix &_Mat;
  Matrix &_PV;
public:
  PVdagMLinearOperator(Matrix& Mat, Matrix& PV) : _Mat(Mat), _PV(PV) {}
  void Op(const Field& in, Field& out) {
    Field tmp(in.Grid()); _Mat.M(in, tmp); _PV.Mdag(tmp, out);
  }
  void AdjOp(const Field& in, Field& out) {
    Field tmp(in.Grid()); _PV.M(in, tmp); _Mat.Mdag(tmp, out);
  }
  void HermOp(const Field& in, Field& out) {
    Field tmp(in.Grid()); Op(in, tmp); AdjOp(tmp, out);
  }
  void OpDiag(const Field&, Field&)                        { assert(0); }
  void OpDir (const Field&, Field&, int, int)              { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)         { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&) { assert(0); }
};

typedef MobiusFermionD  ActionD;
typedef LatticeFermionD FermType;

#ifdef REFF_NONHERMITIAN
typedef PVdagMLinearOperator<ActionD, FermType> ReuseOp;
static const char* kModeName = "NON-HERMITIAN (Mobius PV^dag M, KrylovSchur)";
#else
typedef MdagMLinearOperator<ActionD, FermType>  ReuseOp;
static const char* kModeName = "HERMITIAN (Mobius MdagM, IRL)";
#endif

typedef std::function<void(const FermType&, FermType&)> ApplyOp;

static inline ApplyOp makeApply(LinearOperatorBase<FermType>& A)
{
#ifdef REFF_NONHERMITIAN
  return [&A](const FermType& in, FermType& out) { A.Op(in, out); };
#else
  return [&A](const FermType& in, FermType& out) { A.HermOp(in, out); };
#endif
}

/// Owns the action(s) and the operator for one gauge field.
struct OpHolder {
  ActionD Ddwf;
#ifdef REFF_NONHERMITIAN
  ActionD Dpv;
#endif
  ReuseOp Abare;
  // Everything downstream uses the counted wrapper, so matvecs are attributed
  // to whichever phase is running.  costPerCall = 2: one call of either
  // PV^dag M or MdagM is two applications of the 5D Dirac operator.
  CountingLinearOperator<FermType> A;
  OpHolder(LatticeGaugeField& U,
           GridCartesian* UG, GridRedBlackCartesian* Urb,
           GridCartesian* FG, GridRedBlackCartesian* Frb,
           const ReuseParameters& P)
    : Ddwf(U, *FG, *Frb, *UG, *Urb, P.mass, P.M5, P.mobiusB, P.mobiusC)
#ifdef REFF_NONHERMITIAN
    , Dpv (U, *FG, *Frb, *UG, *Urb, 1.0,    P.M5, P.mobiusB, P.mobiusC)
    , Abare(Ddwf, Dpv)
#else
    , Abare(Ddwf)
#endif
    , A(Abare, 2)
  {}
};


// Chebyshev acceleration only exists in the Hermitian mode, so this check
// is compiled and called there only.
#ifndef REFF_NONHERMITIAN
/**
 * Estimate the upper spectral bound of A.HermOp by power iteration, using
 * the same idiom as the HDCG tests (PowerMethod PM; PM(HermOp, src)).
 *
 * The Chebyshev accelerator suppresses [chebLo, chebHi] and amplifies
 * outside it.  If chebHi < lambda_max the TOP of the spectrum is amplified
 * and the eigensolver silently returns the largest modes instead of the
 * smallest -- a failure that looks like success.  Stop here instead.
 */
static RealD checkSpectralBound(LinearOperatorBase<FermType>& A, GridCartesian* workGrid,
                                GridParallelRNG& RNG, const ReuseParameters& P)
{
  FermType src(workGrid);
  gaussian(RNG, src);
  PowerMethod<FermType> PM;
  RealD lmax = PM(A, src);

  std::cout << GridLogMessage << "spectral bound: lambda_max(HermOp) ~ "
            << lmax << "   Chebyshev window [" << P.chebLo << ", "
            << P.chebHi << "]" << std::endl;

  if (lmax > P.chebHi) {
    std::cout << GridLogError
              << "chebHi = " << P.chebHi << " is BELOW lambda_max = " << lmax
              << ".  The Chebyshev would amplify the top of the spectrum and"
              << " the solver would return the LARGEST modes, not the"
              << " smallest.  Raise chebHi above " << lmax << " and rerun."
              << std::endl;
    assert(0 && "chebHi below lambda_max");
  }
  if (P.chebLo >= lmax) {
    std::cout << GridLogError
              << "chebLo = " << P.chebLo << " is at or above lambda_max = "
              << lmax << "; the whole spectrum would be amplified." << std::endl;
    assert(0 && "chebLo above lambda_max");
  }
  return lmax;
}
#endif


/**
 * Low modes of A, optionally from a supplied starting vector.
 * Returns wall seconds.
 */
static double lowModes(LinearOperatorBase<FermType>& A, GridCartesian* workGrid,
                       GridParallelRNG& RNG, const ReuseParameters& P,
                       const FermType* seed, const std::string& phase,
                       std::vector<FermType>& evec,
                       std::vector<std::complex<double> >& eval)
{
  FermType src(workGrid);
  if (seed) src = *seed;
  else      gaussian(RNG, src);

  double t0 = usecond();

#ifdef REFF_NONHERMITIAN
  KrylovSchur<FermType> KS(A, workGrid, P.lanTol, EvalNormSmall);
  setupCheckpoint(KS, P, phase);
  KS(src, P.lanMaxIter, P.Nm, P.Nk, P.Nstop);
  std::vector<FermType> ev = KS.evecs;
  int keep = std::min((int)P.Nstop, (int)ev.size());
  ev.erase(ev.begin() + keep, ev.end());
  evec.swap(ev);
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
  PlainHermOp<FermType>    hermop_wrap(A);
  Chebyshev<FermType>      Cheb(P.chebLo, P.chebHi, P.chebOrd + 1);
  FunctionHermOp<FermType> Cheb_wrap(Cheb, A);
  ImplicitlyRestartedLanczos<FermType> IRL(Cheb_wrap, hermop_wrap,
                                           P.Nstop, P.Nk, P.Nm,
                                           P.lanTol, P.lanMaxIter);
  setupCheckpoint(IRL, P, phase);
  std::vector<RealD> rval(P.Nm);
  std::vector<FermType> ev(P.Nm, FermType(workGrid));
  int Nconv = 0;
  IRL.calc(rval, ev, src, Nconv);
  int keep = std::min((int)P.Nstop, std::max(Nconv, 0));
  ev.erase(ev.begin() + keep, ev.end());
  evec.swap(ev);
  eval.clear();
  for (int i = 0; i < keep; i++) eval.push_back(std::complex<double>(rval[i], 0.0));
#endif

  double t1 = usecond();
  return (t1 - t0) / 1e6;
}

/// L2 modified Gram-Schmidt in place.
static int orthonormalise(std::vector<FermType>& V)
{
  int keep = 0;
  for (size_t j = 0; j < V.size(); j++) {
    for (int i = 0; i < keep; i++)
      V[j] = V[j] - V[i] * toGridCmplx(toStdCmplx(innerProduct(V[i], V[j])));
    RealD n = std::sqrt(norm2(V[j]));
    if (n > 1e-10) { V[j] = V[j] * (1.0/n); std::swap(V[keep], V[j]); keep++; }
  }
  return keep;
}

/**
 * Compress V to exactly Nreduced vectors by summing consecutive groups.
 * The k inputs are split into Nreduced balanced contiguous groups (sizes
 * differ by at most one), so Nreduced is the output size, not a factor --
 * the same meaning it has in mode 2.
 *
 * Because the groups use disjoint sets of orthonormal vectors, the sums are
 * automatically mutually orthogonal; only normalisation is needed.
 *
 * It is lossy by construction: span is reduced from k to m dimensions, and
 * summing modes with different eigenvalues produces a vector that is not
 * close to an eigenvector of anything.  The eigenvalue spread within each
 * group is reported as a guide to how lossy a given grouping is -- a group
 * spanning a decade of eigenvalue is doing real damage, a group of
 * near-degenerate modes much less.
 *
 * A better but more expensive alternative, not implemented: form R on the
 * new configuration and keep the right singular vectors of R with the
 * SMALLEST singular values.  Those are the directions of span(V) best
 * preserved by the new operator, so the compression is chosen to minimise
 * exactly the residual the block relation has to carry.
 */
static void compressSubspace(const std::vector<FermType>& V,
                             const std::vector<std::complex<double> >& eval,
                             int Nreduced, GridCartesian* grid,
                             std::vector<FermType>& Vc)
{
  int k = V.size();
  int m = std::min(Nreduced, k);
  Vc.clear();
  Vc.reserve(m);
  std::cout << GridLogMessage << "  compressing " << k << " vectors -> " << m
            << " by grouped sum (balanced groups)" << std::endl;

  for (int g = 0; g < m; g++) {
    // balanced contiguous partition into exactly m groups
    int i0 = (int)(( (long)g      * k) / m);
    int i1 = (int)(( (long)(g+1)  * k) / m);
    FermType w(grid); w = Zero();
    double lo = 1e300, hi = 0.0;
    for (int i = i0; i < i1; i++) {
      w = w + V[i];
      if (i < (int)eval.size()) {
        double a = std::abs(eval[i]);
        lo = std::min(lo, a); hi = std::max(hi, a);
      }
    }
    RealD n = std::sqrt(norm2(w));
    if (n > 1e-30) w = w * (1.0/n);
    Vc.push_back(w);
    std::cout << GridLogMessage << "    group " << g << "  modes [" << i0
              << "," << (i1-1) << "]   |lambda| range " << lo << " .. " << hi
              << "   spread " << (lo > 0 ? hi/lo : -1.0) << std::endl;
  }
}

/**
 * Compression by best-preserved directions.
 *
 * Form R = (I - V V^dag) A2 V on the NEW configuration and keep the right
 * singular vectors of R with the SMALLEST singular values.  Those are the
 * combinations of V that the new operator moves least out of span(V), i.e.
 * precisely the directions whose residual a block restart would have to
 * carry.  Compressing this way minimises that residual, rather than
 * discarding span by index as grouped summation does.
 *
 * Uses the new gauge field, which is legitimate here: the k applications of
 * A2 are the same ones phase B needs anyway.
 *
 * Singular values come from the k x k Gram matrix R^dag R, whose
 * eigenvalues are sigma^2 -- far cheaper than an SVD of field vectors.
 * Eigen returns them ascending, so the leading m columns are the ones we
 * want.
 */
static void compressSVD(const std::vector<FermType>& V, ApplyOp applyA2,
                        int m, GridCartesian* grid,
                        std::vector<FermType>& Vc)
{
  int k = V.size();
  m = std::min(m, k);
  std::cout << GridLogMessage << "  compressing " << k << " vectors -> " << m
            << " by smallest-residual directions" << std::endl;

  // R = (I - P) A2 V, built in place
  std::vector<FermType> R(k, FermType(grid));
  Eigen::MatrixXcd H(k, k);
  RealD normAV2 = 0.0;
  for (int j = 0; j < k; j++) {
    applyA2(V[j], R[j]);
    normAV2 += norm2(R[j]);
    for (int i = 0; i < k; i++) H(i, j) = toStdCmplx(innerProduct(V[i], R[j]));
    for (int i = 0; i < k; i++) R[j] = R[j] - V[i] * toGridCmplx(H(i, j));
  }

  Eigen::MatrixXcd G(k, k);
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      G(i, j) = toStdCmplx(innerProduct(R[i], R[j]));

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(G);
  Eigen::VectorXd  ev = es.eigenvalues();      // ascending = smallest sigma first
  Eigen::MatrixXcd Q  = es.eigenvectors();

  RealD scale = std::sqrt(normAV2);
  for (int j = 0; j < k; j++) {
    RealD sig = std::sqrt(std::max(0.0, ev(j))) / scale;
    std::cout << GridLogMessage << "    sigma[" << j << "] = " << sig
              << (j < m ? "   KEPT" : "   dropped") << std::endl;
  }

  // Vc_j = sum_i V_i Q(i,j), j = 0..m-1.  Orthonormal because V is and Q is.
  Vc.clear();
  Vc.reserve(m);
  for (int j = 0; j < m; j++) {
    FermType w(grid); w = Zero();
    for (int i = 0; i < k; i++) w = w + V[i] * toGridCmplx(Q(i, j));
    RealD n = std::sqrt(norm2(w));
    if (n > 1e-30) w = w * (1.0/n);
    Vc.push_back(w);
  }
}

/**
 * Phase B.  Rayleigh-Ritz of V against the new operator, with the EXACT
 * residual rather than an estimate.
 *
 * For a Ritz pair x = V y with H y = lambda y,
 *     A x - lambda x = (A V - V H) y = R y,
 * so the residual is exact once R is formed.  No eigensolver is run.
 */
static void rayleighRitz(std::vector<FermType>& V, ApplyOp apply,
                         GridCartesian* workGrid, RealD tol)
{
  int k = V.size();
  if (k == 0) { std::cout << GridLogMessage << "  empty V" << std::endl; return; }

  std::vector<FermType> R(k, FermType(workGrid));
  Eigen::MatrixXcd H(k, k);
  RealD normAV2 = 0.0;
  for (int j = 0; j < k; j++) {
    apply(V[j], R[j]);
    normAV2 += norm2(R[j]);
    for (int i = 0; i < k; i++) H(i, j) = toStdCmplx(innerProduct(V[i], R[j]));
    for (int i = 0; i < k; i++) R[j] = R[j] - V[i] * toGridCmplx(H(i, j));
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(H);
  Eigen::VectorXcd lam = es.eigenvalues();
  Eigen::MatrixXcd Y   = es.eigenvectors();

  std::cout << GridLogMessage << "  Rayleigh-Ritz of config1 modes on config2"
            << "   ||A2 V||_F = " << std::sqrt(normAV2) << std::endl;

  int ngood = 0;
  for (int m = 0; m < k; m++) {
    Eigen::VectorXcd y = Y.col(m);
    y /= y.norm();
    // exact residual ||R y||
    FermType r(workGrid); r = Zero();
    for (int j = 0; j < k; j++) r = r + R[j] * toGridCmplx(y(j));
    RealD res = std::sqrt(norm2(r));
    if (res < tol) ngood++;
    std::cout << GridLogMessage << "    ritz[" << m << "]  lambda = "
              << lam(m) << "   exact |res| = " << res
              << (res < tol ? "   *" : "") << std::endl;
  }
  std::cout << GridLogMessage << "  modes already below tol=" << tol
            << ":  " << ngood << " / " << k << std::endl;
}

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  ReuseParameters P;
  { XmlReader rd("ReuseParams.xml"); read(rd, "ReuseParameters", P); }
  std::cout << GridLogMessage << P << std::endl;
  std::cout << GridLogMessage << "MODE: " << kModeName << std::endl;

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian*         FGrid   = SpaceTimeGrid::makeFiveDimGrid(P.Ls, UGrid);
  GridRedBlackCartesian* FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(P.Ls, UGrid);

  std::vector<int> seeds5({5,6,7,8});
  GridParallelRNG RNG5(FGrid); RNG5.SeedFixedIntegers(seeds5);

  FieldMetaData h1, h2;
  LatticeGaugeField U1(UGrid), U2(UGrid);
  NerscIO::readConfiguration(U1, h1, P.config1);
  NerscIO::readConfiguration(U2, h2, P.config2);
  {
    LatticeGaugeField dU(UGrid); dU = U2 - U1;
    std::cout << GridLogMessage << "config1 plaq = " << h1.plaquette
              << "   config2 plaq = " << h2.plaquette
              << "   ||U2-U1||/||U1|| = " << std::sqrt(norm2(dU)/norm2(U1))
              << std::endl;
  }

  OpHolder H1(U1, UGrid, UrbGrid, FGrid, FrbGrid, P);
  OpHolder H2(U2, UGrid, UrbGrid, FGrid, FrbGrid, P);

  // ---- Phase A: modes of config1 ----
  std::cout << GridLogMessage << "=== Phase A: eigensolve on config1 ===" << std::endl;
  std::vector<FermType> V;
  std::vector<std::complex<double> > ev1;
  H1.A.reset();
  double tA = lowModes(H1.A, FGrid, RNG5, P, nullptr, "phaseA", V, ev1);
  long mvA = H1.A.matvecs();
  int nv = orthonormalise(V);
  V.erase(V.begin() + nv, V.end());
  std::cout << GridLogMessage << "  config1: " << V.size()
            << " modes in " << tA << " s, " << mvA << " matvecs" << std::endl;
  for (size_t i = 0; i < ev1.size(); i++)
    std::cout << GridLogMessage << "    eval1[" << i << "] = " << ev1[i] << std::endl;

  // ---- optional compression of the reused subspace ----
  std::vector<FermType> Vuse;
  int Nred = (int)P.Nreduced;
  if (P.compressMode != 0 && (int)V.size() > 1 &&
      Nred > 0 && Nred < (int)V.size()) {
    std::cout << GridLogMessage << "=== compressing config1 subspace (mode "
              << P.compressMode << ", Nreduced = " << Nred << ") ===" << std::endl;
    if (P.compressMode == 1) {
      compressSubspace(V, ev1, Nred, FGrid, Vuse);
    } else {
      compressSVD(V, makeApply(H2.A), Nred, FGrid, Vuse);
    }
    int nc = orthonormalise(Vuse);
    Vuse.erase(Vuse.begin() + nc, Vuse.end());
    std::cout << GridLogMessage << "  compressed block size = " << Vuse.size()
              << std::endl;
  } else {
    Vuse = V;
    std::cout << GridLogMessage
              << "=== no compression (compressMode=" << P.compressMode
              << ", Nreduced=" << Nred << ", k=" << V.size() << ") ==="
              << std::endl;
  }

  // ---- Phase B: reuse without solving ----
  std::cout << GridLogMessage << "=== Phase B: reuse V on config2, no solve ==="
            << "   (block size " << Vuse.size() << ")" << std::endl;
  H2.A.reset();
  rayleighRitz(Vuse, makeApply(H2.A), FGrid, P.lanTol);
  long mvB = H2.A.matvecs();
  std::cout << GridLogMessage << "  phase B cost: " << mvB
            << " matvecs (exact residuals, no solve)" << std::endl;

  // ---- Phase C: fresh solve on config2 ----
  double tC = -1.0;
  long   mvC = -1;
  std::vector<std::complex<double> > evC;
  if (P.doFresh) {
    std::cout << GridLogMessage << "=== Phase C: fresh eigensolve on config2 ===" << std::endl;
    std::vector<FermType> Vc;
    H2.A.reset();
    tC = lowModes(H2.A, FGrid, RNG5, P, nullptr, "phaseC", Vc, evC);
    mvC = H2.A.matvecs();
    std::cout << GridLogMessage << "  config2 fresh: " << Vc.size()
              << " modes in " << tC << " s, " << mvC << " matvecs" << std::endl;
    for (size_t i = 0; i < evC.size(); i++)
      std::cout << GridLogMessage << "    evalC[" << i << "] = " << evC[i] << std::endl;
  }

  // ---- Phase D: seeded solve on config2 ----
  double tD = -1.0;
  long   mvD = -1;
  std::vector<std::complex<double> > evD;
  if (P.doSeeded && !Vuse.empty()) {
    std::cout << GridLogMessage << "=== Phase D: seeded eigensolve on config2 ===" << std::endl;
    FermType seed(FGrid); seed = Zero();
    for (size_t i = 0; i < Vuse.size(); i++) seed = seed + Vuse[i];
    RealD n = std::sqrt(norm2(seed));
    if (n > 1e-30) seed = seed * (1.0/n);
    std::vector<FermType> Vd;
    H2.A.reset();
    tD = lowModes(H2.A, FGrid, RNG5, P, &seed, "phaseD", Vd, evD);
    mvD = H2.A.matvecs();
    std::cout << GridLogMessage << "  config2 seeded: " << Vd.size()
              << " modes in " << tD << " s, " << mvD << " matvecs" << std::endl;
    for (size_t i = 0; i < evD.size(); i++)
      std::cout << GridLogMessage << "    evalD[" << i << "] = " << evD[i] << std::endl;
  }

  std::cout << GridLogMessage << "=== summary ===" << std::endl;
  std::cout << GridLogMessage
            << "  phase                     wall_s      matvecs" << std::endl;
  std::cout << GridLogMessage << "  A config1 solve       "
            << std::setw(10) << tA << std::setw(13) << mvA << std::endl;
  std::cout << GridLogMessage << "  B reuse test on cfg2  "
            << std::setw(10) << "-" << std::setw(13) << mvB << std::endl;
  std::cout << GridLogMessage << "  C config2 fresh       "
            << std::setw(10) << tC << std::setw(13) << mvC << std::endl;
  std::cout << GridLogMessage << "  D config2 seeded      "
            << std::setw(10) << tD << std::setw(13) << mvD << std::endl;
  // Matvecs are the figure to trust.  Wall time on a shared GPU has swung by
  // two orders of magnitude in this tree for an unchanged calculation, so a
  // wall-time speedup near 1 is not evidence either way on its own.
  if (mvC > 0 && mvD > 0)
    std::cout << GridLogMessage << "  seeding speedup, matvecs: "
              << (double)mvC/(double)mvD << " x" << std::endl;
  if (tC > 0 && tD > 0)
    std::cout << GridLogMessage << "  seeding speedup, wall   : " << tC/tD << " x" << std::endl;
  std::cout << GridLogMessage
            << "  (phase B is the reuse test itself: " << mvB
            << " matvecs against " << mvC << " to rebuild from scratch)" << std::endl;

  std::cout << GridLogMessage << std::endl << "Done" << std::endl;
  Grid_finalize();
  return 0;
}
