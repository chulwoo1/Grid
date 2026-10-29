/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Non-Hermitian FEAST contour-integration eigensolver.

    References:
      E. Polizzi, "Density-Matrix-Based Algorithms for Solving Eigenvalue
        Problems", Phys. Rev. B 79, 115112 (2009)
      J. Kestyn, E. Polizzi, P. T. P. Tang, "FEAST Eigensolver for
        Non-Hermitian Problems", SIAM J. Sci. Comput. 38(5), S772-S799 (2016)
      FEAST v4.0 User Guide, arXiv:2002.04807  (refs/2002.04807v1.pdf)

    FEAST is subspace iteration in which the usual multiplication by A is
    replaced by an approximate spectral projector built from a contour
    integral.  For a contour C enclosing the wanted eigenvalues,

        P = (1/2 pi i) \oint_C (z - A)^{-1} dz

    is the exact projector onto the invariant subspace of the enclosed
    eigenvalues.  Replacing the integral by an Ne-point quadrature gives the
    rational filter actually applied to the trial block Y:

        Q_R = sum_k w_k (z_k - A)^{-1} Y_R          (right subspace)
        Q_L = sum_k conj(w_k) (conj(z_k) - A^H)^{-1} Y_L   (left subspace)

    For a non-Hermitian A the spectrum is not symmetric about the real axis,
    so a FULL contour is required and BOTH subspaces must be built: the
    Rayleigh-Ritz step is oblique,

        A_q = Q_L^H A Q_R,   B_q = Q_L^H Q_R,   A_q u = lambda B_q u

    and the Ritz vectors are X_R = Q_R u_R, X_L = Q_L u_L.  Eigenvalues
    falling outside the contour are discarded; the surviving ones seed the
    next iteration.  Convergence is linear at a rate set by the decay of the
    rational filter outside C, so it is governed by Ne and by the oversampling
    M0 relative to the true count M inside C.

    The dominant cost is the Ne*M0 shifted solves per iteration.  These are
    delegated through the shiftSolve hook so the caller chooses the solver;
    useBiCGSTAB() installs a reasonable default.

*************************************************************************************/
#ifndef GRID_FEAST_NON_HERMITIAN_H
#define GRID_FEAST_NON_HERMITIAN_H

#include <functional>
#include <numeric>
#include <memory>

NAMESPACE_BEGIN(Grid);

/**
 * Shifted operator adapter:  (z - A)  or  (conj(z) - A^H).
 *
 * Only Op() and AdjOp() are meaningful; BiCGSTAB and GMRES use Op() alone.
 * The remaining LinearOperatorBase virtuals are not implementable for a
 * shifted operator and assert if called.
 */
template<class Field>
class FeastShiftedOperator : public LinearOperatorBase<Field> {
private:
  LinearOperatorBase<Field>& A_;
  LinearOperatorBase<Field>* B_;   // nullptr => B = I (standard problem)
  ComplexD                   z_;
  bool                       useAdjoint_;
  mutable Field*             tmp_;

public:
  FeastShiftedOperator(LinearOperatorBase<Field>& A,
                       LinearOperatorBase<Field>* B = nullptr,
                       Field* tmp = nullptr)
    : A_(A), B_(B), z_(0.0, 0.0), useAdjoint_(false), tmp_(tmp) {}

  /// adj == false : represent (z B - A);  adj == true : represent (conj(z) B^H - A^H)
  void setShift(ComplexD z, bool adj) { z_ = z; useAdjoint_ = adj; }

  void Op(const Field& in, Field& out) {
    ComplexD zc(z_.real(), -z_.imag());
    if (useAdjoint_) {
      A_.AdjOp(in, out);
      if (B_) { B_->AdjOp(in, *tmp_); out = zc * (*tmp_) - out; }
      else    {                       out = zc * in      - out; }
    } else {
      A_.Op(in, out);
      if (B_) { B_->Op(in, *tmp_);    out = z_ * (*tmp_) - out; }
      else    {                       out = z_ * in      - out; }
    }
  }

  void AdjOp(const Field& in, Field& out) {
    // (zB - A)^H = conj(z) B^H - A^H, and (conj(z) B^H - A^H)^H = zB - A
    ComplexD zc(z_.real(), -z_.imag());
    if (useAdjoint_) {
      A_.Op(in, out);
      if (B_) { B_->Op(in, *tmp_);    out = z_ * (*tmp_) - out; }
      else    {                       out = z_ * in      - out; }
    } else {
      A_.AdjOp(in, out);
      if (B_) { B_->AdjOp(in, *tmp_); out = zc * (*tmp_) - out; }
      else    {                       out = zc * in      - out; }
    }
  }

  void OpDiag(const Field&, Field&)                       { assert(0); }
  void OpDir (const Field&, Field&, int, int)             { assert(0); }
  void OpDirAll(const Field&, std::vector<Field>&)        { assert(0); }
  void HermOp(const Field&, Field&)                       { assert(0); }
  void HermOpAndNorm(const Field&, Field&, RealD&, RealD&){ assert(0); }
};

enum FeastQuadrature {
  FeastTrapezoidal,    // uniform in theta; exponentially accurate on a circle
  FeastGaussLegendre   // FEAST's own default
};

template<class Field>
class FeastNonHermitian {
public:
  typedef Eigen::MatrixXcd CMat;
  typedef Eigen::VectorXcd CVec;
  typedef std::complex<double> CplxS;   // Eigen/scalar side
  /// ComplexD is thrust::complex in CUDA builds, so convert explicitly.
  static ComplexD toGrid(const CplxS& c) { return ComplexD(c.real(), c.imag()); }

  /**
   * Solve the shifted system for one right-hand side.
   *   adj == false :  (z        - A  ) out = in
   *   adj == true  :  (conj(z)  - A^H) out = in
   */
  using ShiftSolveFunc =
    std::function<void(ComplexD z, bool adj, const Field& in, Field& out)>;

private:
  LinearOperatorBase<Field>& Linop;
  LinearOperatorBase<Field>* Bop;    // nullptr => standard problem (B = I)
  GridBase*                  Grid_;

  // --- contour (ellipse: centre Emid, horizontal radius r, axis ratio, rotation)
  CplxS    Emid_;
  RealD    radius_;
  RealD    ratio_;      // vertical/horizontal half-axis ratio (1 = circle)
  RealD    angle_;      // rotation of the ellipse, radians
  int      Ne_;         // quadrature nodes on the full contour
  FeastQuadrature quadrature_;

  std::vector<CplxS> nodes_;
  std::vector<CplxS> weights_;

  int   M0_;            // subspace size (must exceed the true count inside C)
  RealD Tolerance_;
  int   maxIter_;

  // outputs
  CVec               evals_;
  std::vector<Field> evecsR_;
  std::vector<Field> evecsL_;
  std::vector<RealD> residuals_;

  std::shared_ptr<FeastShiftedOperator<Field> > shiftedOp_;
  std::shared_ptr<Field>                        shiftScratch_;

public:
  ShiftSolveFunc shiftSolve;
  bool           verbose = true;

  FeastNonHermitian(LinearOperatorBase<Field>& op, GridBase* grid,
                    CplxS Emid, RealD radius, int Ne, int M0,
                    RealD tol = 1e-8, int maxIter = 20,
                    RealD ratio = 1.0, RealD angle = 0.0,
                    FeastQuadrature quad = FeastGaussLegendre)
    : Linop(op), Bop(nullptr), Grid_(grid), Emid_(Emid), radius_(radius), ratio_(ratio),
      angle_(angle), Ne_(Ne), quadrature_(quad), M0_(M0),
      Tolerance_(tol), maxIter_(maxIter)
  {
    assert(Ne_  >= 2  && "FEAST: need at least 2 contour nodes");
    assert(M0_  >= 1  && "FEAST: subspace size must be positive");
    assert(radius_ > 0.0);
    computeContour();
  }

  /**
   * Switch to the generalised problem  A x = lambda B x.
   *
   * The contour quadrature then needs (z B - A)^{-1}, a single shifted solve
   * in the composite operators -- no nested inversion of B.  This is the
   * right way to get the spectrum of B^{-1}A (e.g. PV^{-1}M for domain wall)
   * without ever forming B^{-1} or approximating it by B^dag.
   */
  void setB(LinearOperatorBase<Field>& B)
  {
    // useGMRES()/useBiCGSTAB() capture Bop when they build the shifted
    // operator, so a later setB() would be silently ignored.  Fail loudly
    // instead of quietly solving the wrong problem.
    assert(!shiftedOp_ &&
           "FEAST: setB() must be called before useGMRES()/useBiCGSTAB()");
    Bop = &B;
  }

  CVec               getEvals()     { return evals_;     }
  std::vector<Field> getEvecsR()    { return evecsR_;    }
  std::vector<Field> getEvecsL()    { return evecsL_;    }
  std::vector<RealD> getResiduals() { return residuals_; }
  const std::vector<CplxS>& getNodes()   { return nodes_;   }
  const std::vector<CplxS>& getWeights() { return weights_; }

  /**
   * Install a BiCGSTAB-based default for the shifted solves.  BiCGSTAB is
   * used rather than CG because (z - A) is non-Hermitian for complex z even
   * when A is Hermitian.
   *
   * innerTol can be loose: FEAST tolerates inexact inner solves (the IFEAST
   * observation), the cost being a slower outer rate.
   */
  void useBiCGSTAB(RealD innerTol = 1e-8, int innerMaxIter = 10000)
  {
    shiftScratch_ = std::make_shared<Field>(Grid_);
    shiftedOp_ = std::make_shared<FeastShiftedOperator<Field> >(Linop, Bop, shiftScratch_.get());
    auto op    = shiftedOp_;
    shiftSolve = [op, innerTol, innerMaxIter]
                 (ComplexD z, bool adj, const Field& in, Field& out) {
      op->setShift(z, adj);
      BiCGSTAB<Field> solver(innerTol, innerMaxIter, /*err_on_no_conv=*/false);
      out = Zero();
      solver(*op, in, out);
    };
  }

  /**
   * GMRES default for the shifted solves.
   *
   * Preferred over BiCGSTAB when the contour passes close to the spectrum:
   * (z - A) is then near-singular, and BiCGSTAB can break down (omega -> 0)
   * and return NaN, whereas restarted GMRES cannot break down -- at worst it
   * stagnates.  The price is the restart_length Krylov vectors it must store.
   */
  void useGMRES(RealD innerTol = 1e-8, int innerMaxIter = 10000,
                int restartLength = 50)
  {
    shiftScratch_ = std::make_shared<Field>(Grid_);
    shiftedOp_ = std::make_shared<FeastShiftedOperator<Field> >(Linop, Bop, shiftScratch_.get());
    auto op    = shiftedOp_;
    shiftSolve = [op, innerTol, innerMaxIter, restartLength]
                 (ComplexD z, bool adj, const Field& in, Field& out) {
      op->setShift(z, adj);
      GeneralisedMinimalResidual<Field> solver(innerTol, innerMaxIter,
                                               restartLength,
                                               /*err_on_no_conv=*/false);
      out = Zero();
      solver(*op, in, out);
    };
  }

  /// Is lambda inside the search contour?
  bool insideContour(CplxS lambda) const
  {
    CplxS d = (lambda - Emid_) * CplxS(std::cos(-angle_), std::sin(-angle_));
    RealD u = d.real() / radius_;
    RealD v = d.imag() / (radius_ * ratio_);
    return (u*u + v*v) <= 1.0;
  }

  /**
   * Run FEAST.  Y0 supplies the M0 starting vectors; the same block seeds
   * both the left and right subspaces, which is the usual choice.
   */
  void operator()(const std::vector<Field>& Y0)
  {
    assert((int)Y0.size() >= M0_ && "FEAST: need at least M0 starting vectors");

    std::vector<Field> YR, YL, QR, QL;
    YR.reserve(M0_); YL.reserve(M0_); QR.reserve(M0_); QL.reserve(M0_);
    for (int m = 0; m < M0_; m++) {
      YR.push_back(Y0[m]);
      YL.push_back(Y0[m]);
      QR.push_back(Y0[m]);
      QL.push_back(Y0[m]);
    }

    Field tmp(Grid_), work(Grid_);

    for (int iter = 0; iter < maxIter_; iter++) {

      if (verbose)
        std::cout << GridLogMessage << "FEAST: ---- iteration " << iter
                  << " ----" << std::endl;

      // ---- 1. apply the rational filter: contour integration -------------
      for (int m = 0; m < M0_; m++) { QR[m] = Zero(); QL[m] = Zero(); }

      // A diverging or breaking-down inner solver returns NaN rather than
      // failing, and that NaN then silently poisons Q, Bq and every Ritz
      // value downstream (it shows up only as |det(Bq)| = nan).  Catch it
      // here, where the offending node and column are still known.
      auto guard = [&](const Field& f, int k, int m, const char* which) {
        RealD n = norm2(f);
        if (!std::isfinite(n)) {
          std::cout << GridLogError
                    << "FEAST: inner solve returned a non-finite result at node "
                    << k << " (z = " << nodes_[k] << "), column " << m
                    << ", " << which << " subspace.  The shifted system is"
                    << " likely too ill-conditioned for the chosen solver:"
                    << " try useGMRES(), more contour nodes, or a contour"
                    << " further from the spectrum." << std::endl;
          assert(0 && "FEAST: non-finite inner solve");
        }
      };

      for (int k = 0; k < Ne_; k++) {
        for (int m = 0; m < M0_; m++) {
          shiftSolve(toGrid(nodes_[k]), false, YR[m], tmp);
          guard(tmp, k, m, "right");
          QR[m] = QR[m] + toGrid(weights_[k]) * tmp;

          shiftSolve(toGrid(nodes_[k]), true,  YL[m], tmp);
          guard(tmp, k, m, "left");
          QL[m] = QL[m] + toGrid(std::conj(weights_[k])) * tmp;
        }
        if (verbose)
          std::cout << GridLogMessage << "FEAST:   node " << k << "/" << Ne_
                    << "  z = " << nodes_[k] << std::endl;
      }

      // ---- 2. oblique Rayleigh-Ritz --------------------------------------
      CMat Aq = CMat::Zero(M0_, M0_);
      CMat Bq = CMat::Zero(M0_, M0_);
      Field bwork(Grid_);
      for (int j = 0; j < M0_; j++) {
        Linop.Op(QR[j], work);
        if (Bop) Bop->Op(QR[j], bwork); else bwork = QR[j];
        for (int i = 0; i < M0_; i++) {
          Aq(i, j) = toStdCmplx(innerProduct(QL[i], work ));
          Bq(i, j) = toStdCmplx(innerProduct(QL[i], bwork));
        }
      }

      // Eigen has no complex generalised eigensolver; Bq is well conditioned
      // when M0 oversamples the contour, so reduce to a standard problem.
      Eigen::FullPivLU<CMat> luBq(Bq);
      RealD rcond = std::abs(luBq.determinant());
      if (verbose)
        std::cout << GridLogMessage << "FEAST:   |det(Bq)| = " << rcond << std::endl;
      CMat C = luBq.solve(Aq);

      Eigen::ComplexEigenSolver<CMat> esR(C);
      Eigen::ComplexEigenSolver<CMat> esL(C.adjoint());
      CVec lambdas = esR.eigenvalues();
      CMat U       = esR.eigenvectors();   // right reduced eigenvectors
      CMat W       = esL.eigenvectors();   // left  reduced eigenvectors

      // ---- 3. Ritz vectors, residuals, contour filtering -----------------
      std::vector<int>   keep;
      std::vector<RealD> res;
      std::vector<CplxS> lam;

      evecsR_.clear(); evecsL_.clear();
      std::vector<Field> XR, XL;

      for (int m = 0; m < M0_; m++) {
        CplxS lambda = lambdas(m);

        Field xr(Grid_); xr = Zero();
        Field xl(Grid_); xl = Zero();
        for (int j = 0; j < M0_; j++) {
          xr = xr + QR[j] * toGrid(U(j, m));
          xl = xl + QL[j] * toGrid(W(j, m));
        }
        RealD nr = std::sqrt(norm2(xr));
        RealD nl = std::sqrt(norm2(xl));
        if (nr > 1e-30) xr = xr * (1.0 / nr);
        if (nl > 1e-30) xl = xl * (1.0 / nl);

        XR.push_back(xr);
        XL.push_back(xl);

        Linop.Op(xr, work);
        if (Bop) { Bop->Op(xr, bwork); work = work - toGrid(lambda) * bwork; }
        else     {                     work = work - toGrid(lambda) * xr;    }
        RealD r = std::sqrt(norm2(work));

        lam.push_back(lambda);
        res.push_back(r);
        if (insideContour(lambda)) keep.push_back(m);
      }

      // ---- 4. report & convergence test ----------------------------------
      int nconv = 0;
      for (int idx : keep) if (res[idx] < Tolerance_) nconv++;

      if (verbose) {
        std::cout << GridLogMessage << "FEAST:   inside contour = " << keep.size()
                  << " / " << M0_ << ",  converged = " << nconv << std::endl;
        for (int idx : keep)
          std::cout << GridLogMessage << "FEAST:     lambda = " << lam[idx]
                    << "  |res| = " << res[idx]
                    << (res[idx] < Tolerance_ ? "  *" : "") << std::endl;
      }

      bool done = (!keep.empty() && nconv == (int)keep.size());

      if (done || iter == maxIter_ - 1) {
        evals_.resize(keep.size());
        residuals_.clear();
        evecsR_.clear();
        evecsL_.clear();
        for (size_t i = 0; i < keep.size(); i++) {
          evals_(i) = lam[keep[i]];
          residuals_.push_back(res[keep[i]]);
          evecsR_.push_back(XR[keep[i]]);
          evecsL_.push_back(XL[keep[i]]);
        }
        std::cout << GridLogMessage
                  << (done ? "FEAST: converged with " : "FEAST: maxIter reached with ")
                  << keep.size() << " eigenpairs inside the contour after "
                  << iter + 1 << " iteration(s)." << std::endl;
        return;
      }

      // ---- 5. reseed -----------------------------------------------------
      for (int m = 0; m < M0_; m++) { YR[m] = XR[m]; YL[m] = XL[m]; }
    }
  }

private:
  /// Gauss-Legendre nodes/weights on [-1,1] by Newton iteration.
  static void gaussLegendre(int n, std::vector<RealD>& x, std::vector<RealD>& w)
  {
    x.resize(n); w.resize(n);
    const RealD eps = 1e-15;
    int m = (n + 1) / 2;
    for (int i = 0; i < m; i++) {
      RealD z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), z1, pp;
      int guard = 0;
      do {
        RealD p1 = 1.0, p2 = 0.0;
        for (int j = 0; j < n; j++) {
          RealD p3 = p2; p2 = p1;
          p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
        }
        pp = n * (z * p1 - p2) / (z * z - 1.0);
        z1 = z;
        z  = z1 - p1 / pp;
      } while (std::abs(z - z1) > eps && ++guard < 100);
      x[i]         = -z;
      x[n - 1 - i] =  z;
      w[i]         = 2.0 / ((1.0 - z * z) * pp * pp);
      w[n - 1 - i] = w[i];
    }
  }

  /**
   * Quadrature nodes z_k and weights w_k for
   *     (1/2 pi i) \oint (z-A)^{-1} dz  ~  sum_k w_k (z_k - A)^{-1}
   * on the rotated ellipse  z(theta) = Emid + e^{i angle}(r cos t + i r*ratio sin t).
   *
   * With theta = pi(x+1) mapping x in [-1,1] onto [0,2pi],
   *     w_k = (pi * omega_k) * z'(theta_k) / (2 pi i)
   */
  void computeContour()
  {
    std::vector<RealD> x, w;
    if (quadrature_ == FeastGaussLegendre) {
      gaussLegendre(Ne_, x, w);
    } else {
      x.resize(Ne_); w.resize(Ne_);
      for (int k = 0; k < Ne_; k++) {
        x[k] = -1.0 + (2.0 * k + 1.0) / Ne_;   // midpoint rule in theta
        w[k] = 2.0 / Ne_;
      }
    }

    nodes_.resize(Ne_);
    weights_.resize(Ne_);
    CplxS rot(std::cos(angle_), std::sin(angle_));

    for (int k = 0; k < Ne_; k++) {
      RealD theta = M_PI * (x[k] + 1.0);
      RealD ct = std::cos(theta), st = std::sin(theta);
      CplxS p ( radius_ * ct,          radius_ * ratio_ * st);
      CplxS dp(-radius_ * st,          radius_ * ratio_ * ct);
      nodes_[k]   = Emid_ + rot * p;
      weights_[k] = (M_PI * w[k]) * (rot * dp) / CplxS(0.0, 2.0 * M_PI);
    }

    std::cout << GridLogMessage << "FEAST: contour centre " << Emid_
              << "  radius " << radius_ << "  ratio " << ratio_
              << "  angle " << angle_
              << "  nodes " << Ne_
              << (quadrature_ == FeastGaussLegendre ? " (Gauss-Legendre)"
                                                    : " (trapezoidal)")
              << std::endl;
  }
};

NAMESPACE_END(Grid);

#endif // GRID_FEAST_NON_HERMITIAN_H
