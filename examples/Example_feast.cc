/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./examples/Example_feast.cc

    Non-Hermitian FEAST example for the Wilson Dirac operator.

    Reads a gauge configuration from "config" (NERSC format) and parameters
    from "FeastParams.xml".  Computes the eigenvalues of D_W inside a circular
    (or elliptical) contour in the complex plane.

    Unlike Krylov methods, FEAST targets an explicitly chosen region: the
    contour is an input, not an outcome of a sort criterion.  Pick Emid/radius
    to bracket the modes you want, and M0 comfortably larger than the number
    of eigenvalues you expect inside.

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

struct FeastParameters : Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(FeastParameters,
    RealD,   mass,
    RealD,   EmidRe,      // contour centre, real part
    RealD,   EmidIm,      // contour centre, imaginary part
    RealD,   radius,      // horizontal half-axis
    RealD,   ratio,       // vertical/horizontal ratio (1 = circle)
    RealD,   angle,       // ellipse rotation, radians
    Integer, Ne,          // quadrature nodes on the contour
    Integer, M0,          // subspace size (oversample!)
    Integer, maxIter,
    Integer, gaussLegendre,
    Integer, useGMRES,    // 1 = GMRES (robust near the spectrum), 0 = BiCGSTAB
    Integer, gmresRestart, // GMRES restart length
    RealD,   resid,       // outer (eigenpair) tolerance
    RealD,   innerResid,  // inner solver tolerance
    Integer, innerMaxIter)

  FeastParameters()
    : mass(0.0), EmidRe(1.0), EmidIm(0.0), radius(0.5), ratio(1.0), angle(0.0),
      Ne(8), M0(16), maxIter(20), gaussLegendre(1), useGMRES(1), gmresRestart(50),
      resid(1e-8), innerResid(1e-8), innerMaxIter(10000)
  {}
};

} // namespace Grid

typedef WilsonFermionD WilsonOp;
typedef typename WilsonFermionD::FermionField FermionField;

int main(int argc, char** argv)
{
  Grid_init(&argc, &argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(
      GridDefaultLatt(),
      GridDefaultSimd(Nd, vComplex::Nsimd()),
      GridDefaultMpi());
  GridRedBlackCartesian* UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::vector<int> seeds4({1, 2, 3, 4});
  GridParallelRNG RNG4(UGrid);
  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  FieldMetaData header;
  std::string configFile("config");
  NerscIO::readConfiguration(Umu, header, configFile);
  std::cout << GridLogMessage << "Loaded gauge configuration: " << configFile << std::endl;

  FeastParameters Params;
  {
    XmlReader rd("FeastParams.xml");
    read(rd, "FeastParameters", Params);
  }
  std::cout << GridLogMessage << Params << std::endl;

  std::vector<Complex> boundary = {1, 1, 1, -1};
  WilsonOp::ImplParams WilsonParams(boundary);
  WilsonOp Dwilson(Umu, *UGrid, *UrbGrid, Params.mass, WilsonParams);
  NonHermitianLinearOperator<WilsonOp, FermionField> DLinOp(Dwilson);

  ComplexD Emid(Params.EmidRe, Params.EmidIm);

  FeastNonHermitian<FermionField> feast(
      DLinOp, UGrid, Emid, Params.radius, Params.Ne, Params.M0,
      Params.resid, Params.maxIter, Params.ratio, Params.angle,
      Params.gaussLegendre ? FeastGaussLegendre : FeastTrapezoidal);

  if (Params.useGMRES) {
    std::cout << GridLogMessage << "FEAST inner solver: GMRES (restart "
              << Params.gmresRestart << ")" << std::endl;
    feast.useGMRES(Params.innerResid, Params.innerMaxIter, Params.gmresRestart);
  } else {
    std::cout << GridLogMessage << "FEAST inner solver: BiCGSTAB" << std::endl;
    feast.useBiCGSTAB(Params.innerResid, Params.innerMaxIter);
  }

  // M0 random starting vectors; the same block seeds both subspaces.
  std::vector<FermionField> Y0(Params.M0, FermionField(UGrid));
  for (int m = 0; m < (int)Params.M0; m++) random(RNG4, Y0[m]);

  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " Non-Hermitian FEAST on D_W" << std::endl;
  std::cout << GridLogMessage << "   mass   = " << Params.mass   << std::endl;
  std::cout << GridLogMessage << "   Emid   = " << Emid          << std::endl;
  std::cout << GridLogMessage << "   radius = " << Params.radius << std::endl;
  std::cout << GridLogMessage << "   Ne     = " << Params.Ne     << std::endl;
  std::cout << GridLogMessage << "   M0     = " << Params.M0     << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;

  feast(Y0);

  Eigen::VectorXcd          evals = feast.getEvals();
  std::vector<RealD>        res   = feast.getResiduals();
  std::vector<FermionField> evecs = feast.getEvecsR();

  std::cout << GridLogMessage << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  std::cout << GridLogMessage << " FEAST: " << evals.size()
            << " eigenpairs inside the contour" << std::endl;
  std::cout << GridLogMessage << "*******************************************" << std::endl;
  for (int i = 0; i < (int)evals.size(); i++)
    std::cout << GridLogMessage
              << "  [" << std::setw(3) << i << "]  lambda = " << evals(i)
              << "  |res| = " << res[i] << std::endl;

  std::cout << GridLogMessage << std::endl << "Done" << std::endl;

  Grid_finalize();
  return 0;
}
