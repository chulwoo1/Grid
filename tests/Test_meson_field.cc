#include <Grid/Grid.h>

using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  const int Ls=12;
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});

  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedFixedIntegers(seeds5);

  LatticeFermion src   (FGrid); random(RNG5,src);
  LatticeFermion result(FGrid); result=Zero();
  LatticeFermion result2(FGrid); result2=Zero();

  LatticeGaugeField Umu(UGrid); random(RNG4,Umu);
  std::vector<LatticeColourMatrix> U(4,UGrid);

  Umu=Zero();
  for(int nn=0;nn<4;nn++){
    random(RNG4,U[nn]);
    PokeIndex<LorentzIndex>(Umu,U[nn],nn);
  }

  RealD mass=0.1;
  RealD M5  =1.8;

  MobiusFermionD Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mass, M5, 1.5,0.5);

  Ddwf.M(src,result);
  for (int i=0;i<100;i++) {
    Ddwf.M(src,result2);
    double eps2 = norm2(closure(result - result2)) / norm2(src);
    std::cout << GridLogMessage << "Test M reproducibility: eps[" << i << "]^2 = " << eps2 << std::endl;
  }

  std::cout << "Done!" << std::endl;
  Grid_finalize();

  return 0;
}
