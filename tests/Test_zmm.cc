#include <Grid.h>

using namespace Grid;
using namespace Grid::QCD;


//class DumbComplex{
  

int main(int argc,char **argv)
{
  Grid_init(&argc,&argv);

  // std::vector<int> latt4 = GridDefaultLatt();
  // GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  // GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  // GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  // GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  // std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  // std::vector<int> mpi_layout  = GridDefaultMpi();
  // int threads = GridThread::GetThreads();

  std::vector<int> seeds4({1,2,3,4});
  // std::vector<int> seeds5({5,6,7,8});

  GridSerialRNG sRNG; 
  sRNG.SeedFixedIntegers(seeds4);

  typedef iScalar<iScalar<iMatrix<ComplexD, Nc> > > pooh;

  int iters = 1000;
#pragma acc kernels
  for(int i=0;i<iters;i++){
    //iScalar<iMatrix<float,3> > a;
    //iScalar<iMatrix<float,3> > b;
    //iScalar<iMatrix<float,3> > c = a + b;


    pooh a;
    pooh b;
    pooh c = a+b;


    //ColourMatrixD a;
    //random(sRNG,a);
    
    //ColourMatrixD b;
    //random(sRNG,b);

    //ColourMatrixD c = a+b;
  }

  return 0;
}

