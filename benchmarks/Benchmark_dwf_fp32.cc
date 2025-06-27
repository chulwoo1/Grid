 /*************************************************************************************
    Grid physics library, www.github.com/paboyle/Grid
    Source file: ./benchmarks/Benchmark_dwf.cc
    Copyright (C) 2015

    Author: Peter Boyle <paboyle@ph.ed.ac.uk>
    Author: paboyle <paboyle@ph.ed.ac.uk>

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
#include <Grid/Grid.h>
#ifdef GRID_CUDA
#define CUDA_PROFILE
#endif

#ifdef CUDA_PROFILE
#include <cuda_profiler_api.h>
#endif

using namespace std;
using namespace Grid;

////////////////////////
/// Move to domains ////
////////////////////////

Gamma::Algebra Gmu [] = {
			 Gamma::Algebra::GammaX,
			 Gamma::Algebra::GammaY,
			 Gamma::Algebra::GammaZ,
			 Gamma::Algebra::GammaT
};

//void Benchmark(int Ls, Coordinate Dirichlet);
void Laplace();

#if 0
struct controls {
  int Opt;
  int CommsOverlap;
  Grid::CartesianCommunicator::CommunicatorPolicy_t CommsAsynch;
};

struct time_statistics{
  double mean;
  double err;
  double min;
  double max;

  void statistics(std::vector<double> v){
      double sum = std::accumulate(v.begin(), v.end(), 0.0);
      mean = sum / v.size();

      std::vector<double> diff(v.size());
      std::transform(v.begin(), v.end(), diff.begin(), [=](double x) { return x - mean; });
      double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
      err = std::sqrt(sq_sum / (v.size()*(v.size() - 1)));

      auto result = std::minmax_element(v.begin(), v.end());
      min = *result.first;
      max = *result.second;
}
};
#endif
void Benchmark(int Ls, Coordinate Dirichlet,bool Sloppy);

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);


  int threads = GridThread::GetThreads();

  int Ls=16;
  for(int i=0;i<argc;i++) {
    if(std::string(argv[i]) == "-Ls"){
      std::stringstream ss(argv[i+1]); ss >> Ls;
    }
  }

  //////////////////
  // With comms
  //////////////////
  Coordinate Dirichlet(Nd+1,0);

  std::cout << "\n\n\n\n\n\n" <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  std::cout << GridLogMessage<< " Testing with full communication " <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  
//  Benchmark(Ls,Dirichlet);
//  Laplace();
  Benchmark(Ls,Dirichlet,false);

  std::cout << "\n\n\n\n\n\n" <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  std::cout << GridLogMessage<< " Testing with sloppy communication " <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  
  Benchmark(Ls,Dirichlet,true);

  //////////////////
  // Domain decomposed
  //////////////////
  /*
  Coordinate latt4  = GridDefaultLatt();
  Coordinate mpi    = GridDefaultMpi();
  Coordinate CommDim(Nd);
  Coordinate shm;
  GlobalSharedMemory::GetShmDims(mpi,shm);


  std::cout << "\n\n\n\n\n\n" <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  //  std::cout << GridLogMessage<< " Testing without internode communication " <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;

  for(int d=0;d<Nd;d++) CommDim[d]= (mpi[d]/shm[d])>1 ? 1 : 0;
  Dirichlet[0] = 0;
  Dirichlet[1] = CommDim[0]*latt4[0]/mpi[0] * shm[0];
  Dirichlet[2] = CommDim[1]*latt4[1]/mpi[1] * shm[1];
  Dirichlet[3] = CommDim[2]*latt4[2]/mpi[2] * shm[2];
  Dirichlet[4] = CommDim[3]*latt4[3]/mpi[3] * shm[3];

  Benchmark(Ls,Dirichlet,false);

  std::cout << "\n\n\n\n\n\n" <<std::endl;

  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;
  std::cout << GridLogMessage<< " Testing with sloppy communication " <<std::endl;
  std::cout << GridLogMessage<< "++++++++++++++++++++++++++++++++++++++++++++++++" <<std::endl;

  for(int d=0;d<Nd;d++) CommDim[d]= mpi[d]>1 ? 1 : 0;
  
  Benchmark(Ls,Dirichlet,true);
  */
  
  Grid_finalize();
  exit(0);
}
void Benchmark(int Ls, Coordinate Dirichlet,bool sloppy)
{
  Coordinate latt4 = GridDefaultLatt();
  GridLogLayout();

  long unsigned int single_site_flops = 8*Nc*(7+16*Nc);

  std::vector<int> seeds4({1,2,3,4});
  std::vector<int> seeds5({5,6,7,8});
#define SINGLE
#ifdef SINGLE
  typedef vComplexF          Simd;
  typedef LatticeFermionF    FermionField;
  typedef LatticeGaugeFieldF GaugeField;
  typedef LatticeColourMatrixF ColourMatrixField;
  typedef DomainWallFermionF FermionAction;
#else
  typedef vComplexD          Simd;
  typedef LatticeFermionD    FermionField;
  typedef LatticeGaugeFieldD GaugeField;
  typedef LatticeColourMatrixD ColourMatrixField;
  typedef DomainWallFermionD FermionAction;
#endif
  
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,Simd::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::cout << GridLogMessage << "Initialising 4d RNG" << std::endl;
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedUniqueString(std::string("The 4D RNG"));

  std::cout << GridLogMessage << "Initialising 5d RNG" << std::endl;
  GridParallelRNG          RNG5(FGrid);  RNG5.SeedUniqueString(std::string("The 5D RNG"));

 
  FermionField src   (FGrid); random(RNG5,src);
#if 0
  src = Zero();
  {
    Coordinate origin({0,0,0,latt4[2]-1,0});
    SpinColourVectorF tmp;
    tmp=Zero();
    tmp()(0)(0)=Complex(-2.0,0.0);
    std::cout << " source site 0 " << tmp<<std::endl;
    pokeSite(tmp,src,origin);
  }
#else
  RealD N2 = 1.0/::sqrt(norm2(src));
  src = src*N2;
#endif

  FermionField result(FGrid); result=Zero();
  FermionField    ref(FGrid);    ref=Zero();
  FermionField    tmp(FGrid);
  FermionField    err(FGrid);

  std::cout << GridLogMessage << "Drawing gauge field" << std::endl;
  GaugeField Umu(UGrid);
  GaugeField UmuCopy(UGrid);
  SU<Nc>::HotConfiguration(RNG4,Umu);
  //  SU<Nc>::ColdConfiguration(Umu);
  UmuCopy=Umu;
  std::cout << GridLogMessage << "Random gauge initialised " << std::endl;

  ////////////////////////////////////
  // Apply BCs
  ////////////////////////////////////
  Coordinate Block(4);
  for(int d=0;d<4;d++)  Block[d]= Dirichlet[d+1];

  std::cout << GridLogMessage << "Applying BCs for Dirichlet Block5 " << Dirichlet << std::endl;
  std::cout << GridLogMessage << "Applying BCs for Dirichlet Block4 " << Block << std::endl;

  DirichletFilter<GaugeField> Filter(Block);
  Filter.applyFilter(Umu);
  
  ////////////////////////////////////
  // Naive wilson implementation
  ////////////////////////////////////
  std::vector<ColourMatrixField> U(4,UGrid);
  for(int mu=0;mu<Nd;mu++){
    U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
  }

  std::cout << GridLogMessage << "Setting up Cshift based reference " << std::endl;

  if (1)
  {
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      tmp = Cshift(src,mu+1,1);
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = U_v[ss]*tmp_v[Ls*ss+s];
	  }
	}
      }
      ref=ref + tmp - Gamma(Gmu[mu])*tmp;

      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      tmp =Cshift(tmp,mu+1,-1);
      ref=ref + tmp + Gamma(Gmu[mu])*tmp;
    }
    ref = -0.5*ref;
  }

  RealD mass=0.1;
  RealD M5  =1.8;

  RealD NP = UGrid->_Nprocessors;
  RealD NN = UGrid->NodeCount();

  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Kernel options --dslash-generic, --dslash-unroll, --dslash-asm" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermionR::Dhop                  "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<Simd::Nsimd()<<std::endl;
  std::cout << GridLogMessage<< "* VComplex size is "<<sizeof(Simd)<< " B"<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*****************************************************************" <<std::endl;

  FermionAction::ImplParams p;
  p.dirichlet=Dirichlet;
  FermionAction Dw(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,p);
  Dw.SloppyComms(sloppy);
  Dw.ImportGauge(Umu);
  
  int ncall =300;
  RealD n2e;
  
  if (1) {
    FGrid->Barrier();
    Dw.Dhop(src,result,0);
    std::cout<<GridLogMessage<<"Called warmup"<<std::endl;
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.Dhop(src,result,0);
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=single_site_flops*volume*ncall;

    auto nsimd = Simd::Nsimd();
    auto simdwidth = sizeof(Simd);

    // RF: Nd Wilson * Ls, Nd gauge * Ls, Nc colors
    double data_rf = volume * ((2*Nd+1)*Nd*Nc + 2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    // mem: Nd Wilson * Ls, Nd gauge, Nc colors
    double data_mem = (volume * (2*Nd+1)*Nd*Nc + (volume/Ls) *2*Nd*Nc*Nc) * simdwidth / nsimd * ncall / (1024.*1024.*1024.);

    std::cout<<GridLogMessage << "Called Dw "<<ncall<<" times in "<<t1-t0<<" us"<<std::endl;
    std::cout<<GridLogMessage << "mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per rank =  "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "mflop/s per node =  "<< flops/(t1-t0)/NN<<std::endl;
    err = ref-result;
    n2e = norm2(err);
    std::cout<<GridLogMessage << "norm diff   "<< n2e<< "  Line "<<__LINE__ <<std::endl;

    if(( n2e>1.0e-4) ) {
      std::cout<<GridLogMessage << "WRONG RESULT" << std::endl;
      FGrid->Barrier();
      std::cout<<GridLogMessage << "RESULT" << std::endl;
      //      std::cout << result<<std::endl;
      std::cout << norm2(result)<<std::endl;
      std::cout<<GridLogMessage << "REF" << std::endl;
      std::cout << norm2(ref)<<std::endl;
      std::cout<<GridLogMessage << "ERR" << std::endl;
      std::cout << norm2(err)<<std::endl;
      FGrid->Barrier();
      exit(-1);
    }
    assert (n2e< 1.0e-4 );
  }

  if (1)
  { // Naive wilson dag implementation
    ref = Zero();
    for(int mu=0;mu<Nd;mu++){

      //    ref =  src - Gamma(Gamma::Algebra::GammaX)* src ; // 1+gamma_x
      tmp = Cshift(src,mu+1,1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	autoView( U_v  , U[mu]  , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    int i=s+Ls*ss;
	    ref_v[i]+= U_v[ss]*(tmp_v[i] + Gamma(Gmu[mu])*tmp_v[i]); ;
	  }
	}
      }
      
      {
	autoView( tmp_v  , tmp  , CpuWrite);
	autoView( U_v  , U[mu]  , CpuRead);
	autoView( src_v, src    , CpuRead);
	for(int ss=0;ss<U[mu].Grid()->oSites();ss++){
	  for(int s=0;s<Ls;s++){
	    tmp_v[Ls*ss+s] = adj(U_v[ss])*src_v[Ls*ss+s];
	  }
	}
      }
      //      tmp =adj(U[mu])*src;
      tmp =Cshift(tmp,mu+1,-1);
      {
	autoView( ref_v, ref, CpuWrite);
	autoView( tmp_v, tmp, CpuRead);
	for(int i=0;i<ref_v.size();i++){
	  ref_v[i]+= tmp_v[i] - Gamma(Gmu[mu])*tmp_v[i]; ;
	}
      }
    }
    ref = -0.5*ref;
  }

  Dw.Dhop(src,result,DaggerYes);

  std::cout << GridLogMessage << "----------------------------------------------------------------" << std::endl;
  std::cout << GridLogMessage << "Compare to naive wilson implementation Dag to verify correctness" << std::endl;
  std::cout << GridLogMessage << "----------------------------------------------------------------" << std::endl;

  std::cout<<GridLogMessage << "Called DwDag"<<std::endl;
  std::cout<<GridLogMessage << "norm dag result "<< norm2(result)<<std::endl;
  std::cout<<GridLogMessage << "norm dag ref    "<< norm2(ref)<<std::endl;
  err = ref-result;
  n2e= norm2(err);
  std::cout<<GridLogMessage << "norm dag diff   "<< n2e<< "  Line "<<__LINE__ <<std::endl;

  assert((n2e)<1.0e-4);
  
  FermionField src_e (FrbGrid);
  FermionField src_o (FrbGrid);
  FermionField r_e   (FrbGrid);
  FermionField r_o   (FrbGrid);
  FermionField r_eo  (FGrid);

  std::cout<<GridLogMessage << "Calling Deo and Doe and //assert Deo+Doe == Dunprec"<<std::endl;
  pickCheckerboard(Even,src_e,src);
  pickCheckerboard(Odd,src_o,src);

  std::cout<<GridLogMessage << "src_e"<<norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "src_o"<<norm2(src_o)<<std::endl;


  // S-direction is INNERMOST and takes no part in the parity.
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  std::cout << GridLogMessage<< "* Benchmarking DomainWallFermion::DhopEO                "<<std::endl;
  std::cout << GridLogMessage<< "* Vectorising space-time by "<<Simd::Nsimd()<<std::endl;
#ifdef GRID_OMP
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
  if ( WilsonKernelsStatic::Comms == WilsonKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential comms compute" <<std::endl;
#endif
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using GENERIC Nc WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptHandUnroll) std::cout << GridLogMessage<< "* Using Nc=3       WilsonKernels" <<std::endl;
  if ( WilsonKernelsStatic::Opt == WilsonKernelsStatic::OptInlineAsm ) std::cout << GridLogMessage<< "* Using Asm Nc=3   WilsonKernels" <<std::endl;
  std::cout << GridLogMessage<< "*********************************************************" <<std::endl;
  {
    FGrid->Barrier();
    Dw.DhopEO(src_o,r_e,DaggerNo);
    double t0=usecond();
    for(int i=0;i<ncall;i++){
      Dw.DhopEO(src_o,r_e,DaggerNo);
    }
    double t1=usecond();
    FGrid->Barrier();

    double volume=Ls;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
    double flops=(single_site_flops*volume*ncall)/2.0;

    std::cout<<GridLogMessage << "Deo mflop/s =   "<< flops/(t1-t0)<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per rank   "<< flops/(t1-t0)/NP<<std::endl;
    std::cout<<GridLogMessage << "Deo mflop/s per node   "<< flops/(t1-t0)/NN<<std::endl;
  }
  Dw.DhopEO(src_o,r_e,DaggerNo);
  Dw.DhopOE(src_e,r_o,DaggerNo);
  Dw.Dhop  (src  ,result,DaggerNo);

  std::cout<<GridLogMessage << "r_e"<<norm2(r_e)<<std::endl;
  std::cout<<GridLogMessage << "r_o"<<norm2(r_o)<<std::endl;
  std::cout<<GridLogMessage << "res"<<norm2(result)<<std::endl;

  setCheckerboard(r_eo,r_o);
  setCheckerboard(r_eo,r_e);

  err = r_eo-result;
  n2e= norm2(err);
  std::cout<<GridLogMessage << "norm diff   "<< n2e<<std::endl;
  assert(n2e<1.0e-4);

  pickCheckerboard(Even,src_e,err);
  pickCheckerboard(Odd,src_o,err);
  std::cout<<GridLogMessage << "norm diff even  "<< norm2(src_e)<<std::endl;
  std::cout<<GridLogMessage << "norm diff odd   "<< norm2(src_o)<<std::endl;

  assert(norm2(src_e)<1.0e-4);
  assert(norm2(src_o)<1.0e-4);
}

void Laplace()
  {
    double mflops;
    double mflops_best = 0;
    double mflops_worst= 0;
    std::vector<double> mflops_all;

    ///////////////////////////////////////////////////////
    // Set/Get the layout & grid size
    ///////////////////////////////////////////////////////
    int threads = GridThread::GetThreads();
    Coordinate mpi = GridDefaultMpi(); assert(mpi.size()==4);
//    Coordinate local({L,L,L,L});
//    Coordinate latt4({local[0]*mpi[0],local[1]*mpi[1],local[2]*mpi[2],local[3]*mpi[3]});
    Coordinate latt4 = GridDefaultLatt();
    GridLogLayout();
    
    GridCartesian         * TmpGrid   = SpaceTimeGrid::makeFourDimGrid(latt4,
								       GridDefaultSimd(Nd,vComplex::Nsimd()),
								       GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
//    NN_global=NN;
    uint64_t SHM=NP/NN;


    ///////// Welcome message ////////////
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
//    std::cout<<GridLogMessage << "Benchmark Laplace on "<<L<<"^4 local volume "<<std::endl;
    std::cout<<GridLogMessage << "* Global volume  : "<<GridCmdVectorIntToString(latt4)<<std::endl;
    std::cout<<GridLogMessage << "* ranks          : "<<NP  <<std::endl;
    std::cout<<GridLogMessage << "* nodes          : "<<NN  <<std::endl;
    std::cout<<GridLogMessage << "* ranks/node     : "<<SHM <<std::endl;
    std::cout<<GridLogMessage << "* ranks geom     : "<<GridCmdVectorIntToString(mpi)<<std::endl;
    std::cout<<GridLogMessage << "* Using "<<threads<<" threads"<<std::endl;
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;

    ///////// Lattice Init ////////////
    GridCartesian         * FGrid   = SpaceTimeGrid::makeFourDimGrid(latt4, GridDefaultSimd(Nd,vComplexF::Nsimd()),GridDefaultMpi());
    GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);
    
    ///////// RNG Init ////////////
    std::vector<int> seeds4({1,2,3,4});
    GridParallelRNG          RNG4(FGrid);  RNG4.SeedFixedIntegers(seeds4);
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    RealD mass=0.1;
    RealD c1=9.0/8.0;
    RealD c2=-1.0/24.0;
    RealD u0=1.0;

    typedef LatticeGaugeFieldF Gauge;
    
    Gauge Umu(FGrid);  SU<Nc>::HotConfiguration(RNG4,Umu); 

    typedef typename PeriodicGimplF::LinkField GaugeLinkFieldF;

    ///////// Source preparation ////////////
    GaugeLinkFieldF src   (FGrid); random(RNG4,src);
    GaugeLinkFieldF r_eo  (FGrid);
  
    {

      const int num_cases = 1;
      std::string fmt("G/O/C  ");
      
      controls Cases [] = {
	{  StaggeredKernelsStatic::OptGeneric   ,  StaggeredKernelsStatic::CommsAndCompute  ,CartesianCommunicator::CommunicatorPolicyConcurrent  },
      }; 

      for(int c=0;c<num_cases;c++) {
        CovariantAdjointLaplacianStencil<PeriodicGimplF,typename PeriodicGimplF::LinkField> LapStencilF(FGrid);
        QuadLinearOperator<CovariantAdjointLaplacianStencil<PeriodicGimplF,typename PeriodicGimplF::LinkField>,PeriodicGimplF::LinkField> QuadOpF(LapStencilF,c2,c1,1.);
        LapStencilF.GaugeImport(Umu);
	

	StaggeredKernelsStatic::Comms = Cases[c].CommsOverlap;
	StaggeredKernelsStatic::Opt   = Cases[c].Opt;
	CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);
      
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	if ( StaggeredKernelsStatic::Opt == StaggeredKernelsStatic::OptGeneric   ) std::cout << GridLogMessage<< "* Using Stencil Nc Laplace" <<std::endl;
	if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute ) std::cout << GridLogMessage<< "* Using Overlapped Comms/Compute" <<std::endl;
	if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsThenCompute) std::cout << GridLogMessage<< "* Using sequential Comms/Compute" <<std::endl;
	std::cout << GridLogMessage<< "* SINGLE precision "<<std::endl;
	std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
	
	int nwarm = 10;
	double t0=usecond();
	FGrid->Barrier();
	for(int i=0;i<nwarm;i++){
//	  Ds.DhopEO(src_o,r_e,DaggerNo);
//          QuadOpF.HermOp(src,r_eo);
          LapStencilF.M(src,r_eo);
	}
	FGrid->Barrier();
	double t1=usecond();
	uint64_t ncall = 500;

	FGrid->Broadcast(0,&ncall,sizeof(ncall));

	//	std::cout << GridLogMessage << " Estimate " << ncall << " calls per second"<<std::endl;

	time_statistics timestat;
	std::vector<double> t_time(ncall);
	for(uint64_t i=0;i<ncall;i++){
	  t0=usecond();
//	  Ds.DhopEO(src_o,r_e,DaggerNo);
//          QuadOpF.HermOp(src,r_eo);
          LapStencilF.M(src,r_eo);
	  t1=usecond();
	  t_time[i] = t1-t0;
	}

    GaugeLinkFieldF r_eo2  (FGrid);
    GaugeLinkFieldF err  (FGrid);
	{
		std::vector<GaugeLinkFieldF> U(Nd,FGrid);
//		Gauge  UmuF;
//		precisionChange(UmuF,Umu);
	    GaugeLinkFieldF tmp(FGrid);
    	GaugeLinkFieldF tmp2(FGrid);
      for (int mu = 0; mu < Nd; mu++) 
		U[mu]=PeekIndex<LorentzIndex>(Umu, mu);

//    for (int nu = 0; nu < Nd; nu++) {
      r_eo2 = Zero();
//      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
//      GaugeLinkField out_nu(out.Grid());
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(src, mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * src * U[mu];
        r_eo2 += tmp + Cshift(tmp2, mu, -1) - 2.0 * src;
      }
//      out_nu = (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
//      PokeIndex<LorentzIndex>(out, out_nu, nu);
//    }

	}
	FGrid->Barrier();
	std::cout<<GridLogMessage << "Lap result "<< norm2(r_eo)<<std::endl;
    std::cout<<GridLogMessage << "Lap ref    "<< norm2(r_eo2)<<std::endl;
  err = r_eo-r_eo2;
  auto n2e= norm2(err);
  std::cout<<GridLogMessage << "Lap diff   "<< n2e<< "  Line "<<__LINE__ <<std::endl;
	
	double volume=1;  for(int mu=0;mu<Nd;mu++) volume=volume*latt4[mu];
//Quad
//	double flops=(2*2*8*216.0*volume);
	double flops=(2*8*216.0*volume);
	double mf_hi, mf_lo, mf_err;
	
	timestat.statistics(t_time);
	mf_hi = flops/timestat.min;
	mf_lo = flops/timestat.max;
	mf_err= flops/timestat.min * timestat.err/timestat.mean;

	mflops = flops/timestat.mean;
	mflops_all.push_back(mflops);
	if ( mflops_best == 0   ) mflops_best = mflops;
	if ( mflops_worst== 0   ) mflops_worst= mflops;
	if ( mflops>mflops_best ) mflops_best = mflops;
	if ( mflops<mflops_worst) mflops_worst= mflops;
	
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Laplace mflop/s =   "<< mflops << " ("<<mf_err<<") " << mf_lo<<"-"<<mf_hi <<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Laplace mflop/s per rank   "<< mflops/NP<<std::endl;
	std::cout<<GridLogMessage << std::fixed << std::setprecision(1)<<"Laplace mflop/s per node   "<< mflops/NN<<std::endl;
	FGrid->Barrier();
      
      }

      std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
      std::cout<<GridLogMessage << latt4 <<"  Quad Best  mflop/s        =   "<< mflops_best << " ; " << mflops_best/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage << latt4 <<"  Quad Worst mflop/s        =   "<< mflops_worst<< " ; " << mflops_worst/NN<<" per node " <<std::endl;
      std::cout<<GridLogMessage <<fmt << std::endl;
      std::cout<<GridLogMessage ;
	FGrid->Barrier();

      for(int i=0;i<mflops_all.size();i++){
	std::cout<<mflops_all[i]/NN<<" ; " ;
      }
      std::cout<<std::endl;
    }
    std::cout<<GridLogMessage << "=================================================================================="<<std::endl;
//    return mflops_best;
  }
