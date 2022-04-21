/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2015

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Guido Cossu <guido.cossu@ed.ac.uk>

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

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/Grid.h>

int main(int argc, char **argv) 
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

   // Typedefs to simplify notation
//  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef GenericHMCRunner<ForceGradient> HMCWrapper;  // Uses the default minimum norm

  HMCWrapper TheHMC;

//  TheHMC.Parameters.MD.name    = std::string("MinimumNorm2");
  TheHMC.Parameters.MD.name    = std::string("ForceGradient");
    

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  TopologyObsParameters TopParams;
  TopParams.interval = 5;
  TopParams.do_smearing = true;
  TopParams.Smearing.steps = 200;
  TopParams.Smearing.step_size = 0.01;
  TopParams.Smearing.meas_interval = 50;
  TopParams.Smearing.maxTau = 2.0; 
//  TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
#if 1
  RealD beta = 6.4;
  WilsonGaugeActionD Waction(beta);
  std::cout << "Wilson Gauge beta= " <<beta <<std::endl;
#else
//  RBC c_1 for DBW2
  RealD beta = 1.0038;
  RealD c_1 = -1.4088;
  RBCGaugeActionR Waction(beta,c_1);
  std::cout << "P+R Gauge beta= " <<beta <<"c_1= "<<c_1 <<std::endl;
#endif
  
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  //Level1.push_back(WGMod.getPtr());
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable 
    // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
//    HMCparams.StartingType     =std::string("ColdStart");
//    HMCparams.Kappa=0.0;
//    HMCparams.MD = MD;
  TheHMC.Parameters.StartingType     =std::string("ColdStart");
  TheHMC.Parameters.Kappa=0.0;
  TheHMC.Parameters.MD.MDsteps = 20;
  TheHMC.Parameters.MD.trajL   = 0.1*std::sqrt(2.);

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  if( GridCmdOptionExists(argv,argv+argc,"--trajL") ){
    std::string arg= GridCmdOptionPayload(argv,argv+argc,"--trajL");
    std::vector<int> traj(0);
    GridCmdOptionIntVector(arg,traj);
    assert(traj.size()==1);
    TheHMC.Parameters.MD.trajL *= double(traj[0]);
  }

  std::cout << "trajL= " <<TheHMC.Parameters.MD.trajL <<" steps= "<<TheHMC.Parameters.MD.MDsteps << " integrator= "<<TheHMC.Parameters.MD.name<<std::endl;

  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
