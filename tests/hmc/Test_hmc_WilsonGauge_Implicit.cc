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

  
  std::string arg;

#if 0
  std::vector<int> steps(0);
  if( GridCmdOptionExists(argv,argv+argc,"--MDsteps") ){
    arg= GridCmdOptionPayload(argv,argv+argc,"--MDsteps");
    GridCmdOptionIntVector(arg,steps);
    assert(steps.size()==1);
  }
#endif

   // Typedefs to simplify notation
  typedef GenericHMCRunner<ImplicitMinimumNorm2> HMCWrapper;  // Uses the default minimum norm
//  typedef GenericHMCRunner<ImplicitLeapFrog> HMCWrapper;  // Uses the default minimum norm

  IntegratorParameters MD;
//  MD.name    = std::string("ImplicitLeapFrog");
  MD.name    = std::string("ImplicitMinimumNorm2");
  MD.trajL   = 0.01*std::sqrt(2.);
  MD.MDsteps = 1;
  if( GridCmdOptionExists(argv,argv+argc,"--trajL") ){
    arg= GridCmdOptionPayload(argv,argv+argc,"--trajL");
    std::vector<int> traj(0);
    GridCmdOptionIntVector(arg,traj);
    assert(traj.size()==1);
    MD.trajL *= double(traj[0]);
  }
  MD.RMHMCTol=1e-7;
  MD.RMHMCCGTol=1e-8;
  std::cout << "RMHMCTol= "<<  MD.RMHMCTol<<" RMHMCCGTol= "<<MD.RMHMCCGTol<<std::endl;

  HMCparameters HMCparams;
  HMCparams.StartTrajectory  = 0;
  HMCparams.Trajectories     = 200;
  HMCparams.NoMetropolisUntil=  100;
  // "[HotStart, ColdStart, TepidStart, CheckpointStart]\n";
  HMCparams.StartingType     =std::string("ColdStart");
  HMCparams.Kappa=0.01; //checking against trivial. Pathetic.
  HMCparams.MD = MD;

  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";
  
  HMCWrapper TheHMC(HMCparams);
  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
//  TheHMC.Resources.AddFourDimGridF("gauge_f");
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
  TopParams.interval = 1;
  TopParams.do_smearing = true;
  TopParams.Smearing.steps = 1600;
  TopParams.Smearing.step_size = 0.01;
  TopParams.Smearing.meas_interval = 10;
  TopParams.Smearing.maxTau = 16.0; 
  TheHMC.Resources.AddObservable<QObs>(TopParams);
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
//  RealD beta = 10.0;
  RealD beta = 6.4;
  WilsonGaugeActionR Waction(beta);
  std::cout << "Wilson Gauge beta= " <<beta <<std::endl;
  
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&Waction);
  //Level1.push_back(WGMod.getPtr());
  TheHMC.TheAction.push_back(Level1);
  /////////////////////////////////////////////////////////////

  // HMC parameters are serialisable 
//  TheHMC.Parameters.MD.MDsteps = 20;
//  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  std::cout << "trajL= " <<TheHMC.Parameters.MD.trajL <<" steps= "<<TheHMC.Parameters.MD.MDsteps << " integrator= "<<TheHMC.Parameters.MD.name<<std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();

} // main
