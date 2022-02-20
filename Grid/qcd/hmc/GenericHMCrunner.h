/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/GenericHmcRunner.h

Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef GRID_GENERIC_HMC_RUNNER
#define GRID_GENERIC_HMC_RUNNER

#include <unordered_map>

NAMESPACE_BEGIN(Grid);

// very ugly here but possibly resolved if we had a base Reader class
template < class ReaderClass >
class HMCRunnerBase {
public:
  virtual void Run() = 0;
  virtual void initialize(ReaderClass& ) = 0;
};

template <class Implementation, class ImplementationF,
          template <typename, typename, typename> class Integrator,
          class RepresentationsPolicy = NoHirep, class ReaderClass = XmlReader>
class HMCWrapperTemplate: public HMCRunnerBase<ReaderClass> {
public:
  INHERIT_FIELD_TYPES(Implementation);
  typedef Implementation ImplPolicy;  // visible from outside
  template <typename S = NoSmearing<Implementation> >
  using IntegratorType = Integrator<Implementation, S, RepresentationsPolicy>;

  HMCparameters Parameters;
  std::string ParameterFile;
  HMCResourceManager<Implementation, ImplementationF> Resources;

  // The set of actions (keep here for lower level users, for now)
  ActionSet<Field, RepresentationsPolicy> TheAction;

  HMCWrapperTemplate() = default;

  HMCWrapperTemplate(HMCparameters Par){
    Parameters = Par;
  }

  void initialize(ReaderClass & TheReader){
    std::cout  << "Initialization of the HMC" << std::endl;
    Resources.initialize(TheReader);

    // eventually add smearing

    Resources.GetActionSet(TheAction);    
  }


  void ReadCommandLine(int argc, char **argv) {
    std::string arg;

    if (GridCmdOptionExists(argv, argv + argc, "--StartingType")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--StartingType");

      if (arg != "HotStart" && arg != "ColdStart" && arg != "TepidStart" &&
          arg != "CheckpointStart") {
        std::cout << GridLogError << "Unrecognized option in --StartingType\n";
        std::cout
	  << GridLogError
	  << "Valid [HotStart, ColdStart, TepidStart, CheckpointStart]\n";
        exit(1);
      }
      Parameters.StartingType = arg;
    }

    if (GridCmdOptionExists(argv, argv + argc, "--StartingTrajectory")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--StartingTrajectory");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.StartTrajectory = ivec[0];
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Trajectories")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Trajectories");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.Trajectories = ivec[0];
    }

    if (GridCmdOptionExists(argv, argv + argc, "--Thermalizations")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--Thermalizations");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.NoMetropolisUntil = ivec[0];
    }

    if (GridCmdOptionExists(argv, argv + argc, "--ParameterFile")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--ParameterFile");
      ParameterFile = arg;
    }
    if (GridCmdOptionExists(argv, argv + argc, "--MDsteps")) {
      arg = GridCmdOptionPayload(argv, argv + argc, "--MDsteps");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg, ivec);
      Parameters.MD.MDsteps = ivec[0];
    }
  }


  template <class SmearingPolicy>
  void Run(SmearingPolicy &S) {
    Runner(S);
  }

  void Run(){
    NoSmearing<Implementation> S;
    Runner(S);
  }

  //////////////////////////////////////////////////////////////////

private:
  template <class SmearingPolicy>
  void Runner(SmearingPolicy &Smearing) {
    GridCartesian *UGrid = Resources.GetCartesian();
    Resources.AddRNGs();
    Field U(UGrid);
    std::cout << GridLogMessage << "UGrid= " <<UGrid << std::endl;
    GridCartesian *UGrid_f(NULL);
    UGrid_f = Resources.GetCartesianF();
//    Field U_f(UGrid_f);
    std::cout << GridLogMessage << "UGrid_f="<<UGrid_f << std::endl;
//    exit(-42);

    // Can move this outside?
    typedef IntegratorType<SmearingPolicy> TheIntegrator;
    // Metric
#if 0
    std::cout << GridLogMessage << "Trivial metric" << std::endl;
    TrivialMetric<typename Implementation::Field> Mtr;
    TheIntegrator MDynamics(UGrid, Parameters.MD, TheAction, Smearing, Mtr);
#else
    ConjugateGradient<LatticeGaugeField> CG(1.0e-8,10000);
    LaplacianParams LapPar(0.0001, 1.0, 10000, 1e-8, 12, 64);

// Better to pass the generalised momenta to the integrator
//    RealD Kappa = Parameters.Kappa;
    std::cout << GridLogMessage << "LaplacianRat " << std::endl;
//#include "DBW2_3.h.inc"
//#include "Wilson10_3.h.inc"
#include "g_x3_1.h.inc"
//#include "g_x3_3_sym.h.inc"
//#include "Wilson64.h.inc"
//#include "32I2.8Gev.h.inc"
    gpar.tolerance=Parameters.MD.RMHMCCGTol;
    mpar.tolerance=Parameters.MD.RMHMCCGTol;
    std::cout << GridLogMessage << "gpar a0= " << gpar.a0 <<std::endl;
    std::cout << GridLogMessage << " a1= " << gpar.a1 <<std::endl;
    std::cout << GridLogMessage << " b0= " << gpar.b0 <<std::endl;
    std::cout << GridLogMessage << " b1= " << gpar.b1 <<std::endl;
    std::cout << GridLogMessage << " b2= " << gpar.b2 <<std::endl ;;

    std::cout << GridLogMessage << "mpar a0= " << mpar.a0 <<std::endl;
    std::cout << GridLogMessage << " a1= " << mpar.a1 <<std::endl;
    std::cout << GridLogMessage << " b0= " << mpar.b0 <<std::endl;
    std::cout << GridLogMessage << " b1= " << mpar.b1 <<std::endl;
    std::cout << GridLogMessage << " b2= " << mpar.b2 <<std::endl;
    LaplacianAdjointRat<Implementation, ImplementationF> Laplacian(UGrid, UGrid_f,CG, gpar, mpar);
//    Laplacian.RMHMCCGTol=Parameters.MD.RMHMCCGTol;
//    Laplacian.RMHMCTol=Parameters.MD.RMHMCTol;
    std::cout << GridLogMessage << " RMHMCCGTol= " << Parameters.MD.RMHMCCGTol <<std::endl;
    std::cout << GridLogMessage << " RMHMCTol= " << Parameters.MD.RMHMCTol <<std::endl;
    TheIntegrator MDynamics(UGrid, Parameters.MD, TheAction, Smearing, Laplacian);
#endif

    if (Parameters.StartingType == "HotStart") {
      // Hot start
      Resources.SeedFixedIntegers();
      Implementation::HotConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "ColdStart") {
      // Cold start
      Resources.SeedFixedIntegers();
      Implementation::ColdConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "TepidStart") {
      // Tepid start
      Resources.SeedFixedIntegers();
      Implementation::TepidConfiguration(Resources.GetParallelRNG(), U);
    } else if (Parameters.StartingType == "CheckpointStart") {
      // CheckpointRestart
      Resources.GetCheckPointer()->CheckpointRestore(Parameters.StartTrajectory, U,
						     Resources.GetSerialRNG(),
						     Resources.GetParallelRNG());
    } else {
      // others
      std::cout << GridLogError << "Unrecognized StartingType\n";
      std::cout
	<< GridLogError
	<< "Valid [HotStart, ColdStart, TepidStart, CheckpointStart]\n";
      exit(1);
    }

    Smearing.set_Field(U);

    HybridMonteCarlo<TheIntegrator> HMC(Parameters, MDynamics,
                                        Resources.GetSerialRNG(),
                                        Resources.GetParallelRNG(), 
                                        Resources.GetObservables(), U);

    // Run it
    HMC.evolve();
  }
};

// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunner = HMCWrapperTemplate<PeriodicGimplD, PeriodicGimplF, Integrator>;
//template <template <typename, typename, typename> class Integrator>
//using GenericHMCRunnerF = HMCWrapperTemplate<PeriodicGimplF, Integrator>;
//template <template <typename, typename, typename> class Integrator>
//using GenericHMCRunnerD = HMCWrapperTemplate<PeriodicGimplD, Integrator>;


// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunner = HMCWrapperTemplate<ConjugateGimplD, ConjugateGimplF, Integrator>;
//template <template <typename, typename, typename> class Integrator>
//using ConjugateHMCRunnerF = HMCWrapperTemplate<ConjugateGimplF, Integrator>;
//template <template <typename, typename, typename> class Integrator>
//using ConjugateHMCRunnerD = HMCWrapperTemplate<ConjugateGimplD, Integrator>;



template <class RepresentationsPolicy,
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerHirep =
				     HMCWrapperTemplate<PeriodicGimplD, PeriodicGimplF, Integrator, RepresentationsPolicy>;

template <class Implementation, class ImplementationF, class RepresentationsPolicy, 
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerTemplate = HMCWrapperTemplate<Implementation, ImplementationF, Integrator, RepresentationsPolicy>;

typedef HMCWrapperTemplate<ScalarImplD, ScalarImplF, MinimumNorm2, ScalarFields>
ScalarGenericHMCRunner;

typedef HMCWrapperTemplate<ScalarAdjImplD, ScalarAdjImplF, MinimumNorm2, ScalarMatrixFields>
ScalarAdjGenericHMCRunner;

template <int Colours> 
using ScalarNxNAdjGenericHMCRunner = HMCWrapperTemplate < ScalarNxNAdjImplD<Colours>, ScalarNxNAdjImplF<Colours>, ForceGradient, ScalarNxNMatrixFields<Colours> >;

NAMESPACE_END(Grid);

#endif  // GRID_GENERIC_HMC_RUNNER
