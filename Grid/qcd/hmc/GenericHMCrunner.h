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

template <class Implementation,
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
  HMCResourceManager<Implementation> Resources;

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
    auto UGrid = Resources.GetCartesian();
    Resources.AddRNGs();
    Field U(UGrid);

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
    RealD Kappa = Parameters.Kappa;
#if 0
    std::cout << GridLogMessage << "Kappa = " << Kappa << std::endl;
    LaplacianAdjointField<PeriodicGimplR> Laplacian(UGrid, CG, LapPar, Kappa);
#else
    std::cout << GridLogMessage << "LaplacianRat " << std::endl;
#if 0
//working?
    LaplacianRatParams gpar(1),mpar(1);
    gpar.offset = 1.;
    gpar.a0[0] = 1.+Kappa*0.5;
    gpar.a1[0] = -Kappa*0.5;
    gpar.b0[0] = 1.-Kappa*0.5;
    gpar.b1[0] = +Kappa*0.5;
    gpar.b2=1.;
    mpar.offset = 1.;
    mpar.a0[0] = 1.-Kappa*0.5;
    mpar.a1[0] = +Kappa*0.5;
    mpar.b0[0] = 1.+Kappa*0.5;
    mpar.b1[0] = -Kappa*0.5;
    mpar.b2=1.;
#else
//g_x3_1
    LaplacianRatParams gpar(4,10000),mpar(4,10000);
    gpar. offset = -1.806640805644669e+02;
    gpar. a0[0] = 3.828256737535615e+00;
    gpar. a1[0] = -1.759884800801541e+00;
    gpar. b0[0] = 3.720710594848511e+00;
    gpar. b1[0] = -3.529025938486833e+00;
    gpar. a0[1] = 1.385427197085892e+01;
    gpar. a1[1] = -1.210145337093361e+04;
    gpar. b0[1] = 5.230416808032104e+00;
    gpar. b1[1] = 5.038496433918291e+01;
    gpar. a0[2] = -2.615764258902624e-01;
    gpar. a1[2] = 1.768889411282898e-01;
    gpar. b0[2] = 2.659251295152384e+00;
    gpar. b1[2] = -3.254117569635528e+00;
    gpar. a0[3] = 1.523953680699539e+05;
    gpar. a1[3] = 4.651204348978543e+04;
    gpar. b0[3] = 3.615577944657219e+02;
    gpar. b1[3] = 1.158681982729600e+02;
    gpar.b2=1.;

    mpar. offset = -5.535134581681094e-03;
    mpar. a0[0] = -6.897276105582431e+00;
    mpar. a1[0] = -3.760976465900414e+00;
    mpar. b0[0] = -3.792989299317706e+01;
    mpar. b1[0] = -1.652205896546682e+01;
    mpar. a0[1] = -1.049418884789610e+01;
    mpar. a1[1] = 2.705292434239074e+00;
    mpar. b0[1] = 6.255355656904058e+01;
    mpar. b1[1] = -7.586738463469761e+00;
    mpar. a0[2] = -2.027141884855678e-02;
    mpar. a1[2] = 2.446143264550702e-03;
    mpar. b0[2] = 4.042328199062386e+00;
    mpar. b1[2] = -3.638856523270800e+00;
    mpar. a0[3] = 1.482983537063678e-03;
    mpar. a1[3] = -9.754596728239950e-04;
    mpar. b0[3] = 2.639273996355545e+00;
    mpar. b1[3] = -3.240812191556963e+00;
    mpar.b2=1.;
#endif
    for(int i=0;i<4;i++){
       gpar.a1[i] *=16.;
       gpar.b1[i] *=16.;
       mpar.a1[i] *=16.;
       mpar.b1[i] *=16.;
    }
    gpar.b2 *= 16.*16.;
    mpar.b2 *= 16.*16.;
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
    LaplacianAdjointRat<PeriodicGimplR> Laplacian(UGrid, CG, gpar, mpar);
#endif
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
using GenericHMCRunner = HMCWrapperTemplate<PeriodicGimplR, Integrator>;
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunnerF = HMCWrapperTemplate<PeriodicGimplF, Integrator>;
template <template <typename, typename, typename> class Integrator>
using GenericHMCRunnerD = HMCWrapperTemplate<PeriodicGimplD, Integrator>;


// These are for gauge fields, default integrator MinimumNorm2
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunner = HMCWrapperTemplate<ConjugateGimplR, Integrator>;
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunnerF = HMCWrapperTemplate<ConjugateGimplF, Integrator>;
template <template <typename, typename, typename> class Integrator>
using ConjugateHMCRunnerD = HMCWrapperTemplate<ConjugateGimplD, Integrator>;



template <class RepresentationsPolicy,
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerHirep =
				     HMCWrapperTemplate<PeriodicGimplR, Integrator, RepresentationsPolicy>;

template <class Implementation, class RepresentationsPolicy, 
          template <typename, typename, typename> class Integrator>
using GenericHMCRunnerTemplate = HMCWrapperTemplate<Implementation, Integrator, RepresentationsPolicy>;

typedef HMCWrapperTemplate<ScalarImplR, MinimumNorm2, ScalarFields>
ScalarGenericHMCRunner;

typedef HMCWrapperTemplate<ScalarAdjImplR, MinimumNorm2, ScalarMatrixFields>
ScalarAdjGenericHMCRunner;

template <int Colours> 
using ScalarNxNAdjGenericHMCRunner = HMCWrapperTemplate < ScalarNxNAdjImplR<Colours>, ForceGradient, ScalarNxNMatrixFields<Colours> >;

NAMESPACE_END(Grid);

#endif  // GRID_GENERIC_HMC_RUNNER
