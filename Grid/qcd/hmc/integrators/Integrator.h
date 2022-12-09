/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/hmc/integrators/Integrator.h

Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: Guido Cossu <cossu@post.kek.jp>
Author: Chulwoo Jung <chulwoo@bnl.gov>

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
			   //--------------------------------------------------------------------
#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

#include <memory>

NAMESPACE_BEGIN(Grid);

class IntegratorParameters: Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(IntegratorParameters,
				  std::string, name,      // name of the integrator
				  unsigned int, MDsteps,  // number of outer steps
				  RealD, RMHMCTol,
				  RealD, RMHMCCGTol,
				  RealD, trajL)           // trajectory length

  IntegratorParameters(int MDsteps_ = 10, RealD trajL_ = 1.0)
  : MDsteps(MDsteps_),
    trajL(trajL_),RMHMCTol(1e-8),RMHMCCGTol(1e-8) {};

  template <class ReaderClass, typename std::enable_if<isReader<ReaderClass>::value, int >::type = 0 >
  IntegratorParameters(ReaderClass & Reader)
  {
    std::cout << GridLogMessage << "Reading integrator\n";
    read(Reader, "Integrator", *this);
  }

  void print_parameters() const {
    std::cout << GridLogMessage << "[Integrator] Type               : " << name << std::endl;
    std::cout << GridLogMessage << "[Integrator] Trajectory length  : " << trajL << std::endl;
    std::cout << GridLogMessage << "[Integrator] Number of MD steps : " << MDsteps << std::endl;
    std::cout << GridLogMessage << "[Integrator] Step size          : " << trajL/MDsteps << std::endl;
  }
};

/*! @brief Class for Molecular Dynamics management */
template <class FieldImplementation_, class SmearingPolicy, class RepresentationPolicy>
class Integrator {
protected:
  typedef FieldImplementation_ FieldImplementation;
  typedef typename FieldImplementation::Field MomentaField;  //for readability
  typedef typename FieldImplementation::Field Field;

  int levels;  // number of integration levels
  double t_U;  // Track time passing on each level and for U and for P
  std::vector<double> t_P;  

//  MomentaField P;
  GeneralisedMomenta<FieldImplementation > P;
  SmearingPolicy& Smearer;
  RepresentationPolicy Representations;
  IntegratorParameters Params;

  const ActionSet<Field, RepresentationPolicy> as;

  void update_P(Field& U, int level, double ep) 
  {
    t_P[level] += ep;
    update_P(P.Mom, U, level, ep);

    std::cout << GridLogIntegrator << "[" << level << "] P " << " dt " << ep << " : t_P " << t_P[level] << std::endl;
  }

  void update_P2(Field& U, int level, double ep) 
  {
    t_P[level] += ep;
    update_P2(P.Mom, U, level, ep);

    std::cout << GridLogIntegrator << "[" << level << "] P " << " dt " << ep << " : t_P " << t_P[level] << std::endl;
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _updateP 
  {
    template <class FieldType, class GF, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep,
                    GF& Mom, GF& U, double ep) {
      for (int a = 0; a < repr_set.size(); ++a) {
        FieldType forceR(U.Grid());
        // Implement smearing only for the fundamental representation now
        repr_set.at(a)->deriv(Rep.U, forceR);
        GF force = Rep.RtoFundamentalProject(forceR);  // Ta for the fundamental rep
        Real force_abs = std::sqrt(norm2(force)/(U.Grid()->gSites()));
        std::cout << GridLogIntegrator << "Hirep Force average: " << force_abs << std::endl;
	Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
      }
    }
  } update_P_hireps{};

 
  void update_P(MomentaField& Mom, Field& U, int level, double ep) {
    // input U actually not used in the fundamental case
    // Fundamental updates, include smearing

    for (int a = 0; a < as[level].actions.size(); ++a) {

      double start_full = usecond();
      Field force(U.Grid());
      conformable(U.Grid(), Mom.Grid());

      Field& Us = Smearer.get_U(as[level].actions.at(a)->is_smeared);
      double start_force = usecond();

      std::cout << GridLogMessage << "AuditForce["<<level<<"]["<<a<<"] before"<<std::endl;
      AUDIT();
      
      as[level].actions.at(a)->deriv_timer_start();
      as[level].actions.at(a)->deriv(Us, force);  // deriv should NOT include Ta
      as[level].actions.at(a)->deriv_timer_stop();

      std::cout << GridLogMessage << "AuditForce["<<level<<"]["<<a<<"] after"<<std::endl;
      AUDIT();

      std::cout << GridLogIntegrator << "Smearing (on/off): " << as[level].actions.at(a)->is_smeared << std::endl;
      auto name = as[level].actions.at(a)->action_name();
      if (as[level].actions.at(a)->is_smeared) Smearer.smeared_force(force);

      force = FieldImplementation::projectForce(force); // Ta for gauge fields
      double end_force = usecond();

      //      DumpSliceNorm("force ",force,Nd-1);
      MomFilter->applyFilter(force);
      std::cout << GridLogIntegrator << " update_P : Level [" << level <<"]["<<a <<"] "<<name<<" dt "<<ep<<  std::endl;
      DumpSliceNorm("force filtered ",force,Nd-1);
      
      Real force_abs   = std::sqrt(norm2(force)/U.Grid()->gSites()); //average per-site norm.  nb. norm2(latt) = \sum_x norm2(latt[x]) 
      Real impulse_abs = force_abs * ep * HMC_MOMENTUM_DENOMINATOR;    

      Real force_max   = std::sqrt(maxLocalNorm2(force));
      Real impulse_max = force_max * ep * HMC_MOMENTUM_DENOMINATOR;    

      as[level].actions.at(a)->deriv_log(force_abs,force_max,impulse_abs,impulse_max);
      
      std::cout << GridLogIntegrator<< "["<<level<<"]["<<a<<"] dt           : " << ep <<" "<<name<<std::endl;
      std::cout << GridLogIntegrator<< "["<<level<<"]["<<a<<"] Force average: " << force_abs <<" "<<name<<std::endl;
      std::cout << GridLogIntegrator<< "["<<level<<"]["<<a<<"] Force max    : " << force_max <<" "<<name<<std::endl;
      std::cout << GridLogIntegrator<< "["<<level<<"]["<<a<<"] Fdt average  : " << impulse_abs <<" "<<name<<std::endl;
      std::cout << GridLogIntegrator<< "["<<level<<"]["<<a<<"] Fdt max      : " << impulse_max <<" "<<name<<std::endl;

      Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
      double end_full = usecond();
      double time_full  = (end_full - start_full) / 1e3;
      double time_force = (end_force - start_force) / 1e3;
      std::cout << GridLogMessage << "["<<level<<"]["<<a<<"] P update elapsed time: " << time_full << " ms (force: " << time_force << " ms)"  << std::endl;

    }

    // Force from the other representations
    as[level].apply(update_P_hireps, Representations, Mom, U, ep);
  }

  void update_P2(MomentaField& Mom, Field& U, int level, double ep) {
    // input U actually not used in the fundamental case
    // Fundamental updates, include smearing

    std::cout << GridLogIntegrator << "U before update_P2: " << std::sqrt(norm2(U)) << std::endl;
    // Generalised momenta  
    // Derivative of the kinetic term must be computed before
    // Mom is the momenta and gets updated by the 
    // actions derivatives
    MomentaField MomDer(P.Mom.Grid());
    P.M.ImportGauge(U);
    P.DerivativeU(P.Mom, MomDer);
    std::cout << GridLogIntegrator << "MomDer update_P2: " << std::sqrt(norm2(MomDer)) << std::endl;
//    Mom -= MomDer * ep;
    Mom -= MomDer * ep * HMC_MOMENTUM_DENOMINATOR;
    std::cout << GridLogIntegrator << "Mom update_P2: " << std::sqrt(norm2(Mom)) << std::endl;

    // Auxiliary fields
    P.update_auxiliary_momenta(ep*0.5 );
    P.AuxiliaryFieldsDerivative(MomDer);
    std::cout << GridLogIntegrator << "MomDer(Aux) update_P2: " << std::sqrt(norm2(Mom)) << std::endl;
//    Mom -= MomDer * ep;
    Mom -= MomDer * ep * HMC_MOMENTUM_DENOMINATOR;
    P.update_auxiliary_momenta(ep*0.5 );

    for (int a = 0; a < as[level].actions.size(); ++a) {
      double start_full = usecond();
      Field force(U.Grid());
      conformable(U.Grid(), Mom.Grid());

      Field& Us = Smearer.get_U(as[level].actions.at(a)->is_smeared);
      double start_force = usecond();
      as[level].actions.at(a)->deriv(Us, force);  // deriv should NOT include Ta

      std::cout << GridLogIntegrator << "Smearing (on/off): " << as[level].actions.at(a)->is_smeared << std::endl;
      if (as[level].actions.at(a)->is_smeared) Smearer.smeared_force(force);
      force = FieldImplementation::projectForce(force); // Ta for gauge fields
      double end_force = usecond();
      Real force_abs = std::sqrt(norm2(force)/U.Grid()->gSites());
      std::cout << GridLogIntegrator << "["<<level<<"]["<<a<<"] Force average: " << force_abs << std::endl;
      Mom -= force * ep* HMC_MOMENTUM_DENOMINATOR;; 
      double end_full = usecond();
      double time_full  = (end_full - start_full) / 1e3;
      double time_force = (end_force - start_force) / 1e3;
      std::cout << GridLogMessage << "["<<level<<"]["<<a<<"] P update elapsed time: " << time_full << " ms (force: " << time_force << " ms)"  << std::endl;
    }

    // Force from the other representations
    as[level].apply(update_P_hireps, Representations, Mom, U, ep);
  }

  void implicit_update_P(Field& U, int level, double ep, bool intermediate = false) {
    t_P[level] += ep;

    std::cout << GridLogIntegrator << "[" << level << "] P "
              << " dt " << ep << " : t_P " << t_P[level] << std::endl;
    std::cout << GridLogIntegrator << "U before implicit_update_P: " << std::sqrt(norm2(U)) << std::endl;
    // Fundamental updates, include smearing
    MomentaField Msum(P.Mom.Grid());
    Msum = Zero();
    for (int a = 0; a < as[level].actions.size(); ++a) {
      // Compute the force terms for the lagrangian part
      // We need to compute the derivative of the actions
      // only once
      Field force(U.Grid());
      conformable(U.Grid(), P.Mom.Grid());
      Field& Us = Smearer.get_U(as[level].actions.at(a)->is_smeared);
      as[level].actions.at(a)->deriv(Us, force);  // deriv should NOT include Ta

      std::cout << GridLogIntegrator << "Smearing (on/off): " << as[level].actions.at(a)->is_smeared << std::endl;
      if (as[level].actions.at(a)->is_smeared) Smearer.smeared_force(force);
      force = FieldImplementation::projectForce(force);  // Ta for gauge fields
      Real force_abs = std::sqrt(norm2(force) / U.Grid()->gSites());
      std::cout << GridLogIntegrator << "|Force| site average: " << force_abs
                << std::endl;
      Msum += force;
    }

    MomentaField NewMom = P.Mom;
    MomentaField OldMom = P.Mom;
    double threshold = Params.RMHMCTol;
    P.M.ImportGauge(U);
    MomentaField MomDer(P.Mom.Grid());
    MomentaField MomDer1(P.Mom.Grid());
    MomentaField AuxDer(P.Mom.Grid());
    MomDer1 = Zero();
    MomentaField diff(P.Mom.Grid());
    double factor = 2.0;
    if (intermediate){
      P.DerivativeU(P.Mom, MomDer1);
      factor = 1.0;
    }
//    std::cout << GridLogIntegrator << "MomDer1 implicit_update_P: " << std::sqrt(norm2(MomDer1)) << std::endl;

    // Auxiliary fields
    P.update_auxiliary_momenta(ep*0.5 );
    P.AuxiliaryFieldsDerivative(AuxDer);
    Msum += AuxDer;
    

    // Here run recursively
    int counter = 1;
    RealD RelativeError;
    do {
      std::cout << GridLogIntegrator << "UpdateP implicit step "<< counter << std::endl;

      // Compute the derivative of the kinetic term
      // with respect to the gauge field
      P.DerivativeU(NewMom, MomDer);
      Real force_abs = std::sqrt(norm2(MomDer) / U.Grid()->gSites());
      std::cout << GridLogIntegrator << "|Force| laplacian site average: " << force_abs
                << std::endl;

      NewMom = P.Mom - ep* 0.5 * HMC_MOMENTUM_DENOMINATOR * (2.0*Msum + factor*MomDer + MomDer1);// simplify
      diff = NewMom - OldMom;
      counter++;
      RelativeError = std::sqrt(norm2(diff))/std::sqrt(norm2(NewMom));
      std::cout << GridLogIntegrator << "UpdateP RelativeError: " << RelativeError << std::endl;
      OldMom = NewMom;
    } while (RelativeError > threshold);

    P.Mom = NewMom;
    std::cout << GridLogIntegrator << "NewMom implicit_update_P: " << std::sqrt(norm2(NewMom)) << std::endl;

    // update the auxiliary fields momenta    
    P.update_auxiliary_momenta(ep*0.5 );
  }

  void update_U(Field& U, double ep) 
  {
    update_U(P.Mom, U, ep);

    t_U += ep;
    int fl = levels - 1;
    std::cout << GridLogIntegrator << "   " << "[" << fl << "] U " << " dt " << ep << " : t_U " << t_U << std::endl;
  }
  
  void update_U(MomentaField& Mom, Field& U, double ep) 
  {
    MomentaField MomFiltered(Mom.Grid());
    MomFiltered = Mom;
    MomFilter->applyFilter(MomFiltered);

    // exponential of Mom*U in the gauge fields case
    FieldImplementation::update_field(MomFiltered, U, ep);

    // Update the smeared fields, can be implemented as observer
    Smearer.set_Field(U);

    // Update the higher representations fields
    Representations.update(U);  // void functions if fundamental representation
  }

  void implicit_update_U(Field&U, double ep){
    t_U += ep;
    int fl = levels - 1;
    std::cout << GridLogIntegrator << "   " << "[" << fl << "] U " << " dt " << ep << " : t_U " << t_U << std::endl;
    std::cout << GridLogIntegrator << "U before implicit_update_U: " << std::sqrt(norm2(U)) << std::endl;

    MomentaField Mom1(P.Mom.Grid());
    MomentaField Mom2(P.Mom.Grid());
    RealD RelativeError;
    Field diff(U.Grid());
    Real threshold =  Params.RMHMCTol;
    int counter = 1;
    int MaxCounter = 100;

    Field OldU = U;
    Field NewU = U;

    P.M.ImportGauge(U);
    P.DerivativeP(Mom1); // first term in the derivative 
    std::cout << GridLogIntegrator << "implicit_update_U: Mom1: " << std::sqrt(norm2(Mom1)) << std::endl;

    P.update_auxiliary_fields(ep*0.5);


    MomentaField sum=Mom1;
    do {
      std::cout << GridLogIntegrator << "UpdateU implicit step "<< counter << std::endl;
      
      P.DerivativeP(Mom2); // second term in the derivative, on the updated U
      std::cout << GridLogIntegrator << "implicit_update_U: Mom1: " << std::sqrt(norm2(Mom1)) << std::endl;
      sum = (Mom1 + Mom2);

      for (int mu = 0; mu < Nd; mu++) {
        auto Umu = PeekIndex<LorentzIndex>(U, mu);
        auto Pmu = PeekIndex<LorentzIndex>(sum, mu);
        Umu = expMat(Pmu, ep * 0.5, 12) * Umu;
        PokeIndex<LorentzIndex>(NewU, ProjectOnGroup(Umu), mu);
      }

      diff = NewU - OldU;
      RelativeError = std::sqrt(norm2(diff))/std::sqrt(norm2(NewU));
      std::cout << GridLogIntegrator << "UpdateU RelativeError: " << RelativeError << std::endl;
      
      P.M.ImportGauge(NewU);
      OldU = NewU; // some redundancy to be eliminated
      counter++;
    } while (RelativeError > threshold && counter < MaxCounter);

    U = NewU;
    std::cout << GridLogIntegrator << "NewU implicit_update_U: " << std::sqrt(norm2(U)) << std::endl;
    P.update_auxiliary_fields(ep*0.5);
  }

  virtual void step(Field& U, int level, int first, int last) = 0;

public:
  Integrator(GridBase* grid, IntegratorParameters Par,
             ActionSet<Field, RepresentationPolicy>& Aset,
             SmearingPolicy& Sm, Metric<MomentaField>& M)
    : Params(Par),
      as(Aset),
      P(grid, M),
      levels(Aset.size()),
      Smearer(Sm),
      Representations(grid) 
  {
    t_P.resize(levels, 0.0);
    t_U = 0.0;
    // initialization of smearer delegated outside of Integrator
  };

  virtual ~Integrator() {}

  virtual std::string integrator_name() = 0;

  //Set the momentum filter allowing for manipulation of the conjugate momentum
  void setMomentumFilter(const MomentumFilterBase<MomentaField> &filter){
  MomFilter = &filter;
  }

  //Access the conjugate momentum
  const MomentaField & getMomentum() const{ return P; }



  void reset_timer(void)
  {
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
        as[level].actions.at(actionID)->reset_timer();
      }
    }
  }
  void print_timer(void)
  {
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
    std::cout << GridLogMessage << " Refresh cumulative timings "<<std::endl;
    std::cout << GridLogMessage << "--------------------------- "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	std::cout << GridLogMessage 
		  << as[level].actions.at(actionID)->action_name()
		  <<"["<<level<<"]["<< actionID<<"] "
		  << as[level].actions.at(actionID)->refresh_us*1.0e-6<<" s"<< std::endl;
      }
    }
    std::cout << GridLogMessage << "--------------------------- "<<std::endl;
    std::cout << GridLogMessage << " Action cumulative timings "<<std::endl;
    std::cout << GridLogMessage << "--------------------------- "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	std::cout << GridLogMessage 
		  << as[level].actions.at(actionID)->action_name()
		  <<"["<<level<<"]["<< actionID<<"] "
		  << as[level].actions.at(actionID)->S_us*1.0e-6<<" s"<< std::endl;
      }
    }
    std::cout << GridLogMessage << "--------------------------- "<<std::endl;
    std::cout << GridLogMessage << " Force cumulative timings "<<std::endl;
    std::cout << GridLogMessage << "------------------------- "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	std::cout << GridLogMessage 
		  << as[level].actions.at(actionID)->action_name()
		  <<"["<<level<<"]["<< actionID<<"] "
		  << as[level].actions.at(actionID)->deriv_us*1.0e-6<<" s"<< std::endl;
      }
    }
    std::cout << GridLogMessage << "--------------------------- "<<std::endl;
    std::cout << GridLogMessage << " Force average size "<<std::endl;
    std::cout << GridLogMessage << "------------------------- "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	std::cout << GridLogMessage 
		  << as[level].actions.at(actionID)->action_name()
		  <<"["<<level<<"]["<< actionID<<"] :\n\t\t "
		  <<" force max " << as[level].actions.at(actionID)->deriv_max_average()
		  <<" norm "      << as[level].actions.at(actionID)->deriv_norm_average()
		  <<" Fdt max  "  << as[level].actions.at(actionID)->Fdt_max_average()
		  <<" Fdt norm "  << as[level].actions.at(actionID)->Fdt_norm_average()
		  <<" calls "     << as[level].actions.at(actionID)->deriv_num
		  << std::endl;
      }
    }
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::"<< std::endl;
  }
  
  void print_parameters()
  {
    std::cout << GridLogMessage << "[Integrator] Name : "<< integrator_name() << std::endl;
    Params.print_parameters();
  }

  void print_actions()
  {
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
    std::cout << GridLogMessage << "[Integrator] Action summary: "<<std::endl;
    for (int level = 0; level < as.size(); ++level) {
      std::cout << GridLogMessage << "[Integrator] ---- Level: "<< level << std::endl;
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	std::cout << GridLogMessage << "["<< as[level].actions.at(actionID)->action_name() << "] ID: " << actionID << std::endl;
	std::cout << as[level].actions.at(actionID)->LogParameters();
      }
    }
    std::cout << GridLogMessage << ":::::::::::::::::::::::::::::::::::::::::"<< std::endl;
  }

  void reverse_momenta()
  {
//    P *= -1.0;
    P.Mom *= -1.0;
    P.AuxMom *= -1.0;
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _refresh {
    template <class FieldType, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep, GridSerialRNG & sRNG, GridParallelRNG& pRNG) {
      for (int a = 0; a < repr_set.size(); ++a){
        repr_set.at(a)->refresh(Rep.U, sRNG, pRNG);
      
	std::cout << GridLogDebug << "Hirep refreshing pseudofermions" << std::endl;
      }
    }
  } refresh_hireps{};

  // Initialization of momenta and actions
  void refresh(Field& U, GridSerialRNG & sRNG, GridParallelRNG& pRNG) 
  {
    assert(P.Mom.Grid() == U.Grid());
    std::cout << GridLogIntegrator << "Integrator refresh\n";

    std::cout << GridLogIntegrator << "Generating momentum" << std::endl;
//    FieldImplementation::generate_momenta(P.Mom, pRNG);
    P.M.ImportGauge(U);
    P.MomentaDistribution(sRNG,pRNG);


    // Update the smeared fields, can be implemented as observer
    // necessary to keep the fields updated even after a reject
    // of the Metropolis
    std::cout << GridLogIntegrator << "Updating smeared fields" << std::endl;
    Smearer.set_Field(U);
    // Set the (eventual) representations gauge fields

    std::cout << GridLogIntegrator << "Updating representations" << std::endl;
    Representations.update(U);

    // The Smearer is attached to a pointer of the gauge field
    // automatically gets the correct field
    // whether or not has been accepted in the previous sweep
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
        // get gauge field from the SmearingPolicy and
        // based on the boolean is_smeared in actionID
	auto name = as[level].actions.at(actionID)->action_name();
        std::cout << GridLogMessage << "refresh [" << level << "][" << actionID << "] "<<name << std::endl;

        Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);

	std::cout << GridLogMessage << "AuditRefresh["<<level<<"]["<<actionID<<"] before"<<std::endl;
	AUDIT();
	as[level].actions.at(actionID)->refresh_timer_start();
        as[level].actions.at(actionID)->refresh(Us, sRNG, pRNG);
	as[level].actions.at(actionID)->refresh_timer_stop();
	std::cout << GridLogMessage << "AuditRefresh["<<level<<"]["<<actionID<<"] after"<<std::endl;
	AUDIT();
      }

      // Refresh the higher representation actions
      as[level].apply(refresh_hireps, Representations, pRNG);
    }
  }

  // to be used by the actionlevel class to iterate
  // over the representations
  struct _S {
    template <class FieldType, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep, int level, RealD& H) {
      
      for (int a = 0; a < repr_set.size(); ++a) {
        RealD Hterm = repr_set.at(a)->S(Rep.U);
        std::cout << GridLogMessage << "S Level " << level << " term " << a << " H Hirep = " << Hterm << std::endl;
        H += Hterm;

      }
    }
  } S_hireps{};

  // Calculate action
  RealD S(Field& U) 
  {  // here also U not used

    std::cout.precision(15);
    std::cout << GridLogIntegrator << "Integrator action\n";
    std::cout.precision(15);
    static RealD Saux=0.,Smom=0.,Sg=0.;

    RealD H = - FieldImplementation::FieldSquareNorm(P.Mom)/HMC_MOMENTUM_DENOMINATOR; // - trace (P*P)/denom
    std::cout << GridLogMessage << "S:FieldSquareNorm H_p = " << H << "\n";
    std::cout << GridLogMessage << "S:dSField = " << H-Smom << "\n";
    Smom=H;
    P.M.ImportGauge(U);
    RealD Hterm = - P.MomentaAction();
//    H = - P.MomentaAction()/HMC_MOMENTUM_DENOMINATOR;
    std::cout << GridLogMessage << "S:Momentum action H_p = " << Hterm << "\n";
    std::cout << GridLogMessage << "S:dSMom = " << Hterm-Saux << "\n";
    Saux=Hterm;
    H = Hterm;

    // Actions
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
	AUDIT();
        // get gauge field from the SmearingPolicy and
        // based on the boolean is_smeared in actionID
        Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] action eval " << std::endl;
	        as[level].actions.at(actionID)->S_timer_start();
        Hterm = as[level].actions.at(actionID)->S(Us);
   	        as[level].actions.at(actionID)->S_timer_stop();
//        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] H = " << Hterm << std::endl;
        std::cout << GridLogMessage << "S: action [" << level << "][" << actionID << "] H = " << Hterm << std::endl;
        std::cout << GridLogMessage << "S:dSg = " << Hterm-Sg << "\n";
        Sg=Hterm;
        H += Hterm;
	AUDIT();
      }
      as[level].apply(S_hireps, Representations, level, H);
    }

    std::cout << GridLogMessage << "S:Total  H = " << H << "\n";
    return H;
  }

  struct _Sinitial {
    template <class FieldType, class Repr>
    void operator()(std::vector<Action<FieldType>*> repr_set, Repr& Rep, int level, RealD& H) {
      
      for (int a = 0; a < repr_set.size(); ++a) {
	AUDIT();
        RealD Hterm = repr_set.at(a)->Sinitial(Rep.U);
	AUDIT();
        std::cout << GridLogMessage << "Sinitial Level " << level << " term " << a << " H Hirep = " << Hterm << std::endl;
        H += Hterm;

      }
    }
  } Sinitial_hireps{};

  RealD Sinitial(Field& U) 
  {  // here also U not used

    std::cout << GridLogIntegrator << "Integrator initial action\n";

    RealD H = - FieldImplementation::FieldSquareNorm(P)/HMC_MOMENTUM_DENOMINATOR; // - trace (P*P)/denom

    RealD Hterm;

    // Actions
    for (int level = 0; level < as.size(); ++level) {
      for (int actionID = 0; actionID < as[level].actions.size(); ++actionID) {
        // get gauge field from the SmearingPolicy and
        // based on the boolean is_smeared in actionID
        Field& Us = Smearer.get_U(as[level].actions.at(actionID)->is_smeared);
        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] action eval " << std::endl;
	        as[level].actions.at(actionID)->S_timer_start();
	AUDIT();
        Hterm = as[level].actions.at(actionID)->Sinitial(Us);
   	        as[level].actions.at(actionID)->S_timer_stop();
	AUDIT();
        std::cout << GridLogMessage << "S [" << level << "][" << actionID << "] H = " << Hterm << std::endl;
        H += Hterm;
      }
      as[level].apply(Sinitial_hireps, Representations, level, H);
    }

    return H;
  }

  
  void integrate(Field& U) 
  {
    AUDIT();
    // reset the clocks
    t_U = 0;
    for (int level = 0; level < as.size(); ++level) {
      t_P[level] = 0;
    }

    for (int stp = 0; stp < Params.MDsteps; ++stp) {  // MD step
      int first_step = (stp == 0);
      int last_step = (stp == Params.MDsteps - 1);
      this->step(U, 0, first_step, last_step);
    }

    // Check the clocks all match on all levels
    for (int level = 0; level < as.size(); ++level) {
      assert(fabs(t_U - t_P[level]) < 1.0e-6);  // must be the same
      std::cout << GridLogIntegrator << " times[" << level << "]= " << t_P[level] << " " << t_U << std::endl;
    }
    AUDIT();

    FieldImplementation::Project(U);
    AUDIT();

    // and that we indeed got to the end of the trajectory
    assert(fabs(t_U - Params.trajL) < 1.0e-6);

  }

};

NAMESPACE_END(Grid);

#endif  // INTEGRATOR_INCLUDED

