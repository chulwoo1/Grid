    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourEvenOddRatio.h

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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD_RATIO_H
//#include<Grid/Grid.h>
//#include<Grid/parallelIO/MetaData.h>

NAMESPACE_BEGIN(Grid);

#if 0
// PLEASE FIX 
template<class vobj> static std::string getFormatStringLocal (void)
{
  std::string format;
  typedef typename getPrecision<vobj>::real_scalar_type stype;
  if ( sizeof(stype) == sizeof(float) ) {
    format = std::string("IEEE32BIG");
  }
  if ( sizeof(stype) == sizeof(double) ) {
    format = std::string("IEEE64BIG");
  }
  return format;
}

template <class fobj, class sobj>
struct PFUnmunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

  void operator()(sobj &in, fobj &out) {
    // take word by word and transform accoding to the status
    fobj_stype *out_buffer = (fobj_stype *)&out;
    sobj_stype *in_buffer = (sobj_stype *)&in;
    size_t fobj_words = sizeof(out) / sizeof(fobj_stype);
    size_t sobj_words = sizeof(in) / sizeof(sobj_stype);
    assert(fobj_words == sobj_words);

    for (unsigned int word = 0; word < sobj_words; word++)
      out_buffer[word] = in_buffer[word];  // type conversion on the fly

  }
};

template <class fobj, class sobj>
struct PFMunger {
  typedef typename getPrecision<fobj>::real_scalar_type fobj_stype;
  typedef typename getPrecision<sobj>::real_scalar_type sobj_stype;

  void operator()(fobj &in, sobj &out) {
    // take word by word and transform accoding to the status
    fobj_stype *in_buffer = (fobj_stype *)&in;
    sobj_stype *out_buffer = (sobj_stype *)&out;
    size_t fobj_words = sizeof(in) / sizeof(fobj_stype);
    size_t sobj_words = sizeof(out) / sizeof(sobj_stype);
    assert(fobj_words == sobj_words);

    for (unsigned int word = 0; word < sobj_words; word++)
      out_buffer[word] = in_buffer[word];  // type conversion on the fly

  }
};
#endif



    ///////////////////////////////////////
    // Two flavour ratio
    ///////////////////////////////////////
    template<class Impl>
    class TwoFlavourEvenOddRatioPseudoFermionAction : public Action<typename Impl::GaugeField> {
    public:
      INHERIT_IMPL_TYPES(Impl);
      
    private:
      FermionOperator<Impl> & NumOp;// the basic operator
      FermionOperator<Impl> & DenOp;// the basic operator

      OperatorFunction<FermionField> &DerivativeSolver;
      OperatorFunction<FermionField> &ActionSolver;
      OperatorFunction<FermionField> &HeatbathSolver;

      FermionField PhiOdd;   // the pseudo fermion field for this trajectory
      FermionField PhiEven;  // the pseudo fermion field for this trajectory

      RealD RefreshAction;
      
    public:
      TwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
                                                FermionOperator<Impl>  &_DenOp, 
                                                OperatorFunction<FermionField> & DS,
                                                OperatorFunction<FermionField> & AS ) : 
      TwoFlavourEvenOddRatioPseudoFermionAction(_NumOp,_DenOp, DS,AS,AS) {};

      TwoFlavourEvenOddRatioPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
                                                FermionOperator<Impl>  &_DenOp, 
                                                OperatorFunction<FermionField> & DS,
                                                OperatorFunction<FermionField> & AS, OperatorFunction<FermionField> & HS) :
      NumOp(_NumOp), 
      DenOp(_DenOp), 
      DerivativeSolver(DS), 
      ActionSolver(AS),
      HeatbathSolver(HS),
      PhiEven(_NumOp.FermionRedBlackGrid()),
      PhiOdd(_NumOp.FermionRedBlackGrid()) 
        {
          conformable(_NumOp.FermionGrid(), _DenOp.FermionGrid());
          conformable(_NumOp.FermionRedBlackGrid(), _DenOp.FermionRedBlackGrid());
          conformable(_NumOp.GaugeGrid(), _DenOp.GaugeGrid());
          conformable(_NumOp.GaugeRedBlackGrid(), _DenOp.GaugeRedBlackGrid());
        };

      virtual std::string action_name(){
	std::stringstream sstream;
	sstream<<"TwoFlavourEvenOddRatioPseudoFermionAction det("<<DenOp.Mass()<<") / det("<<NumOp.Mass()<<")";
	return sstream.str();
      }

      virtual std::string LogParameters(){
	std::stringstream sstream;
	sstream<< GridLogMessage << "["<<action_name()<<"] -- No further parameters "<<std::endl;
	return sstream.str();
      } 

      
      const FermionField &getPhiOdd() const{ return PhiOdd; }

      virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) {
	      refresh(U,sRNG, pRNG,0.);
      }

      virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG& pRNG, RealD c1=0.) {
        // P(eta_o) = e^{- eta_o^dag eta_o}
        //
        // e^{x^2/2 sig^2} => sig^2 = 0.5.
        // 
	
	std::cout << GridLogMessage << "["<<action_name()<<"] -- No further parameters "<<std::endl;
        RealD scale = std::sqrt(0.5);

        FermionField eta    (NumOp.FermionGrid());
        gaussian(pRNG,eta); eta = eta * scale;

	refresh(U,eta,c1);
      }

      void refresh(const GaugeField &U, const FermionField &eta,RealD c1=0.) {

	Real c2=sqrt(1.-c1*c1);
	std::cout << GridLogMessage << " c1 "<<c1 <<" "<<c2<<std::endl;
        // P(phi) = e^{- phi^dag Vpc (MpcdagMpc)^-1 Vpcdag phi}
        //
        // NumOp == V
        // DenOp == M
        //
        FermionField etaOdd (NumOp.FermionRedBlackGrid());
        FermionField etaEven(NumOp.FermionRedBlackGrid());
        FermionField tmp    (NumOp.FermionRedBlackGrid());

        pickCheckerboard(Even,etaEven,eta);
        pickCheckerboard(Odd,etaOdd,eta);

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);
	std::cout << " TwoFlavourRefresh:  Imported gauge "<<std::endl;

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

	std::cout << " TwoFlavourRefresh: Diff ops "<<std::endl;

	int fnum=Grid::FieldNum();
	std::string fileO("./PhiOdd."+std::to_string(Grid::traj_num)+"_"+std::to_string(fnum) );
	std::string fileE("./PhiEven."+std::to_string(Grid::traj_num)+"_"+std::to_string(fnum) );
        std::ifstream fsO(fileO);
        std::ifstream fsE(fileE);
	GridBase* grid = PhiEven.Grid();

        if ( fsO.good() && fsE.good() ) {
          fsO.close();fsE.close();
	  // setting checkerboard
          pickCheckerboard(Even,PhiEven,eta);
          pickCheckerboard(Odd,PhiOdd,eta);
	} else {
	
	
          fsO.close();fsE.close();

        // Odd det factors
          Mpc.MpcDag(etaOdd,PhiOdd);
  	std::cout << " TwoFlavourRefresh: MpcDag "<<std::endl;
          tmp=Zero();
  	std::cout << " TwoFlavourRefresh: Zero() guess "<<std::endl;
          HeatbathSolver(Vpc,PhiOdd,tmp);
  	std::cout << " TwoFlavourRefresh: Heatbath solver "<<std::endl;
          Vpc.Mpc(tmp,PhiOdd);            
  	std::cout << " TwoFlavourRefresh: Mpc "<<std::endl;
  
  
          // Even det factors
          DenOp.MooeeDag(etaEven,tmp);
          NumOp.MooeeInvDag(tmp,PhiEven);
 
           uint32_t nersc_csum;
           uint32_t scidac_csuma;
           uint32_t scidac_csumb;
           typedef typename FermionField::vector_object  vobj;
           typedef typename FermionField::scalar_object  sobj;
           PFMunger<sobj,sobj> munge;
           std::string format = getFormatStringLocal<typename FermionField::vector_object>();
	   uint64_t offset; //leave 64 bits for header

	  if(c1 > 0.01){
	 uint32_t hdr_checksum, hdr_size;
         std::string format;
	std::string fileOP("./PhiOdd."+std::to_string(Grid::traj_num-1)+"_"+std::to_string(fnum) );
	std::string fileEP("./PhiEven."+std::to_string(Grid::traj_num-1)+"_"+std::to_string(fnum) );
        std::ifstream fsOP(fileOP);
        std::ifstream fsEP(fileEP);
         FermionField PhiOddP(PhiOdd.Grid());
         FermionField PhiEvenP(PhiEven.Grid());
         uint64_t offset = readHeader(hdr_size, hdr_checksum, format, fileOP);
         std::cout << "Data offset read " << offset << std::endl;
         std::cout << "Data size read " << hdr_size << std::endl;
         BinaryIO::readLatticeObject<vobj,sobj>(PhiOddP,fileOP,munge, offset, format,
                                                   nersc_csum,scidac_csuma,scidac_csumb);
         offset = readHeader(hdr_size, hdr_checksum, format, fileEP);
         std::cout << "Data offset read " << offset << std::endl;
         std::cout << "Data size read " << hdr_size << std::endl;
         BinaryIO::readLatticeObject<vobj,sobj>(PhiEvenP,fileEP,munge, offset, format,
                                                   nersc_csum,scidac_csuma,scidac_csumb);
          fsOP.close();fsEP.close();

	  PhiEven *=c2;
	  PhiOdd *=c2;
	  PhiEvenP *=c1;
	  PhiOddP *=c1;

	  PhiEven += PhiEvenP;
	  PhiOdd += PhiOddP;

	  }
  //         emptyUserRecord record;

           if ( grid->IsBoss() ) {
	     std::ofstream fout(fileO,std::ios::out);
             offset = writeHeader(1, 0, format, fileO);
           }
           grid->Broadcast(0,(void *)&offset,sizeof(offset)); //use as a barrier

  
           BinaryIO::writeLatticeObject<vobj,sobj>(PhiOdd,fileO,munge, offset, format,
                                                     nersc_csum,scidac_csuma,scidac_csumb);
           std::cout << GridLogMessage << " PhiOdd "<<fileO <<" checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::dec<< std::endl;

           if ( grid->IsBoss() ) {
	     std::ofstream fout(fileE,std::ios::out);
             offset = writeHeader(1, 0, format, fileE);
           }
           grid->Broadcast(0,(void *)&offset,sizeof(offset)); //use as a barrier
  
           BinaryIO::writeLatticeObject<vobj,sobj>(PhiEven,fileE,munge, offset, format,
                                                     nersc_csum,scidac_csuma,scidac_csumb);
           std::cout << GridLogMessage << " PhiEven "<<fileE <<" checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::dec<< std::endl;
	   std::cout << std::dec <<std::endl;


	}

	{
         uint32_t nersc_csum;
         uint32_t scidac_csuma;
         uint32_t scidac_csumb;
         typedef typename FermionField::vector_object  vobj;
         typedef typename FermionField::scalar_object  sobj;

         PFUnmunger<sobj,sobj> munge;

	 uint32_t hdr_checksum, hdr_size;
         std::string format;
         uint64_t offset = readHeader(hdr_size, hdr_checksum, format, fileO);
         std::cout << "Data offset read " << offset << std::endl;
         std::cout << "Data size read " << hdr_size << std::endl;
//         assert(data.size() == hdr_size);
//         std::string format = getFormatStringLocal<typename FermionField::vector_object>();

         BinaryIO::readLatticeObject<vobj,sobj>(PhiOdd,fileO,munge, offset, format,
                                                   nersc_csum,scidac_csuma,scidac_csumb);
         std::cout << GridLogMessage << " PhiOdd "<<fileO <<" checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::dec<< std::endl;

         offset = readHeader(hdr_size, hdr_checksum, format, fileE);
         BinaryIO::readLatticeObject<vobj,sobj>(PhiEven,fileE,munge, offset, format,
                                                   nersc_csum,scidac_csuma,scidac_csumb);
         std::cout << GridLogMessage << " PhiEven "<<fileE <<" checksums "<<std::hex << scidac_csuma << " "<<scidac_csumb<<std::dec<< std::endl;
	 std::cout << std::dec <<std::endl;


	}
	std::cout << " TwoFlavourRefresh: Mee "<<std::endl;

	RefreshAction = norm2(etaEven)+norm2(etaOdd);
	std::cout  << std::setprecision(14) << " refresh " <<action_name()<< " action "<<RefreshAction<<std::endl;
      };

      //////////////////////////////////////////////////////
      // S = phi^dag V (Mdag M)^-1 Vdag phi
      //////////////////////////////////////////////////////
      virtual RealD Sinitial(const GaugeField &U) {
	std::cout << GridLogMessage << "Returning stored two flavour refresh action "<<RefreshAction<<std::endl;
	return RefreshAction;
      }
      virtual RealD S(const GaugeField &U) {

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

        FermionField X(NumOp.FermionRedBlackGrid());
        FermionField Y(NumOp.FermionRedBlackGrid());

	std::cout  << std::setprecision(14) << " PhiOdd " <<norm2(PhiOdd)<<std::endl;
        Vpc.MpcDag(PhiOdd,Y);           // Y= Vdag phi
	std::cout  << std::setprecision(14) << " Y " <<norm2(Y)<<std::endl;
        X=Zero();
        ActionSolver(Mpc,Y,X);          // X= (MdagM)^-1 Vdag phi
	std::cout  << std::setprecision(14) << " X " <<norm2(X)<<std::endl;
        //Mpc.Mpc(X,Y);                   // Y=  Mdag^-1 Vdag phi
        // Multiply by Ydag
        RealD action = real(innerProduct(Y,X));
	std::cout  << std::setprecision(14) << " X Y " << action <<std::endl;

        //RealD action = norm2(Y);

        // The EE factorised block; normally can replace with zero if det is constant (gauge field indept)
        // Only really clover term that creates this. Leave the EE portion as a future to do to make most
        // rapid progresss on DWF for now.
        //
	std::cout  << std::setprecision(14) << " PhiEven " <<norm2(PhiEven)<<std::endl;
        NumOp.MooeeDag(PhiEven,X);
	std::cout  << std::setprecision(14) << " X " <<norm2(X)<<std::endl;
        DenOp.MooeeInvDag(X,Y);
	std::cout  << std::setprecision(14) << " Y " <<norm2(Y)<<std::endl;
        action = action + norm2(Y);

	std::cout  << std::setprecision(14) << " S " <<action_name()<< " action "<<action<<std::endl;
        return action;
      };

      //////////////////////////////////////////////////////
      // dS/du = phi^dag dV (Mdag M)^-1 V^dag  phi
      //       - phi^dag V (Mdag M)^-1 [ Mdag dM + dMdag M ]  (Mdag M)^-1 V^dag  phi
      //       + phi^dag V (Mdag M)^-1 dV^dag  phi
      //////////////////////////////////////////////////////
      virtual void deriv(const GaugeField &U,GaugeField & dSdU) {

        NumOp.ImportGauge(U);
        DenOp.ImportGauge(U);

        SchurDifferentiableOperator<Impl> Mpc(DenOp);
        SchurDifferentiableOperator<Impl> Vpc(NumOp);

        FermionField  X(NumOp.FermionRedBlackGrid());
        FermionField  Y(NumOp.FermionRedBlackGrid());

        // This assignment is necessary to be compliant with the HMC grids
	GaugeField force(dSdU.Grid());

        //Y=Vdag phi
        //X = (Mdag M)^-1 V^dag phi
        //Y = (Mdag)^-1 V^dag  phi
        Vpc.MpcDag(PhiOdd,Y);          // Y= Vdag phi
	std::cout << GridLogMessage <<" Y "<<norm2(Y)<<std::endl;
        X=Zero();
        DerivativeSolver(Mpc,Y,X);     // X= (MdagM)^-1 Vdag phi
	std::cout << GridLogMessage <<" X "<<norm2(X)<<std::endl;
        Mpc.Mpc(X,Y);                  // Y=  Mdag^-1 Vdag phi
	std::cout << GridLogMessage <<" Y "<<norm2(Y)<<std::endl;

        // phi^dag V (Mdag M)^-1 dV^dag  phi
        Vpc.MpcDagDeriv(force , X, PhiOdd );   dSdU = force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;
  
        // phi^dag dV (Mdag M)^-1 V^dag  phi
        Vpc.MpcDeriv(force , PhiOdd, X );      dSdU = dSdU+force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;

        //    -    phi^dag V (Mdag M)^-1 Mdag dM   (Mdag M)^-1 V^dag  phi
        //    -    phi^dag V (Mdag M)^-1 dMdag M   (Mdag M)^-1 V^dag  phi
        Mpc.MpcDeriv(force,Y,X);              dSdU = dSdU-force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;
        Mpc.MpcDagDeriv(force,X,Y);           dSdU = dSdU-force;
	std::cout << GridLogMessage <<" deriv "<<norm2(force)<<std::endl;

        // FIXME No force contribution from EvenEven assumed here
        // Needs a fix for clover.
        assert(NumOp.ConstEE() == 1);
        assert(DenOp.ConstEE() == 1);

        dSdU = -dSdU;
        
      };
    };
NAMESPACE_END(Grid);
#endif
