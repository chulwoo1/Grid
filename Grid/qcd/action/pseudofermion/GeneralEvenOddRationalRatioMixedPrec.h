    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/GeneralEvenOddRationalRatioMixedPrec.h

    Copyright (C) 2015

    Author: Christopher Kelly <ckelly@bnl.gov>
    Author: Peter Boyle <paboyle@ph.ed.ac.uk>

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
#ifndef QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_MIXED_PREC_H
#define QCD_PSEUDOFERMION_GENERAL_EVEN_ODD_RATIONAL_RATIO_MIXED_PREC_H

#include <Grid/algorithms/iterative/ConjugateGradientMultiShiftCleanup.h>

NAMESPACE_BEGIN(Grid);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generic rational approximation for ratios of operators utilizing the mixed precision multishift algorithm
    // cf. GeneralEvenOddRational.h for details
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
    template<class ImplD, class ImplF>
    class GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction : public GeneralEvenOddRatioRationalPseudoFermionAction<ImplD> {
    private:
      typedef typename ImplD::FermionField FermionFieldD;
      typedef typename ImplF::FermionField FermionFieldF;

      FermionOperator<ImplD> & NumOpD;
      FermionOperator<ImplD> & DenOpD;

      FermionOperator<ImplF> & NumOpF;
      FermionOperator<ImplF> & DenOpF;

      Integer ReliableUpdateFreq;
    protected:

      //Action evaluation
      //Allow derived classes to override the multishift CG
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, FermionFieldD &out){
	std::cout << "NoEXT multiShiftInverse (b,MSF,i,V,V)"<<std::endl;
         
#if 1
	SchurDifferentiableOperator<ImplD> schurOp(numerator ? NumOpD : DenOpD);
	ConjugateGradientMultiShift<FermionFieldD> msCG(MaxIter, approx);
	msCG(schurOp,in, out);
#else
	SchurDifferentiableOperator<ImplD> schurOpD(numerator ? NumOpD : DenOpD);
	SchurDifferentiableOperator<ImplF> schurOpF(numerator ? NumOpF : DenOpF);
	FermionFieldD inD(NumOpD.FermionRedBlackGrid());
	FermionFieldD outD(NumOpD.FermionRedBlackGrid());

	// Action better with higher precision?
	ConjugateGradientMultiShiftMixedPrec<FermionFieldD, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	msCG(schurOpD, in, out);
#endif
	std::cout << "NoEXT multiShiftInverse (b,MSF,i,V,V) done"<<std::endl;
      }
      //Force evaluation
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, std::vector<FermionFieldD> &out_elems, FermionFieldD &out){
	std::cout << "NoEXT multiShiftInverse (b,MSF,i,V,V*,V) order="<<approx.order<<" poly= "<<approx.poly.size()<<std::endl;
	SchurDifferentiableOperator<ImplD> schurOpD(numerator ? NumOpD : DenOpD);
	SchurDifferentiableOperator<ImplF>  schurOpF(numerator ? NumOpF  : DenOpF);

	FermionFieldD inD(NumOpD.FermionRedBlackGrid());
	FermionFieldD outD(NumOpD.FermionRedBlackGrid());
//	std::vector<FermionFieldD> out_elemsD(out_elems.size(),NumOpD.FermionRedBlackGrid());
	ConjugateGradientMultiShiftMixedPrecCleanup<FermionFieldD, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	out.Checkerboard()=in.Checkerboard();
	msCG(schurOpD, in, out_elems, out);
	std::cout << "NoEXT multiShiftInverse (b,MSF,i,V,V*,V) done "<<std::endl;
      }
      //Allow derived classes to override the gauge import
      virtual void ImportGauge(const typename ImplD::GaugeField &Ud){

	typename ImplF::GaugeField Uf(NumOpF.GaugeGrid());
	precisionChange(Uf, Ud);

	std::cout << "Importing "<<norm2(Ud)<<" "<< norm2(Uf)<<" " <<std::endl;
	
	NumOpD.ImportGauge(Ud);
	DenOpD.ImportGauge(Ud);

	NumOpF.ImportGauge(Uf);
	DenOpF.ImportGauge(Uf);
      }
      
    public:
      GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction(FermionOperator<ImplD>  &_NumOpD, FermionOperator<ImplD>  &_DenOpD, 
							      FermionOperator<ImplF>  &_NumOpF, FermionOperator<ImplF>  &_DenOpF, 
							      const RationalActionParams & p, Integer _ReliableUpdateFreq
							      ) : GeneralEvenOddRatioRationalPseudoFermionAction<ImplD>(_NumOpD, _DenOpD, p),
								  ReliableUpdateFreq(_ReliableUpdateFreq),
								  NumOpD(_NumOpD), DenOpD(_DenOpD),
								  NumOpF(_NumOpF), DenOpF(_DenOpF)
      {}
      
      virtual std::string action_name(){return "GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction";}
    };

    template<class ImplD, class ImplF>
    class GeneralEvenOddRatioRationalMixedPrecPseudoFermionActionExt : public GeneralEvenOddRatioRationalPseudoFermionActionExt<ImplD> {
    private:
      typedef typename ImplD::FermionField FermionFieldD;
      typedef typename ImplF::FermionField FermionFieldF;

      FermionOperator<ImplD> & NumOpD;
      FermionOperator<ImplD> & DenOpD;

      FermionOperator<ImplF> & NumOpF;
      FermionOperator<ImplF> & DenOpF;

      Integer ReliableUpdateFreq;
    protected:

      //Action evaluation
      //Allow derived classes to override the multishift CG
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, FermionFieldD &out){
	std::cout << "EXT multiShiftInverse (b,MSF,i,V,V) order="<<approx.order <<
		" poly= "<<approx.poly.size()<<std::endl;
	out.Checkerboard()== in.Checkerboard();
#if 1
	SchurDifferentiableOperator<ImplD> schurOp(numerator ? NumOpD : DenOpD);
	ConjugateGradientMultiShift<FermionFieldD> msCG(MaxIter, approx);
	msCG(schurOp,in, out);
#else
	SchurDifferentiableOperator<ImplD> schurOpD(numerator ? NumOpD : DenOpD);
	SchurDifferentiableOperator<ImplF> schurOpF(numerator ? NumOpF : DenOpF);
	FermionFieldD inD(NumOpD.FermionRedBlackGrid());
	FermionFieldD outD(NumOpD.FermionRedBlackGrid());

	// Action better with higher precision?
	ConjugateGradientMultiShiftMixedPrec<FermionFieldD, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	msCG(schurOpD, in, out);
#endif
	std::cout << "EXT multiShiftInverse (b,MSF,i,V,V) done"<<std::endl;
      }
      //Force evaluation
      virtual void multiShiftInverse(bool numerator, const MultiShiftFunction &approx, const Integer MaxIter, const FermionFieldD &in, std::vector<FermionFieldD> &out_elems, FermionFieldD &out){

	std::cout << "EXT multiShiftInverse (b,MSF,i,V,V*,V) order="<<approx.order <<
		" out_elems.size()= "<< out_elems.size() <<
		" poly= "<<approx.poly.size()<<std::endl;
	out.Checkerboard()== in.Checkerboard();
	for(int i=0;i<out_elems.size();i++)
	    out_elems[i].Checkerboard() == in.Checkerboard();

	SchurDifferentiableOperator<ImplD> schurOpD(numerator ? NumOpD : DenOpD);
	SchurDifferentiableOperator<ImplF>  schurOpF(numerator ? NumOpF  : DenOpF);

//	FermionFieldD inD(NumOpD.FermionRedBlackGrid());
//	FermionFieldD outD(NumOpD.FermionRedBlackGrid());
//	std::vector<FermionFieldD> out_elemsD(out_elems.size(),NumOpD.FermionRedBlackGrid());
	ConjugateGradientMultiShiftMixedPrecCleanup<FermionFieldD, FermionFieldF> msCG(MaxIter, approx, NumOpF.FermionRedBlackGrid(), schurOpF, ReliableUpdateFreq);
	std::cout << "msCG "<<std::endl;
	msCG(schurOpD, in, out_elems, out);
	std::cout << "msCG done "<<std::endl;
	std::cout << "EXT multiShiftInverse (b,MSF,i,V,V*,V) done"<<std::endl;
      }
      //Allow derived classes to override the gauge import
      virtual void ImportGauge(const typename ImplD::GaugeField &Ud){

	typename ImplF::GaugeField Uf(NumOpF.GaugeGrid());
	precisionChange(Uf, Ud);

	std::cout << "Importing "<<norm2(Ud)<<" "<< norm2(Uf)<<" " <<std::endl;
	
	NumOpD.ImportGauge(Ud);
	DenOpD.ImportGauge(Ud);

	NumOpF.ImportGauge(Uf);
	DenOpF.ImportGauge(Uf);
      }
      
    public:
      GeneralEvenOddRatioRationalMixedPrecPseudoFermionActionExt(FermionOperator<ImplD>  &_NumOpD, FermionOperator<ImplD>  &_DenOpD, 
							      FermionOperator<ImplF>  &_NumOpF, FermionOperator<ImplF>  &_DenOpF, 
							      const RationalActionParams & p, Integer _ReliableUpdateFreq
							      ) : GeneralEvenOddRatioRationalPseudoFermionActionExt<ImplD>(_NumOpD, _DenOpD, p),
								  ReliableUpdateFreq(_ReliableUpdateFreq),
								  NumOpD(_NumOpD), DenOpD(_DenOpD),
								  NumOpF(_NumOpF), DenOpF(_DenOpF)
      {}
      
      virtual std::string action_name(){return "GeneralEvenOddRatioRationalMixedPrecPseudoFermionActionExt";}
    };

NAMESPACE_END(Grid);

#endif
