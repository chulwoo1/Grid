    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourEvenOddRatioShifted.h

    Copyright (C) 2015

    Author: Chulwoo Jung (chulwoo@bnl.gov)

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
#ifndef QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD__RATIO_SHIFTED_H
#define QCD_PSEUDOFERMION_TWO_FLAVOUR_EVEN_ODD__RATIO_SHIFTED_H

NAMESPACE_BEGIN(Grid);

    ///////////////////////////////////////
    // One flavour rational
    ///////////////////////////////////////

    // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
    //
    // Here P/Q \sim R_{1/4}  ~ (V^dagV)^{1/4}  
    // Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}  
  
    template<class Impl>
    class TwoFlavourEvenOddRatioShiftedPseudoFermionAction : public GeneralEvenOddRatioRationalPseudoFermionActionExt<Impl> {
    public:
      typedef RationalActionParams Params;
#if 0
    private:
      static RationalActionParams transcribe(const Params &in){
	RationalActionParams out;
	out.inv_pow = 2;
	out.lo = in.lo;
	out.hi = in.hi;
	out.MaxIter = in.MaxIter;
	out.action_tolerance = out.md_tolerance = in.tolerance;
	out.action_degree = out.md_degree = in.degree;
	out.precision = in.precision;
	out.BoundsCheckFreq = in.BoundsCheckFreq;
	return out;
      }
#endif

    public:
      TwoFlavourEvenOddRatioShiftedPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
							FermionOperator<Impl>  &_DenOp, 
							const Params & p
							) : 
	GeneralEvenOddRatioRationalPseudoFermionActionExt<Impl>(_NumOp, _DenOp, p){}

      virtual std::string action_name(){return "TwoFlavourEvenOddRatioShiftedPseudoFermionAction";}      
    };

    template<class Impl,class ImplF>
    class TwoFlavourEvenOddRatioShiftedMixedPrecPseudoFermionAction
      : public GeneralEvenOddRatioRationalMixedPrecPseudoFermionActionExt<Impl,ImplF> {
    public:
      typedef RationalActionParams Params;
#if 0
    private:
      static RationalActionParams transcribe(const Params &in){
	RationalActionParams out;
	out.inv_pow = 2;
	out.lo = in.lo;
	out.hi = in.hi;
	out.MaxIter = in.MaxIter;
	out.action_tolerance = out.md_tolerance = in.tolerance;
	out.action_degree = out.md_degree = in.degree;
	out.precision = in.precision;
	out.BoundsCheckFreq = in.BoundsCheckFreq;
	return out;
      }
#endif

    public:
      TwoFlavourEvenOddRatioShiftedMixedPrecPseudoFermionAction(FermionOperator<Impl>  &_NumOp, 
								 FermionOperator<Impl>  &_DenOp, 
								 FermionOperator<ImplF>  &_NumOpF, 
								 FermionOperator<ImplF>  &_DenOpF, 
								 const Params & p, Integer ReliableUpdateFreq
							) : 
	GeneralEvenOddRatioRationalMixedPrecPseudoFermionActionExt<Impl,ImplF>(_NumOp, _DenOp,_NumOpF, _DenOpF, p,ReliableUpdateFreq){
		std::cout << GridLogMessage << action_name() <<std::endl;
	}

      virtual std::string action_name(){return "TwoFlavourEvenOddRatioShiftedMixedPrecPseudoFermionAction";}      
    };

NAMESPACE_END(Grid);

#endif
