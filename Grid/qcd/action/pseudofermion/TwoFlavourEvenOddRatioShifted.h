    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/action/pseudofermion/TwoFlavourEvenOddRatioShifted.h

    Copyright (C) 2015

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
#ifndef QCD_PSEUDOFERMION_ONE_FLAVOUR_EVEN_ODD_RATIONAL_RATIO_H
#define QCD_PSEUDOFERMION_ONE_FLAVOUR_EVEN_ODD_RATIONAL_RATIO_H

NAMESPACE_BEGIN(Grid);

    ///////////////////////////////////////
    // One flavour rational
    ///////////////////////////////////////

    // S_f = chi^dag* P(V^dag*V)/Q(V^dag*V)* N(M^dag*M)/D(M^dag*M)* P(V^dag*V)/Q(V^dag*V)* chi       
    //
    // Here P/Q \sim R_{1/4}  ~ (V^dagV)^{1/4}  
    // Here N/D \sim R_{-1/2} ~ (M^dagM)^{-1/2}  


    template<class Impl>
    class TwoFlavourEvenOddRatioShiftedPseudoFermionAction : public GeneralEvenOddRatioRationalPseudoFermionAction<Impl> {
    public:
      typedef RationalActionParams Params;

    public:
      TwoFlavourEvenOddRatioShiftedPseudoFermionAction(FermionOperator<Impl>  &_NumOp,
                                                        FermionOperator<Impl>  &_DenOp,
                                                        const Params & p
                                                        ) :
        GeneralEvenOddRatioRationalPseudoFermionAction<Impl>(_NumOp, _DenOp, p){}

      virtual std::string action_name(){return "OneFlavourEvenOddRatioRationalPseudoFermionAction";}
    };

    template<class Impl,class ImplF>
    class TwoFlavourEvenOddRatioShiftedMixedPrecPseudoFermionAction
      : public GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction<Impl,ImplF> {
    public:
//      typedef OneFlavourRationalParams Params;
      typedef RationalActionParams Params;


    public:
      TwoFlavourEvenOddRatioShiftedMixedPrecPseudoFermionAction (FermionOperator<Impl>  &_NumOp,
                                                                 FermionOperator<Impl>  &_DenOp,
                                                                 FermionOperator<ImplF>  &_NumOpF,
                                                                 FermionOperator<ImplF>  &_DenOpF,
                                                                 const Params & p, Integer ReliableUpdateFreq
                                                        ) :
        GeneralEvenOddRatioRationalMixedPrecPseudoFermionAction<Impl,ImplF>(_NumOp, _DenOp,_NumOpF, _DenOpF, p,ReliableUpdateFreq){}

      virtual std::string action_name(){return "OneFlavourEvenOddRatioRationalPseudoFermionAction";}
    };

NAMESPACE_END(Grid);

#endif


