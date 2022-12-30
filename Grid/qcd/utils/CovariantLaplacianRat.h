/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/scalar/CovariantLaplacianRat.h

Copyright (C) 2016

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
#pragma once 

NAMESPACE_BEGIN(Grid);

struct LaplacianRatParams {

  RealD offset;
  int order;
  std::vector<RealD> a0;
  std::vector<RealD> a1;
  std::vector<RealD> b0;
  std::vector<RealD> b1;
  RealD b2; //for debugging
  int   MaxIter;
  RealD tolerance;
  int   precision;
  
  // constructor 
  LaplacianRatParams(int ord = 1,
                  int maxit     = 1000,
                  RealD tol     = 1.0e-8, 
                  int precision = 64)
    : offset(1.), order(ord),b2(1.),
      MaxIter(maxit),
      tolerance(tol),
      precision(precision){ 
      a0.resize(ord,0.);
      a1.resize(ord,0.);
      b0.resize(ord,0.);
      b1.resize(ord,0.);
      };
};



////////////////////////////////////////////////////////////
// Laplacian operator L on adjoint fields
//
// phi: adjoint field
// L: D_mu^dag D_mu
//
// L phi(x) = Sum_mu [ U_mu(x)phi(x+mu)U_mu(x)^dag + 
//                     U_mu(x-mu)^dag phi(x-mu)U_mu(x-mu)
//                     -2phi(x)]
//
// Operator designed to be encapsulated by
// an HermitianLinearOperator<.. , ..>
////////////////////////////////////////////////////////////

template <class Impl, class ImplF>
class LaplacianAdjointRat: public Metric<typename Impl::Field> {
  OperatorFunction<typename Impl::Field> &Solver;
  LaplacianRatParams Gparam;
  LaplacianRatParams Mparam;
  GridBase *grid;
  GridBase *grid_f;
public:
  INHERIT_GIMPL_TYPES(Impl);
  typedef typename ImplF::Field GaugeFieldF;
  GaugeField Usav;
  GaugeFieldF UsavF;

	  LaplacianAdjointRat(GridBase* _grid, GridBase* _grid_f, OperatorFunction<GaugeField>& S, LaplacianRatParams& gpar, LaplacianRatParams& mpar)
    : grid(_grid),grid_f(_grid_f), U(Nd, _grid), Solver(S), Gparam(gpar), Mparam(mpar),Usav(_grid), UsavF(_grid_f)  {
//    std::cout<<GridLogMessage << "Generating degree "<<param.degree<<" for x^(1/2)"<<std::endl;
    this->triv=0;
        

  };
  LaplacianAdjointRat(){this->triv=0; printf("triv=%d\n",this->Trivial());}
  void Mdir(const GaugeField&, GaugeField&, int, int){ assert(0);}
  void MdirAll(const GaugeField&, std::vector<GaugeField> &){ assert(0);}
  void Mdiag(const GaugeField&, GaugeField&){ assert(0);}

  void ImportGauge(const GaugeField& _U) {
    RealD total=0.;
    for (int mu = 0; mu < Nd; mu++) {
      U[mu] = PeekIndex<LorentzIndex>(_U, mu);
      total += norm2(U[mu]);
    }
    Usav = _U;
    precisionChange(UsavF,Usav);
    std::cout <<GridLogDebug << "ImportGauge:norm2(_U) = "<<" "<<total<<std::endl;
  }

#if 0
  void Lap(const GaugeField& in, GaugeField& out) {
    // in is an antihermitian matrix
    // test
    //GaugeField herm = in + adj(in);
    //std::cout << "AHermiticity: " << norm2(herm) << std::endl;

   
    GaugeLinkField tmp(in.Grid());
    GaugeLinkField tmp2(in.Grid());
    GaugeLinkField sum(in.Grid());

    RealD kappa=1.;
    for (int nu = 0; nu < Nd; nu++) {
      sum = Zero();
      GaugeLinkField in_nu = PeekIndex<LorentzIndex>(in, nu);
      GaugeLinkField out_nu(out.Grid());
      for (int mu = 0; mu < Nd; mu++) {
        tmp = U[mu] * Cshift(in_nu, mu, +1) * adj(U[mu]);
        tmp2 = adj(U[mu]) * in_nu * U[mu];
        sum += tmp + Cshift(tmp2, mu, -1) - 2.0 * in_nu;
      }
      out_nu = (1.0 - kappa) * in_nu - kappa / (double(4 * Nd)) * sum;
      PokeIndex<LorentzIndex>(out, out_nu, nu);
    }
  }

  // separating this temporarily
  void LapDeriv(const GaugeField& left, const GaugeField& right,
              GaugeField& der) {
    // in is anti-hermitian
    RealD kappa=1.;
    RealD factor = -kappa / (double(4 * Nd));

    for (int mu = 0; mu < Nd; mu++) {
      GaugeLinkField der_mu(der.Grid());
      der_mu = Zero();
      for (int nu = 0; nu < Nd; nu++) {
        GaugeLinkField left_nu = PeekIndex<LorentzIndex>(left, nu);
        GaugeLinkField right_nu = PeekIndex<LorentzIndex>(right, nu);
        der_mu += U[mu] * Cshift(left_nu, mu, 1) * adj(U[mu]) * right_nu;
        der_mu += U[mu] * Cshift(right_nu, mu, 1) * adj(U[mu]) * left_nu;
      }
      PokeIndex<LorentzIndex>(der, -factor * der_mu, mu);
    }
    std::cout <<GridLogDebug << "MDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }
#endif

  // separating this temporarily
  void MDerivInt(LaplacianRatParams &par, const GaugeField& left, const GaugeField& right,
              GaugeField& der) {
    GaugeField LMinvMom(left.Grid());

    GaugeField GMom(left.Grid());
    GaugeField MinvGMom(left.Grid());
    GaugeField LMinvGMom(left.Grid());

    GaugeField AGMom(left.Grid());
    GaugeField MinvAGMom(left.Grid());
    GaugeField LMinvAGMom(left.Grid());

    GaugeField AMinvMom(left.Grid());
    GaugeField LMinvAMom(left.Grid());
//    GaugeField MinvAMom(left.Grid());
    GaugeField temp(left.Grid());
    GaugeField temp2(left.Grid());
    std::vector<GaugeField> MinvMom(par.order,left.Grid());


    ConjugateGradient<GaugeField> CG(1.0e-8,10000,false);
    ConjugateGradient<GaugeFieldF> CG_f(1.0e-8,10000,false);
    LaplacianParams LapPar(0.0001, 1.0, 10000, 1e-8, 12, 64);
    LaplacianAdjointField<Impl> Laplacian(left.Grid(), CG, LapPar, 1.,false);
    LaplacianAdjointField<ImplF> LaplacianF(grid_f, CG_f, LapPar, 1.,false);
    Laplacian.ImportGauge(Usav);
    LaplacianF.ImportGauge(UsavF);
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(Laplacian);
//    HermitianLinearOperator<LaplacianAdjointField<ImplF>,GaugeFieldF> HermOpF(Laplacian);
    

    GMom = par.offset * right;
    for(int i =0;i<par.order;i++){
    QuadLinearOperator<LaplacianAdjointField<Impl>,GaugeField> QuadOp(Laplacian,par.b0[i],par.b1[i],par.b2);
    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeFieldF> QuadOpF(LaplacianF,par.b0[i],par.b1[i],par.b2);
//    MixedPrecisionConjugateGradient<GaugeField,GaugeFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
//    MixedCG.InnerTolerance=par.tolerance;
    GaugeField Gtemp2(left.Grid());
//    MixedCG(right,MinvMom[i]);
    CG(QuadOp,right,MinvMom[i]);
    
    GMom += par.a0[i]*MinvMom[i]; 
    HermOp.HermOp(MinvMom[i],Gtemp2);
    GMom += par.a1[i]*Gtemp2; 
    }
    for(int i =0;i<par.order;i++){
    QuadLinearOperator<LaplacianAdjointField<Impl>,GaugeField> QuadOp(Laplacian,par.b0[i],par.b1[i],par.b2);
    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeFieldF> QuadOpF(LaplacianF,par.b0[i],par.b1[i],par.b2);
//    MixedPrecisionConjugateGradient<GaugeField,GaugeFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
//    MixedCG.InnerTolerance=par.tolerance;
    GaugeField Gtemp(left.Grid());
    GaugeField Gtemp2(left.Grid());

//    Solver(QuadOp,GMom,MinvGMom);
//    MixedCG(GMom,MinvGMom);
    CG(QuadOp,GMom,MinvGMom);
    Laplacian.M(MinvGMom, LMinvGMom);
//    MixedCG(right,MinvMom[i]);
    CG(QuadOp,right,MinvMom[i]);

    Laplacian.M(MinvMom[i], LMinvMom);
    Laplacian.M(MinvMom[i], AMinvMom);
    AMinvMom = par.a1[i]*LMinvMom;
    AMinvMom += par.a0[i]*MinvMom[i];

    Laplacian.M(AMinvMom,LMinvAMom);
    Laplacian.M(MinvGMom,temp);
    MinvAGMom = par.a1[i]*temp;
    MinvAGMom += par.a0[i]*MinvGMom;
    Laplacian.M(MinvAGMom,LMinvAGMom);


    RealD coef=0.5;
//    RealD coef=1;
    std::cout<<GridLogMessage << "coef =  "<< coef <<std::endl;
    Laplacian.MDeriv(GMom,MinvMom[i],temp); der += coef*2*par.a1[i]*temp;
    Laplacian.MDeriv(left,MinvGMom,temp); der += coef*2*par.a1[i]*temp;
    Laplacian.MDeriv(LMinvAGMom,MinvMom[i],temp); der += coef*-2.*par.b2*temp;
    Laplacian.MDeriv(LMinvAMom,MinvGMom,temp); der += coef*-2.*par.b2*temp;
    Laplacian.MDeriv(MinvAGMom,LMinvMom,temp); der += coef*-2.*par.b2*temp;
    Laplacian.MDeriv(AMinvMom,LMinvGMom,temp); der += coef*-2.*par.b2*temp;
    Laplacian.MDeriv(MinvAGMom,MinvMom[i],temp); der += coef*-2.*par.b1[i]*temp;
    Laplacian.MDeriv(AMinvMom,MinvGMom,temp); der += coef*-2.*par.b1[i]*temp;

    }
  }

  void MDeriv(const GaugeField& in, GaugeField& der) {
    MDeriv(in,in, der);
  }

  void MDeriv(const GaugeField& left, const GaugeField& right,
              GaugeField& der) {
    der=Zero();
    MDerivInt(Mparam, left, right, der);
    std::cout <<GridLogDebug << "MDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }

  void MinvDeriv(const GaugeField& in, GaugeField& der) {
    der=Zero();
    MDerivInt(Gparam, in, in, der);
    std::cout <<GridLogDebug << "MinvDeriv:norm2(der) = "<<norm2(der)<<std::endl;
  }


  void MSquareRootInt(LaplacianRatParams &par, GaugeField& P){
    GaugeField Gp(P.Grid());
//    GaugeField Gp_f(grid_f);
    Gp = par.offset * P;
    ConjugateGradient<GaugeField> CG(par.tolerance,10000);
    ConjugateGradient<GaugeFieldF> CG_f(1.0e-8,10000);
    LaplacianParams LapPar(0.0001, 1.0, 10000, 1e-8, 12, 64);
    LaplacianAdjointField<Impl> Laplacian(P.Grid(), CG, LapPar, 1.,false);
    LaplacianAdjointField<ImplF> LaplacianF(grid_f, CG_f, LapPar, 1.,false);
    Laplacian.ImportGauge(Usav);
    LaplacianF.ImportGauge(UsavF);
    HermitianLinearOperator<LaplacianAdjointField<Impl>,GaugeField> HermOp(Laplacian);

    for(int i =0;i<par.order;i++){
    QuadLinearOperator<LaplacianAdjointField<Impl>,GaugeField> QuadOp(Laplacian,par.b0[i],par.b1[i],par.b2);
    QuadLinearOperator<LaplacianAdjointField<ImplF>,GaugeFieldF> QuadOpF(LaplacianF,par.b0[i],par.b1[i],par.b2);
//    MixedPrecisionConjugateGradient<GaugeField,GaugeFieldF> MixedCG(par.tolerance,10000,10000,grid_f,QuadOpF,QuadOp);
//    MixedCG.InnerTolerance=par.tolerance;
    GaugeField Gtemp(P.Grid());
    GaugeField Gtemp2(P.Grid());
//    MixedCG(P,Gtemp);
    CG(QuadOp,P,Gtemp);
    Gp += par.a0[i]*Gtemp; 
    HermOp.HermOp(Gtemp,Gtemp2);
    Gp += par.a1[i]*Gtemp2; 
    }
    P = Gp;
  }

  void MSquareRoot(GaugeField& P){
    MSquareRootInt(Mparam,P);
    std::cout <<GridLogDebug << "MSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }

  void MInvSquareRoot(GaugeField& P){
    MSquareRootInt(Gparam,P);
    std::cout <<GridLogDebug << "MInvSquareRoot:norm2(P) = "<<norm2(P)<<std::endl;
  }

  void M(const GaugeField& in, GaugeField& out) {
      out = in;
      MSquareRoot(out);
      MSquareRoot(out);
      std::cout <<GridLogDebug << "M:norm2(out) = "<<norm2(out)<<std::endl;
  }

  void Minv(const GaugeField& in, GaugeField& inverted){
//    HermitianLinearOperator<LaplacianAdjointRat<Impl>,GaugeField> HermOp(*this);
//    Solver(HermOp, in, inverted);
      inverted = in;
      MInvSquareRoot(inverted);
      MInvSquareRoot(inverted);
      std::cout <<GridLogDebug << "Minv:norm2(inverted) = "<<norm2(inverted)<<std::endl;
  }



private:
//  RealD kappa;
  std::vector<GaugeLinkField> U;
};

NAMESPACE_END(Grid);
