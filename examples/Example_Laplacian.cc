#include <Grid/Grid.h>
using namespace Grid;

/*
/////////////////////////////////////////////////////////////////////////////////////////////
// Grid/algorithms/SparseMatrix.h: Interface defining what I expect of a general sparse matrix, such as a Fermion action
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Field> class SparseMatrixBase {
public:
  virtual GridBase *Grid(void) =0;

  virtual void  M    (const Field &in, Field &out)=0;
  virtual void  Mdag (const Field &in, Field &out)=0;
  virtual void  MdagM(const Field &in, Field &out) {
    Field tmp (in.Grid());
    M(in,tmp);
    Mdag(tmp,out);
  }
  virtual  void Mdiag    (const Field &in, Field &out)=0;
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp)=0;
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)=0;
};
*/

const std::vector<int> directions   ({Xdir,Ydir,Zdir,Xdir,Ydir,Zdir});
const std::vector<int> displacements({1,1,1,-1,-1,-1});
const std::vector<int> directions4D   ({Xdir,Ydir,Zdir,Tdir,Xdir,Ydir,Zdir,Tdir});
const std::vector<int> displacements4D({1,1,1,1,-1,-1,-1,-1});

template<class Field, int Dim > class FreeLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  GridBase *grid;
  FreeLaplacianCshift(GridBase *_grid)
  {
    grid=_grid;
  };
  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &in, Field &out)
  {
    out = Zero();
    for(int mu=0;mu<Dim;mu++) {
      out = out + Cshift(in,mu,1) + Cshift(in,mu,-1) - 2.0 * in;
    }
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

template<class Gimpl,class Field> class CovariantLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  GridBase *grid;
  GaugeField U;
  
  CovariantLaplacianCshift(GaugeField &_U)    :
    grid(_U.Grid()),
    U(_U) {  };

  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &in, Field &out)
  {
    out=Zero();
    for(int mu=0;mu<Nd-1;mu++) {
      GaugeLinkField Umu = PeekIndex<LorentzIndex>(U, mu); // NB: Inefficent
      out = out + Gimpl::CovShiftForward(Umu,mu,in);    
      out = out + Gimpl::CovShiftBackward(Umu,mu,in);    
      out = out - 2.0*in;
    }
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

template<class Gimpl,class Field> class CovariantAdjointLaplacianCshift : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  GridBase *grid;
  GaugeField U;
//  RealD kappa;
  
  CovariantAdjointLaplacianCshift(GaugeField &_U )    :
    grid(_U.Grid()),
    U(_U) {  };

  virtual GridBase *Grid(void) { return grid; };

  void M(const Field& in, Field& out) {

    Field tmp(in.Grid());
    Field tmp2(in.Grid());
    Field tmp3(in.Grid());
    Field sum(in.Grid());

      sum = Zero();
      for (int mu = 0; mu < Nd; mu++) {
        GaugeLinkField Umu = PeekIndex<LorentzIndex>(U, mu); // NB: Inefficent
        GaugeLinkField Udag = adj(Umu);
//        tmp = Umu * Cshift(in, mu, +1) * adj(Umu);
        tmp3 = Umu* Cshift(in, mu, +1);
        tmp = adj(tmp3);
        tmp3 = Umu*tmp;
	tmp = adj(tmp3);
//        tmp2 = adj(Umu) * in* Umu;
        tmp3 = Udag *in ;
	tmp2 = adj(tmp3);
	tmp3 = Udag*tmp2;
	tmp2 = adj(tmp3);
        sum = sum + tmp + Cshift(tmp2, mu, -1) - 2.0 * in;
      }
//      out= (1.0 - kappa) * in+ kappa / (double(4 * Nd)) * sum;
      out= sum;
  }
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};


#define LEG_LOAD(Dir)						 \
  SE = st.GetEntry(ptype, Dir, ss);				 \
  if (SE->_is_local ) {						 \
    int perm= SE->_permute;					 \
    chi = coalescedReadPermute(in[SE->_offset],ptype,perm,lane); \
  } else {							 \
    chi = coalescedRead(buf[SE->_offset],lane);			 \
  }								 \
  acceleratorSynchronise();

template<class Field> class FreeLaplacianStencil : public SparseMatrixBase<Field>
{
public:
  typedef typename Field::vector_object siteObject;
  typedef CartesianStencil<siteObject, siteObject, int> StencilImpl;

  GridBase *grid;
  StencilImpl Stencil;
  SimpleCompressor<siteObject> Compressor;
  
  FreeLaplacianStencil(GridBase *_grid)
    : Stencil    (_grid,6,Even,directions,displacements,0), grid(_grid)
  {  };
  
  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &_in, Field &_out)
  {

    ///////////////////////////////////////////////
    // Halo exchange for this geometry of stencil
    ///////////////////////////////////////////////
    Stencil.HaloExchange(_in, Compressor);

    ///////////////////////////////////
    // Arithmetic expressions
    ///////////////////////////////////

    // Views; device friendly/accessible pointers
    auto st = Stencil.View(AcceleratorRead);
    auto buf = st.CommBuf();
    autoView( in     , _in    , AcceleratorRead);
    autoView( out    , _out   , AcceleratorWrite);

    typedef typename Field::vector_object        vobj;
    typedef decltype(coalescedRead(in[0])) calcObj;

    const int      Nsimd = vobj::Nsimd();
    const uint64_t NN = grid->oSites();

    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	int ptype;

	res                 = coalescedRead(in[ss])*(-6.0);
	LEG_LOAD(0);	res = res + chi;
	LEG_LOAD(1);	res = res + chi;
	LEG_LOAD(2);	res = res + chi;
	LEG_LOAD(3);	res = res + chi;
	LEG_LOAD(4);	res = res + chi;
	LEG_LOAD(5);	res = res + chi;

	coalescedWrite(out[ss], res,lane);
	
    });
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

template<class Gimpl,class Field> class CovariantLaplacianStencil : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);

  typedef typename Field::vector_object siteObject;

  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Nc> >, Nds>;
  typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;

  typedef CartesianStencil<siteObject, siteObject, int> StencilImpl;

  GridBase *grid;
  StencilImpl Stencil;
  SimpleCompressor<siteObject> Compressor;
  DoubledGaugeField Uds;
  CovariantLaplacianStencil(GaugeField &Umu)
    :
      grid(Umu.Grid()),
      Stencil    (grid,6,Even,directions,displacements,0),
      Uds(grid)
  {
    for (int mu = 0; mu < Nd; mu++) {
      auto U = PeekIndex<LorentzIndex>(Umu, mu);
      PokeIndex<LorentzIndex>(Uds, U, mu );
      U = adj(Cshift(U, mu, -1));
      PokeIndex<LorentzIndex>(Uds, U, mu + 4);
    }
  };
  
  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &_in, Field &_out)
  {
    ///////////////////////////////////////////////
    // Halo exchange for this geometry of stencil
    ///////////////////////////////////////////////
    Stencil.HaloExchange(_in, Compressor);

    ///////////////////////////////////
    // Arithmetic expressions
    ///////////////////////////////////
    auto st = Stencil.View(AcceleratorRead);
    auto buf = st.CommBuf();

    autoView( in     , _in    , AcceleratorRead);
    autoView( out    , _out   , AcceleratorWrite);
    autoView( U     , Uds    , AcceleratorRead);

    typedef typename Field::vector_object        vobj;
    typedef decltype(coalescedRead(in[0]))    calcObj;
    typedef decltype(coalescedRead(U[0](0))) calcLink;

    const int      Nsimd = vobj::Nsimd();
    const uint64_t NN = grid->oSites();

    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	calcObj Uchi;
	calcLink UU;
	int ptype;

	res                 = coalescedRead(in[ss])*(-6.0);

#define LEG_LOAD_MULT(leg,polarisation)			\
	UU = coalescedRead(U[ss](polarisation));	\
	LEG_LOAD(leg);					\
	mult(&Uchi(), &UU, &chi());			\
	res = res + Uchi;
	
	LEG_LOAD_MULT(0,Xp);
	LEG_LOAD_MULT(1,Yp);
	LEG_LOAD_MULT(2,Zp);
	LEG_LOAD_MULT(3,Xm);
	LEG_LOAD_MULT(4,Ym);
	LEG_LOAD_MULT(5,Zm);

	coalescedWrite(out[ss], res,lane);
    });
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

#undef LEG_LOAD_MULT

template<class Gimpl,class Field> class CovariantAdjointLaplacianStencil : public SparseMatrixBase<Field>
{
public:
  INHERIT_GIMPL_TYPES(Gimpl);
//  RealD kappa;

  typedef typename Field::vector_object siteObject;

  template <typename vtype> using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Nc> >, Nds>;
  typedef iImplDoubledGaugeField<Simd> SiteDoubledGaugeField;
  typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;

  typedef CartesianStencil<siteObject, siteObject, int> StencilImpl;

  GridBase *grid;
  StencilImpl Stencil;
  SimpleCompressor<siteObject> Compressor;
  DoubledGaugeField Uds;
  CovariantAdjointLaplacianStencil(GaugeField &Umu)
    :
      grid(Umu.Grid()),
      Stencil    (grid,8,Even,directions4D,displacements4D,0),
      Uds(grid)
  {
    for (int mu = 0; mu < Nd; mu++) {
      auto U = PeekIndex<LorentzIndex>(Umu, mu);
      PokeIndex<LorentzIndex>(Uds, U, mu );
      U = adj(Cshift(U, mu, -1));
      PokeIndex<LorentzIndex>(Uds, U, mu + 4);
    }
  };
  
  virtual GridBase *Grid(void) { return grid; };

  virtual void  M    (const Field &_in, Field &_out)
  {
    ///////////////////////////////////////////////
    // Halo exchange for this geometry of stencil
    ///////////////////////////////////////////////
    Stencil.HaloExchange(_in, Compressor);

    ///////////////////////////////////
    // Arithmetic expressions
    ///////////////////////////////////
    auto st = Stencil.View(AcceleratorRead);
    auto buf = st.CommBuf();

    autoView( in     , _in    , AcceleratorRead);
    autoView( out    , _out   , AcceleratorWrite);
    autoView( U     , Uds    , AcceleratorRead);

    typedef typename Field::vector_object        vobj;
    typedef decltype(coalescedRead(in[0]))    calcObj;
    typedef decltype(coalescedRead(U[0](0))) calcLink;

    const int      Nsimd = vobj::Nsimd();
    const uint64_t NN = grid->oSites();

    accelerator_for( ss, NN, Nsimd, {

	StencilEntry *SE;
	
	const int lane=acceleratorSIMTlane(Nsimd);

	calcObj chi;
	calcObj res;
	calcObj Uchi;
	calcObj Utmp;
	calcObj Utmp2;
	calcLink UU;
	calcLink Udag;
	int ptype;

	res                 = coalescedRead(in[ss])*(-8.0);

#define LEG_LOAD_MULT(leg,polarisation)			\
	UU = coalescedRead(U[ss](polarisation));	\
	Udag = adj(UU);					\
	LEG_LOAD(leg);					\
	mult(&Utmp(), &UU, &chi());			\
	Utmp2 = adj(Utmp);				\
	mult(&Utmp(), &UU, &Utmp2());			\
	Uchi = adj(Utmp);				\
	res = res + Uchi;
	
	LEG_LOAD_MULT(0,Xp);
	LEG_LOAD_MULT(1,Yp);
	LEG_LOAD_MULT(2,Zp);
	LEG_LOAD_MULT(3,Tp);
	LEG_LOAD_MULT(4,Xm);
	LEG_LOAD_MULT(5,Ym);
	LEG_LOAD_MULT(6,Zm);
	LEG_LOAD_MULT(7,Tm);

	coalescedWrite(out[ss], res,lane);
    });
  };
  virtual void  Mdag (const Field &in, Field &out) { M(in,out);}; // Laplacian is hermitian
  virtual  void Mdiag    (const Field &in, Field &out)                  {assert(0);}; // Unimplemented need only for multigrid
  virtual  void Mdir     (const Field &in, Field &out,int dir, int disp){assert(0);}; // Unimplemented need only for multigrid
  virtual  void MdirAll  (const Field &in, std::vector<Field> &out)     {assert(0);}; // Unimplemented need only for multigrid
};

#undef LEG_LOAD_MULT
#undef LEG_LOAD
	
int main(int argc, char ** argv)
{
  Grid_init(&argc, &argv);

  typedef LatticeColourVector Field;

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();

  GridCartesian    Grid(latt_size,simd_layout,mpi_layout);
  GridParallelRNG  RNG(&Grid);  RNG.SeedFixedIntegers(std::vector<int>({45,12,81,9}));

  FreeLaplacianCshift<Field,3>  FLcs(&Grid);
  FreeLaplacianStencil<Field> FLst(&Grid);
// For adjoint testing  
  FreeLaplacianCshift<Field,4>  FLcs4(&Grid);

  LatticeGaugeField U(&Grid);

  SU<Nc>::ColdConfiguration(RNG,U);

  std::cout << " Gauge field has norm " <<norm2(U)<<std::endl;

  CovariantLaplacianCshift <PeriodicGimplR,Field> CLcs(U);
  CovariantLaplacianStencil<PeriodicGimplR,Field> CLst(U);
  CovariantAdjointLaplacianCshift <PeriodicGimplR,Field> CLcsA(U);
  CovariantAdjointLaplacianStencil<PeriodicGimplR,Field> CLstA(U);

  Field in(&Grid); gaussian(RNG,in);
  Field out_FLcs(&Grid);
  Field out_FLst(&Grid);
  Field out_CLcs(&Grid);
  Field out_CLst(&Grid);
  Field out_CLcsA(&Grid);
  Field out_CLstA(&Grid);
  Field out_FLcs4(&Grid);
  Field diff(&Grid);

  ////////////////////////////////////////////////////////
  // First test: in free field these should all agree
  ////////////////////////////////////////////////////////
  int Nwarm=10;
  int Napp=100;
  FLcs.M(in,out_FLcs);
  FLcs4.M(in,out_FLcs4);

  std:: cout << "******************************************************************" <<std::endl;
  std:: cout << " Test A: consistency of four different Laplacian implementations " <<std::endl;
  std:: cout << "******************************************************************" <<std::endl;
  std:: cout << " Input test vector " <<norm2(in)<<std::endl;
  std:: cout << "--------------------------------------------------------" <<std::endl;
  std:: cout << " Free cshift  output vector " <<norm2(out_FLcs)<<std::endl;
  std:: cout << " Free cshift  4d output vector " <<norm2(out_FLcs4)<<std::endl;

  for (int i=0;i<Nwarm;i++) 
  FLst.M(in,out_FLst);
  double dtime = -usecond();
  for (int i=0;i<Napp;i++) 
  FLst.M(in,out_FLst);
  dtime += usecond();
  std:: cout << " Free stencil output vector " <<norm2(out_FLst)<<" time= "<<dtime/Napp<<" us " <<std::endl;

  for (int i=0;i<Nwarm;i++) 
  CLcs.M(in,out_CLcs);
  dtime = -usecond();
  for (int i=0;i<Napp;i++) 
  CLcs.M(in,out_CLcs);
  dtime += usecond();
  std:: cout << " Cov  cshift  output vector " <<norm2(out_CLcs)<<" time= "<<dtime/Napp<<" us " <<std::endl;
  for (int i=0;i<Nwarm;i++) 
  CLst.M(in,out_CLst);
  dtime = -usecond();
  for (int i=0;i<Napp;i++) 
  CLst.M(in,out_CLst);
  dtime += usecond();
  std:: cout << " Cov  stencil output vector " <<norm2(out_CLst)<<" time= "<<dtime/Napp<<" us " <<std::endl;

  for (int i=0;i<Nwarm;i++) 
  CLcsA.M(in,out_CLcsA);
  dtime = -usecond();
  for (int i=0;i<Napp;i++) 
  CLcsA.M(in,out_CLcsA);
  dtime += usecond();
  std:: cout << " Adjoint Cov  cshift  output vector " <<norm2(out_CLcsA)<<" time= "<<dtime/Napp<<" us " <<std::endl;
  for (int i=0;i<Nwarm;i++) 
  CLstA.M(in,out_CLstA);
  dtime = -usecond();
  for (int i=0;i<Napp;i++) 
  CLstA.M(in,out_CLstA);
  dtime += usecond();
  std:: cout << " Adjoint Cov  stencil output vector " <<norm2(out_CLstA)<<" time= "<<dtime/Napp<<" us " <<std::endl;
  std:: cout << "--------------------------------------------------------" <<std::endl;
  

  diff = out_FLcs - out_FLst;
  std:: cout << " Difference between free Cshift Laplacian and free Stencil Laplacian      = " 
	  <<norm2(diff)<< std::endl;

  diff = out_FLcs - out_CLcs;
  std:: cout << " Difference between free Cshift Laplacian and covariant Cshift Laplacian  = " 
	  <<norm2(diff)<< std::endl;

  diff = out_FLcs - out_CLcs;
  std:: cout << " Difference between free Cshift Laplacian and covariant Cshift Laplacian = " <<norm2(diff)<<std::endl;
  diff = out_CLcs - out_CLst;
  std:: cout << " Difference between covariant Cshift Laplacian and covariant Stencil Laplacian = " <<norm2(diff)<<std::endl;

  diff = out_FLcs4 - out_CLcsA;
  std:: cout << " Difference between 4D free Cshift Laplacian and Adjoint covariant Stencil Laplacian = " <<norm2(diff)<<std::endl;
  diff = out_CLcsA - out_CLstA;
  std:: cout << " Difference between Adjoint covariant Cshift Laplacian and Adjoint covariant Stencil Laplacian = " <<norm2(diff)<<std::endl;

  std:: cout << "--------------------------------------------------------" <<std::endl;
  

  std:: cout << "******************************************************************" <<std::endl;
  std:: cout << " Test B: gauge covariance  " <<std::endl;
  std:: cout << "******************************************************************" <<std::endl;

  LatticeGaugeField     U_GT(&Grid); // Gauge transformed field
  LatticeColourMatrix   g(&Grid);   // local Gauge xform matrix

  U_GT = U;
  // Make a random xform to teh gauge field
  SU<Nc>::RandomGaugeTransform(RNG,U_GT,g); // Unit gauge

  Field in_GT(&Grid); 
  Field out_GT(&Grid);

  Field out_CLcs_GT(&Grid);
  Field out_CLst_GT(&Grid);
  Field out_CLcs_GTA(&Grid);
  Field out_CLst_GTA(&Grid);

  CovariantLaplacianCshift <PeriodicGimplR,Field> CLcs_GT(U_GT);
  CovariantLaplacianStencil<PeriodicGimplR,Field> CLst_GT(U_GT);
  CovariantAdjointLaplacianCshift <PeriodicGimplR,Field> CLcs_GTA(U_GT);
  CovariantAdjointLaplacianStencil<PeriodicGimplR,Field> CLst_GTA(U_GT);

  in_GT  = g*in;
  out_GT = g*out_FLcs;

  // Check M^GT_xy in_GT = g(x) M_xy g^dag(y) g(y) in = g(x) out(x)
  CLcs_GT.M(in_GT,out_CLcs_GT);
  CLst_GT.M(in_GT,out_CLst_GT);
  CLcs_GTA.M(in_GT,out_CLcs_GTA);
  CLst_GTA.M(in_GT,out_CLst_GTA);

  diff = out_CLcs_GT - out_GT;
  std:: cout << " Difference between Gauge xformed result and covariant Cshift Laplacian in xformed gauge  = " <<norm2(diff)<<std::endl;

  diff = out_CLst_GT - out_GT;
  std:: cout << " Difference between Gauge xformed result and covariant Stencil Laplacian in xformed gauge  = " <<norm2(diff)<<std::endl;

  diff = out_CLst_GTA - out_CLcs_GTA ;
  std:: cout << " Difference between Gauge xformed result and covariant Stencil Laplacian in xformed gauge  (Adjoint) = " <<norm2(diff)<<std::endl;

  std:: cout << "--------------------------------------------------------" <<std::endl;


  std:: cout << "******************************************************************" <<std::endl;
  std:: cout << " Test C: compare in free Field to \"Feynman rule\"  " <<std::endl;
  std:: cout << "******************************************************************" <<std::endl;

  std::vector<int> dim_mask({1,1,1,0}); // 3d FFT
  FFT theFFT(&Grid);
  Field out(&Grid);
  Field F_out(&Grid);
  Field F_in(&Grid);

  // FFT the random input vector
  theFFT.FFT_dim_mask(F_in,in,dim_mask,FFT::forward);

  // Convolution theorem: multiply by Fourier representation of (discrete) Laplacian to apply diff op
  LatticeComplexD    lap(&Grid); lap = Zero();
  LatticeComplexD    kmu(&Grid); 
  ComplexD ci(0.0,1.0);
  for(int mu=0;mu<3;mu++) {

    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];

    LatticeCoordinate(kmu,mu);
    kmu = TwoPiL * kmu;

    // (e^ik_mu + e^-ik_mu - 2) = 2( cos kmu - 1) ~ 2 (1 - k_mu^2/2 -1 ) =  - k_mu^2 + O(k^4)
    lap = lap + 2.0*cos(kmu) - 2.0;
    
  }  
  F_out = lap * F_in;
  
  // Inverse FFT the result
  theFFT.FFT_dim_mask(out,F_out,dim_mask,FFT::backward);
  
  std::cout<<"Fourier xformed (in)             "<<norm2(F_in)<<std::endl;
  std::cout<<"Fourier xformed Laplacian x (in) "<<norm2(F_out)<<std::endl;
  
  std::cout<<"Momentum space Laplacian application  "<< norm2(out)<<std::endl;
  std::cout<<"Stencil Laplacian application         "<< norm2(out_CLcs)<<std::endl;
    
  diff = out_CLcs - out;
  std::cout<<"diff "<< norm2(diff)<<std::endl;

  Grid_finalize();
}
