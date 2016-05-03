    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/simd/Grid_empty.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

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
//----------------------------------------------------------------------
/*! @file Grid_sse4.h
  @brief Empty Optimization libraries for debugging

  Using intrinsics
*/
// Time-stamp: <2015-06-09 14:28:02 neo>
//----------------------------------------------------------------------

namespace Grid {
namespace Optimization {
  typedef std::complex<double> dtype;
  typedef std::complex<float> ftype;

  struct Vsplat{
    //Complex float
    inline ftype operator()(float a, float b){
      return ftype(a,b);
    }
    // Real float
    inline ftype operator()(float a){
      return ftype(a,a);
    }
    //Complex double
    inline dtype operator()(double a, double b){
      return dtype(a,b);
    }
    //Real double
    inline dtype operator()(double a){
      return dtype(a,a);
    }
    //Integer
    inline int operator()(Integer a){
      return a;
    }
  };

  struct Vstore{
    //Float 
    inline void operator()(ftype a, float* F){
      memcpy(F,reinterpret_cast<float*>(&a),2*sizeof(float));
    }
    //Double
    inline void operator()(dtype a, double* D){
      memcpy(D,reinterpret_cast<double*>(&a),2*sizeof(double));
    }
    //Integer
    inline void operator()(int a, Integer* I){
      I[0] = a;
    }

  };

  struct Vstream{
    //Float
    inline void operator()(float * a, ftype b){
      memcpy(a,reinterpret_cast<float*>(&b),2*sizeof(float));
    }
    //Double
    inline void operator()(double * a, dtype b){
      memcpy(a,reinterpret_cast<double*>(&b),2*sizeof(double));
    }
  };

  struct Vset{
    // Complex float 
    inline ftype operator()(Grid::ComplexF *a){
      return ftype(a[0].real(),a[0].imag());
    }
    // Complex double 
    inline dtype operator()(Grid::ComplexD *a){
      return dtype(a[0].real(),a[0].imag());
    }
    // Real float 
    inline ftype operator()(float *a){
      return ftype(a[0],a[1]);
    }
    // Real double
    inline dtype operator()(double *a){
      return dtype(a[0],a[1]);
    }
    // Integer
    inline int operator()(Integer *a){
      return a[0];
    }


  };

  template <typename Out_type, typename In_type>
  struct Reduce{
    //Need templated class to overload output type
    //General form must generate error if compiled
    inline Out_type operator()(In_type in){
      printf("Error, using wrong Reduce function\n");
      exit(1);
      return 0;
    }
  };

  /////////////////////////////////////////////////////
  // Arithmetic operations
  /////////////////////////////////////////////////////
  struct Sum{
    //Complex/Real float
    inline ftype operator()(ftype a, ftype b){
      return a+b;
    }
    //Complex/Real double
    inline dtype operator()(dtype a, dtype b){
      return a+b;
    }
    //Integer
    inline int operator()(int a, int b){
      return a + b;
    }
  };

  struct Sub{
    //Complex/Real float
    inline ftype operator()(ftype a, ftype b){
      return a-b;
    }
    //Complex/Real double
    inline dtype operator()(dtype a, dtype b){
      return a-b;
    }
    //Integer
    inline int operator()(int a, int b){
      return a-b;
    }
  };

  struct MultComplex{
    // Complex float
    inline ftype operator()(ftype a, ftype b){
      return a*b;
    }
    // Complex double
    inline dtype operator()(dtype a, dtype b){
      return a*b;
    }
  };

  struct Mult{
    // Real float
    inline ftype operator()(ftype a, ftype b){
      return ftype(a.real()*b.real(),a.imag()*b.imag());
    }
    // Real double
    inline dtype operator()(dtype a, dtype b){
      return dtype(a.real()*b.real(),a.imag()*b.imag());
    }
    // Integer
    inline int operator()(int a, int b){
      return a*b;
    }
  };

  struct Conj{
    // Complex single
    inline ftype operator()(ftype in){
      return std::conj(in);
    }
    // Complex double
    inline dtype operator()(dtype in){
      return std::conj(in);
    }
    // do not define for integer input
  };

  struct TimesMinusI{
    //Complex single
    inline ftype operator()(ftype in, ftype ret){ //note ret is ignored
      return ftype(in.imag(),-in.real());
    }
    //Complex double
    inline dtype operator()(dtype in, dtype ret){ //note ret is ignored
      return dtype(in.imag(),-in.real());
    }
  };

  struct TimesI{
   //Complex single
    inline ftype operator()(ftype in, ftype ret){ //note ret is ignored
      return ftype(-in.imag(),in.real());
    }
    //Complex double
    inline dtype operator()(dtype in, dtype ret){ //note ret is ignored
      return dtype(-in.imag(),in.real());
    }
  };

  //////////////////////////////////////////////
  // Some Template specialization
  struct Permute{
    static inline ftype Permute0(ftype in){ //AB -> BA
      return ftype(in.imag(),in.real());
    };
    static inline ftype Permute1(ftype in){
      return in;
    };
    static inline ftype Permute2(ftype in){
      return in;
    };
    static inline ftype Permute3(ftype in){
      return in;
    };

    static inline dtype Permute0(dtype in){ //AB -> BA
      return dtype(in.imag(),in.real());
    };
    static inline dtype Permute1(dtype in){
      return in;
    };
    static inline dtype Permute2(dtype in){
      return in;
    };
    static inline dtype Permute3(dtype in){
      return in;
    };

  };
  
  template < typename vtype > 
    void permute(vtype &a, vtype b, int perm) {
   }; 

  //Complex float Reduce
  template<>
  inline Grid::ComplexF Reduce<Grid::ComplexF, ftype>::operator()(ftype in){ //1 complex
    return Grid::ComplexF(in.real(),in.imag());
  }
  //Real float Reduce
  template<>
  inline Grid::RealF Reduce<Grid::RealF, ftype>::operator()(ftype in){ //2 floats
    return in.real() + in.imag();
  }
  
  
  //Complex double Reduce
  template<>
  inline Grid::ComplexD Reduce<Grid::ComplexD, dtype>::operator()(dtype in){ //1 complex
    return Grid::ComplexD(in.real(),in.imag());
  }
  
  //Real double Reduce
  template<>
  inline Grid::RealD Reduce<Grid::RealD, dtype>::operator()(dtype in){ //2 doubles
    return in.real() + in.imag();
  }

  //Integer Reduce
  template<>
  inline Integer Reduce<Integer, int>::operator()(int in){
    // FIXME unimplemented
   printf("Reduce : Missing integer implementation -> FIX\n");
    assert(0);
  }
}

//////////////////////////////////////////////////////////////////////////////////////
// Here assign types 

  typedef Optimization::ftype SIMD_Ftype;  // Single precision type
  typedef Optimization::dtype SIMD_Dtype; // Double precision type
  typedef int SIMD_Itype; // Integer type

  // prefetch utilities
  inline void v_prefetch0(int size, const char *ptr){};
  inline void prefetch_HINT_T0(const char *ptr){};



  // Gpermute function
  template < typename VectorSIMD > 
    inline void Gpermute(VectorSIMD &y,const VectorSIMD &b, int perm ) {
    Optimization::permute(y.v,b.v,perm);
  }


  // Function name aliases
  typedef Optimization::Vsplat   VsplatSIMD;
  typedef Optimization::Vstore   VstoreSIMD;
  typedef Optimization::Vset     VsetSIMD;
  typedef Optimization::Vstream  VstreamSIMD;
  template <typename S, typename T> using ReduceSIMD = Optimization::Reduce<S,T>;

 


  // Arithmetic operations
  typedef Optimization::Sum         SumSIMD;
  typedef Optimization::Sub         SubSIMD;
  typedef Optimization::Mult        MultSIMD;
  typedef Optimization::MultComplex MultComplexSIMD;
  typedef Optimization::Conj        ConjSIMD;
  typedef Optimization::TimesMinusI TimesMinusISIMD;
  typedef Optimization::TimesI      TimesISIMD;

}
