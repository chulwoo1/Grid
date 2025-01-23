/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/action/ActionBase.h

Copyright (C) 2015-2016

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
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

#ifndef ACTION_BASE_H
#define ACTION_BASE_H

NAMESPACE_BEGIN(Grid);

///////////////////////////////////
// Smart configuration base class
///////////////////////////////////
template< class Field >
class ConfigurationBase
{
public:
  ConfigurationBase() {}
  virtual ~ConfigurationBase() {}
  virtual void set_Field(Field& U) =0;
  virtual void smeared_force(Field&) = 0;
  virtual Field& get_SmearedU() =0;
  virtual Field &get_U(bool smeared = false) = 0;
};

template <class GaugeField >
class Action 
{
public:
  bool is_smeared = false;
  bool keep_mom = false;
  RealD deriv_norm_sum;
  RealD deriv_max_sum;
  RealD Fdt_norm_sum;
  RealD Fdt_max_sum;
  int   deriv_num;
  RealD deriv_us;
  RealD S_us;
  RealD refresh_us;
  void  reset_timer(void)        {
    deriv_us = S_us = refresh_us = 0.0;
    deriv_norm_sum = deriv_max_sum=0.0;
    Fdt_max_sum =  Fdt_norm_sum = 0.0;
    deriv_num=0;
  }
  void  deriv_log(RealD nrm, RealD max,RealD Fdt_nrm,RealD Fdt_max) {
    if ( max > deriv_max_sum ) {
      deriv_max_sum=max;
    }
    deriv_norm_sum+=nrm;
    if ( Fdt_max > Fdt_max_sum ) {
      Fdt_max_sum=Fdt_max;
    }
    Fdt_norm_sum+=Fdt_nrm; deriv_num++;
  }
  RealD deriv_max_average(void)       { return deriv_max_sum; };
  RealD deriv_norm_average(void)      { return deriv_norm_sum/deriv_num; };
  RealD Fdt_max_average(void)         { return Fdt_max_sum; };
  RealD Fdt_norm_average(void)        { return Fdt_norm_sum/deriv_num; };
  RealD deriv_timer(void)        { return deriv_us; };
  RealD S_timer(void)            { return S_us; };
  RealD refresh_timer(void)      { return refresh_us; };
  void deriv_timer_start(void)   { deriv_us-=usecond(); }
  void deriv_timer_stop(void)    { deriv_us+=usecond(); }
  void refresh_timer_start(void) { refresh_us-=usecond(); }
  void refresh_timer_stop(void)  { refresh_us+=usecond(); }
  void S_timer_start(void)       { S_us-=usecond(); }
  void S_timer_stop(void)        { S_us+=usecond(); }
  /////////////////////////////
  // Heatbath?
  /////////////////////////////
  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG, RealD c1) = 0; // refresh pseudofermions
  virtual RealD S(const GaugeField& U) = 0;                             // evaluate the action
  virtual RealD Sinitial(const GaugeField& U) { return this->S(U); } ;  // if the refresh computes the action, can cache it. Alternately refreshAndAction() ?
  virtual void deriv(const GaugeField& U, GaugeField& dSdU) = 0;        // evaluate the action derivative
 
  /////////////////////////////////////////////////////////////
  // virtual smeared interface through configuration container
  /////////////////////////////////////////////////////////////
  virtual void refresh(ConfigurationBase<GaugeField> & U, GridSerialRNG &sRNG, GridParallelRNG& pRNG, RealD c1=0.)
  {
    refresh(U.get_U(is_smeared),sRNG,pRNG,c1);
  }
  virtual RealD S(ConfigurationBase<GaugeField>& U)
  {
    return S(U.get_U(is_smeared));
  }
  virtual RealD Sinitial(ConfigurationBase<GaugeField>& U) 
  {
    return Sinitial(U.get_U(is_smeared));
  }
  virtual void deriv(ConfigurationBase<GaugeField>& U, GaugeField& dSdU)
  {
    deriv(U.get_U(is_smeared),dSdU); 
    if ( is_smeared ) {
      U.smeared_force(dSdU);
    }
  }
  ///////////////////////////////
  // Logging
  ///////////////////////////////
  virtual std::string action_name()    = 0;                             // return the action name
  virtual std::string LogParameters()  = 0;                             // prints action parameters
  virtual ~Action(){}
//  static int traj_num;
//  static int field_num;
//  static int FieldNum(){return field_num++;}

};

template <class GaugeField >
class EmptyAction : public Action <GaugeField>
{
  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG, RealD c1=0.) { assert(0);}; // refresh pseudofermions
  using Action<GaugeField>::refresh;
  using Action<GaugeField>::Sinitial;
  using Action<GaugeField>::deriv;

//  virtual void refresh(const GaugeField& U, GridSerialRNG &sRNG, GridParallelRNG& pRNG) { assert(0);}; // refresh pseudofermions
  virtual RealD S(const GaugeField& U) { return 0.0;};                             // evaluate the action
  virtual void deriv(const GaugeField& U, GaugeField& dSdU) { assert(0); };        // evaluate the action derivative

  ///////////////////////////////
  // Logging
  ///////////////////////////////
  virtual std::string action_name()    { return std::string("Level Force Log"); };
  virtual std::string LogParameters()  { return std::string("No parameters");};
};

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

static uint64_t writeHeader(const uint32_t size, const uint32_t checksum, const std::string &format, const std::string &file){
  std::ofstream fout(file,std::ios::out|std::ios::in);
  fout.seekp(0,std::ios::beg);
  fout << std::setw(10) << size << std::endl;
  fout << std::hex << std::setw(10) << checksum << std::endl;
  fout << format << std::endl;
  return fout.tellp();
}

static uint64_t readHeader(uint32_t &size, uint32_t &checksum, std::string &format, const std::string &file){
  std::ifstream fin(file);
  std::string line;
  getline(fin,line);
  {
    std::stringstream ss; ss <<line ; ss >> size;
  }
  getline(fin,line);
  {
    std::stringstream ss; ss <<line ; ss >> std::hex >> checksum;
  }
  getline(fin,format);
  removeWhitespace(format);

  return fin.tellg();
}




NAMESPACE_END(Grid);

#endif // ACTION_BASE_H
