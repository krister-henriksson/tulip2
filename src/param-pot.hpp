


#ifndef PARAM_POT_HPP
#define PARAM_POT_HPP


#include "utils.hpp"
#include "utils-math.hpp"
#include "utils-vector.hpp"

#include "param.hpp"

#include "potinfo.hpp"


  
class ParamPot
  : public Param
{
public:
  PotentialInformationFit * p_potinfo;

  ParamPot();
  virtual ~ParamPot();

  ParamPot(const ParamPot & sv);
  ParamPot & operator=(const ParamPot & sv);

  // Other constructors:
  ParamPot(PotentialInformationFit * p_pi);

  virtual void update_pot();
  virtual void update_par();
  virtual utils::Vector<double> Xupdate(const utils::Vector<double> & xi);


} ;





#endif

