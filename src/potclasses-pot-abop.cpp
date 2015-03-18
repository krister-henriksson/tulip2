


#include <string>
#include <new>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "potclasses.hpp"
#include "param.hpp"


using namespace utils;



Potential_ABOP::Potential_ABOP()
  :
  elemname1("none"), elemname2("none"),
  parval(10,0), parname(10)
{
  parname[0] = "D0";
  parname[1] = "r0";
  parname[2] = "beta";
  parname[3] = "S";
  parname[4] = "gamma";
  parname[5] = "c";
  parname[6] = "d";
  parname[7] = "h";
  parname[8] = "R";
  parname[9] = "D";
}


void Potential_ABOP::init_lims(void){
  partype.resize(10);
  parmin.resize(10);
  parmax.resize(10);
  for (int i=0; i<10; ++i){
    partype[i] = PARAM_FIXED;
    parmin[i] = 1;
    parmax[i] = 1;
  }
}



