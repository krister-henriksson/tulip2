


#include <string>
#include <new>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "potclasses.hpp"

using namespace utils;





Potential_Reppot::Potential_Reppot()
  : elemname1("none"), elemname2("none"), Nr_rep(0)
{}

double Potential_Reppot::rcut(void){
  return -1.0;
}

