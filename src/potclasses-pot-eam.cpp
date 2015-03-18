


#include <string>
#include <new>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"
#include "potclasses.hpp"


using namespace utils;






Potential_EAM::Potential_EAM()
  : elemname1("none"), elemname2("none"),
    Nr(0), Nrho_s(0), Nrho_p(0), Nrho_d(0)
{}

double Potential_EAM::rcut(void){
  return r[Nr-1];
}

