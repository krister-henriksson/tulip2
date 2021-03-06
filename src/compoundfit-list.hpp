


#ifndef COMPOUNDFIT_LIST_HPP
#define COMPOUNDFIT_LIST_HPP



#include <iostream>
#include <string>

#include <cstdio>

#include "utils.hpp"
#include "utils-vector.hpp"
#include "utils-matrix.hpp"
#include "utils-matrix3.hpp"

//#include "elem-iacs.hpp"
//#include "specs-fit-prop-pot.hpp"
#include "mdsettings.hpp"
//#include "mdsystem.hpp"
//#include "param-pot.hpp"


#include "compoundfit.hpp"



#include "utils-vector3.hpp"
#include "utils-matrixsq3.hpp"


using utils::Vector3;
using utils::MatrixSq3;


using namespace utils;





// #####################################################################
// #####################################################################




class CompoundListFit {
public:
  Elements elem;
  int ncompounds;
  Vector<CompoundStructureFit> compounds;

  // Read structures:
  CompoundListFit(Elements & el,
		  MDSettings     & mds_specs,
		  std::string         filename);
  int NData();
} ;




#endif

